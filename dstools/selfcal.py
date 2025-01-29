import logging
import os
import re

import numpy as np
from casaconfig import config

config.logfile = "/dev/null"

from pathlib import Path

import astropy.units as u
import dask.array as da
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from astropy.time import Time
from casatools import table
from dstools.casa import applycal, cvel, flagdata, gaincal, mstransform, split
from dstools.utils import DataError, column_exists, prompt
from matplotlib.gridspec import GridSpec
from numpy.typing import ArrayLike

logger = logging.getLogger(__name__)


def split_multi_spw(ms, nspws):

    if nspws == 1:
        return ms, None

    logger.info(f"Transforming from 1 to {nspws} spectral windows")

    multi_spw_ms = ms.with_suffix(".multi-spw.ms")

    mstransform(
        vis=str(ms),
        outputvis=str(multi_spw_ms),
        regridms=True,
        nspw=nspws,
        mode="channel_b",
        datacolumn="all",
        combinespws=False,
        nchan=-1,
        start=0,
        width=1,
        chanbin=1,
        createmms=False,
    )

    # Replace original MS with multi-SPW copy
    return multi_spw_ms, ms


def combine_multi_spw(ms, multi_spw_ms, nspws):

    if nspws == 1:
        return ms

    logger.info(f"Transforming from {nspws} to 1 spectral windows")

    one_spw_ms = ms.with_suffix(".one-spw.ms")

    cvel(
        vis=str(ms),
        outputvis=str(one_spw_ms),
        mode="channel_b",
        nchan=-1,
        start=0,
        width=1,
    )

    # The FEED table is corrupted by mstransform / gaincal / applycal / cvel loop
    # growing in size by a factor of nspws and driving up run-time, so we copy
    # the original here and overwrite the output FEED table after each loop.
    feed_table = table(f"{multi_spw_ms}/FEED")
    feed_table.copy(f"{one_spw_ms}/FEED")
    feed_table.close()

    # Replace multi-SPW ms with combined copy
    os.system(f"rm -r {ms} {ms.with_suffix('.ms.flagversions')}")
    os.system(f"rm -r {multi_spw_ms}")
    os.system(f"mv {one_spw_ms} {multi_spw_ms}")
    ms = multi_spw_ms

    return ms


def increment_selfcal_round(ms: Path) -> Path:

    # Insert selfcal1 before suffix if first round
    if not re.match(r"\S*.selfcal\d*.ms", str(ms)):
        return ms.with_suffix(".selfcal1.ms")

    # Otherwise increment the round
    r = int(re.sub(r"\S*.selfcal(\d*).ms", r"\1", ms.name))
    round_name = ms.name.replace(f"selfcal{r}", f"selfcal{r+1}")

    return ms.with_name(round_name)


def plot_gain_solutions(
    caltable: Path,
    calmode: str,
    nspws: int,
) -> None:

    # Read time and gains from calibration table
    t = table(str(caltable))

    time = t.getcol("TIME")
    spw_ids = t.getcol("SPECTRAL_WINDOW_ID")
    gains = t.getcol("CPARAM")

    npols = gains.shape[0]
    nants = len(np.unique(t.getcol("ANTENNA1")))

    t.close()

    # Calculate phase / amplitude from complex gain solutions
    if calmode == "p":
        gains = xr.apply_ufunc(
            da.angle,
            gains,
            dask="allowed",
            kwargs=dict(deg=True),
        )
    elif calmode == "a":
        gains = da.absolute(gains)
    else:
        raise ValueError("Parameter 'calmode' must be 'p' or 'a'.")

    # Use figures with 2x3 subplots each
    num_figures = int(nants // 6)

    for subfig in range(num_figures):

        # Plot solutions against time
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(2, 3)

        # Color based on number of instrumental pols
        colors = ("k", "r")
        pol = ("X", "Y") if npols == 2 else ("X+Y",)

        for subplot in range(6):
            antaxis = subfig * 6 + subplot
            row = subplot // 3
            col = subplot % 3

            ax = fig.add_subplot(gs[row, col])

            for spw in range(nspws):

                for polaxis in range(npols):

                    t = time[np.where(spw_ids == spw)]
                    time_start = Time(
                        t[0] * u.s.to(u.day),
                        format="mjd",
                        scale="utc",
                    ).iso
                    t -= t[0]

                    # Select polarisation and current SPW
                    g = gains[polaxis, 0, np.where(spw_ids == spw)]

                    # Select current antenna
                    g = g.reshape(-1, nants)[:, antaxis]
                    t = t.reshape(-1, nants)[:, antaxis]

                    color = None if nspws > 1 else colors[polaxis]
                    label = None if nspws > 1 else pol[polaxis]

                    ax.scatter(
                        t / 3600,
                        g,
                        color=color,
                        s=1,
                        alpha=0.2,
                        label=label,
                    )

                    if nspws == 1:
                        ax.legend()

                    ax.set_xlabel(f"Hours from UTC {time_start}")
                    if calmode == "p":
                        maxval = np.abs(gains).max()
                        ax.set_ylabel("Phase [deg]")
                        ax.set_ylim(-maxval, maxval)
                    else:
                        maxval = gains.max() - 1
                        ax.set_ylabel("Amplitude")
                        ax.set_ylim(1 - 2 * maxval, 1 + 2 * maxval)

        fig.tight_layout()

        subfig = "" if subfig == 0 else subfig
        savefile = caltable.with_suffix(f".{calmode}.cal{subfig}.png")
        fig.savefig(savefile, format="png")

    return


def run_selfcal(
    ms: Path,
    calmode: str,
    gaintype: str,
    interval: str,
    split_data: bool,
    interactive: bool,
    refant: str = None,
    nspws: int = 1,
) -> Path:
    """Perform self-calibration on MS with field model in the MODEL_DATA column."""

    try:
        unit = "min" if "min" in interval else "s" if "s" in interval else ""
        int(interval.replace(unit, ""))
    except ValueError:
        raise ValueError(
            "Argument 'interval' must have format <int>[min/s] (e.g. 10s, 1min, 5min)."
        )

    if not column_exists(ms, "MODEL_DATA"):
        raise DataError(f"{ms} does not contain a MODEL_DATA column.")

    # Select reference antenna
    if refant is None:
        if interactive:

            # Calculate antenna flagging statistics
            flagstats = flagdata(vis=str(ms), mode="summary")
            df = pd.DataFrame(flagstats["antenna"]).T.reset_index(names="antenna")
            df["percentage"] = (100 * df.flagged / df.total).round(1)

            # Hide index to avoid confusing prompt
            df.index = [""] * len(df)

            print(f"Antenna flagging statistics:\n{df}")

            refant = input("Select reference antenna: ")
            while refant not in flagstats["antenna"]:
                print(f"Reference antenna must be in: {set(df.antenna)}")
                refant = input("Select reference antenna: ")
        else:
            refant = df.sort_values("percentage", ascending=True).iloc[0].antenna

    cal_table = ms.with_suffix(".cal")

    # Produce MS with multiple spectral windows
    ms, one_spw_ms = split_multi_spw(ms, nspws)

    gains = "phase" if calmode == "p" else "amp + phase"
    logger.info(f"Solving for {gains} over {nspws} spws and {interval} intervals")

    # Solve for self calibration solutions
    gaincal(
        vis=str(ms),
        caltable=str(cal_table),
        solint=interval,
        minblperant=3,
        refant=refant,
        calmode=calmode,
        gaintype=gaintype,
    )

    # Generate phase and amplitude calibration plots
    for mode in calmode:
        plot_gain_solutions(
            caltable=cal_table,
            calmode=mode,
            nspws=nspws,
        )

    if interactive:
        plt.show(block=False)

    # Confirm solution is good before applying
    cal_good = prompt(
        msg="Apply gain solutions?",
        bypass=not interactive,
        default_response=True,
    )

    # If unacceptable, remove calibration tables and multi-spw MS and return
    if not cal_good:
        if nspws > 1:
            os.system(f"rm -r {ms} {ms.with_suffix('ms.flagversions')}")
            ms = one_spw_ms

        os.system(f"rm -r {cal_table}")

        return ms

    # Otherwise proceed with applying calibration solutions
    applycal(
        vis=str(ms),
        gaintable=[str(cal_table)],
        interp="linear",
    )

    # Transform back to single SPW MS
    ms = combine_multi_spw(ms, one_spw_ms, nspws)

    # Split out calibrated MS
    if split_data:

        outms = increment_selfcal_round(ms)

        split(
            vis=str(ms),
            outputvis=str(outms),
            datacolumn="corrected",
        )
        ms = outms

    plt.close("all")

    return ms
