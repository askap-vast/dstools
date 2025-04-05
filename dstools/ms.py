import logging
import os
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import astropy.units as u
import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from astropy.coordinates import SkyCoord
from astropy.time import Time
from casatools.table import table
from matplotlib.gridspec import GridSpec

from dstools.casa import (
    applycal,
    cvel,
    flagdata,
    gaincal,
    mstransform,
    phaseshift,
    split,
    uvsub,
)
from dstools.utils import DataError

logger = logging.getLogger(__name__)


@dataclass
class Table:
    """Base class for interacting with calibration tables and measurement sets."""

    path: Path | str

    def __post_init__(self):
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"MeasurementSet {self.path} not found")

        if isinstance(self.path, str):
            self.path = Path(self.path)

    @contextmanager
    def open_table(
        self,
        subtable: Optional[str] = None,
        nomodify: bool = True,
        query: Optional[str] = None,
    ):
        path = f"{self.path}/{subtable}" if subtable else str(self.path)
        t = table()

        try:
            t.open(path, nomodify=nomodify)
            if query is not None:
                t = t.query(query)
            yield t
        finally:
            t.close()

    def getcolumn(self, column: str, subtable: Optional[str] = None):
        with self.open_table(subtable=subtable) as t:
            col = t.getcol(column)
        return col

    @property
    def antennas(self):
        # Throw away autocorrelations and get list of all antennas
        with self.open_table(query="ANTENNA1 != ANTENNA2") as t:
            ant1 = t.getcol("ANTENNA1")
            ant2 = t.getcol("ANTENNA2")

        antennas = np.unique(
            np.append(ant1, ant2),
        )

        # Remove flagged antennas with -1 index
        antennas = antennas[np.where(antennas != -1)]

        return antennas

    @property
    def nantennas(self):
        return len(self.antennas)

    @property
    def times(self):
        return self.getcolumn("TIME")


class CalTable(Table):
    """An interface to CASA calibration tables."""

    @property
    def gains(self):
        return self.getcolumn("CPARAM")

    @property
    def spw_ids(self):
        return self.getcolumn("SPECTRAL_WINDOW_ID")

    @property
    def nspws(self):
        return len(np.unique(self.spw_ids.flatten()))

    @property
    def npols(self):
        return self.gains.shape[0]

    @property
    def phase_solutions(self):
        phases = xr.apply_ufunc(
            da.angle,
            self.gains,
            dask="allowed",
            kwargs=dict(deg=True),
        )

        return phases

    @property
    def amp_solutions(self):
        return da.absolute(self.gains)

    def plot_solutions(self, calmode: str) -> None:
        # Calculate phase / amplitude from complex gain solutions
        if calmode == "p":
            gains = self.phase_solutions
        elif calmode == "a":
            gains = self.amp_solutions
        else:
            raise ValueError("Parameter 'calmode' must be 'p' or 'a'.")

        # Use figures with 2x3 subplots each
        num_figures = int(self.nantennas // 6)

        for subfig in range(num_figures):
            # Plot solutions against time
            fig = plt.figure(figsize=(12, 8))
            gs = GridSpec(2, 3)

            # Color based on number of instrumental pols
            colors = ("k", "r")
            pol = ("X", "Y") if self.npols == 2 else ("X+Y",)

            for subplot in range(6):
                antaxis = subfig * 6 + subplot
                row = subplot // 3
                col = subplot % 3

                ax = fig.add_subplot(gs[row, col])

                for spw in range(self.nspws):
                    for polaxis in range(self.npols):
                        t = self.times[np.where(self.spw_ids == spw)]
                        time_start = Time(
                            t[0] * u.s.to(u.day),
                            format="mjd",
                            scale="utc",
                        ).iso
                        t -= t[0]

                        # Select polarisation and current SPW
                        g = gains[polaxis, 0, np.where(self.spw_ids == spw)]

                        # Select current antenna
                        g = g.reshape(-1, self.nantennas)[:, antaxis]
                        t = t.reshape(-1, self.nantennas)[:, antaxis]

                        color = None if self.nspws > 1 else colors[polaxis]
                        label = None if self.nspws > 1 else pol[polaxis]

                        ax.scatter(
                            t / 3600,
                            g,
                            color=color,
                            s=1,
                            alpha=0.2,
                            label=label,
                        )

                        if self.nspws == 1:
                            ax.legend()

                        ax.set_xlabel(f"Hours from UTC {time_start}")
                        if calmode == "p":
                            maxval = np.abs(gains).max()
                            ax.set_ylabel("Phase [deg]")
                            ax.set_ylim(-maxval, maxval)
                        else:
                            maxval = np.abs(gains - 1).max()
                            ax.set_ylabel("Amplitude")
                            ax.set_ylim(1 - 2 * maxval, 1 + 2 * maxval)

            fig.tight_layout()

            subfig = "" if subfig == 0 else subfig
            savefile = self.path.with_suffix(f".{calmode}.cal{subfig}.png")
            fig.savefig(savefile, format="png")

        return


class MeasurementSet(Table):
    def __post_init__(self):
        super().__post_init__()
        self.original_path = None

    def __str__(self):
        return f"{self.path}"

    @property
    def nspws(self):
        spw_channels = self.getcolumn("NUM_CHAN", subtable="SPECTRAL_WINDOW")
        return len(spw_channels)

    @property
    def nbaselines(self):
        nants = self.nantennas
        return nants * (nants - 1) // 2

    @property
    def integrations(self):
        return len(self.times)

    @property
    def channels(self):
        chans = self.getcolumn("CHAN_FREQ", subtable="SPECTRAL_WINDOW")
        return chans.reshape(-1)

    @property
    def nchannels(self):
        return len(self.channels)

    @property
    def npols(self):
        return 4

    @property
    def dimensions(self):
        return (self.nbaselines, self.integrations, self.nchannels, self.npols)

    @property
    def ncorrelations(self):
        return np.prod(self.dimensions)

    @property
    def telescope(self):
        return str(self.getcolumn("TELESCOPE_NAME", subtable="OBSERVATION")[0])

    @property
    def feedtype(self):
        feedtype = self.getcolumn("POLARIZATION_TYPE", subtable="FEED")[0, 0]

        feedtype = {
            "X": "linear",
            "Y": "linear",
            "R": "circular",
            "L": "circular",
        }.get(feedtype)

        if feedtype is None:
            raise ValueError(
                f"Feed has polarisation type {feedtype} which cannot be recognised."
            )

        return feedtype

    @property
    def phasecentre(self):
        phasecentre_coords = self.getcolumn("PHASE_DIR", subtable="FIELD")[:, 0, 0]
        phasecentre = SkyCoord(
            ra=phasecentre_coords[0],
            dec=phasecentre_coords[1],
            unit="rad",
        )
        return phasecentre

    @property
    def colnames(self):
        with self.open_table() as t:
            columns = t.colnames()
        return columns

    def column_exists(self, column) -> bool:
        return column in self.colnames

    def header(self, datacolumn: str, pb_scale: float) -> dict:
        header = {
            "telescope": self.telescope,
            "datacolumn": datacolumn,
            "feeds": self.feedtype,
            "antennas": self.nantennas,
            "baselines": self.nbaselines,
            "integrations": self.integrations,
            "channels": self.channels,
            "polarisations": self.npols,
            "correlations": self.ncorrelations,
            "phasecentre": self.phasecentre.to_string("hmsdms"),
            "pb_scale": pb_scale,
        }

        return header

    def _combine_multi_spw(self) -> None:
        if self.original_path is None:
            raise DataError("This MS has not been split with .to_nspws()")

        logger.info(f"Transforming from {self.nspws} to 1 spectral windows")

        # Our current path is the multi-SPW copy of the data
        multi_spw_ms = self.path

        one_spw_ms = self.path.with_suffix(".1spw.ms")

        cvel(
            vis=str(multi_spw_ms),
            outputvis=str(one_spw_ms),
            mode="channel_b",
            nchan=-1,
            start=0,
            width=1,
        )

        # The FEED and SOURCE tables are corrupted by mstransform / gaincal / cvel loop
        # growing in size by a factor of nspws and driving up run-time, so we overwrite
        # the final FEED / SOURCE tables with those from the original MS
        for ms_table in ("FEED", "SOURCE"):
            original_table = table(f"{self.original_path}/{ms_table}")
            original_table.copy(newtablename=f"{one_spw_ms}/{ms_table}")
            original_table.close()

        # Remove original and multi-SPW copies
        os.system(f"rm -r {self.original_path}")
        os.system(f"rm -r {self.path}")

        # Replace original path with recombined copy and reset original_path state
        self.path = self.original_path
        self.original_path = None
        os.system(f"mv {one_spw_ms} {self.path}")

        return

    def _split_multi_spw(self, nspws):
        logger.info(f"Transforming from 1 to {nspws} spectral windows")

        multi_spw_ms = self.path.with_suffix(f".{nspws}spw.ms")

        mstransform(
            vis=str(self.path),
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
        self.original_path = self.path
        self.path = multi_spw_ms

        return

    def to_nspws(self, nspws: int) -> None:
        if nspws == 1 and self.nspws == 1:
            return

        if nspws == 1 and self.nspws > 1:
            self._combine_multi_spw()
            return

        if nspws > 1 and self.nspws == 1:
            self._split_multi_spw(nspws)
            return

        raise ValueError(
            f"Cannot transform from {self.nspws} to {nspws} spectral windows."
        )

    def rotate_phasecentre(self, ra, dec, inplace=False):
        if inplace:
            raise NotImplementedError(
                "Have not yet implemented inplace phasecentre rotation"
            )

        logger.debug(f"Rotating phasecentre to {ra} {dec}")

        # Apply phasecentre rotation
        rotated_ms = self.path.with_suffix(f".dstools-temp.rotated{self.path.suffix}")

        phaseshift(
            vis=str(self.path),
            outputvis=str(rotated_ms),
            phasecenter=f"J2000 {ra} {dec}",
        )

        return MeasurementSet(path=rotated_ms)

    def average_baselines(self, minuvdist: float = 0):
        logger.debug(f"Averaging over baseline axis with uvdist > {minuvdist}m")
        outputvis = self.path.with_suffix(f".dstools-temp.baseavg{self.path.suffix}")

        # Set antenna pairs equal to prepare for baseline averaging
        with self.open_table(nomodify=False) as t:
            ant1 = t.getcol("ANTENNA1")
            ant2 = t.getcol("ANTENNA2")

            # Set all antenna pairs equal for baseline averaging
            nrows = t.nrows()
            t.putcol("ANTENNA1", np.zeros(nrows))
            t.putcol("ANTENNA2", np.ones(nrows))

            # Average over baselines by setting timeaverage interval to less than one scan cycle
            interval = t.getcol("INTERVAL")
            timebin = "{}s".format(min(interval) * 1e-2)

            mstransform(
                vis=str(self.path),
                outputvis=str(outputvis),
                datacolumn="all",
                uvrange=f">{minuvdist}m",
                timeaverage=True,
                timebin=timebin,
                keepflags=False,
            )

            # Replace original antenna names
            # with self.open_table(nomodify=False) as t:
            t.putcol("ANTENNA1", ant1)
            t.putcol("ANTENNA2", ant2)

        return MeasurementSet(path=outputvis)

    def subtract_model(self, split_ms: bool = False):
        if not self.column_exists("MODEL_DATA"):
            raise DataError(
                f"{self.path} does not contain a MODEL_DATA column. Create or insert a model first!"
            )

        uvsub(str(self.path))

        if not split_ms:
            return self

        subtracted_ms = self.path.with_suffix(".subbed.ms")
        split(
            str(self.path),
            outputvis=str(subtracted_ms),
            datacolumn="corrected",
        )

        return MeasurementSet(subtracted_ms)

    def increment_selfcal_round(self) -> Path:
        # Insert selfcal1 before suffix if first round
        if not re.match(r"\S*.selfcal\d*.ms", str(self.path)):
            return self.path.with_suffix(".selfcal1.ms")

        # Otherwise increment the round
        r = int(re.sub(r"\S*.selfcal(\d*).ms", r"\1", self.path.name))
        round_name = self.path.name.replace(f"selfcal{r}", f"selfcal{r + 1}")

        return self.path.with_name(round_name)

    def calc_flag_statistics(self) -> pd.DataFrame:
        # Calculate antenna flagging statistics
        flagstats = flagdata(vis=str(self.path), mode="summary")
        df = pd.DataFrame(flagstats["antenna"]).T.reset_index(names="antenna")
        df["percentage"] = (100 * df.flagged / df.total).round(1)

        # Hide index to avoid confusing prompt
        df.index = [""] * len(df)

        print(f"Antenna flagging statistics:\n{df}")

        return df

    def get_reference_antenna(self, interactive: bool = False) -> str:
        flagstats = self.calc_flag_statistics()

        if interactive:
            refant = input("Select reference antenna: ")
            ant_options = sorted(set(flagstats.antenna))
            while refant not in ant_options:
                print(f"Reference antenna must be in: {ant_options}")
                refant = input("Select reference antenna: ")
        else:
            refant = flagstats.sort_values("percentage", ascending=True).iloc[0].antenna

        return refant

    def solve_gains(
        self,
        interval: str,
        calmode: str,
        gaintype: str,
        minblperant: int = 3,
        refant: Optional[str] = None,
        interactive: bool = False,
    ):
        # Select reference antenna
        if refant is None:
            refant = self.get_reference_antenna(interactive=interactive)

        gains = "phase" if calmode == "p" else "amp + phase"
        logger.info(
            f"Solving for {gains} over {self.nspws} spws and {interval} intervals"
        )

        caltable_path = self.path.with_suffix(".cal")
        gaincal(
            vis=str(self.path),
            caltable=str(caltable_path),
            solint=interval,
            calmode=calmode,
            gaintype=gaintype,
            minblperant=minblperant,
            refant=refant,
        )

        self.caltable = CalTable(path=caltable_path)

        return

    def split_selfcal_round(self):
        outms = self.increment_selfcal_round()

        split(
            vis=str(self.path),
            outputvis=str(outms),
            datacolumn="corrected",
        )

        return MeasurementSet(outms)

    def applycal(
        self,
        # gaintable: Path,
        # unapply: bool = False,
    ):
        # TODO: Implement applycal myself along with an unapply option
        applycal(
            vis=str(self.path),
            gaintable=[str(self.caltable.path)],
            interp="linear",
        )

        return


def combine_spws(ms: MeasurementSet) -> MeasurementSet:
    outvis = ms.path.with_suffix(f".dstools-temp.comb{ms.path.suffix}")

    # Combine spectral windows if more than 1
    logger.debug("Combining spectral windows")
    combine = ms.nspws > 1
    mstransform(
        vis=str(ms.path),
        combinespws=combine,
        datacolumn="all",
        outputvis=str(outvis),
    )

    return MeasurementSet(path=outvis)
