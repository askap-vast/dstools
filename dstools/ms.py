import itertools as it
import logging
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.time import Time
from casacore.tables import table, tablecopy
from matplotlib.gridspec import GridSpec
from numba import njit

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
from dstools.utils import (
    LOCATIONS,
    DataError,
    chunk_iterator,
    get_available_cpus,
    prompt,
)

logger = logging.getLogger(__name__)


def invert_caltable(caltable, inv_caltable):
    if os.path.exists(inv_caltable):
        raise FileExistsError(f"{inv_caltable} already exists")

    # Copy calibration table
    shutil.copytree(caltable.as_posix(), inv_caltable.as_posix())

    # Calculate inverse of gain solutions
    cal = Table(caltable.as_posix())
    C = cal.getcolumn("CPARAM")
    Cinv = 1.0 / C

    # Write to inverted calibration table
    inv_cal = Table(inv_caltable.as_posix())
    with inv_cal.open_table(readonly=False) as t:
        t.putcol("CPARAM", Cinv)

    return inv_caltable.as_posix()


@njit
def swap_xy_feeds(data: np.ndarray) -> np.ndarray:
    """Correct for mislabelling of the X and Y feeds (e.g. MeerKAT).

    Given the correlation Stokes parameter definitions:
    I = (XX + YY) / 2
    Q = (XX - YY) / 2
    U = (XY + YX) / 2
    V = 1j * (YX - XY) / 2

    swapping the X and Y feeds on all antennas has the effect:
    I ->  I
    Q -> -Q
    U ->  U
    V -> -V

    Here we swap the XX <-> YY and XY <-> YX correlations to correct for this.
    """

    nvis, nchan, _ = data.shape
    corrected = np.empty_like(data)

    for i in range(nvis):
        for f in range(nchan):
            corrected[i, f, 0] = data[i, f, 3]
            corrected[i, f, 1] = data[i, f, 2]
            corrected[i, f, 2] = data[i, f, 1]
            corrected[i, f, 3] = data[i, f, 0]

    return corrected


@njit
def rotate_circular_feeds(data: np.ndarray, chi_array: np.ndarray) -> np.ndarray:
    """Rotate polarisation basis from circularly polarised feed frame to sky frame.

    Circular feeds need a complex phase rotation by 2*chi (RL) and -2*chi (LR)
    in the cross-hand polarisations.

    Input shape: (nvis, nchan, 4), with pol order: (RR, RL, LR, LL)
    """
    nvis, nchan, _ = data.shape
    rotated = np.empty_like(data)
    for i in range(nvis):
        rot = np.exp(2j * chi_array[i])
        inv_rot = np.exp(-2j * chi_array[i])

        for f in range(nchan):
            rotated[i, f, 0] = data[i, f, 0]
            rotated[i, f, 1] = data[i, f, 1] * rot
            rotated[i, f, 2] = data[i, f, 2] * inv_rot
            rotated[i, f, 3] = data[i, f, 3]

    return rotated


@njit
def rotate_linear_feeds(
    data: np.ndarray,
    chi_array: np.ndarray,
) -> np.ndarray:
    """Rotate polarisation basis from linearly polarised feed frame to sky frame.

    Linear feeds require a real-valued rotation of the visibility matrix by 2*chi.

    Input shape: (nvis, nchan, 4), with pol order: (XX, XY, YX, YY)
    """
    nvis, nchan, _ = data.shape
    rotated = np.empty_like(data)

    for i in range(nvis):
        cos_chi = np.cos(chi_array[i])
        sin_chi = np.sin(chi_array[i])

        rot = np.array(
            [
                [cos_chi, sin_chi],
                [-sin_chi, cos_chi],
            ],
            dtype=np.complex128,
        )
        inv_rot = rot.T.conj()

        for f in range(nchan):
            V = np.empty((2, 2), dtype=np.complex128)
            V[0, 0] = data[i, f, 0]
            V[0, 1] = data[i, f, 1]
            V[1, 0] = data[i, f, 2]
            V[1, 1] = data[i, f, 3]

            V_rot = inv_rot @ V @ rot

            rotated[i, f, 0] = V_rot[0, 0]
            rotated[i, f, 1] = V_rot[0, 1]
            rotated[i, f, 2] = V_rot[1, 0]
            rotated[i, f, 3] = V_rot[1, 1]

    return rotated


@dataclass
class Table:
    """Base class for interacting with calibration tables and measurement sets."""

    path: Path | str

    def __post_init__(self):
        if isinstance(self.path, str):
            self.path = Path(self.path)

    @contextmanager
    def open_table(
        self,
        subtable: Optional[str] = None,
        readonly: bool = True,
        query: Optional[str] = None,
    ):
        path = self.path / subtable if subtable else self.path
        t = None

        try:
            t = table(
                path.as_posix(),
                readonly=readonly,
                ack=False,
            )
            if query is not None:
                t = t.query(query)
            yield t
        finally:
            if t is not None:
                t.unlock()
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

        antennas = np.unique(np.append(ant1, ant2))

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
        return self.gains.shape[2]

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
                        g = gains[np.where(self.spw_ids == spw), 0, polaxis]

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

        return fig


class MeasurementSet(Table):
    def __post_init__(self):
        super().__post_init__()

        if not os.path.exists(self.path):
            raise FileNotFoundError(f"MeasurementSet {self.path} not found")

        if not self.path.suffix == ".ms":
            raise ValueError(f"Path {self.path} is not a MeasurementSet!")

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
        return self.getcolumn("NUM_CORR", subtable="POLARIZATION")[0]

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
    def location(self):
        return LOCATIONS.get(self.telescope)

    @property
    def receptor_angle(self):
        with self.open_table(subtable="FEED") as t:
            angle = t.getcol("RECEPTOR_ANGLE")
        return angle[0, 0] * u.radian

    @property
    def feedtype(self):
        poltype_col = self.getcolumn("POLARIZATION_TYPE", subtable="FEED")
        poltype = poltype_col.get("array")[0]

        feedtype = {
            "X": "linear",
            "Y": "linear",
            "R": "circular",
            "L": "circular",
        }.get(poltype)

        if feedtype is None:
            raise ValueError(
                f"Feed has polarisation type {poltype} which cannot be recognised."
            )

        return feedtype

    @property
    def row_size_bytes(self):
        row_size = 16 * self.npols * self.nbaselines
        return row_size

    @property
    def phasecentre(self):
        phasecentre_coords = self.getcolumn("PHASE_DIR", subtable="FIELD")[0, 0, :]
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
            tablecopy(
                tablename=f"{self.original_path}/{ms_table}",
                newtablename=f"{one_spw_ms}/{ms_table}",
            )

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

    def swap_xy_feeds(self, datacolumn: str = "CORRECTED_DATA"):
        """Correct for mislabelling of the X and Y feeds (e.g. MeerKAT)."""

        logger.info("Correcting X/Y feed orientations")
        with self.open_table(readonly=False) as t:
            nrows = t.nrows()
            for startrow, chunk_size in chunk_iterator(nrows, self.row_size_bytes):
                data = t.getcol(datacolumn, startrow=startrow, nrow=chunk_size)

                corrected = swap_xy_feeds(data)

                t.putcol(datacolumn, corrected, startrow=startrow, nrow=chunk_size)

                del data, corrected

        # Set FEED reference angle to 0 now that we have fixed the positions
        with self.open_table(subtable="FEED", readonly=False) as t:
            receptor = t.getcol("RECEPTOR_ANGLE")
            receptor *= 0
            t.putcol("RECEPTOR_ANGLE", receptor)

        return

    def correct_feed_rotation(self, datacolumn: str = "CORRECTED_DATA"):
        """Apply corrections for parallactic angle rotation of feeds."""

        # Disable correction for ASKAP which has fixed sky-frame due to roll axis
        if self.telescope == "ASKAP":
            logger.warning(
                f"Correction for feed rotation not required for {self.telescope}. Will not apply."
            )
            return

        # Determine receptor angle offset from IAU standard (X - North, Y - East at zenith)
        logger.info("Correcting feed rotation by parallactic angle")

        with self.open_table(readonly=False) as t:
            observer = Observer(location=self.location)

            # Iterate through table in chunks to limit memory footprint
            nrows = t.nrows()
            for startrow, chunk_size in chunk_iterator(nrows, self.row_size_bytes):
                data = t.getcol(datacolumn, startrow=startrow, nrow=chunk_size)

                # Compute parallactic angle of each timesample
                mjd_sec = t.getcol("TIME", startrow=startrow, nrow=chunk_size)
                time = Time(mjd_sec * u.s.to(u.day), format="mjd", scale="utc")
                chi = observer.parallactic_angle(time, self.phasecentre).to(u.radian)

                # Offset parallactic angle by the receptor angle
                recep = self.receptor_angle
                logger.info(
                    f"Accounting for feed angle offset of {recep.to(u.deg):.0f}"
                )
                chi = (chi - recep).to_value(u.radian)

                # Apply parallactic angle corrections
                if self.feedtype == "linear":
                    data_rot = rotate_linear_feeds(data, chi)
                elif self.feedtype == "circular":
                    data_rot = rotate_circular_feeds(data, chi)

                # Write chunk back to table
                t.putcol(datacolumn, data_rot, startrow=startrow, nrow=chunk_size)

                # Free up memory
                del data, data_rot, mjd_sec, time, chi

        return

    def rotate_phasecentre(
        self,
        position: SkyCoord,
        inplace=False,
        threshold=0.1 * u.arcsec,
    ):
        if inplace:
            raise NotImplementedError(
                "Have not yet implemented inplace phasecentre rotation"
            )

        ra, dec = position.to_string(style="hmsdms").split()

        # Ensure new phasecentre differs from current phasecentre to avoid wasted processing
        if position.separation(self.phasecentre) < threshold:
            current_ra, current_dec = self.phasecentre.to_string(style="hmsdms").split()
            logger.debug(f"Phasecentre already set to {current_ra} {current_dec}")
            return self

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
        with self.open_table(readonly=False) as t:
            ant1 = t.getcol("ANTENNA1")
            ant2 = t.getcol("ANTENNA2")

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
        with self.open_table(readonly=False) as t:
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

        subtracted_ms = self.path.with_suffix(".subtracted.ms")
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

    def applycal(self, unapply: bool = False):
        if self.caltable is None:
            raise ValueError("Gain solutions not found, run .solve_gains() first!")

        if unapply:
            caltable = invert_caltable(
                caltable=self.caltable.path,
                inv_caltable=self.caltable.path.with_suffix(".inv.cal"),
            )
        else:
            caltable = self.caltable.path.as_posix()

        applycal(
            vis=self.path.as_posix(),
            gaintable=[caltable],
            interp="linear",
            applymode="calonly",
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


def run_selfcal(
    ms: MeasurementSet,
    calmode: str,
    gaintype: str,
    interval: str,
    split_data: bool,
    interactive: bool,
    refant: Optional[str] = None,
    nspws: int = 1,
) -> Path:
    """Perform self-calibration on MS with field model in the MODEL_DATA column."""

    if not re.fullmatch(r"(\d+)(s|min)", interval):
        raise ValueError(
            "Argument 'interval' must have format <int>[min/s] (e.g. 10s, 5min)."
        )

    if not ms.column_exists("MODEL_DATA"):
        raise DataError(f"{ms} does not contain a MODEL_DATA column.")

    # Produce MS with multiple spectral windows
    ms.to_nspws(nspws)

    # Solve for self calibration solutions
    ms.solve_gains(
        interval,
        calmode,
        gaintype,
        refant=refant,
        interactive=interactive,
    )

    # Generate phase and amplitude calibration plots
    for mode in calmode:
        ms.caltable.plot_solutions(mode)

    if interactive:
        plt.show(block=False)

    # Confirm solution is good before applying
    cal_good = prompt(
        msg="Apply gain solutions?",
        bypass=not interactive,
        default_response=True,
    )

    # If unacceptable, remove calibration tables, plots, and multi-spw MS and return
    if not cal_good:
        os.system(f"rm -r {ms.caltable.path}")
        os.system(f"rm {ms.path.stem}*.png")

        if nspws > 1:
            os.system(f"rm -r {ms.path} ")
            ms.path = ms.original_path
            ms.original_path = None

        return ms

    # Otherwise proceed with applying calibration solutions
    ms.applycal()

    # Transform back to single SPW MS
    ms.to_nspws(1)

    # Split out calibrated MS
    if split_data:
        ms = ms.split_selfcal_round()

    plt.close("all")

    return ms


def extract_baseline(
    ms: MeasurementSet,
    baseline: tuple[int, tuple[str, str]],
    datacolumn: str,
) -> dict:
    i, (ant1, ant2) = baseline

    with ms.open_table(query=f"(ANTENNA1=={ant1}) && (ANTENNA2=={ant2})") as bl_tab:
        # Identify missing integrations on this baseline (e.g. caused by correlator dropouts)
        bl_time = bl_tab.getcol("TIME")
        missing_times = [t for t in ms.times if t not in bl_time]

        # Add back to time column and identify indices of good integrations
        bl_time = np.sort(np.append(bl_time, missing_times))
        data_idx = np.argwhere(~np.isin(bl_time, missing_times)).ravel()

        # Calculate UVrange for each baseline
        bl_uvw = bl_tab.getcol("UVW")
        bl_uvdist = np.sqrt(np.sum(np.square(bl_uvw), axis=1))

        data_col = bl_tab.getcol(datacolumn)

        data = {
            "baseline": i,
            "data_idx": data_idx,
            "data": data_col,
            "flags": bl_tab.getcol("FLAG"),
            "uvdist": np.nanmean(bl_uvdist, axis=0),
        }

    return data


def get_polslice(pol_products: np.ndarray[str]) -> slice:
    """
    Get a slice object representing the polarisation axis indices that should contain the data.
    """

    # CORR_TYPE definitions from casacore StokesTypes enum
    CORRTYPES = {
        1: "I",
        2: "Q",
        3: "U",
        4: "V",
        5: "RR",
        6: "RL",
        7: "LR",
        8: "LL",
        9: "XX",
        10: "XY",
        11: "YX",
        12: "YY",
    }
    products = [CORRTYPES.get(product) for product in pol_products]

    # In DynamicSpectrum we expect the polarisation axis to be either
    # (RR, RL, LR, LL) or (XX, XY, YX, YY)
    # We also allow just the parallel hands or single polarisation,
    # and return the appropriate slice to avoid copying the data
    # across multiple polarisations
    match products:
        case ["XX", "XY", "YX", "YY"] | ["RR", "RL", "LR", "LL"]:
            polslice = slice(0, 4)
        case ["XX", "YY"] | ["RR", "LL"]:
            polslice = slice(0, 4, 3)
        case ["XX"] | ["RR"]:
            polslice = slice(0, 1)
        case _:
            raise DataError(f"Feed correlation types not supported: {products}")

    return polslice


def extract_baselines(
    ms: MeasurementSet,
    datacolumn: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # We always initialise all four instrumental pols whether or not they are
    # present, to ensure downstream processing can rely upon their presence
    dimensions = (ms.nbaselines, ms.integrations, ms.nchannels, 4)
    baselines = list(it.combinations(ms.antennas, 2))

    # Initialise output arrays
    visibilities = np.full(dimensions, np.nan + np.nan * 1j, dtype=complex)
    flags = np.full(dimensions, np.nan, dtype=bool)
    uvdist = np.full(ms.nbaselines, np.nan)

    # If more than 1 CPU available, use multiple processes to extract baselines in parallel
    ncpus = get_available_cpus()

    if ncpus > 1 and ms.nbaselines > 1:
        with ProcessPoolExecutor(max_workers=ncpus) as executor:
            processes = executor.map(
                extract_baseline,
                [ms] * ms.nbaselines,
                enumerate(baselines),
                [datacolumn] * ms.nbaselines,
            )
            results = [p for p in as_completed(processes)]
    else:
        results = [
            extract_baseline(ms, baseline, datacolumn)
            for baseline in enumerate(baselines)
        ]

    pol_products = ms.getcolumn("CORR_TYPE", subtable="POLARIZATION")[0]
    polslice = get_polslice(pol_products)
    for baseline in results:
        baseline_idx, data_idx = baseline["baseline"], baseline["data_idx"]
        visibilities[baseline_idx, data_idx, :, polslice] = baseline["data"]
        flags[baseline_idx, data_idx, :, polslice] = baseline["flags"]
        uvdist[baseline_idx] = baseline["uvdist"]

    return visibilities, flags, uvdist
