import itertools as it
import logging
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, wait
from pathlib import Path

import click
import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import FITSFixedWarning
from casaconfig import config
from numpy.typing import ArrayLike

config.logfile = "/dev/null"
from importlib.metadata import version

from casatools import table
from dstools.casa import mstransform, phaseshift
from dstools.imaging import get_pb_correction
from dstools.logger import setupLogger
from dstools.utils import column_exists, parse_coordinates

warnings.filterwarnings("ignore", category=FITSFixedWarning, append=True)

logger = logging.getLogger(__name__)


def average_baselines(ms: Path, minuvdist: float = 0) -> Path:
    logger.debug(f"Averaging over baseline axis with uvdist > {minuvdist}m")
    outputvis = ms.with_suffix(f".dstools-temp.baseavg{ms.suffix}")

    tab = table(str(ms), nomodify=False)

    ant1 = tab.getcol("ANTENNA1")
    ant2 = tab.getcol("ANTENNA2")

    # Set all antenna pairs equal for baseline averaging
    nrows = tab.nrows()
    tab.putcol("ANTENNA1", np.zeros(nrows))
    tab.putcol("ANTENNA2", np.ones(nrows))

    # Average over baselines by setting timeaverage interval to less than one scan cycle
    interval = tab.getcol("INTERVAL")
    timebin = "{}s".format(min(interval) * 1e-2)

    mstransform(
        vis=str(ms),
        outputvis=str(outputvis),
        datacolumn="all",
        uvrange=f">{minuvdist}m",
        timeaverage=True,
        timebin=timebin,
        keepflags=False,
    )

    # Replace original antenna names
    tab.putcol("ANTENNA1", ant1)
    tab.putcol("ANTENNA2", ant2)

    tab.unlock()
    tab.close()

    return outputvis


def get_header_properties(ms: Path, datacolumn: str, pb_scale: float) -> dict:

    # Calculate data dimensions
    times, freqs, antennas, nbaselines = get_data_dimensions(ms)

    # Infer polarisation of feeds
    feedtype = get_feed_polarisation(ms)

    # Get telescope name
    ta = table(f"{ms}/OBSERVATION", ack=False)
    telescope = str(ta.getcol("TELESCOPE_NAME")[0])
    ta.close()

    # Parse phasecentre direction
    ta = table(f"{ms}/FIELD", ack=False)
    phasecentre_coords = ta.getcol("PHASE_DIR")[:, 0, 0]
    phasecentre = SkyCoord(
        ra=phasecentre_coords[0],
        dec=phasecentre_coords[1],
        unit="rad",
    )
    ta.close()

    # Create header
    ncorrelations = nbaselines * len(freqs) * len(times) * 4
    header = {
        "telescope": telescope,
        "datacolumn": datacolumn,
        "feeds": feedtype,
        "antennas": len(antennas),
        "baselines": nbaselines,
        "integrations": len(times),
        "channels": len(freqs),
        "polarisations": 4,
        "correlations": ncorrelations,
        "phasecentre": phasecentre.to_string("hmsdms"),
        "pb_scale": pb_scale,
    }

    return header


def get_feed_polarisation(ms: Path) -> str:

    tf = table(f"{ms}/FEED", ack=False)
    feedtype = tf.getcol("POLARIZATION_TYPE")[0, 0]
    tf.close()

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


def combine_spws(ms: Path) -> Path:

    outvis = ms.with_suffix(f".dstools-temp.comb{ms.suffix}")

    # Determine number of spectral windows
    tab = table(f"{ms}/SPECTRAL_WINDOW")

    nspws = len(tab.getcol("FREQ_GROUP_NAME"))

    tab.unlock()
    tab.close()

    # Combine spectral windows if more than 1
    combine = nspws > 1
    mstransform(
        vis=str(ms),
        combinespws=combine,
        datacolumn="all",
        outputvis=str(outvis),
    )

    return outvis


def get_data_dimensions(ms: Path) -> tuple[ArrayLike, ArrayLike, ArrayLike, int]:

    # Get antenna count and time / frequency arrays
    tab = table(str(ms), ack=False)

    # Throw away autocorrelations
    tab = tab.query("ANTENNA1 != ANTENNA2")

    # Get time, frequency, and baseline axes
    times = np.unique(tab.getcol("TIME"))
    tf = table(f"{ms}/SPECTRAL_WINDOW", ack=False)
    freqs = tf.getcol("CHAN_FREQ").reshape(-1)

    antennas = np.unique(
        np.append(
            tab.getcol("ANTENNA1"),
            tab.getcol("ANTENNA2"),
        ),
    )

    tf.close()
    tab.close()

    # Calculate number of baselines
    nbaselines = len(antennas) * (len(antennas) - 1) // 2

    return times, freqs, antennas, nbaselines


def rotate_phasecentre(ms: Path, ra, dec) -> Path:
    logger.debug(f"Rotating phasecentre to {ra} {dec}")

    # Apply phasecentre rotation
    rotated_ms = ms.with_suffix(f".dstools-temp.rotated{ms.suffix}")

    phaseshift(
        vis=str(ms),
        outputvis=str(rotated_ms),
        phasecenter=f"J2000 {ra} {dec}",
    )

    return rotated_ms


def process_baseline(
    ms: Path,
    times: ArrayLike,
    baseline: tuple[str, str],
    datacolumn: str,
) -> dict:

    i, (ant1, ant2) = baseline

    tab = table(str(ms), ack=False)
    bl_tab = tab.query(f"(ANTENNA1=={ant1}) && (ANTENNA2=={ant2})")

    # Identify missing integrations on this baseline
    bl_time = bl_tab.getcol("TIME")
    missing_times = [t for t in times if t not in bl_time]

    # Add back to time column and identify indices of good integrations
    bl_time = np.sort(np.append(bl_time, missing_times))
    data_idx = np.argwhere(~np.in1d(bl_time, missing_times)).ravel()

    # Calculate UVrange for each baseline
    bl_uvw = bl_tab.getcol("UVW").T
    bl_uvdist = np.sqrt(np.sum(np.square(bl_uvw), axis=1))

    data = {
        "baseline": i,
        "data_idx": data_idx,
        "data": bl_tab.getcol(datacolumn).T,
        "flags": bl_tab.getcol("FLAG").T,
        "uvdist": np.nanmean(bl_uvdist, axis=0),
    }

    tab.close()
    bl_tab.close()

    return data


@click.command(context_settings={"show_default": True})
@click.option(
    "-d",
    "--datacolumn",
    type=click.Choice(["data", "corrected", "model"]),
    default="data",
    help="Selection of DATA, CORRECTED_DATA, or MODEL column.",
)
@click.option(
    "-p",
    "--phasecentre",
    type=str,
    nargs=2,
    default=None,
    help="Coordinates of phasecentre at which to extract DS (provide as separate values, e.g. -p <RA> <DEC>).",
)
@click.option(
    "-P",
    "--primary-beam",
    type=Path,
    default=None,
    help="Path to primary beam image with which to correct flux scale. Must also provide phasecentre.",
)
@click.option(
    "-F",
    "--noflag",
    is_flag=True,
    default=False,
    help="Remove flagging mask.",
)
@click.option(
    "-B",
    "--baseline-average",
    is_flag=True,
    default=True,
    help="Average over baseline axis.",
)
@click.option(
    "-u",
    "--minuvdist",
    type=float,
    default=0,
    help="Minimum UV distance in meters to retain if averaging over baseline axis.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms", type=Path)
@click.argument("outfile", type=Path)
def main(
    ms,
    outfile,
    datacolumn,
    phasecentre,
    primary_beam,
    noflag,
    baseline_average,
    minuvdist,
    verbose,
):

    setupLogger(verbose=verbose)

    columns = {
        "data": "DATA",
        "corrected": "CORRECTED_DATA",
        "model": "MODEL_DATA",
    }
    datacolumn = columns[datacolumn]

    # Check that selected column exists in MS
    if not column_exists(ms, datacolumn):
        logger.error(f"{ms} does not contain {datacolumn} column.")
        exit(1)

    # Combine multiple spectral windows (e.g. VLA)
    # This also appears to fix an MS corrupted by model insertion
    # which has otherwise been very difficult to debug
    ms = combine_spws(ms)

    # Optionally rotate phasecentre to new coordinates
    if phasecentre is not None:
        ra, dec = parse_coordinates(phasecentre)
        ms = rotate_phasecentre(ms, ra, dec)

    # Get primary beam correction
    pb_scale = get_pb_correction(primary_beam, ra, dec) if primary_beam else 1

    # Construct header with observation properties
    header = get_header_properties(ms, datacolumn, pb_scale)

    # Optionally average over baselines
    if baseline_average:
        ms = average_baselines(ms, minuvdist)

    # Calculate final dimensions of DS
    times, freqs, antennas, nbaselines = get_data_dimensions(ms)
    data_shape = (nbaselines, len(times), len(freqs), 4)

    # Initialise output arrays
    waterfall = np.full(data_shape, np.nan, dtype=complex)
    flags = np.full(data_shape, np.nan, dtype=bool)
    uvdist = np.full(nbaselines, np.nan)

    # Construct 4D data and flag cubes on each baseline separately
    # to verify indices of missing data (e.g. due to correlator dropouts)
    with ProcessPoolExecutor(max_workers=15) as executor:
        processes = [
            executor.submit(
                process_baseline,
                ms,
                np.copy(times),
                baseline,
                datacolumn,
            )
            for baseline in enumerate(it.combinations(antennas, 2))
        ]
        wait(processes)

    # Insert data into 4D data / flag cubes
    results = [p.result() for p in processes]
    for baseline in results:
        baseline_idx, data_idx = baseline["baseline"], baseline["data_idx"]
        waterfall[baseline_idx, data_idx] = baseline["data"]
        flags[baseline_idx, data_idx] = baseline["flags"]
        uvdist[baseline_idx] = baseline["uvdist"]

    # Apply flags
    if not noflag:
        waterfall[flags] = np.nan

    # Apply primary beam correction
    waterfall /= header["pb_scale"]

    # Write all data to file
    with h5py.File(outfile, "w", track_order=True) as f:

        f.attrs["dstools_version"] = version("radio-dstools")
        for attr in header:
            f.attrs[attr] = header[attr]

        f.create_dataset("time", data=times)
        f.create_dataset("frequency", data=freqs)
        f.create_dataset("uvdist", data=uvdist)
        f.create_dataset("flux", data=waterfall)

    # Clean up intermediate files
    ms_dir = ms.parent
    os.system(f"rm -r {ms_dir}/*dstools-temp*.ms 2>/dev/null")
    os.system("rm *.pre *.last 2>/dev/null")


if __name__ == "__main__":
    main()
