import logging
import os
import warnings
from importlib.metadata import version
from pathlib import Path

import click
import h5py
import numpy as np
from astropy.wcs import FITSFixedWarning

from dstools.logger import setupLogger
from dstools.utils import parse_coordinates

try:
    from dstools.imaging import get_pb_correction
    from dstools.ms import MeasurementSet, combine_spws, extract_baselines

    HAS_CASA_SUPPORT = True
except ImportError:
    HAS_CASA_SUPPORT = False

warnings.filterwarnings("ignore", category=FITSFixedWarning, append=True)

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-d",
    "--datacolumn",
    type=click.Choice(["data", "corrected", "model"]),
    default="data",
    help="Selection of DATA, CORRECTED_DATA, or MODEL_DATA column.",
)
@click.option(
    "-p",
    "--phasecentre",
    type=str,
    nargs=2,
    default=None,
    help=(
        "Coordinates of phasecentre at which to extract DS "
        "(provide as separate values, e.g. -p <RA> <DEC>)."
    ),
)
@click.option(
    "-P",
    "--primary-beam",
    type=Path,
    default=None,
    help=(
        "Path to primary beam image with which to correct flux scale. "
        "Must also provide phasecentre. Provide non-existent path to compute PB separately."
    ),
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
    "--baseline-average/--no-baseline-average",
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

    if not HAS_CASA_SUPPORT:
        logger.error("Dynamic spectrum extraction is not supported on this system.")
        raise SystemExit(1)

    ms = MeasurementSet(ms)

    columns = {
        "data": "DATA",
        "corrected": "CORRECTED_DATA",
        "model": "MODEL_DATA",
    }
    datacolumn = columns[datacolumn]

    # Check that selected column exists in MS
    if not ms.column_exists(datacolumn):
        logger.error(f"{ms} does not contain {datacolumn} column.")
        raise SystemExit(1)

    # Combine multiple spectral windows (e.g. VLA)
    # This also appears to fix an MS corrupted by model insertion
    # which has otherwise been very difficult to debug
    ms = combine_spws(ms)

    # Optionally rotate phasecentre to new coordinates
    if phasecentre is not None:
        position = parse_coordinates(phasecentre)
        ms = ms.rotate_phasecentre(position)

    # Get primary beam correction
    if primary_beam is not None and phasecentre is not None:
        pb_scale = get_pb_correction(ms, position, primary_beam)
    else:
        pb_scale = 1

    # Construct header with observation properties
    header = ms.header(
        datacolumn=datacolumn,
        pb_scale=pb_scale,
    )

    # Optionally average over baselines
    if baseline_average:
        ms = ms.average_baselines(minuvdist)

    # Construct 4D data and flag cubes on each baseline separately
    # to verify indices of missing data (e.g. due to correlator dropouts)
    visibilities, flags, uvdist = extract_baselines(ms, datacolumn)

    # Apply flags
    if not noflag:
        visibilities[flags] = np.nan

    # Apply primary beam correction
    visibilities /= header["pb_scale"]

    # Write all data to file
    with h5py.File(outfile, "w", track_order=True) as f:
        f.attrs["dstools_version"] = version("radio-dstools")
        for attr in header:
            f.attrs[attr] = header[attr]

        f.create_dataset("time", data=ms.times)
        f.create_dataset("frequency", data=ms.channels)
        f.create_dataset("uvdist", data=uvdist)
        f.create_dataset("flux", data=visibilities)

    # Clean up intermediate files
    os.system(f"rm -r {ms.path.parent}/*dstools-temp*.*ms 2>/dev/null")


if __name__ == "__main__":
    main()
