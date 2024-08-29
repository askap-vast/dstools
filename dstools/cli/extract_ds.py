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
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from casacore.tables import table

import dstools
from dstools.logger import setupLogger
from dstools.utils import parse_coordinates

warnings.filterwarnings("ignore", category=FITSFixedWarning, append=True)

DSTOOLS_PATH = dstools.__path__[0]

logger = logging.getLogger(__name__)


def get_header_properties(ms, datacolumn, pb_scale):

    # Calculate data dimensions
    times, freqs, antennas, nbaselines = get_data_dimensions(ms)

    # Infer polarisation of feeds
    feedtype = get_feed_polarisation(ms)

    # Get telescope name
    ta = table(f"{ms}::OBSERVATION", ack=False)
    telescope = ta.getcol("TELESCOPE_NAME")[0]
    ta.close()

    # Parse phasecentre direction
    ta = table(f"{ms}::FIELD", ack=False)
    phasecentre_coords = ta.getcol("PHASE_DIR")[0][0]
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


def get_feed_polarisation(ms):

    tf = table(f"{ms}::FEED", ack=False)
    feedtype = tf.getcol("POLARIZATION_TYPE")["array"][0]
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


def combine_spws(ms):

    outvis = ms.replace(".ms", ".dstools-temp.comb.ms")
    os.system(f"_dstools-combine-spws {ms} {outvis} 1>/dev/null")

    return outvis


def get_data_dimensions(ms):

    # Get antenna count and time / frequency arrays
    tab = table(ms, ack=False, lockoptions="autonoread")

    # Throw away autocorrelations
    tab = tab.query("ANTENNA1 != ANTENNA2")

    # Get time, frequency, and baseline axes
    times = np.unique(tab.getcol("TIME"))
    tf = table(f"{ms}::SPECTRAL_WINDOW", ack=False)
    freqs = tf[0]["CHAN_FREQ"]

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


def validate_datacolumn(ms, datacolumn):

    tab = table(ms, ack=False, lockoptions="autonoread")
    col_exists = datacolumn in tab.colnames()
    tab.close()

    return col_exists


def rotate_phasecentre(ms, ra, dec):
    logger.debug(f"Rotating phasecentre to {ra} {dec}")

    # Apply phasecentre rotation
    os.system(f"_dstools-rotate -- {ms} {ra} {dec} 1>/dev/null")

    # Use rotated MeasurementSet for subsequent processing
    rotated_ms = ms.replace(".ms", ".dstools-temp.rotated.ms")

    return rotated_ms


def get_pb_correction(primary_beam, ra, dec):

    if primary_beam is None:
        return 1

    position = SkyCoord(ra=ra, dec=dec, unit=("hourangle", "deg"))

    if ".fits" not in primary_beam:
        pbfits = primary_beam + ".dstools.fits"
        cmd = f'exportfits(imagename="{primary_beam}", fitsimage="{pbfits}", overwrite=True)'
        os.system(f"casa --nologger --nologfile -c '{cmd}' 1>/dev/null")
        primary_beam = pbfits

    with fits.open(primary_beam) as hdul:
        header, data = hdul[0].header, hdul[0].data
        data = data[0, 0, :, :]

    if "dstools" in primary_beam:
        os.system(f"rm {pbfits} 2>/dev/null")

    wcs = WCS(header, naxis=2)
    x, y = wcs.wcs_world2pix(position.ra, position.dec, 1)
    x, y = int(x // 1), int(y // 1)
    xmax, ymax = data.shape

    # Check position is within limits of supplied PB image
    im_outside_limit = [
        x < 0,
        x > xmax,
        y < 0,
        y > ymax,
    ]
    if any(im_outside_limit):
        logger.warning(
            f"Position {ra} {dec} outside of supplied PB image, disabling PB correction."
        )
        return 1

    scale = data[x, y]

    logger.debug(
        f"PB correction scale {scale:.4f} measured at pixel {x},{y} in image of size {xmax},{ymax}"
    )

    return scale


def process_baseline(ms, times, baseline, datacolumn):

    i, (ant1, ant2) = baseline

    tab = table(ms, ack=False, lockoptions="autonoread")
    bl_tab = tab.query("(ANTENNA1==$ant1) && (ANTENNA2==$ant2)")

    # Identify missing integrations on this baseline
    bl_time = bl_tab.getcol("TIME")
    missing_times = [t for t in times if t not in bl_time]

    # Add back to time column and identify indices of good integrations
    bl_time = np.sort(np.append(bl_time, missing_times))
    data_idx = np.argwhere(~np.in1d(bl_time, missing_times)).ravel()

    # Calculate UVrange for each baseline
    bl_uvw = bl_tab.getcol("UVW")
    bl_uvdist = np.sqrt(np.sum(np.square(bl_uvw), axis=1))

    data = {
        "baseline": i,
        "data_idx": data_idx,
        "data": bl_tab.getcol(datacolumn),
        "flags": bl_tab.getcol("FLAG"),
        "uvdist": np.nanmean(bl_uvdist, axis=0),
    }

    tab.close()
    bl_tab.close()

    return data


@click.command()
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
    help="RA and Dec of phasecentre at which to extract DS.",
)
@click.option(
    "-P",
    "--primary-beam",
    type=click.Path(),
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
    help="Disable averaging over baseline axis.",
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
@click.argument("ms")
@click.argument("outfile")
def main(
    datacolumn,
    phasecentre,
    primary_beam,
    noflag,
    baseline_average,
    minuvdist,
    verbose,
    ms,
    outfile,
):

    setupLogger(verbose=verbose)

    columns = {
        "data": "DATA",
        "corrected": "CORRECTED_DATA",
        "model": "MODEL_DATA",
    }
    datacolumn = columns[datacolumn]

    # Combine multiple spectral windows (e.g. VLA)
    # This also appears to fix an MS corrupted by model insertion
    # which has otherwise been very difficult to debug
    ms = combine_spws(ms)

    # Check that selected column exists in MS
    if not validate_datacolumn(ms, datacolumn):
        logger.error(f"{datacolumn} column does not exist in {ms}")
        exit(1)

    # Optionally rotate phasecentre to new coordinates
    pb_scale = 1
    if phasecentre is not None:
        ra, dec = parse_coordinates(phasecentre)
        ms = rotate_phasecentre(ms, ra, dec)
        pb_scale = get_pb_correction(primary_beam, ra, dec)

    # Construct header with observation properties
    header = get_header_properties(ms, datacolumn, pb_scale)

    # Optionally average over baselines
    if baseline_average:
        logger.debug(f"Averaging over baseline axis with uvdist > {minuvdist}m")
        os.system(f"_dstools-avg-baselines -u {minuvdist} {ms} 1>/dev/null")
        ms = ms.replace(".ms", ".dstools-temp.baseavg.ms")

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
        for attr in header:
            f.attrs[attr] = header[attr]
        f.create_dataset("time", data=times)
        f.create_dataset("frequency", data=freqs)
        f.create_dataset("uvdist", data=uvdist)
        f.create_dataset("flux", data=waterfall)

    # Clean up intermediate files
    ms_dir = Path(ms).parent
    os.system(f"rm -r {ms_dir}/*dstools-temp*.ms 2>/dev/null")
    os.system("rm *.pre *.last 2>/dev/null")


if __name__ == "__main__":
    main()
