import logging
import os
from pathlib import Path

import click
from dstools.imaging import WSClean
from dstools.logger import setupLogger
from dstools.utils import BANDS, CONFIGS, Array, parse_coordinates

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-I",
    "--imsize",
    default=None,
    type=int,
    help="Image size in pixels.",
)
@click.option(
    "-c",
    "--cell",
    default=None,
    type=float,
    help="Cell / pixel size in arcseconds.",
)
@click.option(
    "-B",
    "--band",
    default="AT_L",
    type=click.Choice(BANDS),
    help="Observing band, used to calculate image and pixel sizes if unspecified.",
)
@click.option(
    "-C",
    "--config",
    type=click.Choice(CONFIGS),
    default="6km",
    help="Array configuration, used to calculate image and pixel sizes if unspecified. ASKAP is equivalent to 6km.",
)
@click.option(
    "-N",
    "--iterations",
    default=30000,
    help="Maximum number of clean iterations.",
)
@click.option(
    "-g",
    "--gain",
    default=0.8,
    help="Deconvolution minor cycle loop gain.",
)
@click.option(
    "-t",
    "--threshold",
    default=3,
    help="Clean threshold in multiples of RMS.",
)
@click.option(
    "-f",
    "--channels-out",
    default=8,
    help="Number of sub-band images to produce.",
)
@click.option(
    "-d",
    "--deconvolution-channels",
    default=8,
    help="Number of sub-bands over which to run deconvolution.",
)
@click.option(
    "-s",
    "--subimages",
    default=1,
    help="Number of subimage planes in each axis over which to run parallel deconvolution.",
)
@click.option(
    "-n",
    "--nterms",
    default=3,
    help="Number of Taylor terms to use in deconvolution.",
)
@click.option(
    "-u",
    "--minuvw",
    default=0,
    help="Minimum uv distance in meters.",
)
@click.option(
    "-m",
    "--mask-threshold",
    default=5,
    help="Automask threshold in multiples of RMS.",
)
@click.option(
    "-r",
    "--robust",
    default=0.5,
    help="Briggs weighting robust parameter.",
)
@click.option(
    "-p",
    "--phasecentre",
    type=str,
    nargs=2,
    default=None,
    help="Coordinates of imaging phasecentre (provide as separate values, e.g. -p <RA> <DEC>).",
)
@click.option(
    "-M",
    "--name",
    type=str,
    default="wsclean",
    help="Image prefix to provide to wsclean.",
)
@click.option(
    "-o",
    "--out_directory",
    type=Path,
    default="wsclean_model",
    help="Directory path in which to store WSclean image and model products.",
)
@click.option(
    "-v",
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms", type=Path)
def main(
    ms,
    imsize,
    cell,
    config,
    band,
    iterations,
    gain,
    minuvw,
    threshold,
    mask_threshold,
    channels_out,
    deconvolution_channels,
    subimages,
    nterms,
    robust,
    phasecentre,
    name,
    out_directory,
    verbose,
):

    setupLogger(verbose=verbose)

    # Set up and create working directories
    # --------------------------------------

    proj_dir = ms.parent.absolute()

    field_model_path = proj_dir / out_directory
    field_model_path.mkdir(exist_ok=True)

    # Set imaging parameters:
    # If cellsize / imsize not specified they will be estimated from array config / observing band
    # -----------------------

    array = Array(band, config)

    cell = cell if cell is not None else array.cell
    imsize = imsize if imsize is not None else array.imsize

    subimage_size = imsize // subimages
    cellsize = f"{cell}asec"

    if phasecentre:
        ra, dec = parse_coordinates(phasecentre)
        phasecentre = f"-shift {ra} {dec}"
    else:
        phasecentre = ""

    # Run WSclean
    # -----------------

    # Move into working directory to store imaging products
    ms = ms.absolute()
    os.chdir(field_model_path)

    wsc = WSClean(
        ms=ms,
        name=name,
        imsize=imsize,
        cellsize=cellsize,
        nterms=nterms,
        minuvw=minuvw,
        iterations=iterations,
        channels_out=channels_out,
        deconvolution_channels=deconvolution_channels,
        robust=robust,
        gain=gain,
        threshold=threshold,
        mask_threshold=mask_threshold,
        phasecentre=phasecentre,
        subimage_size=subimage_size,
        verbose=verbose,
    )

    wsc.run()

    # Return to start directory
    os.chdir(proj_dir)

    return


if __name__ == "__main__":
    main()
