import logging
import os
from pathlib import Path

import click

from dstools.imaging import WSClean
from dstools.logger import setupLogger
from dstools.ms import MeasurementSet
from dstools.utils import BANDS, CONFIGS, Array

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
    "--mgain",
    default=0.8,
    help="Deconvolution major cycle loop gain.",
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
    "--spectral-pol-terms",
    default=3,
    help="Number of polynomial terms used to model spectral structure in MFS deconvolution.",
)
@click.option(
    "--minuvw_m",
    default=None,
    help="Minimum uv distance in meters.",
)
@click.option(
    "--minuvw_l",
    default=None,
    help="Minimum uv distance in wavelengths.",
)
@click.option(
    "-S",
    "--multiscale",
    is_flag=True,
    default=False,
    help="Enable multiscale deconvolution",
)
@click.option(
    "--multiscale-scale-bias",
    type=float,
    default=0.7,
    help="Deconvolution scale bias. Higher values give more weight to large scales.",
)
@click.option(
    "--multiscale-max-scales",
    type=int,
    default=8,
    help="Maximum number of multiscale scales to use in deconvolution.",
)
@click.option(
    "-M",
    "--fits-mask",
    default=None,
    type=Path,
    help="FITS image to use as deconvolution mask.",
)
@click.option(
    "-T",
    "--target-mask",
    default=None,
    type=Path,
    help="FITS image containing target emission pixels to exclude from mask.",
)
@click.option(
    "-G",
    "--galvin-clip-mask",
    default=None,
    type=Path,
    help="FITS image to use in minimum absolute clip mask.",
)
@click.option(
    "--erode-beam-shape",
    is_flag=True,
    default=False,
    help="Enable binary erosion of deconvolution mask with synthesised beam kernel.",
)
@click.option(
    "--local-rms-window",
    type=int,
    default=None,
    help="Window size in multiples of PSF over which to calculate local RMS.",
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
    "--name",
    type=str,
    default="wsclean",
    help="Image prefix to provide to wsclean.",
)
@click.option(
    "-o",
    "--out-dir",
    type=Path,
    default=Path("wsclean_model"),
    help="Directory path in which to store WSclean image and model products.",
)
@click.option(
    "--temp-dir",
    type=Path,
    help="Directory path in which to store temporary WSclean products.",
)
@click.option(
    "-L",
    "--savelogs",
    is_flag=True,
    default=False,
    help="Store processing logs.",
)
@click.option(
    "-j",
    "--threads",
    type=int,
    default=None,
    help="Number of CPU cores to use.",
)
@click.option(
    "--abs-mem",
    type=int,
    default=None,
    help="Maximum memory limit in gigabytes.",
)
@click.option(
    "--parallel-gridding",
    type=int,
    default=None,
    help="Execute this number of gridders in parallel.",
)
@click.option(
    "--parallel-reordering",
    type=int,
    default=None,
    help="Execute reordering with this number of threads.",
)
@click.option(
    "-v",
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms", type=MeasurementSet)
def main(
    ms,
    imsize,
    cell,
    config,
    band,
    iterations,
    mgain,
    minuvw_m,
    minuvw_l,
    multiscale,
    multiscale_scale_bias,
    multiscale_max_scales,
    fits_mask,
    target_mask,
    galvin_clip_mask,
    erode_beam_shape,
    threshold,
    mask_threshold,
    local_rms_window,
    channels_out,
    deconvolution_channels,
    subimages,
    spectral_pol_terms,
    robust,
    phasecentre,
    name,
    temp_dir,
    out_dir,
    savelogs,
    threads,
    abs_mem,
    parallel_gridding,
    parallel_reordering,
    verbose,
):
    os.system(f"mkdir -p {ms.path.parent.absolute() / out_dir}")
    logfile = (
        ms.path.parent.absolute() / out_dir / "create-model.log" if savelogs else None
    )
    setupLogger(verbose=verbose, filename=logfile)

    # Set imaging parameters:
    # If cellsize / imsize not specified they will be estimated from array config / observing band
    # -----------------------

    array = Array(band, config)

    cell = cell if cell is not None else array.cell
    imsize = imsize if imsize is not None else array.imsize

    parallel_deconvolution = imsize // subimages if subimages > 1 else None
    cellsize = f"{cell}asec"

    fits_mask = fits_mask.absolute() if fits_mask else None
    galvin_clip_mask = galvin_clip_mask.absolute() if galvin_clip_mask else None

    # Run WSclean
    # -----------------

    wsc = WSClean(
        imsize=imsize,
        cellsize=cellsize,
        spectral_pol_terms=spectral_pol_terms,
        minuvw_m=minuvw_m,
        minuvw_l=minuvw_l,
        multiscale=multiscale,
        multiscale_scale_bias=multiscale_scale_bias,
        multiscale_max_scales=multiscale_max_scales,
        iterations=iterations,
        channels_out=channels_out,
        deconvolution_channels=deconvolution_channels,
        robust=robust,
        mgain=mgain,
        auto_threshold=threshold,
        mask_threshold=mask_threshold,
        local_rms_window=local_rms_window,
        fits_mask=fits_mask,
        target_mask=target_mask,
        galvin_clip_mask=galvin_clip_mask,
        erode_beam_shape=erode_beam_shape,
        phasecentre=phasecentre,
        threads=threads,
        abs_mem=abs_mem,
        parallel_deconvolution=parallel_deconvolution,
        parallel_reordering=parallel_reordering,
        parallel_gridding=parallel_gridding,
        out_dir=out_dir,
        temp_dir=temp_dir,
        verbose=verbose,
    )

    wsc.run(ms=ms, name=name)

    return


if __name__ == "__main__":
    main()
