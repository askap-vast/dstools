import glob
import os
import subprocess
from pathlib import Path

import click

from dstools.utils import BANDS, CONFIGS, Array, parse_coordinates


@click.command()
@click.option(
    "-C",
    "--config",
    type=click.Choice(CONFIGS),
    default="6km",
    help="Array configuration, used to calculate image and pixel sizes. ASKAP is equivalent to 6km.",
)
@click.option(
    "-B",
    "--band",
    default="AT_L",
    type=click.Choice(BANDS),
    help="Observing band, used to calculate image and pixel sizes.",
)
@click.option(
    "-N",
    "--iterations",
    default=10000,
    help="Maximum number of clean iterations.",
)
@click.option(
    "-t",
    "--threshold",
    default=8e-5,
    help="Clean threshold in Jy.",
)
@click.option(
    "-n",
    "--nterms",
    default=3,
    help="Number of Taylor terms to use in imaging.",
)
@click.option(
    "-S",
    "--scale",
    default=[0, 4, 8],
    type=int,
    multiple=True,
    help="Multi-scale clean pixel smoothing scales.",
)
@click.option(
    "-w",
    "--widefield/--no-widefield",
    default=False,
    help="Enable w-projection to correct for widefield artefacts in off-axis sources.",
)
@click.option(
    "-r",
    "--robust",
    default=0.5,
    help="Briggs weighting robust parameter.",
)
@click.option(
    "-l",
    "--pblim",
    default=-0.1,
    help="Primary beam fractional power cutoff, default is no cutoff.",
)
@click.option(
    "-I",
    "--interactive/--no-interactive",
    default=True,
    help="Run tclean in interactive mode",
)
@click.option(
    "-k",
    "--killcoord",
    type=str,
    nargs=2,
    help="Coordinates of source to model and subtract.",
)
@click.option(
    "-i",
    "--imsize",
    default=400,
    help="Image size in pixels.",
)
@click.option(
    "-m",
    "--mpi",
    is_flag=True,
    default=False,
    help="Use parallel processing with mpi.",
)
@click.argument("data")
# @click.argument("killcoord", nargs=2)
def main(
    data,
    config,
    band,
    iterations,
    threshold,
    nterms,
    scale,
    widefield,
    robust,
    pblim,
    interactive,
    killcoord,
    imsize,
    mpi,
):

    # Set up directories
    # ------------------

    datapath = Path(data)
    proj_dir = str(datapath.parent.absolute()) + "/"
    offaxis_impath = f"{proj_dir}/field_model/offaxis_subtraction/"

    os.system(f"mkdir -p {offaxis_impath}")

    # Parameter Settings
    # ------------------

    array = Array(band, config)
    freq = array.frequency
    cell = array.cell
    imsize = int(imsize)

    field = "0"
    cellsize = f"{cell}arcsec"
    reffreq = f"{freq}MHz"
    clean_scales = list(scale)

    if widefield:
        gridder = "widefield"
        wprojplanes = -1
    else:
        gridder = "standard"
        wprojplanes = 1

    ra, dec = parse_coordinates(killcoord)
    killcoord = f"J2000 {ra} {dec}"

    killcount = len(glob.glob(f"{offaxis_impath}/*image.tt0"))

    # Generate offaxis source model
    # -----------------------------

    tclean(
        vis=data,
        field=field,
        cell=[cellsize],
        imsize=[imsize],
        imagename=f"{offaxis_impath}/offaxis{killcount+1}",
        nterms=nterms,
        niter=2000,
        deconvolver="mtmfs",
        scales=clean_scales,
        reffreq=reffreq,
        weighting="briggs",
        stokes="IQUV",
        robust=robust,
        pbmask=0.0,
        gridder=gridder,
        wprojplanes=wprojplanes,
        phasecenter=killcoord,
        interactive=interactive,
        pblimit=pblim,
        parallel=mpi,
    )

    # Subtract model
    # --------------

    scaleargs = " ".join([f"-S {scale}" for scale in clean_scales])
    mpi = "" if mpi else ""
    interactive = " --no-interactive" if not interactive else ""
    models = " ".join(glob.glob(f"{offaxis_impath}/*.model.tt*"))
    call = f"dstools-subtract-model -v{mpi}p {ra} {dec} {scaleargs}{interactive} {data} {models}".split(
        " "
    )
    subprocess.run(call)

    # Split to new MS
    # ---------------

    split(
        vis=data,
        outputvis=data.replace(".ms", f".offaxis{killcount+1}_subbed.ms"),
        datacolumn="corrected",
    )


if __name__ == "__main__":
    main()
