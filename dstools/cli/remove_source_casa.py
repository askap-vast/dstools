import click
import glob
import os
import subprocess
import astropy.units as u
from pathlib import Path

from dstools.utils import BANDS, CONFIGS, Array, prompt, update_param


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
    "-p",
    "--phasecenter",
    default="",
    help="Imaging phasecentre (in format J2000 hms dms).",
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
    "-a",
    "--automask/--no-automask",
    default=False,
    help="Use auto-multithresh mode to automatically generate clean mask.",
)
@click.option(
    "-k",
    "--killcoord",
    type=str,
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
    phasecenter,
    pblim,
    interactive,
    automask,
    killcoord,
    imsize,
    mpi,
):

    datapath = Path(data)
    proj_dir = str(datapath.parent.absolute()) + "/"
    source = str(data.split(".")[0])

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

    # Set cycleniter=niter if using both automasking and interactive mode
    # to ensure major cycles trigger at expected times
    usemask = "auto-multithresh" if automask else "user"
    cycleniter = iterations if automask and interactive else -1

    field_model_path = f"{proj_dir}/field_model/{band}/"
    kill_impath = f"{field_model_path}/offaxis_kill/"

    killcount = len(glob.glob(f"{field_model_path}/offaxis_kill/*image.tt0"))

    os.system(f"mkdir -p {kill_impath}")

    tclean(
        vis=data,
        field=field,
        cell=[cellsize],
        imsize=[imsize],
        threshold=threshold,
        niter=iterations,
        imagename=f"{kill_impath}/offaxis{killcount+1}",
        nterms=nterms,
        deconvolver="mtmfs",
        scales=clean_scales,
        reffreq=reffreq,
        weighting="briggs",
        stokes="IQUV",
        robust=robust,
        usemask=usemask,
        pbmask=0.0,
        cycleniter=cycleniter,
        gridder=gridder,
        wprojplanes=wprojplanes,
        phasecenter=killcoord,
        interactive=interactive,
        pblimit=pblim,
        parallel=mpi,
    )
    bgmodel = f"{kill_impath}/offaxis{killcount+1}.model"

    # Insert masked background model into visibilities and subtract
    tclean(
        vis=data,
        field=field,
        cell=[cellsize],
        imsize=[imsize],
        startmodel=[f"{bgmodel}.tt{i}" for i in range(nterms)],
        savemodel="modelcolumn",
        niter=0,
        imagename=f"{kill_impath}/offaxis{killcount+1}.im_presub",
        nterms=nterms,
        deconvolver="mtmfs",
        scales=clean_scales,
        reffreq=reffreq,
        weighting="briggs",
        stokes="IQUV",
        robust=robust,
        gridder=gridder,
        wprojplanes=wprojplanes,
        phasecenter=killcoord,
        pblimit=pblim,
        parallel=mpi,
    )
    os.system(f"rm -r {kill_impath}/offaxis{killcount+1}.im_presub* >/dev/null 2>&1")

    # Perform UV-subtraction.
    uvsub(vis=data)

    split(
        vis=data,
        outputvis=data.replace(".ms", f".offaxis{killcount+1}_subbed.ms"),
        datacolumn="corrected",
    )


if __name__ == "__main__":
    main()
