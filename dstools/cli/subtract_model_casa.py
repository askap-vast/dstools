import logging
import os
from pathlib import Path

import astropy.units as u
import click
import numpy as np

from dstools.logger import setupLogger
from dstools.utils import parse_coordinates

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "-p",
    "--phasecenter",
    type=str,
    nargs=2,
    default=None,
    help="Imaging phasecenter",
)
@click.option(
    "-r",
    "--robust",
    default=0.5,
    help="Briggs weighting robust parameter.",
)
@click.option(
    "-t",
    "--target_position",
    type=str,
    nargs=2,
    default=None,
    help="Coordinates around which to mask model images before subtraction.",
)
@click.option(
    "-I",
    "--interactive/--no-interactive",
    is_flag=True,
    default=True,
    help="Mask out target source in interactive tclean.",
)
@click.option(
    "-w",
    "--widefield/--no-widefield",
    default=False,
    help="Enable w-projection to correct for widefield artefacts in off-axis sources.",
)
@click.option(
    "-m",
    "--mpi",
    is_flag=True,
    default=False,
    help="Use parallel processing with mpi.",
)
@click.option(
    "-v",
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms")
@click.argument("model", nargs=-1)
def main(
    phasecenter,
    robust,
    target_position,
    interactive,
    widefield,
    mpi,
    verbose,
    ms,
    model,
):

    setupLogger(verbose=verbose)

    # Read image and cell sizes from model_base image properties
    # ------------------

    model_dir = str(Path(model[0]).parent)
    model_base = model[0].replace(".tt0", "")
    image_props = imhead(imagename=model[0])

    imsize = image_props["shape"][0]
    cell = round(abs(image_props["incr"][0]) * u.rad.to(u.arcsec), 3)
    freqaxis = np.argwhere(image_props["axisnames"] == "Frequency")
    freq = image_props["refval"][freqaxis][0, 0] * u.Hz.to(u.MHz)

    # Parameter Settings
    # ------------------

    field = "0"
    cellsize = f"{cell}arcsec"
    reffreq = f"{freq}MHz"
    pblim = -0.1
    nterms = len(model)

    if widefield:
        gridder = "widefield"
        wprojplanes = -1
    else:
        gridder = "standard"
        wprojplanes = 1

    if phasecenter:
        ra, dec = parse_coordinates(phasecenter)
        phasecenter = f"J2000 {ra} {dec}"
    else:
        phasecenter = ""

    # Create mask for target source
    # --------------------

    bgmodel = model_base.replace(".model", ".bgmodel")
    maskgen = model_base.replace(".model", ".maskgen")

    if not interactive and not target_position:
        source_mask = f"circle[[{imsize-1}pix, {imsize-1}pix], {1}pix]"
        logger.debug("No mask provided for target position. Subtracting full model.")
    elif target_position:
        ra, dec = parse_coordinates(target_position)
        source_mask = f"circle[[{ra}, {dec}], {10*cell}arcsec]"
        logger.debug(f"Masking target source from field model at {ra} {dec}")
    else:
        source_mask = model_base.replace(".model", ".source.mask")
        logger.debug("Masking target source from field model interactively")
        tclean(
            vis=ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            niter=1,
            imagename=maskgen,
            deconvolver="mtmfs",
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            interactive=True,
            phasecenter=phasecenter,
            pblimit=pblim,
            parallel=mpi,
        )
        os.system(f"mv {maskgen}.mask {source_mask} >/dev/null 2>&1")
        os.system(f"rm -r {model_dir}/*maskgen* >/dev/null 2>&1")

    os.system(f"rm -r {bgmodel}* >/dev/null 2>&1")

    # Apply mask to field model
    # --------------------------

    for tt in [f"tt{i}" for i in range(nterms)]:
        mask = f"{maskgen}.modelmask.{tt}"
        modelfile = f"{model_base}.{tt}"
        bgmodelfile = f"{bgmodel}.{tt}"

        makemask(
            mode="copy",
            inpimage=modelfile,
            output=mask,
            inpmask=source_mask,
            overwrite=True,
        )
        immath(
            imagename=[modelfile, mask],
            outfile=bgmodelfile,
            expr="IM0*(1-IM1)",
        )

    # Insert masked background model into visibilities and subtract
    # --------------------------------------------------------------

    tclean(
        vis=ms,
        field=field,
        cell=[cellsize],
        imsize=[imsize],
        startmodel=[f"{model_base}.tt{i}" for i in range(nterms)],
        savemodel="modelcolumn",
        niter=0,
        imagename=f"{maskgen}.im_presub",
        nterms=nterms,
        deconvolver="mtmfs",
        reffreq=reffreq,
        weighting="briggs",
        stokes="IQUV",
        robust=robust,
        gridder=gridder,
        wprojplanes=wprojplanes,
        phasecenter=phasecenter,
        pblimit=pblim,
        parallel=mpi,
    )
    os.system(f"rm -r {maskgen}.im_presub* >/dev/null 2>&1")

    # Perform field model subtraction
    # -----------------------

    uvsub(vis=ms)

    # Reimage to confirm field subtraction
    # -------------------------------------

    subbed_image = model_base.replace(".model", ".subbed")
    os.system(f"rm -r {subbed_image}* >/dev/null 2>&1")
    tclean(
        vis=ms,
        field=field,
        datacolumn="corrected",
        cell=[cellsize],
        imsize=[imsize],
        niter=0,
        imagename=subbed_image,
        nterms=nterms,
        deconvolver="mtmfs",
        reffreq=reffreq,
        weighting="briggs",
        stokes="IQUV",
        robust=robust,
        gridder=gridder,
        wprojplanes=wprojplanes,
        phasecenter=phasecenter,
        pblimit=pblim,
        pbcor=True,
        parallel=mpi,
    )

    # Export to FITS format
    # ----------------------

    for stokes in "IQUV":
        subimagename = f"{subbed_image}.image.{stokes}.tt0"
        os.system(f"rm -r {subimagename} >/dev/null 2>&1")

        imsubimage(
            imagename=f"{subbed_image}.image.tt0",
            outfile=subimagename,
            stokes=stokes,
        )
        exportfits(
            imagename=subimagename,
            fitsimage=f"{subimagename}.fits",
            overwrite=True,
        )

    return


if __name__ == "__main__":
    main()
