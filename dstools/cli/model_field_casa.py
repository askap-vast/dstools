import glob
import os
from pathlib import Path

import astropy.units as u
import click

from dstools.utils import BANDS, CONFIGS, Array, parse_coordinates, prompt, update_param


def validate_parameters(cell, imsize, nterms):
    while True:
        imradius = imsize * cell * u.arcsec.to(u.deg) / 2
        accept_params = prompt(
            f"Imaging with properties:\n  cell: {cell} arcsec\n  pixels: {imsize}\n  imradius: {imradius:.2f} deg\n  nterms: {nterms}\n\nProceed?"
        )

        if accept_params:
            break

        # Update parameters
        cell = update_param("cell", cell, float)
        imsize = update_param("imsize", imsize, int)
        nterms = update_param("nterms", nterms, int)

    return cell, imsize, nterms


def run_selfcal(ms):
    """Perform self-calibration on MS with field model in the MODEL_DATA column."""

    # TODO: Write check for existence of MODEL_DATA column

    path = str(Path(ms).parent.absolute())
    selfcal_round = len(glob.glob(f"{path}/*.cal")) + 1

    while True:
        # Select self calibration mode
        calmode_prompt = "Perform phase only (p) or phase and amplitude (ap) selfcal: "
        calmode = input(calmode_prompt)
        while calmode not in ["p", "ap"]:
            print("Selection must be one of either p or ap.")
            calmode = input(calmode_prompt)

        # Select gaintype
        gaintype_prompt = prompt("Combine polarisations in each solution interval?")
        gaintype = "T" if gaintype_prompt else "G"

        # Select refant
        flagstats = flagdata(vis=ms, mode="summary")
        refant = input("Select reference antenna: ")
        while refant not in flagstats["antenna"]:
            print(f"Reference antenna must be in: {flagstats['antenna']}")
            refant = input("Select reference antenna: ")

        # Select solution interval
        interval = input("Select solution interval (in min/s): ")
        try:
            unit = "min" if "min" in interval else "s" if "s" in interval else ""
            int(interval.replace(unit, ""))
        except ValueError:
            print("Solution interval must be in format <int>[min/s].")
            continue

        # Solve for self calibration solutions
        cal_table = ms.replace(".ms", f".round{selfcal_round}.cal")
        gaincal(
            vis=ms,
            caltable=cal_table,
            solint=interval,
            minblperant=3,
            refant=refant,
            calmode=calmode,
            gaintype=gaintype,
        )

        # Generate phase and amplitude calibration plots
        plotms(
            vis=cal_table,
            xaxis="time",
            yaxis="phase",
            iteraxis="antenna",
            coloraxis="corr",
            plotrange=[0, 0, -30, 30],
            showgui=True,
        )

        # Confirm solution is good before applying, else trial another solution
        cal_good = prompt("Is selfcal a good solution?")

        if cal_good:
            selfcal_ms = ms.replace(".ms", ".selfcal.ms")
            applycal(
                vis=ms,
                gaintable=[cal_table],
                interp="linear",
            )
            split(
                vis=ms,
                outputvis=selfcal_ms,
                datacolumn="corrected",
            )

            # Replace original MS with self-calibrated copy
            os.system(f"rm -r {ms}")
            os.system(f"mv {selfcal_ms} {ms}")
            break
        else:
            os.system(f"rm -r {cal_table}")

    return


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
    "--phasecentre",
    type=str,
    nargs=2,
    default=None,
    help="RA and Dec of phasecentre at which to extract DS.",
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
    "-m",
    "--mpi",
    is_flag=True,
    default=False,
    help="Use parallel processing with mpi.",
)
@click.option(
    "-c",
    "--copy/--no-copy",
    is_flag=True,
    default=True,
    help="Work on copy of MeasurementSet, preserving original. Disable to save on disk space.",
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
    phasecentre,
    pblim,
    interactive,
    automask,
    mpi,
    copy,
):

    # Set up and create working directories
    # --------------------------------------

    datapath = Path(data)
    proj_dir = str(datapath.parent.absolute()) + "/"
    source = str(datapath.parts[-1].split(".")[0])

    work_ms = f"{proj_dir}/{source}.{band}.ms"
    if copy:
        os.system(f"cp -r {data} {work_ms} >/dev/null 2>&1")
    else:
        os.system(f"mv {data} {work_ms} >/dev/null 2>&1")

    field_model_path = f"{proj_dir}/field_model/"
    os.system(f"mkdir -p {field_model_path}")

    # Set imaging parameters
    # -----------------------

    array = Array(band, config)
    freq = array.frequency
    cell = array.cell
    imsize = array.imsize
    imradius = array.imradius

    cell, imsize, nterms = validate_parameters(cell, imsize, nterms)

    field = "0"
    cellsize = f"{cell}arcsec"
    reffreq = f"{freq}MHz"
    clean_scales = list(scale)

    if phasecentre:
        ra, dec = parse_coordinates(phasecentre)
        phasecentre = f"J2000 {ra} {dec}"
    else:
        phasecentre = ""

    if widefield:
        gridder = "widefield"
        wprojplanes = -1
    else:
        gridder = "standard"
        wprojplanes = 1

    # Check for existing clean mask
    clean_mask = f"{proj_dir}/{source}.{band}.mask"
    clean_mask = clean_mask if os.path.exists(clean_mask) else ""

    # Set cycleniter=niter if using both automasking and interactive mode
    # to ensure major cycles trigger at expected times
    cycleniter = iterations if automask and interactive else -1
    usemask = "auto-multithresh" if automask else "user"
    lownoisethreshold = 2
    sidelobethreshold = 1.25

    # Flagging
    # ---------

    autoflag = prompt("Run automatic flagging?")
    if autoflag:
        flagdata(vis=work_ms, mode="rflag")
        flagdata(vis=work_ms, mode="tfcrop")

    # Clean / self-calibrate field
    # -----------------

    files = glob.glob(f"{field_model_path}/{source}.{band}.im[0-9].image.tt0")
    i = len(files)

    while True:

        i += 1

        # Check for a backed up mask in main directory
        imagename = f"{field_model_path}/{source}.{band}.im{i}"

        tclean(
            vis=work_ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            savemodel="modelcolumn",
            threshold=threshold,
            niter=10000,
            imagename=imagename,
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            robust=robust,
            stokes="IQUV",
            interactive=interactive,
            gridder=gridder,
            wprojplanes=wprojplanes,
            phasecenter=phasecentre,
            mask=clean_mask,
            usemask=usemask,
            sidelobethreshold=sidelobethreshold,
            lownoisethreshold=lownoisethreshold,
            cycleniter=cycleniter,
            pblimit=pblim,
            pbmask=0.0,
            pbcor=True,
            parallel=mpi,
        )

        # Self-calibrate with field model
        if prompt("Run self calibration on latest field model?"):
            run_selfcal(work_ms)

        if not prompt("Continue with another round of cleaning?"):
            break
        else:
            clean_mask = ""

    # Export to FITS format
    # ----------------------

    images = glob.glob(f"{field_model_path}/{source}.{band}.im[0-9].image.tt0")
    num_images = len(images)
    for stokes in "IQUV":
        subimagename = f"{field_model_path}/{source}.{band}.image.{stokes}.tt0"
        os.system(f"rm -r {subimagename} >/dev/null 2>&1")

        imsubimage(
            imagename=f"{field_model_path}/{source}.{band}.im{num_images}.image.tt0",
            outfile=subimagename,
            stokes=stokes,
        )
        exportfits(
            imagename=subimagename,
            fitsimage=f"{subimagename}.fits",
            overwrite=True,
        )

    os.system("rm *.last")

    return


if __name__ == "__main__":
    main()
