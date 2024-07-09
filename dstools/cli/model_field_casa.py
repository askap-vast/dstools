import glob
import os
import subprocess
from functools import partial

import astropy.units as u
import click
from dstools.utils import BANDS, CONFIGS, Array, prompt, update_param


def import_data(input_file, proj_dir, msname):
    os.system(f"mv {proj_dir}/{msname} {proj_dir}/{msname}.bak >/dev/null 2>&1")
    os.system(f"rm -r {proj_dir}/{msname}.flagversions >/dev/null 2>&1")

    # Convert Miriad visibilities to UVFITS format for reliable import
    if input_file.endswith(".cal"):
        os.system(f"fits in={input_file} out={input_file}.fits op=uvout")
        input_file = f"{input_file}.fits"

    if input_file.endswith(".fits"):
        importuvfits(
            fitsfile=input_file,
            vis=f"{proj_dir}/{msname}",
        )
    elif input_file.endswith(".ms"):
        os.system(f"cp -r {input_file} {proj_dir}/{msname} >/dev/null 2>&1")

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
    "-N", "--iterations", default=10000, help="Maximum number of clean iterations."
)
@click.option("-t", "--threshold", default=8e-5, help="Clean threshold in Jy.")
@click.option(
    "-n", "--nterms", default=3, help="Number of Taylor terms to use in imaging."
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
@click.option("-r", "--robust", default=0.5, help="Briggs weighting robust parameter.")
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
    "-P",
    "--proj_dir",
    type=click.Path(),
    default=None,
    help="Project namespace / sub-directory, default is to infer from data path",
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
    "-m", "--mpi", is_flag=True, default=False, help="Use parallel processing with mpi."
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
    proj_dir,
    interactive,
    automask,
    mpi,
):
    # PARAMETER SETTINGS
    # ------------------

    # If no project directory provided, use directory containing supplied data
    if proj_dir is None:
        proj_dir = "/".join(data.split("/")[:-1]) + "/"
        source = proj_dir.split("/")[-2].lower()
    else:
        source = "source"

    array = Array(band, config)
    freq = array.frequency
    cell = array.cell
    imsize = array.imsize
    imradius = array.imradius

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

    clean_mask = f"{proj_dir}/{source}.{band}.mask"

    # Set cycleniter=niter if using both automasking and interactive mode
    # to ensure major cycles trigger at expected times
    usemask = "auto-multithresh" if automask else "user"
    cycleniter = iterations if automask and interactive else -1

    # -------------------------------------------------

    while True:
        accept_params = prompt(
            f"Imaging with properties:\n  cell: {cell} arcsec\n  pixels: {imsize}\n  imradius: {imradius:.2f} deg\n  nterms: {nterms}\n\nProceed?"
        )

        if accept_params:
            break

        # Update parameters
        cell = update_param("cell", cell, float)
        cellsize = f"{cell}arcsec"

        imsize = update_param("imsize", imsize, int)
        imradius = imsize * cell * u.arcsec.to(u.deg) / 2
        nterms = update_param("nterms", nterms, int)

    os.system(f"mkdir -p {proj_dir}")

    # Data Import.
    # ------------
    msname = f"{source}.{band}.ms"
    calibrated_ms = f"{proj_dir}/{source}.{band}.selfcal.ms"

    if os.path.exists(f"{proj_dir}/{msname}"):
        reimport = prompt("Redo data import?")
    else:
        reimport = True

    if reimport:
        import_data(
            data,
            proj_dir,
            msname,
        )

    # Imaging Parameter Experimentation.
    # ----------------------------------

    test_params = prompt("Experiment with imsize / cellsize / nterms / flagging?")

    if test_params:
        os.system(f"rm -r {proj_dir}/test_params >/dev/null 2>&1")
        os.system(f"mkdir {proj_dir}/test_params")
        testms = f"{proj_dir}/test_params/testms.ms"
        os.system(f"cp -r {proj_dir}/{msname} {testms} >/dev/null 2>&1")

        while True:
            test_image = prompt("Make test image?")
            if test_image:
                # Optionally specify low-res and short timerange for quicker test imaging
                lowres = prompt("Test with low resolution?")
                testcell = f"{cell*5:2f}arcsec" if lowres else cellsize
                testimsize = imsize // 5 if lowres else imsize

                # Optionally specify short timerange for quicker test imaging
                listobs(vis=testms)
                timerange = input(
                    "Enter test imaging timerange (empty for full observation): "
                )

                tclean(
                    vis=testms,
                    field=field,
                    cell=[testcell],
                    imsize=[testimsize],
                    threshold=threshold,
                    niter=2000,
                    imagename=f"{proj_dir}/test_params/test_im",
                    nterms=nterms,
                    deconvolver="mtmfs",
                    scales=clean_scales,
                    reffreq=reffreq,
                    weighting="briggs",
                    robust=robust,
                    stokes="IQUV",
                    timerange=timerange,
                    interactive=interactive,
                    phasecenter=phasecenter,
                    pblimit=pblim,
                    parallel=mpi,
                )
                os.system(f"rm -r {proj_dir}/test_params/test_im*")

            # Run some additional flagging
            reflag = prompt("Perform flagging?")
            if reflag:
                manualflag = prompt("Run manual flagging?")
                if manualflag:
                    plotms(
                        vis=proj_dir + msname,
                        xaxis="channel",
                        yaxis="amp",
                        avgbaseline=False,
                    )

                    flagchan = prompt("Flag channel ranges?")
                    if flagchan:
                        flagvals = input(
                            "Specify channels to flag (; delimited, ~ range):"
                        )
                        flagdata(
                            vis=proj_dir + msname,
                            mode="manual",
                            spw="0:" + flagvals,
                        )

                    flagtime = prompt("Flag time ranges?")

                autoflag = prompt("Run automatic flagging?")
                if autoflag:
                    flagdata(vis=proj_dir + msname, mode="rflag")
                    flagdata(vis=proj_dir + msname, mode="tfcrop")

            # Update parameters
            imsize = update_param("imsize", imsize, int)
            pblim = update_param("pblim", pblim, float)
            cellsize = update_param("cellsize", cellsize, str)
            nterms = update_param("nterms", nterms, int)

            accept_params = prompt("Finished experimenting?")
            if accept_params:
                break

    os.system(f"rm -rf {proj_dir}/test_params >/dev/null 2>&1")

    # Self calibration.
    # -----------------

    selfcal = prompt("Run selfcal?")

    if selfcal:
        os.system(f"mkdir -p {proj_dir}/selfcal/{band}")

        files = glob.glob(f"{proj_dir}/selfcal/{band}/{source}.{band}.selfcal.[0-9].ms")
        i = len(files)
        cont = "y"
        while cont:
            i += 1
            selfcal_template = "{}/selfcal/{}/{}.selfcal_im{}"
            selfcal_im = selfcal_template.format(proj_dir, band, source, i)

            # Use freshly split ms each iteration after first loop
            selfcal_ms = (
                proj_dir
                + f"selfcal/{band}/"
                + msname.replace(".ms", f".selfcal.{i-1}.ms")
            )
            selfcal_ms = selfcal_ms if os.path.exists(selfcal_ms) else proj_dir + msname

            # Set mask from previous iteration as input to next tclean loop
            selfcal_mask = (
                selfcal_template.format(proj_dir, band, source, i - 1) + ".mask"
            )

            # If no selfcal mask exists, check for a backed up mask in main directory
            if os.path.exists(selfcal_mask):
                init_mask = selfcal_mask
            else:
                init_mask = clean_mask if os.path.exists(clean_mask) else ""

            tclean(
                vis=selfcal_ms,
                field=field,
                cell=[cellsize],
                imsize=[imsize],
                savemodel="modelcolumn",
                threshold=threshold,
                niter=2000,
                imagename=selfcal_im,
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
                phasecenter=phasecenter,
                mask=init_mask,
                usemask=usemask,
                cycleniter=cycleniter,
                pblimit=pblim,
                parallel=mpi,
            )

            # Trial selfcal solution intervals
            while True:
                # Select self calibration mode
                calmode_prompt = (
                    "Perform phase only (p) or phase and amplitude (ap) selfcal: "
                )
                calmode = input(calmode_prompt)
                while calmode not in ["p", "ap"]:
                    print("Selection must be one of either p or ap.")
                    calmode = input(calmode_prompt)

                # Select gaintype
                gaintype_prompt = prompt(
                    "Combine polarisations in each solution interval?"
                )
                gaintype = "T" if gaintype_prompt else "G"

                # Select refant
                flagstats = flagdata(vis=selfcal_ms, mode="summary")
                refant = input("Select reference antenna: ")
                while refant not in flagstats["antenna"]:
                    print(f"Reference antenna must be in: {flagstats['antenna']}")
                    refant = input("Select reference antenna: ")

                # Select solution interval
                interval = input("Select solution interval (in min/s): ")
                try:
                    unit = (
                        "min" if "min" in interval else "s" if "s" in interval else ""
                    )
                    int(interval.replace(unit, ""))
                except ValueError:
                    print("Solution interval must be in format <int>[min/s].")
                    continue

                # Solve for self calibration solutions
                cal_file = f"{proj_dir}/selfcal/{band}/{source}.{calmode}_selfcal_{interval}.{i}"
                cal_table = cal_file + ".cal"
                gaincal(
                    vis=selfcal_ms,
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
                    subprocess.run(
                        f"ragavi-gains -t {cal_table} -o {cal_file}", shell=True
                    )
                    applycal(vis=selfcal_ms, gaintable=[cal_table], interp="linear")
                    split(
                        vis=selfcal_ms,
                        outputvis=proj_dir
                        + f"selfcal/{band}/"
                        + msname.replace(".ms", f".selfcal.{i}.ms"),
                        datacolumn="corrected",
                    )
                    break
                else:
                    os.system(f"rm -r {cal_table}")

            # Break from loop once solution is sufficient
            cont = prompt("Proceed with more selfcal?")

        # Backup clean mask to main project directory
        backup_mask = selfcal_template.format(proj_dir, band, source, i) + ".mask"
        if os.path.exists(clean_mask):
            replace_mask = prompt("Update the existing clean mask with current mask?")

            if replace_mask:
                os.system(f"rm -r {clean_mask}")
                os.system(f"cp -r {backup_mask} {clean_mask}")
        else:
            os.system(f"cp -r {backup_mask} {clean_mask}")

        # Copy calibrated MS to root folder
        os.system(f"rm -r {calibrated_ms} >/dev/null 2>&1")
        os.system(f"cp -r {selfcal_ms} {calibrated_ms}")

    else:
        # If restarting a run with selfcal already complete, use existing calibrated MS.
        # Otherwise use the most up-to-date selfcal version
        if not os.path.exists(calibrated_ms):
            selfcal_path = f"{proj_dir}/selfcal/{band}/{source}.{band}.selfcal.[0-9].ms"
            files = glob.glob(selfcal_path)
            i = len(files)
            if i == 0:
                # If selfcal was skipped, use the original MS
                selfcal_ms = proj_dir + msname
            else:
                selfcal_ms = sorted(files)[-1]

            os.system(f"cp -r {selfcal_ms} {calibrated_ms}")

            # Check for existing clean mask
            if not os.path.exists(clean_mask):
                clean_mask = ""

    # Preparation for Deep Clean
    # --------------------------

    field_model_path = f"{proj_dir}/field_model/{band}/"
    deep_mask = clean_mask

    if os.path.exists(field_model_path):
        # If continuing with an existing model, we should remove
        # the mask parameter to ensure the existing clean mask is used
        clear_model = prompt("Start from fresh field model?")
        if clear_model:
            os.system(f"mv {field_model_path} {field_model_path[:-1]}_backup")

        keep_mask = prompt("Continue with existing clean mask?")
        if keep_mask:
            deep_mask = ""

    os.system(f"mkdir -p {field_model_path}")

    # Deep clean to produce field model.
    # ----------------------------------
    # This is within an endless loop to allow iterative cleaning without babysitting.
    # Set a conservative clean mask and let it run until stopping criteria met,
    # then assess whether further cleaning is required with an updated mask.

    deep_clean = prompt("Perform deep clean?")

    work_ms = calibrated_ms.replace(".ms", ".subbed.ms")
    if not os.path.exists(work_ms):
        os.system(f"cp -r {calibrated_ms} {work_ms}")

    while deep_clean:
        clean_round = len(glob.glob(f"{field_model_path}/*im_deep*.image.tt0"))
        imname = f"im_deep{clean_round}"

        tclean(
            vis=work_ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            threshold=threshold,
            niter=iterations,
            imagename=f"{field_model_path}/{source}.{imname}",
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            mask=deep_mask,
            usemask=usemask,
            pbmask=0.0,
            cycleniter=cycleniter,
            gridder=gridder,
            wprojplanes=wprojplanes,
            phasecenter=phasecenter,
            interactive=interactive,
            pblimit=pblim,
            parallel=mpi,
        )

        cont = prompt("Continue with further cleaning?")
        if not cont:
            deep_mask = f"{field_model_path}/{source}.{imname}.mask"
            os.system(f"rm -r {clean_mask} >/dev/null 2>&1")
            os.system(f"cp -r {deep_mask} {clean_mask} >/dev/null 2>&1")

            break
        else:
            deep_mask = ""

    # Use latest image / model for remaining processing
    clean_round = len(glob.glob(f"{field_model_path}/*im_deep*.image.tt0")) - 1
    imname = f"im_deep{clean_round}"

    # Mask out the source.
    # --------------------
    # Either use a generic circular mask at image center or specify via a custom tclean mask

    model = f"{field_model_path}/{source}.{imname}.model"
    bgmodel = f"{field_model_path}/{source}.{imname}.bgmodel"

    if not prompt("Mask out region of field model?"):
        source_mask = "circle[[{imsize-1}pix, {imsize-1}pix], {1}pix]"
    else:
        tclean(
            vis=work_ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            threshold=threshold,
            niter=1,
            imagename=f"{field_model_path}/{source}.maskgen",
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            interactive=True,
            phasecenter=phasecenter,
            pblimit=pblim,
            parallel=mpi,
        )
        source_mask = f"{field_model_path}/{source}.source.mask"
        os.system(
            f"mv {field_model_path}/{source}.maskgen.mask {source_mask} >/dev/null 2>&1"
        )
        os.system(f"rm -r {field_model_path}/*maskgen* >/dev/null 2>&1")

    os.system(f"rm -r {bgmodel}* >/dev/null 2>&1")

    for tt in [f"tt{i}" for i in range(nterms)]:
        mask = f"{field_model_path}/{source}.{imname}.mask.{tt}"
        modelfile = f"{model}.{tt}"
        bgmodelfile = f"{bgmodel}.{tt}"

        # Remove target source from field model
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
    tclean(
        vis=work_ms,
        field=field,
        cell=[cellsize],
        imsize=[imsize],
        startmodel=[f"{bgmodel}.tt{i}" for i in range(nterms)],
        savemodel="modelcolumn",
        niter=0,
        imagename=f"{field_model_path}/{source}.im_presub",
        nterms=nterms,
        deconvolver="mtmfs",
        scales=clean_scales,
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
    os.system(f"rm -r {field_model_path}/{source}.im_presub* >/dev/null 2>&1")

    # Perform UV-subtraction.
    uvsub(vis=work_ms)

    # Reimage to confirm field subtraction
    os.system(f"rm -r {field_model_path}/{source}.im_subbed* >/dev/null 2>&1")
    tclean(
        vis=work_ms,
        field=field,
        datacolumn="corrected",
        cell=[cellsize],
        imsize=[imsize],
        threshold=threshold,
        niter=iterations // 4,
        imagename=f"{field_model_path}/{source}.im_subbed",
        nterms=nterms,
        deconvolver="mtmfs",
        scales=clean_scales,
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

    # Export to FITS format.
    # ----------------------

    imtypes = [f"deep{i}" for i in range(clean_round)] + ["subbed"]
    taylorterms = [f"tt{i}" for i in range(nterms)]
    stokes_params = ["I", "V"]

    for tt in taylorterms:
        for stokes in stokes_params:
            for imtype in imtypes:
                os.system(
                    f"rm -r {field_model_path}/{source}.im_{imtype}.{stokes}.image* >/dev/null 2>&1"
                )
                imagename = f"{field_model_path}/{source}.im_{imtype}.image.{tt}/"
                subimagename = (
                    f"{field_model_path}/{source}.im_{imtype}.{stokes}.image.{tt}"
                )
                fitsname = f"{field_model_path}/{source}.im_{stokes}.{imtype}.{tt}.fits"

                imsubimage(
                    imagename=imagename,
                    outfile=subimagename,
                    stokes=stokes,
                )
                exportfits(
                    imagename=subimagename,
                    fitsimage=fitsname,
                    overwrite=True,
                )


if __name__ == "__main__":
    main()
