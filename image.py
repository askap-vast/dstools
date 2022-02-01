import os
import sys
import glob


def colored(msg):

    return "\033[91m{}\033[0m".format(msg)


def prompt(msg):

    msg = msg[:-1] + " (y/n)?\n"

    resp = input(colored(msg))
    if resp not in ["y", "n"]:
        resp = input(colored(msg))

    return True if resp == "y" else False


def nearest_power(number, powerlim=7):
    original = number
    exps = range(powerlim, 14)
    pows = [2 ** i for i in exps]

    keep = []
    # Find highest power of 2 that is lower than number,
    # iteratively subtract this and find next highest power
    # until residual is below threshold (default 128)
    while number > 2 ** powerlim:
        gtpow = [num for num, power in zip(exps, pows) if number > power][-1]
        keep.append(gtpow)

        number -= 2 ** gtpow

    # Add threshold to get next power of 2 above original number
    keep.append(powerlim)

    return sum([2 ** i for i in keep])


def update_param(name, val, dtype):
    while True:
        newval = input("Enter new {} (currently {}): ".format(name, val))

        # Accept old value if nothing entered
        if not newval:
            break

        # Test for correct input type
        try:
            val = dtype(newval)
            break
        except ValueError:
            print("{} must be of type {}".format(name, type(val)))

    return val

def resolve_array_config(band, config):
    """Determine reffreq, primary beam and cell sizes from array paramters."""

    wavelengths = {
        'L': 0.0967,
        'C': 0.0461,
        'X': 0.0200,
    }
    frequencies = {
        'L': '2100',
        'C': '5500',
        'X': '9000',
    }
    primary_beams = {
        'L': 0.75,
        'C': 0.25,
        'X': 0.25,
    }
    max_baselines = {
        '6km': 6000,
        '750_no6': 750,
        '750_6': 5020,
        'H168': 185,
    }

    freq = frequencies[band]
    
    wavelength = wavelengths[band]
    baseline = max_baselines[config]

    # Convert to resolution in arcsec
    resolution = wavelength / baseline * 206264.8
    
    imradius = primary_beams[band]
    cell = round(resolution / 5, 2)
         
    return freq, imradius, cell


# CLI ARGUMENT PARSING (arg 0, 1, and 2 are consumed by casa, -c, and image.py)
# --------------------

args = sys.argv[3:]

if len(args) == 1:
    band = 'L'
    config = '6km'
elif len(args) == 2:
    band = args[1]
    config = '6km'
else:
    band = args[1]
    config = args[2]

proj_dir = args[0]
if proj_dir[-1] != "/":
    proj_dir += "/"

# Interpret 'obs' in a source name as one of multiple observations in their own directory
# Interpret 'epoch' in a source name as a directory containing multiple sources (e.g. C3431)
if "obs" in proj_dir:
    source = proj_dir.split("_")[0].lower()
elif "epoch" in proj_dir or "reduced" in proj_dir or "merged" in proj_dir:
    source = proj_dir.split("/")[-2].lower()
else:
    source = proj_dir.split("/")[0].lower()

# PARAMETER SETTINGS
# ------------------

freq, imradius, cell = resolve_array_config(band, config)

cellsize = "{}arcsec".format(cell)
imsize = nearest_power(imradius * 3600 / cell)
pblim = -0.1

masksize = 20
iterations = 2000
threshold = 8e-5

nterms = 3
clean_scales = [0, 3, 8]
robust = 0.5

field = "0"
reffreq = "{}MHz".format(freq)
int_freqs = [".2"] if reffreq == "2100MHz" else [""]
interactive = True

clean_mask = "{}/{}.{}.mask".format(proj_dir, source, band)
# -------------------------------------------------

accept_params = prompt(
    "Imaging with radius of {} deg, using {} pixels with size {}. Proceed?".format(
        imradius, imsize, cellsize
    )
)
if accept_params:

    for int_freq in int_freqs:
        mirfile = "{}.{}{}.cal".format(source, reffreq[:-3], int_freq)
        msname = "{}.{}.ms".format(source, band)
        calibrated_ms = "{}/{}_{}_cal.ms".format(proj_dir, source, band)

        # Check if version exists with field sources removed
        if os.path.exists(proj_dir + mirfile + ".zapped"):
            mirfile += ".zapped"

        # DATA IMPORT / FLAGGING
        if os.path.exists(proj_dir + msname):
            reimport = prompt("Redo data import?")

            if reimport:
                os.system("rm -r {}/{}".format(proj_dir, msname))
                os.system("rm -r {}/{}.flagversions".format(proj_dir, msname))

                importmiriad(
                    mirfile="{}/{}".format(proj_dir, mirfile),
                    vis="{}/{}".format(proj_dir, msname),
                )
        else:
            importmiriad(
                mirfile="{}/{}".format(proj_dir, mirfile),
                vis="{}/{}".format(proj_dir, msname),
            )

        test_params = prompt("Experiment with imsize / cellsize / nterms / flagging?")

        # Check if default imsize is acceptable (e.g. no strong out of field sources)
        if test_params:
            os.system("rm -r {}/test_params".format(proj_dir))
            os.system("mkdir {}/test_params".format(proj_dir))
            testms = "{}/test_params/testms.ms".format(proj_dir)
            os.system("cp -r {} {}".format(proj_dir + msname, testms))

            while True:

                tclean(
                    vis=testms,
                    field=field,
                    cell=[cellsize],
                    imsize=[imsize],
                    savemodel="modelcolumn",
                    threshold=threshold,
                    niter=2000,
                    imagename="{}/test_params/test_im".format(proj_dir),
                    nterms=nterms,
                    deconvolver="mtmfs",
                    scales=clean_scales,
                    reffreq=reffreq,
                    weighting="briggs",
                    robust=robust,
                    stokes="IQUV",
                    interactive=interactive,
                    pblimit=pblim,
                )

                accept_params = prompt("Finished experimenting?")
                if accept_params:
                    break
                else:
                    os.system("rm -r {}/test_params/test_im*".format(proj_dir))

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

                        rflag = prompt("Run rflag?")
                        if rflag:
                            flagdata(vis=proj_dir + msname, mode="rflag")

                        tfcrop = prompt("Run tfcrop?")
                        if tfcrop:
                            flagdata(vis=proj_dir + msname, mode="tfcrop")

                        flag45 = prompt("Flag baseline 4-5?")
                        if flag45:
                            flagdata(
                                vis=proj_dir + msname, mode="manual", antenna="3&4"
                            )

                    # Update parameters
                    param = update_param("imsize", imsize, int)
                    pblim = update_param("pblim", pblim, float)
                    cellsize = update_param("cellsize", cellsize, str)
                    nterms = update_param("nterms", nterms, int)

        os.system("rm -rf {}/test_params".format(proj_dir))

        selfcal = prompt("Run selfcal?")

        if selfcal:

            os.system("mkdir -p {}/selfcal/{}".format(proj_dir, band))

            # Run selfcal
            files = glob.glob("{}/selfcal/{}/{}.{}.selfcal.[0-9].ms".format(proj_dir, band, source, band))
            i = len(files)
            cont = "y"
            while cont:
                i += 1
                selfcal_template = "{}/selfcal/{}/{}{}.selfcal_im{}"
                selfcal_im = selfcal_template.format(proj_dir, band, source, int_freq, i)

                # Use freshly split ms each iteration after first loop
                selfcal_ms = (
                    proj_dir
                    + "selfcal/{}/".format(band)
                    + msname.replace(".ms", ".selfcal.{}.ms".format(i - 1))
                )
                selfcal_ms = (
                    selfcal_ms if os.path.exists(selfcal_ms) else proj_dir + msname
                )

                # Set mask from previous iteration as input to next tclean loop
                selfcal_mask = (
                    selfcal_template.format(proj_dir, band, source, int_freq, i - 1) + ".mask"
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
                    mask=init_mask,
                    pblimit=pblim,
                )

                # Trial self-cal solution intervals
                while True:
                    interval = input("Select solution interval (in min/s): ")
                    try:
                        unit = (
                            "min"
                            if "min" in interval
                            else "s"
                            if "s" in interval
                            else ""
                        )
                        num = int(interval.replace(unit, ""))
                    except ValueError:
                        print(
                            "Invalid solution interval entered, must be format <int>[min/s]."
                        )
                        continue

                    # Save self-cal plots
                    cal_file = "{}/selfcal/{}/{}{}.phase_selfcal_{}.{}".format(
                        proj_dir, band, source, int_freq, interval, i
                    )
                    cal_table = cal_file + ".cal"
                    cal_plot = cal_file + ".png"
                    gaincal(
                        vis=selfcal_ms,
                        caltable=cal_table,
                        solint=interval,
                        calmode="p",
                        gaintype="G",
                    )
                    plotms(
                        vis=cal_table,
                        xaxis="time",
                        yaxis="phase",
                        plotrange=[0, 0, -30, 30],
                        plotfile=cal_plot,
                        showgui=True,
                    )

                    # Confirm solution is good before applying, else trial another solution
                    cal_good = prompt("Is self-cal a good solution?")

                    if cal_good:
                        applycal(vis=selfcal_ms, gaintable=[cal_table], interp="linear")
                        split(
                            vis=selfcal_ms,
                            outputvis=proj_dir
                            + "selfcal/{}/".format(band)
                            + msname.replace(".ms", ".selfcal.{}.ms".format(i)),
                            datacolumn="corrected",
                        )
                        break
                    else:
                        os.system("rm -r {} {}".format(cal_table, cal_plot))

                # Break from loop once solution is sufficient
                cont = prompt("Proceed with more selfcal?")

            # Backup clean mask to main project directory
            backup_mask = (
                selfcal_template.format(proj_dir, band, source, int_freq, i) + ".mask"
            )
            if os.path.exists(clean_mask):
                replace_mask = prompt(
                    "Update the existing clean mask with current mask?"
                )

                if replace_mask:
                    os.system("cp -r {} {}".format(backup_mask, clean_mask))
            else:
                os.system("cp -r {} {}".format(backup_mask, clean_mask))

            # Copy calibrated MS to root folder
            os.system("cp -r {} {}".format(selfcal_ms, calibrated_ms))

        else:
            selfcal_path = "{}/selfcal/{}/{}.{}.selfcal.[0-9].ms".format(proj_dir, band, source, band)
            files = glob.glob(selfcal_path)
            i = len(files)
            if i == 0:
                # If selfcal was skipped, use the original MS
                selfcal_ms = proj_dir + msname
            else:
                selfcal_ms = sorted(files)[-1]

            os.system("cp -r {} {}".format(selfcal_ms, calibrated_ms))

            # Check for existing clean mask
            if not os.path.exists(clean_mask):
                clean_mask = ""

        # -----------------
        #    DEEP CLEAN
        # -----------------
        
        field_model_path = "{}/field_model/{}/".format(proj_dir, band)
        if os.path.exists(field_model_path) and prompt("Start from fresh field model?"):
            os.system("rm -r {}".format(field_model_path))
            os.system("mkdir -p {}".format(field_model_path))

        # Deep clean to produce field model
        tclean(
            vis=calibrated_ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            threshold=threshold,
            niter=iterations * 5,
            imagename="{}/{}{}.im_deep".format(field_model_path, source, int_freq),
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            mask=clean_mask,
            interactive=interactive,
            pblimit=pblim,
        )

        # Mask out the source
        # Either use a generic circular mask at image center or specify via a custom tclean mask
        automask = prompt("Automatically generate source mask?")
        if automask:
            source_mask = "circle[[{}pix, {}pix], {}pix]".format(
                imsize // 2, imsize // 2, masksize * 3
            )
        else:
            tclean(
                vis=calibrated_ms,
                field=field,
                cell=[cellsize],
                imsize=[imsize],
                threshold=threshold,
                niter=1,
                imagename="{}/{}{}.maskgen".format(field_model_path, source, int_freq),
                nterms=nterms,
                deconvolver="mtmfs",
                scales=clean_scales,
                reffreq=reffreq,
                weighting="briggs",
                stokes="IQUV",
                robust=robust,
                interactive=interactive,
                pblimit=pblim,
            )
            source_mask = "{}/{}{}.source.mask".format(
                field_model_path, source, int_freq
            )
            os.system(
                "mv {}/{}{}.maskgen.mask {}".format(
                    field_model_path, source, int_freq, source_mask
                )
            )
            os.system(
                "rm -r {}/{}{}.maskgen*".format(field_model_path, source, int_freq)
            )

        model = "{}/{}{}.im_deep.model".format(field_model_path, source, int_freq)
        bgmodel = "{}/{}{}.im_deep.bgmodel".format(field_model_path, source, int_freq)

        # os.system('rm -r {}*'.format(bgmodel))
        for tt in ["tt{}".format(i) for i in range(nterms)]:
            mask = "{}/{}.mask.{}".format(field_model_path, source, tt)
            modelfile = "{}.{}".format(model, tt)
            makemask(
                mode="copy",
                inpimage=modelfile,
                output=mask,
                inpmask=source_mask,
                overwrite=True,
            )
            immath(
                imagename=[modelfile, mask],
                outfile=bgmodel + "." + tt,
                expr="IM0*(1-IM1)",
            )

        # Insert masked background model into visibilities and subtract
        os.system("rm -r {}/{}{}.im_presub*".format(field_model_path, source, int_freq))
        tclean(
            vis=calibrated_ms,
            field=field,
            cell=[cellsize],
            imsize=[imsize],
            startmodel=[bgmodel + ".tt{}".format(i) for i in range(nterms)],
            savemodel="modelcolumn",
            niter=0,
            imagename="{}/{}{}.im_presub".format(field_model_path, source, int_freq),
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            pblimit=pblim,
        )

        # Perform UV-subtraction and move to new file
        subbed_ms = calibrated_ms.replace(".ms", ".subbed.ms")
        os.system(
            "cp -r {} {}".format(calibrated_ms, calibrated_ms.replace(".ms", ".bak.ms"))
        )
        uvsub(vis=calibrated_ms)
        os.system("mv {} {}".format(calibrated_ms, subbed_ms))
        os.system(
            "mv {} {}".format(calibrated_ms.replace(".ms", ".bak.ms"), calibrated_ms)
        )

        # Reimage to confirm field subtraction
        os.system("rm -r {}/{}{}.im_subbed*".format(field_model_path, source, int_freq))
        tclean(
            vis=subbed_ms,
            field=field,
            datacolumn="corrected",
            cell=[cellsize],
            imsize=[imsize],
            threshold=threshold,
            niter=iterations // 4,
            imagename="{}/{}{}.im_subbed".format(field_model_path, source, int_freq),
            nterms=nterms,
            deconvolver="mtmfs",
            scales=clean_scales,
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=robust,
            pblimit=pblim,
        )

        # Export to FITS format
        for tt in ["tt{}".format(i) for i in range(nterms)]:
            for stokes in ["I", "V"]:
                for imtype in ["deep", "subbed"]:

                    os.system(
                        "rm -r {}/{}{}.im_{}.{}.image*".format(
                            field_model_path, source, int_freq, imtype, stokes
                        )
                    )

                    imsubimage(
                        imagename="{}/{}{}.im_{}.image.{}/".format(
                            field_model_path, source, int_freq, imtype, tt
                        ),
                        outfile="{}/{}{}.im_{}.{}.image.{}".format(
                            field_model_path, source, int_freq, imtype, stokes, tt
                        ),
                        stokes=stokes,
                    )
                    exportfits(
                        imagename="{}/{}{}.im_{}.{}.image.{}".format(
                            field_model_path, source, int_freq, imtype, stokes, tt
                        ),
                        fitsimage="{}/{}{}.im_{}.{}.{}.fits".format(
                            field_model_path, source, int_freq, stokes, imtype, tt
                        ),
                        overwrite=True,
                    )
