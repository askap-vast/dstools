## Installation / Setup ##

Simply clone this repo to a working directory which will contain your data. 

These tools rely on your RPFITS files residing in a `data/` subdirectory. Any extra files with the same target name (e.g. backup 1934 scans) can be placed in a subdirectory of `data/` to avoid being auto-ingested.

### NOTE: If Skipping ATCA Calibration ###

If beginning from ATCA RPFITS files, the other necessary directory structure will be created automatically. If you already have calibrated data, however, you will want to create the following sub-directories:
```
reduced/<PROJECT_NAME>/<TARGET_NAME>/
reduced/<PROJECT_NAME>/miriad/
```
where `<PROJECT_NAME>` can be whatever you like and `<TARGET_NAME>` is the name of your target source, as written in your observing schedule. Then make sure to place your calibrated visibilities (`.uv` format) in the `miriad` subdirectory.

## Dependencies ##

* CASA 6
* python-casacore
* click
* numpy
* astropy
* pandas

## ATCA Calibration ##

`atca_cal.sh` is a miriad script to calibrate and flag ATCA data prior to use with CASA in the other scripts. I have designed the flow of this to match my needs which may not necessarily meet yours! If not, just produce flagged and calibrated visibilities in `.uv` format as you otherwise would, and make sure to drop them in the `miriad` subdirectory as describe in Installation / Setup

You can optionally edit the first few lines of the script to change the `refant`, `mfinterval` and `bpinterval` parameters that will be used for each miriad task.

Run the script with
```
atca_cal.sh <PROJECT_NAME> <PROJECT_CODE>
```
where `<PROJECT_CODE>` is your observation's project code and `<PROJECT_NAME>` is the name of a sub-directory to store processing output, call this whatever you like. 

The script then steps through prompts to:
* load the RPFITS files, 
* select the IF to be processed, 
* declare which scan targets to use for primary and secondary calibrators and the target source,
* perform primary calibrator flagging and bandpass / flux calibration,
* perform secondary calibrator flagging and gain calibration,
* perform science target flagging,
* and finally transfer all calibration solutions to the science target.

## Imaging, Field Modeling \& Subtraction ##

`image.py` is a CASA script that imports calibrated visibilities and performs self-calibration and imaging tasks. The deconvolution stage also stores the model for all cleaned field sources and produces a uv-subtracted MeasurementSet with (ideally) only the target source remaining---to be further used in constructing lightcurves and dynamic spectra.

Run the script with
```
casa -c image.py reduced/<PROJECT_NAME>/<TARGET_NAME> <BAND> <CONFIG>
```
NOTE: the command interface for this script will be updated in future to avoid writing the full directory path
The `<BAND>` and `<CONFIG>` arguments specify the frequency band (`L`, `C`, or `X`) and array configuration (`6km`, `750_no6`, `750_6`, or `H168`), which will be used to set imaging parameters. If not provided, `<BAND>` will default to a value of `L` and `<CONFIG>` will default to a value of `6km`.

You can further adjust parameters in this script, such the:
* robust parameter,
* primary beam limit,
* number of Taylor terms used in mtmfs deconvolution,
* number and size of clean scales,
* and whether to use interactive `tclean`.

The script then progresses through prompts to load the data into a MeasurementSet and optionally test the imaging parameters are appropriate and perform additional flagging before performing further processing. 

The next stage is an optional self-calibration loop to improve phase calibration, with prompts to adjust solution intervals and perform additional rounds until happy with the results. Intermediate self-calibration products are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/selfcal/<BAND>/`.

Next is a deep clean to produce the field source model. By default this is set to interactive mode to place a clean mask carefully over each field source. Once happy that the mask covers all sources of interest you can set the threshold and iterations to desired levels before running the remainder of the task automatically. 

Once the field model has been produced, the next step is to mask the target source from this model so that it is untouched by uv subtraction. If your source is at the phase / image centre you can use automatic masking to cut out a circle of radius 10 pixels. The better approach is to opt out of this which will trigger another interactive `tclean` run so that you can draw a clean mask manually around your target source. Click the green arrow to finish `tclean` and this will be used to remove your source from the field model. The masked field model will then be inserted into the MeasurementSet and subtracted from the visibilities. The resulting subtracted MeasurementSet is stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/` with a `subbed.ms` suffix.

`tclean` will run once more to produce a field-subtracted image so that you can confirm no field sources remain. Results of both the deep and subtracted clean are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/field_model/<BAND>/`.

## Baseline Averaging / Dynamic Spectrum Forming ##

These are a few small CASA / Python scripts used to collapse a MeasurementSet (either the subtracted version from `image.py` or another MeasurementSet entirely) into a time vs frequency array ready for plotting.

`fix_phasecentre.py` is a CASA script that rotates a MeasurementSet's phase centre to provided coordinates. 

Run the script with
```
casa -c fix_phasecentre.py reduced/<PROJECT_NAME><MS> <COORDS>
```
where `<MS>` is the MeasurementSet you want to work with and `<COORDS>` are the target coordinates in the format `"J2000 hms dms"` format (e.g. `"J2000 12h30m41.2s -45d11m04.1s"`). `:` delimiters will result in the declination being treated in hourangle units.

`avg_baselines.py` is a CASA script that averages the data over all baselines.

Run the script with
```
casa -c avg_baselines.py reduced/<PROJECT_NAME>/<MS> <MIN_UVDIST>
```
where `<MIN_UVDIST>` is an optional lower cutoff to uv distance, used to improve results with RFI affected short baselines. 

`make_dspec.py` is a python script that converts the flux, frequency, and time values in a MeasurementSet data column into numpy arrays for each instrumental polarisation. These are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/dynamic_spectra/<BAND>/`. If a corrected data column exists it will be used, though the subtracted products from `image.py` will be fresh MeasurementSets with no corrected column.

Run the script with
```
python make_dspec.py -B <BAND> reduced/<PROJECT_NAME>/<MS>
```

## Plotting ##

`plot_dspec.py` is a convenience script for plotting the dynamic spectra produced by `make_dspec.py`.

Run the script with
```
python plot_dspec.py -f <FREQ_AVG> -t <TIME_AVG> -B <BAND> reduced/<PROJECT_NAME>
```
where `<FREQ_AVG>` and `<TIME_AVG>` are integer factors to average / rebin the data by. See other optional arguments / flags with
```
python plot_dspec.py --help
```
