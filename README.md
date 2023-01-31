# DStools #

A processing pipeline to produce dynamic spectra from ATCA and ASKAP visibilities.

## Installation / Setup ##

Clone this repo and install with `pip` or `poetry`, preferably into a virtual environment. In order to make the modules from this package available within CASA, create or edit `~/.casa/config.py` to add the following lines
```
import sys

dstools_path = '/path/to/dstools/installation'

sys.path.append('')
sys.path.append(dstools_path)
```
where `dstools_path` should point to the install path. If you are unsure where this is, open a python REPL after installing `dstools` and run
```
import dstools

print(dstools.__path__[0])
```

The command line scripts should now be accessible from anywhere.

These tools rely on your RPFITS files residing in a `data/` subdirectory of your working directory. Any extra files with the same target name (e.g. backup 1934 scans) can be placed in a subdirectory of `data/` to avoid being auto-ingested.

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
* colorlog
* numpy
* astropy
* pandas

## ATCA Calibration ##

`dstools-cal` is a miriad script to calibrate and flag ATCA data prior to use with CASA in the other scripts. I have designed the flow of this to match my needs which may not necessarily meet yours! If not, just produce flagged and calibrated visibilities in `.uv` format as you otherwise would, and make sure to drop them in the `miriad` subdirectory as describe in Installation / Setup

```
dstools-cal <PROJECT_NAME> <PROJECT_CODE>
```
where `<PROJECT_CODE>` is your observation's project code and `<PROJECT_NAME>` is the name of a sub-directory to store processing output, call this whatever you like. You can supply further options to modify gain and bandpass solution intervals and the reference antenna, see details with `dstools-cal --help`.

The script steps through prompts to:
* load the RPFITS files, 
* select the IF to be processed, 
* declare which scan targets to use for primary and secondary calibrators and the target source,
* perform primary calibrator flagging and bandpass / flux calibration,
* perform secondary calibrator flagging and gain calibration,
* perform science target flagging,
* and finally transfer all calibration solutions to the science target.

## Imaging, Field Modeling \& Subtraction ##

`dstools-model-field` is a CASA script that imports calibrated visibilities and performs self-calibration and imaging tasks. The deconvolution stage also stores the model for all cleaned field sources and produces a uv-subtracted MeasurementSet with (ideally) only the target source remaining---to be further used in constructing lightcurves and dynamic spectra.

Run the script with
```
dstools-model-field reduced/<PROJECT_NAME>/<TARGET_NAME>/<DATA>
```
where `<DATA>` are calibrated visibilities either in miriad (ending in .cal), UV FITS (ending in .calfits), or CASA MeasurementSet (ending in .ms) format.

You can supply further options such as:
* array configuration (`6km`, `750_no6`, `750_6`, or `H168`),
* frequency band (`low`, `mid`, `L`, `C`, or `X`),
* robust parameter,
* clean threshold,
* maximum clean iterations,
* primary beam limit,
* number of Taylor terms used in mtmfs deconvolution,
* number and size of clean scales,
* use of widefield gridder with w-projection,
* and whether to use interactive `tclean`.

The script then progresses through prompts to load the data into a MeasurementSet and optionally test the imaging parameters are appropriate and perform additional flagging before performing further processing. 

The next stage is an optional self-calibration loop to improve phase calibration, with prompts to adjust solution intervals and perform additional rounds until happy with the results. Intermediate self-calibration products are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/selfcal/<BAND>/`.

Next is a deep clean to produce the field source model. By default this is set to interactive mode to place a clean mask carefully over each field source. Once happy that the mask covers all sources of interest you can set the threshold and iterations to desired levels before running the remainder of the task automatically. 

Once the field model has been produced, the next step is to mask the target source from this model so that it is untouched by uv subtraction. If your source is at the phase / image centre you can use automatic masking to cut out a circle of radius 30 pixels. The better approach is to opt out of this which will trigger another interactive `tclean` run so that you can draw a clean mask manually around your target source. Click the green arrow to finish `tclean` and this will be used to remove your source from the field model. The masked field model will then be inserted into the MeasurementSet and subtracted from the visibilities. The resulting subtracted MeasurementSet is stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/` with a `subbed.ms` suffix.

`tclean` will run once more to produce a field-subtracted image so that you can confirm no field sources remain. Results of both the deep and subtracted clean are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/field_model/<BAND>/`.

## Baseline Averaging / Dynamic Spectrum Forming ##

These are a few small CASA / Python scripts used to collapse a MeasurementSet (either the subtracted version from `dstools-model-field` or another MeasurementSet entirely) into a time vs frequency array ready for plotting.

`dstools-rotate` is a script that rotates a MeasurementSet's phase centre to provided coordinates. 

Run the script with
```
dstools-rotate reduced/<PROJECT_NAME>/<MS> <COORDS>
```
where `<MS>` is the MeasurementSet you want to work with and `<COORDS>` are the target coordinates in the format `"J2000 hms dms"` format (e.g. `"J2000 12h30m41.2s -45d11m04.1s"`). `:` delimiters will result in the declination being treated in hourangle units.

`dstools-avg-baselines` is a script that averages the data over all baselines. Options are available to select either the DATA or CORRECTED_DATA column with `-c`, and include a lower cutoff to uv distance used to improve results with RFI affected short baselines with `-u`. 

Run the script with
```
dstools-avg-baselines reduced/<PROJECT_NAME>/<MS>
```
and see option details with
```
dstools-avg-baselines --help
```

`dstools-make-dspec` is a python script that converts the flux, frequency, and time values in a MeasurementSet data column into numpy arrays for each instrumental polarisation. These are stored in `reduced/<PROJECT_NAME>/<TARGET_NAME>/dynamic_spectra/<BAND>/`. Options include selection of data column, frequency band, and removal of flagging mask.

Run the script supplying the path to the baseline averaged MeasurementSet. If using products from the field modeling stage, this is
```
dstools-make-dspec reduced/<PROJECT_NAME>/<TARGET_NAME>/<MS>
```
See option details with
```
dstools-make-dspec --help
```

## Plotting ##

`dstools-plot-dspec` is a convenience script for plotting the dynamic spectra produced by `dstools-make-dspec`, as well as functionality to produce 1D lightcurves and spectra, 2D auto-correlation functions, and fold the data to a specified period before plotting.

To produce a basic dynamic spectrum, run the script with
```
dstools-plot-dspec -d -f <FREQ_AVG> -t <TIME_AVG> reduced/<PROJECT_NAME>/<TARGET_NAME>
```
where `<FREQ_AVG>` and `<TIME_AVG>` are integer factors to average / rebin the data by. See other option details with
```
dstools-plot-dspec --help
```

`dstools` can also be imported into your own scripts/notebooks as a package for more customised plotting. The main object is `DynamicSpectrum` which can be created as follows:
```
import matplotlib.pyplot as plt
from dstools.dynamic_spectrum import DynamicSpectrum

# Create DS object
ds = DynamicSpectrum(
    project='reduced/<PROJECT_NAME>/<TARGET_NAME>',
    tavg=<TIME_AVG>,
    favg=<FREQ_AVG>,
)

# Plot dynamic spectrum with real visibilities and color-scale clipped at 20 mJy/beam
fig, ax = ds.plot_ds(stokes='I', cmax=20, imag=False)

# Add or modify custom plot elements here using the fig and ax objects
...

plt.show()
```
