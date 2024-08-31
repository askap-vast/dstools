# DStools #

`DStools` is a processing pipeline to produce and post-process dynamic spectrum data products from radio interferometer visibilities. `DStools` currently directly supports extraction of dynamic spectra from the following telescopes:
* ATCA
* ASKAP
* MeerKAT
* VLA

## Table of Contents
- [Installation](#installation)
- [Command Line Scripts](#cli)
    - [ATCA Calibration](#atca-cal)
    - [ASKAP Pre-processing](#askap-preprocess)
    - [Field Modeling](#model-field)
    - [Model Subtraction](#model-subtraction)
    - [Dynamic Spectrum Extraction](#ds-extraction)
    - [Plotting](#ds-plotting)
- [DStools Library](#dstools-library)

<a name="installation"></a>
## Installation

### Dependencies ###

`DStools` is built on top of `CASA 6.6` for the majority of tasks, and also uses `miriad` for pre-processing and calibration of ATCA observations. Make sure these tools are installed on your system:

* `CASA 6.6`
* `miriad`

### Installation / Configuration ###

Install `DStools` using `pip` or your preferred package manager:
```
pip install radio-dstools
```

To make `DStools` available within your CASA environment, run the following setup script:
```
dstools-setup
```
You will only need to run this once.

<a name="cli"></a>
## Command Line Scripts ##

Generating dynamic spectra with `DStools` generally involves some level of pre-processing (e.g. calibration, flagging, field source subtraction), extraction of the visibilities into a custom HDF5 data structure, and some level of post-processing (e.g. averaging in time/frequency, period folding, time/frequency/baseline filtering, polarisation processing, etc.). The following commands are provided to perform these steps:
| Command                    | Description                                                                    |
| -------------------------- | ------------------------------------------------------------------------------ |
| `dstools-cal`              | convenience script for flagging/calibration of raw ATCA visibilities           |
| `dstools-askap-preprocess` | script to set the correct flux scale and pointing centre of ASKAP beams        |
| `dstools-model-field`      | script to produce a model of field sources with optional self-calibration loop |
| `dstools-subtract-model`   | script to subtract a (multi-term) field model image from visibilities          |
| `dstools-extract-ds`       | script to extract visibilities for use with the `DStools` library              |
| `dstools-plot-ds`          | convenience script to post-process and plot dynamic spectra in various ways    |

The following scripts are used in the above commands, but are also available for more modular processing needs:

| Command                  | Description                                                    |
| -------------------------| -------------------------------------------------------------- |
| `_dstools-combine-spws`  | script to combine multiple spectral windows                    |
| `_dstools-avg-baselines` | script to average visibilities over all baselines              |
| `_dstools-rotate`        | script to rotate the phasecentre to a given set of coordinates |

Below some common workflows are described. For further details on usage, run any of these commands with the `--help` flag.

<a name="atca-cal"></a>
### ATCA Calibration ###

`dstools-cal` is a convenience script to calibrate raw ATCA data prior to use with CASA in the other scripts.

```
dstools-cal -d <DATA_DIR> <PROJECT_DIR> <PROJECT_CODE>
```
where `<DATA_DIR>` is a directory containing your raw ATCA RPFITS files, `<PROJECT_DIR>` is the name of a sub-directory to store processing output, and `<PROJECT_CODE>` is your observation's project code. You can supply further options to modify gain and bandpass solution intervals, set the reference antenna, and run automatic calibration/flagging. 

In manual mode the script steps through prompts to:
* load the RPFITS files, 
* select the IF to be processed, 
* declare which scan targets to use for primary and secondary calibrators and the target source,
* perform primary calibrator flagging and bandpass / flux calibration,
* perform secondary calibrator flagging and gain calibration,
* and finally transfer all calibration solutions to the science target and export to MeasurementSet format.

In auto mode (use flag `-A`) these steps will be performed non-interactively.

<a name="askap-preprocess"></a>
### ASKAP Pre-processing ###

ASKAP data requires extra pre-processing to 
1) set the instrumental polarisation flux scale to agree with CASA conventions (e.g. `I = (XX + YY)/2`), and 
2) set the reference frame of the beam phase centre to the correct coordinates.

These corrections should be applied before any further imaging or dynamic spectrum tasks. You can run both steps with:
```
dstools-askapsoft-preprocess <MS>
```
where `<MS>` is the path to your data in MeasurementSet format.

<a name="model-field"></a>
### Imaging and Field Modeling ###

`dstools-model-field` is a CASA script to image calibrated visibilities, optionally run a self-calibration loop, and produce a final model of selected sources in the field.

Run the script with
```
dstools-model-field <MS>
```

You can supply further options (see details with `dstools-model-field --help`) such as:
* array configuration and frequency band (to help choose imaging parameters),
* imaging phasecentre,
* robust parameter,
* clean threshold,
* maximum clean iterations,
* primary beam fractional power cutoff,
* number of Taylor terms used in mtmfs deconvolution,
* number and size of clean scales,
* use of automatic clean masking,
* use of the widefield gridder with w-projection,
* and whether to use interactive mode in `tclean`.

The script then progresses through prompts to:
* run autoflaggers,
* run an optional self-calibration loop, with intermediate self-calibration products stored in a `selfcal/<BAND>` directory,
* store the field-model imaging products in a `field_model/<BAND>` directory.

<a name="model-subtraction"></a>
### Model Subtraction ###

`dstools-subtract-model` is a CASA script to subtract a field model (either produced by `dstools-field-model` or elsewhere) from the visibilities. Model images can be supplied for multiple Taylor terms  You can optionally mask the model components of your target from the field model so that it is untouched by model subtraction.

Run the script with
```
dstools-subtract-model <MS> <MODEL>.tt0 <MODEL>.tt1 ...
```
where `<MODEL>.tt0` (etc.) are the model image path for each Taylor term used in the production of your model.

You can supply further options (see details with `dstools-subtract-model --help` such as:
* phasecentre coordinates at which the model images are centred,
* robust parameter,
* use of the widefield gridder with w-projection,
* whether to mask your target interactively,
* or alternatively the coordinates of your target for automatic masking.

The masked field model will then be subtracted from the visibilities with the result stored in the `CORRECTED_DATA` column of the MeasurementSet, and a model-subtracted image stored in the `field_model` directory.

<a name="ds-extraction"></a>
### Dynamic Spectrum Extraction ###

`dstoools-extract-ds` is the main task to extract dynamic spectra from your MeasurementSet. This command stores dynamic spectra as a 4D cube with dimensions of `baselines x integrations x channels x instrumental polarisations` as an HDF5 data structure, which can be read and post-processed with the `DynamicSpectrum` class provided with the `DStools` library. By default the script will average the visibilities over the baseline axis to save memory and disk space.

Run the script with:
```
dstools-extract-ds <MS> <DS>
```
where `<MS>` is the path to your data and `<DS>` is the path to store your output dynamic spectrum.

You can supply further options (see details with `dstools-extract-ds --help`) to:
* set the phasecentre at which to extract the dynamic spectrum with `-p <RA> <DEC>` (coordinates can be in sexagesimal or decimal degree formats),
* select extraction from either the `DATA`, `CORRECTED_DATA`, or `MODEL_DATA` column,
* throw away baselines shorter than some threshold in meters with (for example) `-u 500`
* disable averaging over the baseline axis with `-B`,
* correct for primary beam attenuation by supplying a primary beam map (e.g. from tclean) with `-P <PB_PATH>.pb.tt0`,
* disable masking of flagged data with `-F`.

<a name="ds-plotting"></a>
### Plotting ###

`dstools-plot-ds` is a convenience script to plot the dynamic spectra produced by `dstools-extract-ds`, as well as perform post-processing to produce 1D lightcurves and spectra, average the data in time and frequency, fold the data to a specified period,

To produce a basic dynamic spectrum, run the script with
```
dstools-plot-ds -d <DS>
```
where `<DS>` is your HDF5 dynamic spectrum file. 

Some other simple options (see details with `dstools-plot-ds --help` include:
* choose which Stokes parameters to plot with a subset of `{I, Q, U, V, L}` (e.g. `-s IQUV`)
* plot a channel-averaged lightcurve with `-l`,
* plot a time-averaged spectrum with `-p`,
* produce a summary plot including a lightcurve, spectrum, and dynamic spectra in all polarisations with `-Y`,
* average in time (`-t`) or frequency (`-f`) by an integer factor (e.g. `-t 5 -f 10` to average every five integrations and 10 channels),
* perform RM synthesis and correct for Faraday rotation with `-E`,
* plot the Faraday dispersion function with `-R`,
* perform 2D auto-correlation of the dynamic spectra with `-a` to highlight periodic features,
* fold the data to a specified period with `-FT <PERIOD>`.

<a name="dstools-library"></a>
## DStools Library ##

`dstools` can also be imported into your own scripts/notebooks as a package for more customised plotting. The main object is `DynamicSpectrum` which can be created as follows:
```
import matplotlib.pyplot as plt
from dstools.dynamic_spectrum import DynamicSpectrum

# Create DS object
ds = DynamicSpectrum(ds_path='path/to/dynamic_spectrum.hdf5')

# Plot Stokes I dynamic spectrum with real visibilities and color-scale clipped at 20 mJy
fig, ax = ds.plot_ds(stokes='I', cmax=20, imag=False)

# Add or modify custom plot elements here using the fig and ax objects
...

plt.show()
```

The `DynamicSpectrum` class takes the following keyword arguments:
| Parameter                 | Type             | Default | Description                                                   |
| ------------------------- | -----------------|-------- | ------------------------------------------------------------- |
| `tavg`                    | int              | 1       | factor by which to average the data across time               |
| `favg`                    | int              | 1       | factor by which to average the data across frequency channels |
| `mintime` / `maxtime`     | float            | None    | min and max cuts on frequency in units of `tunit`             |
| `minfreq` / `maxfreq`     | float            | None    | min and max cuts on frequency in units of MHz                 |
| `minuvdist` / `maxuvdist` | float            | None    | min and max cuts on baseline distance in units of meters      |
| `minuvwave` / `maxuvwave` | float            | None    | min and max cuts on baseline distance in units of wavelengths |
| `tunit`                   | astropy Quantity | u.hour  | time unit to use for selection and plotting                   |
| `corr_dumptime`           | astropy Quantity | 10*u.s  | correlator dumptime, used to detect calibrator scan breaks    |
| `derotate`                | bool             | False   | Apply Faraday de-rotation to linear polarisations             |
| `fold`                    | bool             | False   | enable folding, must also provide `period` keyword            |
| `period`                  | float            | None    | period on which to fold the data in units of `tunit`          |
| `period_offset`           | float            | 0.0     | period phase offset in units of `period`                      |
| `fold_periods`            | float            | 2       | number of folded periods to display for visualisation         |
| `calscans`                | bool             | True    | insert breaks during off-source time                          |
| `trim`                    | bool             | True    | remove flagged channel ranges at band edges                   |

Note: selection on baseline distance requires DS extraction without averaging over baselines (see `dstools-extract-ds`)
