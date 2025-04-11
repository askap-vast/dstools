# DStools # 
[![CI](https://github.com/askap-vast/dstools/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/askap-vast/dstools/actions)
[![codecov](https://codecov.io/gh/askap-vast/dstools/branch/main/graph/badge.svg)](https://codecov.io/gh/askap-vast/dstools)

[![pypi](https://img.shields.io/pypi/v/radio-dstools)](https://pypi.org/project/radio-dstools/)
[![versions](https://img.shields.io/pypi/pyversions/radio-dstools)](https://pypi.org/project/radio-dstools/)

[![DOI](https://zenodo.org/badge/306534672.svg)](https://zenodo.org/doi/10.5281/zenodo.13626182)


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
  - [Field Modeling / Insertion](#model-field)
  - [Model Subtraction](#model-subtraction)
  - [Dynamic Spectrum Extraction](#ds-extraction)
  - [Plotting](#ds-plotting)
- [DStools Library](#dstools-library)

<a name="installation"></a>
## Installation

### Dependencies ###

`DStools` uses `miriad` for pre-processing and calibration of ATCA observations, and `WSclean` for imaging and model insertion. Make sure these tools are installed on your system if you wish to use these parts of the package:

* `miriad`
* `WSclean 3.5`

`DStools` is built on top of modular `CASA` which is presently only available on linux.

### Installation / Configuration ###

Install `DStools` using `pip` or your preferred package manager:
```
pip install radio-dstools
```

<a name="cli"></a>
## Command Line Scripts ##

| Command                                                       | Description                                                                 |
|---------------------------------------------------------------|-----------------------------------------------------------------------------|
| `dstools-atca-cal`<span style="color:red">&dagger;</span>     | convenience script for flagging/calibration of raw ATCA visibilities        |
| `dstools-askap-preprocess`                                    | script to set the correct flux scale and reference frame of ASKAP beams     |
| `dstools-create-model`<span style="color:red">&Dagger;</span> | script to image and model the field with WSclean                            |
| `dstools-insert-model`<span style="color:red">&Dagger;</span> | script to predict model images into visibilities                            |
| `dstools-selfcal`                                             | script to run a self-calibration loop using inserted model visibilities     |
| `dstools-subtract-model`                                      | script to subtract an inserted model from the visibilities                  |
| `dstools-extract-ds`                                          | script to extract visibilities for use with the `DStools` library           |
| `dstools-plot-ds`                                             | convenience script to post-process and plot dynamic spectra in various ways |

*NOTE:* \
<span style="color:red">&dagger;</span> `miriad` is a requirement to run `dstools-atca-cal` \
<span style="color:red">&Dagger;</span> `WSclean` is a requirement to run `dstools-create-model` and `dstools-insert-model`

Below some common workflows are described. For further details on usage, run any of these commands with the `--help` flag.

<a name="atca-cal"></a>
### ATCA Calibration ###

`dstools-atca-cal` is a convenience script to calibrate ATCA (CABB) data prior to use in the other scripts.

```
dstools-atca-cal -d <DATA_DIR> <OUT_DIR> <PROJECT_CODE>
```
where `<DATA_DIR>` is a directory containing your raw ATCA RPFITS files, `<OUT_DIR>` is the name of a sub-directory to store processing output, and `<PROJECT_CODE>` is your observation's project code. The script will by default look for PKS 1934-638 as a bandpass and flux calibrator, and attempt to identify your intended target and the corresponding gain calibrator. If you have multiple targets, an alternate primary calibrator, or want to include a polarisation leakage calibrator, you can provide their names directly with the `--target`, `--primary-cal`, `--gain-cal`, and `--leakage-cal` options. 

You can provide further options (see details with `dstools-atca-cal --help` to:
* select the band and intermediate frequency to process,
* select the reference antenna,
* shift offsets between pointing and phase centre,
* tune calibration solution intervals,
* solve for absolute XY phase and leakage (requires a good leakage calibrator),
* produce diagnostic plots of calibration solutions,
* use the pipeline interactively, which may be necessary with heavily RFI affected observations.

The output from this script will be a calibrated target dataset in MeasurementSet format, ready for use with other scripts in this package.

<a name="askap-preprocess"></a>
### ASKAP Pre-processing ###

ASKAP data requires extra pre-processing to 
1) set the instrumental polarisation flux scale to agree with CASA conventions (e.g. `I = (XX + YY)/2`), 
2) set the instrumental polarisation axis to align with IAU convention (`X` and `Y` receivers aligned N-S and E-W respectively), and 
2) set the reference frame of the beam phase centre to the correct coordinates (by default the phasecentre is oriented to the mosaicked field centre coordinates).

These corrections should be applied before any further imaging or dynamic spectrum tasks, and are implemented using the [FixMS](https://github.com/AlecThomson/fix-ms) package. You can run both steps with:
```
dstools-askapsoft-preprocess <MS>
```
where `<MS>` is the path to your data in MeasurementSet format.

<a name="model-field"></a>
### Imaging and Self Calibration ###

`dstools-create-field` is a script to image calibrated visibilities with WSclean and produce a model of sources in the field. See the [WSclean documentation](https://wsclean.readthedocs.io/en/latest/) for further information on imaging modes and parameters.

Run the script with
```
dstools-create-model <MS>
```

You can supply further options (see details with `dstools-create-model --help`) including:
* array configuration and frequency band (to help choose imaging parameters),
* imaging phasecentre,
* mask threshold,
* clean threshold,
* robust parameter,
* maximum clean major cycles,
* maximum clean iterations,
* number of frequency sub-bands to seperately deconvolve,
* and number of polynomial terms used to model spectral structure in MFS deconvolution.

Image and model products will be stored in a directory determined by the `--out-dir` parameter.

`dstools-selfcal` is a script to run self-calibration using model visibilities stored in the `MODEL_DATA` column of your MeasurementSet.

Run the script with
```
dstools-selfcal <MS>
```

You can provide further options (see details with `dstools-selfcal --help`) to:
* solve for phase-only `--calmode p` or phase and amplitude `--calmode ap` calibration solutions,
* adjust the solution interval with `-i`,
* solve for gains independently over separate spectral windows with `-n`
* solve for either polarisation-dependent (`--no-combine-pols`) or polarisation-indepdendent (`--combine-pols`) gains.

### Field Modeling / Subtraction ###

`dstools-insert-model` is a script to predict model images into visibilities stored in the `MODEL_DATA` column of your MeasurementSet.

```
dstools-insert-model <MODEL_DIR> <MS>
```

An interactive viewer will be launched allowing you to draw masks around any portions of the image that you want removed from the model prior to insertion. This is generally where you would remove any model components that are associated with your target of interest. Click to create polygon vertices, then either press `x` to mask the enclosed region or `c` to unmask enclosed pixels that have previously been masked. When you are satisfied with the mask close the viewer and the model images will be masked and predicted into visibilities.

<a name="model-subtraction"></a>

`dstools-subtract-model` is a script that will subtract the model visibilities from the data, storing the result in the `CORRECTED_DATA` column of your MeasurementSet (`CORRECTED_DATA = DATA - MODEL_DATA`).

Run the script with
```
dstools-subtract-model <MS>
```

You can optionally use the `-S / --split-data` option, which will split off a new MeasurementSet (name with a `.subtracted.ms` suffix) with the `DATA` column containing the model subtracted visibilities. This can be useful for iterative model subtraction (e.g. removing bright sources far outside the primary beam prior to modelling the primary beam field).

<a name="ds-extraction"></a>
### Dynamic Spectrum Extraction ###

`dstoools-extract-ds` is the main task to extract dynamic spectra from your MeasurementSet. This command extracts dynamic spectra as a 4D cube with dimensions of `baselines x integrations x channels x instrumental polarisations`. Along with observation metadata this is stored as an HDF5 data structure which can be read and post-processed with the `DynamicSpectrum` class provided with the `DStools` library. By default the script will average the visibilities over the baseline axis to save memory and disk space.

Run the script with:
```
dstools-extract-ds <MS> <DS>
```
where `<MS>` is the path to your data and `<DS>` is the path to store your output dynamic spectrum.

You can supply further options (see details with `dstools-extract-ds --help`) to:
* set the phasecentre at which to extract the dynamic spectrum with `-p <RA> <DEC>` (coordinates can be in sexagesimal or decimal degree formats),
* select extraction from either the `DATA`, `CORRECTED_DATA`, or `MODEL_DATA` column with `-d`,
* throw away baselines shorter than some threshold in meters with (for example) `-u 500`
* disable averaging over the baseline axis with `-B`,
* correct for primary beam attenuation by supplying a primary beam map with `-P <PB_PATH>.fits`,
* disable masking of flagged data with `-F`.

<a name="ds-plotting"></a>
### Plotting ###

`dstools-plot-ds` is a convenience script to plot the dynamic spectra produced by `dstools-extract-ds`, as well as perform post-processing to produce 1D lightcurves and spectra, average the data in time and frequency, fold the data to a specified period, and perform polarisation analysis.

To produce a basic dynamic spectrum, run the script with
```
dstools-plot-ds -d <DS>
```
where `<DS>` is your HDF5 dynamic spectrum file. 

Some other simple options (see details with `dstools-plot-ds --help`) include:
* choose which Stokes parameters to plot with a subset of `{I, Q, U, V, L}` (e.g. `-s IQUV`)
* plot a channel-averaged lightcurve with `-l`,
* plot a time-averaged spectrum with `-p`,
* include polarisation fraction / angle / ellipticity in lightcurve / spectrum plots with `-P`,
* produce a summary plot including a lightcurve, spectrum, and dynamic spectra in all polarisations with `-Y`,
* average in time (`-t`) or frequency (`-f`) by an integer factor (e.g. `-t 5 -f 10` to average every five integrations and 10 channels),
* perform RM synthesis (using the [RM-lite](https://github.com/AlecThomson/rm-lite) package) and correct for Faraday rotation with `--derotate`,
* provide a known RM with `--rm`,
* perform 2D auto-correlation of the dynamic spectra with `-a` to highlight periodic features,
* fold the data to a specified period with `--fold -T <PERIOD>`.

<a name="dstools-library"></a>
## DStools Library ##

`dstools` can also be imported into your own scripts/notebooks as a package for more customised plotting. The main object is `DynamicSpectrum` which can be used as follows, supposing you extracted the dynamic spectrum to a file located at `dynamic_spectrum.ds`:
```
import matplotlib.pyplot as plt
from dstools.dynamic_spectrum import DynamicSpectrum
from dstools.plotting import plot_ds

# Create DS object
ds = DynamicSpectrum(ds_path='dynamic_spectrum.ds')

# Plot Stokes I dynamic spectrum with real visibilities and color-scale clipped at 20 mJy
fig, ax = plot_ds(ds, stokes='I', cmax=20, imag=False)

# Add or modify custom plot elements here using the fig and ax objects
...

plt.show()
```

The `DynamicSpectrum` class takes the following keyword arguments:
| Parameter                 | Type             | Default | Description                                                   |
|---------------------------|------------------|---------|---------------------------------------------------------------|
| `tavg`                    | int              | 1       | factor by which to average the data across time               |
| `favg`                    | int              | 1       | factor by which to average the data across frequency channels |
| `mintime` / `maxtime`     | float            | None    | min and max cuts on time in units of `tunit`                  |
| `minfreq` / `maxfreq`     | float            | None    | min and max cuts on frequency in units of MHz                 |
| `minuvdist` / `maxuvdist` | float            | None    | min and max cuts on baseline distance in units of meters      |
| `minuvwave` / `maxuvwave` | float            | None    | min and max cuts on baseline distance in units of wavelengths |
| `tunit`                   | astropy Quantity | u.hour  | time unit to use for selection and plotting                   |
| `corr_dumptime`           | astropy Quantity | 10*u.s  | correlator dumptime, used to detect calibrator scan breaks    |
| `derotate`                | bool             | False   | Apply Faraday de-rotation to linear polarisations             |
| `RM`                      | float            | None    | User provided rotation measure in units of rad / m^2          |
| `fold`                    | bool             | False   | enable folding, must also provide `period` keyword            |
| `period`                  | float            | None    | period on which to fold the data in units of `tunit`          |
| `period_offset`           | float            | 0.0     | period phase offset in units of `period`                      |
| `fold_periods`            | float            | 2       | number of folded periods to display for visualisation         |
| `calscans`                | bool             | True    | insert breaks during off-source time                          |
| `trim`                    | bool             | True    | remove fully flagged channel ranges at band edges             |

Note: selection on baseline distance requires DS extraction without averaging over baselines (see `dstools-extract-ds`)
