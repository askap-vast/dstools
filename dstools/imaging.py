import glob
import logging
import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from casaconfig import config

config.logfile = "/dev/null"

import dask.array as da
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from astropy.time import Time
from astropy.wcs import WCS
from casatasks import (applycal, exportfits, flagdata, gaincal, imsubimage,
                       split, tclean)
from casatools import table
from dstools.utils import prompt
from numpy.typing import ArrayLike

logger = logging.getLogger(__name__)


@dataclass
class Model(ABC):

    model_dir: Path

    def __post_init__(self):
        self._load()
        self._validate()

    @abstractmethod
    def _load(self):
        """Load model images setting self.image property."""

    @abstractmethod
    def _validate(self):
        """Validate existence of model images."""

    @abstractmethod
    def insert_into(self):
        """Predict model visibilities into MODEL_DATA column of a Measurementset."""

    @abstractmethod
    def apply_mask(self, mask: ArrayLike):
        """Back up each model image and apply mask."""


@dataclass
class CASAModel(Model):

    def __post_init__(self):
        super().__post_init__()

        self.name = str(self.image).replace(".tt0.fits", "")

        self.nterms = len(self.model_images)

    def _load(self):
        image_path = list(self.model_dir.glob("*.model.tt0"))[0]

        self.model = str(image_path).replace(".model.tt0", ".model.I.tt0.fits")
        self.image = self.model.replace(".model", ".image")
        self.residual = self.model.replace(".model", ".residual")

        for imtype in [".model", ".residual", ".image"]:

            im = str(image_path).replace(".model.tt0", f"{imtype}.tt0")
            stokesI_im = im.replace(f"{imtype}.tt0", f"{imtype}.I.tt0")
            subim = im.replace(f"{imtype}.tt0", f"{imtype}.I.tt0.fits")

            imsubimage(
                imagename=im,
                outfile=stokesI_im,
                stokes="I",
                overwrite=True,
            )
            exportfits(
                imagename=stokesI_im,
                fitsimage=subim,
                overwrite=True,
            )
            os.system(f"rm -r {stokesI_im}")

        self.model_images = list(self.model_dir.glob("*.model.tt*"))

        return

    def _validate(self):

        # Check for existence of tt0 Stokes I FITS image
        if len(self.image) == 0:
            msg = f"Path {self.model_dir} does not contain any CASA model images with pattern '*.model.I.tt0.fits'."
            raise ValueError(msg)

        # Check for existence of TT model images
        tt_images = [im for im in self.model_images]
        if len(tt_images) == 0:
            msg = f"Path {self.model_dir} does not contain any CASA model images with pattern '*.model.tt*'."
            raise ValueError(msg)

        return

    def insert_into(self, ms):

        image_props = imhead(imagename=self.model)
        model_base = model[0].replace(".tt0", "")

        imsize = image_props["shape"][0]
        cell = round(abs(image_props["incr"][0]) * u.rad.to(u.arcsec), 3)
        freqaxis = np.argwhere(image_props["axisnames"] == "Frequency")
        freq = image_props["refval"][freqaxis][0, 0] * u.Hz.to(u.MHz)

        cellsize = f"{cell}arcsec"
        reffreq = f"{freq}MHz"

        tclean(
            vis=ms,
            cell=[cellsize],
            imsize=[imsize],
            startmodel=[f"{model_base}.tt{i}" for i in range(self.nterms)],
            savemodel="modelcolumn",
            niter=0,
            imagename=f"{maskgen}.im_presub",
            nterms=self.nterms,
            deconvolver="mtmfs",
            reffreq=reffreq,
            weighting="briggs",
            stokes="IQUV",
            robust=0.5,
            gridder="widefield",
            wprojplanes=-1,
            pblimit=-1,
        )
        os.system(f"rm -r {maskgen}.im_presub* >/dev/null 2>&1")

        return

    def apply_mask(self, mask: ArrayLike):

        for image in self.model_images:
            backup = str(image).replace(".model", ".model.premask")
            if not os.path.exists(backup):
                os.system(f"cp -r {image} {backup}")

            t = table()
            t.open(str(image), nomodify=False)

            data = t.getcol("map")
            data[~mask.T, :, :, :] = 0
            t.putcol("map", data)

            t.close()

        return


@dataclass
class WSCleanModel(Model):

    def __post_init__(self):
        super().__post_init__()

        self.name = str(self.image).replace("-MFS-I-image.fits", "")

        chan_images = [
            im for im in self.model_images if "MFS" not in str(im) and "-I-" in str(im)
        ]
        self.channels_out = len(chan_images)

    def _load(self):
        self.model = [p for p in self.model_dir.glob("*-MFS-I-model.fits")][0]
        self.image = str(self.model).replace("-model", "-image")
        self.residual = str(self.model).replace("-model", "-residual")

        self.model_images = [p for p in self.model_dir.glob("*-model.fits")]

    def _validate(self):

        # Check for existence of sub-channel model images
        if len(self.image) == 0:
            msg = f"Path {self.model_dir} does not contain any wsclean model images with pattern '*-MFS-I-model.fits'."
            raise ValueError(msg)

        # Check for existence of MFS model image
        chan_images = [im for im in self.model_images if "MFS" not in str(im)]
        if len(chan_images) == 0:
            msg = f"Path {self.model_dir} does not contain any wsclean model images with pattern '*-model.fits'."
            raise ValueError(msg)

        return

    def insert_into(self, ms: str):

        wsclean_cmd = [
            "wsclean",
            "-multiscale",
            "-pol iquv",
            f"-name {self.name}",
            f"-channels-out {self.channels_out}",
            "-predict",
            # "-quiet",
            f"{ms}",
        ]
        wsclean_cmd = " ".join(wsclean_cmd)
        os.system(wsclean_cmd)

        return

    def apply_mask(self, mask: ArrayLike):

        error_msg = "Cannot mask image of shape {} with mask of shape {}"

        for image in self.model_images:
            backup = str(image).replace(".fits", ".premask.fits")
            os.system(f"cp {image} {backup}")

            with fits.open(image, mode="update") as hdul:
                data = hdul[0].data
                if data[0, 0, :, :].shape != mask.shape:
                    raise ValueError(error_msg.format(data.shape, mask.shape))

                data[0, 0, ~mask] = 0
                hdul[0].data = data

        return


@dataclass
class WSClean:

    ms: str
    name: str
    imsize: int
    cellsize: str
    nterms: int
    minuvw: float = 0
    iterations: int = 30000
    channels_out: int = 8
    deconvolution_channels: int = 8
    robust: float = 0.5
    gain: float = 0.8
    threshold: float = 3
    mask_threshold: float = 5
    phasecentre: str = ""
    subimage_size: int = 1
    intervals: int = 1
    verbose: bool = False

    def run(self):

        logger.debug(
            f"Imaging {self.ms} with {self.imsize}x{self.imsize} pixels, {self.cellsize} cellsize, and {self.nterms} Taylor terms."
        )

        verbose = "-quiet" if not self.verbose else ""
        intervals = f"-intervals-out {self.intervals}" if self.intervals > 1 else ""

        if self.channels_out == 1:
            chan_str = ""
        else:
            chan_str = f"-join-channels -channels-out {self.channels_out} -deconvolution-channels {self.deconvolution_channels}"

        # Run WSclean
        wsclean_cmd = [
            "wsclean",
            f"--size {self.imsize} {self.imsize}",
            f"-scale {self.cellsize}",
            "-multiscale",
            f"-niter {self.iterations}",
            f"-mgain {self.gain}",
            f"-auto-threshold {self.threshold}",
            f"-auto-mask {self.mask_threshold}",
            "-pol iquv",
            chan_str,
            f"-minuvw-m {self.minuvw}",
            f"-weight briggs {self.robust}",
            f"-fit-spectral-pol {self.nterms}",
            f"-parallel-deconvolution {self.subimage_size} ",
            intervals,
            self.phasecentre,
            f"-name {self.name}",
            verbose,
            f"{self.ms}",
        ]
        wsclean_cmd = " ".join(wsclean_cmd)
        os.system(wsclean_cmd)

        return


def plot_gain_solutions(caltable: str, calmode: str) -> None:

    # Read time and gains from calibration table
    t = table(caltable)

    time = t.getcol("TIME")
    time_start = Time(time[0] / 3600 / 24, format="mjd", scale="utc").iso
    time -= time[0]
    time /= 86400

    gains = t.getcol("CPARAM")

    t.close()

    # Calculate phase / amplitude from complex gain solutions
    if calmode == "p":
        gains = xr.apply_ufunc(
            da.angle,
            gains,
            dask="allowed",
            kwargs=dict(deg=True),
        )
    elif calmode == "a":
        gains = da.absolute(gains)
    else:
        raise ValueError("Parameter 'calmode' must be 'p' or 'a'.")

    # Plot solutions against time
    fig, ax = plt.subplots()

    colors = ("k", "r")
    for polaxis in range(gains.shape[0]):
        g = gains[polaxis, 0, :]

        ax.scatter(
            time,
            g,
            color=colors[polaxis],
            s=1,
            alpha=0.2,
        )

    ax.set_xlabel(f"Hours from UTC {time_start}")
    if calmode == "p":
        ax.set_ylabel("Phase [deg]")
        ax.set_ylim(-50, 50)
    else:
        ax.set_ylabel("Amplitude")

    fig.tight_layout()

    savefile = caltable.replace(".cal", ".cal.png")
    fig.savefig(savefile, format="png")

    return


def run_selfcal(
    ms: Path,
    calmode: str,
    gaintype: str,
    interval: str,
    split_data: bool,
    interactive: bool,
) -> None:
    """Perform self-calibration on MS with field model in the MODEL_DATA column."""

    try:
        unit = "min" if "min" in interval else "s" if "s" in interval else ""
        int(interval.replace(unit, ""))
    except ValueError:
        raise ValueError(
            "Argument 'interval' must have format <int>[min/s] (e.g. 10s, 1min, 5min)."
        )

    # TODO: Write check for existence of MODEL_DATA column

    path = ms.parent.absolute()
    selfcal_round = len(list(path.glob("*.cal"))) + 1

    ms = str(ms)

    # Select refant
    flagstats = flagdata(vis=ms, mode="summary")
    df = pd.DataFrame(flagstats["antenna"]).T.reset_index(names="antenna")
    df["percentage"] = (100 * df.flagged / df.total).round(1)
    logger.info(f"Antenna flagging statistics:\n{df}")

    if interactive:
        refant = input("Select reference antenna: ")
        while refant not in flagstats["antenna"]:
            print(f"Reference antenna must be in: {flagstats['antenna']}")
            refant = input("Select reference antenna: ")
    else:
        refant = df.sort_values("percentage", ascending=True).iloc[0].antenna

    logger.info(f"Solving for gains using antenna {refant} as reference antenna...")

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
    for mode in calmode:
        plot_gain_solutions(
            caltable=cal_table,
            calmode=mode,
        )

    if interactive:
        plt.show(block=False)

    # Confirm solution is good before applying
    cal_good = prompt(
        msg="Is selfcal a good solution?",
        bypass=not interactive,
        default_response=True,
    )

    if not cal_good:
        os.system(f"rm -r {cal_table}")
        return

    applycal(
        vis=ms,
        gaintable=[cal_table],
        interp="linear",
    )

    if split_data:
        split(
            vis=ms,
            outputvis=ms.replace(".ms", ".selfcal.ms"),
            datacolumn="corrected",
        )

    return


def get_pb_correction(primary_beam, ra, dec):

    if primary_beam is None:
        return 1

    position = SkyCoord(ra=ra, dec=dec, unit=("hourangle", "deg"))

    if ".fits" not in primary_beam:
        pbfits = primary_beam + ".dstools.fits"
        exportfits(
            imagename=primary_beam,
            fitsimage=pbfits,
            overwrite=True,
        )
        primary_beam = pbfits

    with fits.open(primary_beam) as hdul:
        header, data = hdul[0].header, hdul[0].data
        data = data[0, 0, :, :]

    if "dstools" in primary_beam:
        os.system(f"rm {pbfits} 2>/dev/null")

    wcs = WCS(header, naxis=2)
    x, y = wcs.wcs_world2pix(position.ra, position.dec, 1)
    x, y = int(x // 1), int(y // 1)
    xmax, ymax = data.shape

    # Check position is within limits of supplied PB image
    im_outside_limit = [
        x < 0,
        x > xmax,
        y < 0,
        y > ymax,
    ]
    if any(im_outside_limit):
        logger.warning(
            f"Position {ra} {dec} outside of supplied PB image, disabling PB correction."
        )
        return 1

    scale = data[x, y]

    logger.debug(
        f"PB correction scale {scale:.4f} measured at pixel {x},{y} in image of size {xmax},{ymax}"
    )

    return scale
