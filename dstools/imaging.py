import glob
import logging
import os
import subprocess
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from casaconfig import config

config.logfile = "/dev/null"

from astropy.wcs import WCS
from casatasks import exportfits, imsubimage, tclean
from casatools import table
from dstools.logger import parse_stdout_stderr
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

    imsize: int
    cellsize: str

    # deconvolution
    iterations: int = 30000
    major_cycles: int = 20
    mgain: float = 0.85
    channels_out: int = 8
    deconvolution_channels: int = 8
    spectral_pol_terms: int = 3

    # masking / thresholds
    fits_mask: Path | None = None
    mask_threshold: float = 5
    auto_threshold: float = 3
    local_rms_window: int = 25

    # weight / gridding
    robust: float = 0.5
    parallel_deconvolution: int | None = None
    phasecentre: tuple[str] | None = None

    # data selection
    pol: str = "iquv"
    data_column: str | None = None
    minuvw_m: float | None = None
    minuvw_l: float | None = None
    intervals_out: int = 1

    # multiscale
    multiscale: bool = False
    multiscale_scale_bias: float = 0.7
    multiscale_max_scales: int = 8

    # I/O
    reuse_psf: Path | None = None
    reuse_dirty: Path | None = None
    no_dirty: bool = False
    save_source_list: bool = False
    save_reordered: bool = False
    reuse_reordered: bool = False
    temp_dir: Path | None = None
    out_dir: Path = Path(".")

    verbose: bool = False

    def __post_init__(self):
        self.optional_args = (
            "fits_mask",
            "parallel_deconvolution",
            "data_column",
            "minuvw_m",
            "minuvw_l",
            "intervals_out",
            "reuse_psf",
            "reuse_dirty",
            "no_dirty",
            "save_source_list",
            "save_reordered",
            "reuse_reordered",
            "temp_dir",
        )

        if self.temp_dir is None:
            self.temp_dir = self.out_dir
        self.temp_dir = Path(self.temp_dir).absolute()

    @property
    def _multiscale_args(self):

        if self.multiscale is None:
            return ""

        return (
            "-multiscale "
            f"-multiscale-scale-bias {self.multiscale_scale_bias} "
            f"-multiscale-max-scales {self.multiscale_max_scales}"
        )

    @property
    def _phasecentre_args(self):

        if self.phasecentre is None:
            return ""

        ra, dec = self.phasecentre
        return f"-shift {ra} {dec}"

    @property
    def _spectral_args(self):

        if self.channels_out == 1:
            return ""

        return (
            f"-join-channels "
            f"-channels-out {self.channels_out} "
            f"-deconvolution-channels {self.deconvolution_channels} "
            f"-fit-spectral-pol {self.spectral_pol_terms}"
        )

    @property
    def _verbosity(self):
        return "" if self.verbose else "-quiet"

    def _format_optional_argument(self, arg: str):

        val = getattr(self, arg)
        if val is None:
            return ""

        # Convert to WSclean argument format
        arg = arg.replace("_", "-")

        # Strip boolean flags of value
        if isinstance(val, bool):
            return f"-{arg}" if val else ""

        return f"-{arg} {val}"

    def run(self, ms: Path, name: str):

        logger.debug(
            f"Imaging {ms} with {self.imsize}x{self.imsize} pixels, {self.cellsize} cellsize, and {self.spectral_pol_terms} spectral terms."
        )

        # Add all essential arguments
        wsclean_cmd = [
            "wsclean",
            f"-name {name}",
            f"-size {self.imsize} {self.imsize}",
            f"-scale {self.cellsize}",
            f"-niter {self.iterations}",
            f"-nmiter {self.major_cycles}",
            f"-mgain {self.mgain}",
            f"-pol {self.pol}",
            f"-weight briggs {self.robust}",
            f"-auto-threshold {self.auto_threshold}",
            f"-auto-mask {self.mask_threshold}",
            f"-local-rms-window {self.local_rms_window}",
            self._phasecentre_args,
            self._multiscale_args,
            self._spectral_args,
            self._verbosity,
        ]

        # Add provided optional keyword arguments
        for arg in self.optional_args:
            argstr = self._format_optional_argument(arg)
            wsclean_cmd.append(argstr)

        # Add MS positional argument
        ms = Path(ms).absolute()
        wsclean_cmd.append(str(ms))
        wsclean_cmd = " ".join(wsclean_cmd)

        # Create output directory
        model_path = ms.parent.absolute() / self.out_dir
        model_path.mkdir(exist_ok=True)

        # Move into working directory to store imaging products
        cwd = Path(".").absolute()
        os.chdir(model_path)

        # Run WSclean
        p = subprocess.Popen(
            wsclean_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )

        parse_stdout_stderr(p, logger, print_stdout=False)

        # Return to start directory
        os.chdir(cwd)

        return


def get_pb_correction(primary_beam: Path, ra: str, dec: str) -> float:

    position = SkyCoord(ra=ra, dec=dec, unit=("hourangle", "deg"))

    # If PB image is CASA format, create a temporary FITS copy
    if primary_beam.suffix != ".fits":
        pbfits = primary_beam.with_suffix(".dstools.fits")
        exportfits(
            imagename=primary_beam,
            fitsimage=pbfits,
            overwrite=True,
        )
        primary_beam = pbfits

    with fits.open(primary_beam) as hdul:
        header, data = hdul[0].header, hdul[0].data
        data = data[0, 0, :, :]

    # Remove temporary FITS file
    if "dstools" in primary_beam.name:
        os.system(f"rm {primary_beam} 2>/dev/null")

    # Find pixel coordinates of position in FITS image
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
