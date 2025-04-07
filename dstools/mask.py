import logging
from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.io import fits
from radio_beam import Beam
from scipy.ndimage import binary_erosion, minimum_filter
from scipy.signal import fftconvolve

logger = logging.getLogger(__name__)


class SkewResult(NamedTuple):
    positive_pixel_frac: np.ndarray
    """The fraction of positive pixels in a boxcar function"""
    skew_mask: np.ndarray
    """Mask of pixel positions indicating which positions failed the skew test"""
    box_size: int
    """Size of the boxcar window applies"""
    skew_delta: float
    """The test threshold for skew"""


def create_boxcar_skew_mask(
    image: np.ndarray,
    skew_delta: float,
    box_size: int,
) -> np.ndarray:
    assert 0.0 < skew_delta < 0.5, f"{skew_delta=}, but should be 0.0 to 0.5"
    assert len(image.shape) == 2, (
        f"Expected two dimensions, got image shape of {image.shape}"
    )
    logger.debug(f"Computing boxcar skew with {box_size=} and {skew_delta=}")
    positive_pixels = (image > 0.0).astype(np.float32)

    # Counting positive pixel fraction here. The su
    window_shape = (box_size, box_size)
    positive_pixel_fraction = fftconvolve(
        in1=positive_pixels,
        in2=np.ones(window_shape, dtype=np.float32),
        mode="same",
    ) / np.prod(window_shape)
    positive_pixel_fraction = np.clip(
        positive_pixel_fraction,
        0.0,
        1.0,
    )  # trust nothing

    skew_mask = positive_pixel_fraction > (0.5 + skew_delta)
    logger.debug(f"{np.sum(skew_mask)} pixels above {skew_delta=} with {box_size=}")

    return SkewResult(
        positive_pixel_frac=positive_pixel_fraction,
        skew_mask=skew_mask,
        skew_delta=skew_delta,
        box_size=box_size,
    )


def _minimum_absolute_clip(
    image: np.ndarray,
    increase_factor: float = 2.0,
    box_size: int = 100,
) -> np.ndarray:
    """Given an input image or signal array, construct a simple image mask by applying a
    rolling boxcar minimum filter, and then selecting pixels above a cut of
    the absolute value value scaled by `increase_factor`. This is a pixel-wise operation.

    Args:
        image (np.ndarray): The input array to consider
        increase_factor (float, optional): How large to scale the absolute minimum by.
        box_size (int, optional): Size of the rolling boxcar minimum filter

    Returns:
        np.ndarray: The mask of pixels above the locally varying threshold
    """

    logger.debug(f"Minimum absolute clip, {increase_factor=} {box_size=}")
    rolling_box_min = minimum_filter(image, box_size)

    image_mask = image > (increase_factor * np.abs(rolling_box_min))

    return image_mask


def _adaptive_minimum_absolute_clip(
    image: np.ndarray,
    increase_factor: float = 2.0,
    box_size: int = 100,
    adaptive_max_depth: int = 3,
    adaptive_box_step: float = 2.0,
    adaptive_skew_delta: float = 0.2,
) -> np.ndarray:
    logger.debug(
        f"Using adaptive minimum absolute clip with {box_size=} {adaptive_skew_delta=}"
    )
    min_value = minimum_filter(image, size=box_size)

    for box_round in range(adaptive_max_depth, 0, -1):
        skew_results = create_boxcar_skew_mask(
            image=image,
            skew_delta=adaptive_skew_delta,
            box_size=box_size,
        )
        if np.all(~skew_results.skew_mask):
            logger.info("No skewed islands detected")
            break
        if any([box_size > dim for dim in image.shape]):
            logger.info(f"{box_size=} larger than a dimension in {image.shape=}")
            break

        logger.debug(f"({box_round}) Growing {box_size=} {adaptive_box_step=}")
        box_size = int(box_size * adaptive_box_step)
        minval = minimum_filter(image, box_size)
        logger.debug("Slicing minimum values into place")

        min_value[skew_results.skew_mask] = minval[skew_results.skew_mask]

    mask = image > (np.abs(min_value) * increase_factor)

    return mask


def minimum_absolute_clip(
    image: np.ndarray,
    increase_factor: float = 2.0,
    box_size: int = 100,
    adaptive_max_depth: int | None = None,
    adaptive_box_step: float = 2.0,
    adaptive_skew_delta: float = 0.2,
) -> np.ndarray:
    """Adaptive minimum absolute clip (author: Tim Galvin).

    Implements minimum absolute clip method. A minimum filter of a particular
    boxc size is applied to the input image. The absolute of the output is taken
    and increased by a guard factor, which forms the clipping level used to construct
    a clean mask:

    >>> image > (absolute(minimum_filter(image, box)) * factor)

    The idea is only valid for zero mean and normally distributed pixels, with
    positive definite flux, making it appropriate for Stokes I.

    Larger box sizes and guard factors will make the mask more conservative. Should
    the boxcar be too small relative to some feature it is aligned it is possible
    that an excess of positive pixels will produce an less than optimal clipping
    level. An adaptive box size mode, if activated, attempts to use a larger box
    around these regions.

    The basic idea being detecting regions where the boxcar is too small is around
    the idea that there should be a similar number of positive to negative pixels.
    Should there be too many positive pixels in a region it is likely there is an

    Args:
        image (np.ndarray): Image to create a mask for
        increase_factor (float, optional):
            The guard factor used to inflate the absolute of the minimum filter.
        box_size (int, optional):
            Size of the box car of the minimum filter.
        adaptive_max_depth (Optional[int], optional):
            The maximum number of rounds that the adaptive mode is allowed to perform
            when rescaling boxcar results in certain directions.
        adaptive_box_step (float, optional):
            A multiplicative factor to increase the boxcar size by each round.
        adaptive_skew_delta (float, optional):
            Minimum deviation from 0.5 that needs to be met to classify a region as skewed.

    Returns:
        np.ndarray: Final mask
    """

    if adaptive_max_depth is None:
        return _minimum_absolute_clip(
            image=image,
            box_size=box_size,
            increase_factor=increase_factor,
        )

    adaptive_max_depth = int(adaptive_max_depth)

    return _adaptive_minimum_absolute_clip(
        image=image,
        increase_factor=increase_factor,
        box_size=box_size,
        adaptive_max_depth=adaptive_max_depth,
        adaptive_box_step=adaptive_box_step,
        adaptive_skew_delta=adaptive_skew_delta,
    )


def create_beam_mask_kernel(
    fits_header: fits.Header,
    kernel_size=100,
    minimum_response: float = 0.6,
) -> np.ndarray:
    """Make a mask using the shape of a beam in a FITS Header object. The
    beam properties in the header are used to generate the two-dimensional
    Gaussian main lobe, from which a cut is made based on the minimum
    power.

    Args:
        fits_header (fits.Header):
            The FITS header to create beam from
        kernel_size (int, optional):
            Size of the output kernel in pixels. Will be a square.
        minimum_response (float, optional):
            Minimum response of the beam shape for the mask to be constructed from.

    Raises:
        KeyError: Raised if CDELT1 and CDELT2 missing

    Returns:
        np.ndarray: Boolean mask of the kernel shape
    """

    assert 0.0 < minimum_response < 1.0, (
        f"{minimum_response=}, should be between 0 to 1 (exclusive)"
    )

    POSITION_KEYS = ("CDELT1", "CDELT2")
    if not all([key in fits_header for key in POSITION_KEYS]):
        raise KeyError(f"{POSITION_KEYS=}  all need to be present")

    beam = Beam.from_fits_header(fits_header)
    assert isinstance(beam, Beam)

    cdelt1, cdelt2 = np.abs(fits_header["CDELT1"]), np.abs(fits_header["CDELT2"])  # type: ignore
    assert np.isclose(cdelt1, cdelt2), (
        f"Pixel scales {cdelt1=} {cdelt2=}, but must be equal"
    )

    k = beam.as_kernel(
        pixscale=cdelt1 * u.Unit("deg"),
        x_size=kernel_size,
        y_size=kernel_size,
    )

    return k.array > (np.max(k.array) * minimum_response)


def beam_shape_erode(
    mask: np.ndarray,
    fits_header: fits.Header,
    minimum_response: float = 0.6,
) -> np.ndarray:
    """Construct a kernel representing the shape of the restoring beam at
    a particular level, and use it as the basis of a binary erosion of the
    input mask.

    The ``fits_header`` is used to construct the beam shape that matches the
    same pixel size

    Args:
        mask (np.ndarray):
            The current mask that will be eroded based on the beam shape
        fits_header (fits.Header):
            The fits header of the mask used to generate the beam kernel shape
        minimum_response (float, optional):
            The minimum response of the main restoring beam to craft the shape from.

    Returns:
        np.ndarray: The eroded beam shape
    """

    if not all([key in fits_header for key in ["BMAJ", "BMIN", "BPA"]]):
        logger.warning(
            "Beam parameters missing. Not performing the beam shape erosion. "
        )
        return mask

    logger.debug(f"Eroding the mask using the beam shape with {minimum_response=}")
    beam_mask_kernel = create_beam_mask_kernel(
        fits_header=fits_header,
        minimum_response=minimum_response,
    )

    # This handles any unsqueezed dimensions
    beam_mask_kernel = beam_mask_kernel.reshape(
        mask.shape[:-2] + beam_mask_kernel.shape
    )

    erode_mask = binary_erosion(
        input=mask,
        iterations=1,
        structure=beam_mask_kernel,
    )

    return erode_mask.astype(mask.dtype)
