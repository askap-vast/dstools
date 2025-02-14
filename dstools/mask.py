import logging
from typing import NamedTuple

import numpy as np
from astropy.io import fits
from scipy.ndimage import minimum_filter

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
        in1=positive_pixels, in2=np.ones(window_shape, dtype=np.float32), mode="same"
    ) / np.prod(window_shape)
    positive_pixel_fraction = np.clip(
        positive_pixel_fraction, 0.0, 1.0
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
        increase_factor (float, optional): How large to scale the absolute minimum by. Defaults to 2.0.
        box_size (int, optional): Size of the rolling boxcar minimum filtr. Defaults to 100.

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
            logger.debug("No skewed islands detected")
            break
        if any([box_size > dim for dim in image.shape]):
            logger.debug(f"{box_size=} larger than a dimension in {image.shape=}")
            break

        logger.debug(f"({box_round}) Growing {box_size=} {adaptive_box_step=}")
        box_size = int(box_size * adaptive_box_step)
        _min_value = minimum_filter(image, box_size)
        logger.debug("Slicing minimum values into place")

        min_value[skew_results.skew_mask] = _min_value[skew_results.skew_mask]

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
        increase_factor (float, optional): The guard factor used to inflate the absolute of the minimum filter. Defaults to 2.0.
        box_size (int, optional): Size of the box car of the minimum filter. Defaults to 100.
        adaptive_max_depth (Optional[int], optional): The maximum number of rounds that the adaptive mode is allowed to perform when rescaling boxcar results in certain directions. Defaults to None.
        adaptive_box_step (float, optional): A multiplicative factor to increase the boxcar size by each round. Defaults to 2.0.
        adaptive_skew_delta (float, optional): Minimum deviation from 0.5 that needs to be met to classify a region as skewed. Defaults to 0.2.

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
