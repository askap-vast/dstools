import numpy as np
import pytest
from astropy.io import fits

from dstools.mask import (
    _adaptive_minimum_absolute_clip,
    _minimum_absolute_clip,
    beam_shape_erode,
    create_beam_mask_kernel,
    create_boxcar_skew_mask,
    minimum_absolute_clip,
)


@pytest.fixture
def test_beam_header():
    header = fits.Header()
    header["BMAJ"] = 0.01
    header["BMIN"] = 0.01
    header["BPA"] = 0.0
    header["CDELT1"] = -0.001
    header["CDELT2"] = 0.001

    return header


@pytest.fixture
def simple_image():
    """Simple test image with known values and separations."""

    size = 500
    image = np.zeros((size, size), dtype=float)

    # Place 300 mJy source at centre
    pos_source = size // 2
    image[pos_source, pos_source] = 300

    # Place -100 mJy artefact diagonally off-centre
    pos_artefact = size // 2 - size // 4
    image[pos_artefact, pos_artefact] = -100

    # Place 100 mJy artefact diagonally off-centre
    pos_artefact = size // 2 + size // 4
    image[pos_artefact, pos_artefact] = 100

    return image


@pytest.fixture
def skewed_image():
    """Test image with known values in central square."""

    rng = np.random.default_rng(42)
    size = 500
    scale = 50

    # Begin with random Gaussian noise at 50 mJy level
    image = rng.random((size, size), dtype=float) * 50

    # Place -100 mJy artefact off-centre
    pos_artefact = size // 2 - size // 8
    image[
        pos_artefact - scale : pos_artefact + scale,
        pos_artefact - scale : pos_artefact + scale,
    ] = -100

    # Place 100 mJy artefact off-centre
    pos_artefact = size // 2 + size // 8
    image[
        pos_artefact - scale : pos_artefact + scale,
        pos_artefact - scale : pos_artefact + scale,
    ] = 100

    return image


def test_boxcar_skew_mask_output_shape_and_types(simple_image):
    result = create_boxcar_skew_mask(simple_image, skew_delta=0.2, box_size=3)
    assert isinstance(result.positive_pixel_frac, np.ndarray)
    assert isinstance(result.skew_mask, np.ndarray)
    assert result.positive_pixel_frac.shape == simple_image.shape
    assert result.skew_mask.shape == simple_image.shape
    assert result.skew_delta == 0.2
    assert result.box_size == 3


@pytest.mark.parametrize(
    "box_size, total_flux",
    [
        # 300 mJy + 100 mJy separated by 250 pixels (outside 200 pixel box)
        (200, 400),
        # 300 mJy + 100 mJy separated by 250 pixels (within 500 pixel box)
        (500, 300),
    ],
)
def test_minimum_absolute_clip_box_size(box_size, total_flux, simple_image):
    mask = _minimum_absolute_clip(
        simple_image,
        increase_factor=2.0,
        box_size=box_size,
    )
    assert isinstance(mask, np.ndarray)
    assert mask.shape == simple_image.shape
    assert simple_image[mask].sum() == total_flux


@pytest.mark.parametrize(
    "increase_factor, total_flux",
    [
        (0.5, 400),
        (2.0, 300),
    ],
)
def test_minimum_absolute_clip_increase_factor(
    increase_factor,
    total_flux,
    simple_image,
):
    mask = _minimum_absolute_clip(
        simple_image,
        increase_factor=increase_factor,
        box_size=500,
    )
    assert isinstance(mask, np.ndarray)
    assert mask.shape == simple_image.shape
    assert simple_image[mask].sum() == total_flux


@pytest.mark.parametrize(
    "depth, expected_flux",
    [
        # With three box expansions we reach a scale enclosing both components
        # and successfully filter the positive artefact
        (3, 100),
        # With one box expansion we are still within each separate component
        # and fail to filter the positive artefact
        (1, 101),
    ],
)
def test_adaptive_clip_reduces_skew(depth, expected_flux, skewed_image):
    mask = _adaptive_minimum_absolute_clip(
        skewed_image,
        box_size=60,
        adaptive_max_depth=depth,
    )
    assert skewed_image[mask].max() < expected_flux


def test_adaptive_clip_box_too_large_breaks_cleanly(skewed_image):
    image = skewed_image[250 - 50 : 250 + 50, 250 - 50 : 250 + 50]
    mask = _adaptive_minimum_absolute_clip(
        image,
        increase_factor=2,
        box_size=101,
        adaptive_max_depth=3,
        adaptive_box_step=2.0,
        adaptive_skew_delta=0.2,
    )
    assert image[mask].max() == 100


def test_adaptive_clip_no_skew_breaks_cleanly(skewed_image):
    image = skewed_image[250 - 50 : 250 + 50, 250 - 50 : 250 + 50]
    mask = _adaptive_minimum_absolute_clip(
        image,
        increase_factor=2,
        box_size=50,
        adaptive_max_depth=3,
        adaptive_box_step=2.0,
        adaptive_skew_delta=0.2,
    )
    assert image[mask].max() == 100


def test_minimum_absolute_clip_no_adaptive_mode(skewed_image):
    mask = minimum_absolute_clip(
        skewed_image,
        adaptive_max_depth=None,
    )

    assert skewed_image[mask].max() == 100


def test_minimum_absolute_clip_adaptive_mode(skewed_image):
    mask = minimum_absolute_clip(
        skewed_image,
        adaptive_max_depth=3,
    )

    assert skewed_image[mask].max() < 100


def test_create_beam_mask_kernel_generates_boolean_array(test_beam_header):
    mask = create_beam_mask_kernel(test_beam_header, kernel_size=50)

    assert mask.dtype == bool
    assert mask.shape == (50, 50)
    assert np.any(mask)


def test_create_beam_mask_invalid_header_raises_error(test_beam_header):
    test_beam_header.pop("CDELT1")
    with pytest.raises(KeyError):
        create_beam_mask_kernel(test_beam_header, kernel_size=50)


def test_create_beam_mask_unveven_cdelt_raises_error(test_beam_header):
    test_beam_header.update({"CDELT1": -0.002})
    with pytest.raises(AssertionError):
        create_beam_mask_kernel(test_beam_header, kernel_size=50)


def test_beam_shape_erode(test_beam_header):
    # Create a beam shaped initial mask
    mask = create_beam_mask_kernel(test_beam_header, kernel_size=50) * 10

    # Erode using the same beam header,
    # we should get only one pixel in the eroded mask
    eroded = beam_shape_erode(mask, test_beam_header)

    assert mask.shape == eroded.shape
    assert np.sum(mask) == 600
    assert np.sum(eroded) == 1


def test_beam_shape_erode_invalid_header_returns_original_mask(test_beam_header):
    mask = create_beam_mask_kernel(test_beam_header, kernel_size=50) * 10

    test_beam_header.pop("BMAJ")
    eroded = beam_shape_erode(mask, test_beam_header)

    assert np.all(mask == eroded)
