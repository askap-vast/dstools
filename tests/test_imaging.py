import numpy as np
import pytest

from dstools.imaging import Image, Model, WSClean, WSCleanModel, get_pb_correction


def test_image_construction(im_paths):
    im_path = im_paths["image"]
    image = Image(path=im_path, name="image")

    assert "DATE-OBS" in image.header.keys()
    assert image.data.shape == (500, 500)


def test_model_image_normalisation(im_paths):
    im_path = im_paths["image"]
    image = Image(path=im_path, name="image")

    assert image.norm.vmin < 0

    model_path = im_paths["model"]
    model_image = Image(path=model_path, name="model")

    assert model_image.norm.vmin == 0


def test_wsclean_model(im_paths):
    model_path = im_paths["model"]
    model = WSCleanModel(model_dir=model_path.parent)

    assert model.channels_out == 1


def test_wsclean_model_applymask_alltrue(temp_environment):
    model_path = temp_environment["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((500, 500), True)
    model.apply_mask(mask)

    model_image = Image(model_path / "test-MFS-I-model.fits", name="model")

    assert np.any(model_image.data > 0)


def test_wsclean_model_applymask_allfalse(temp_environment):
    model_path = temp_environment["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((500, 500), False)
    model.apply_mask(mask)

    model_image = Image(model_path / "test-MFS-I-model.fits", name="model")

    assert np.all(model_image.data == 0)


def test_wsclean_model_applymask_invalid_shape_raises_error(temp_environment):
    model_path = temp_environment["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((499, 499), True)

    with pytest.raises(ValueError):
        model.apply_mask(mask)


# def test_wsclean_fits_mask_non_boolean(mocker, ms_path, im_paths):
#     image_path = im_paths["image"]
#     wsclean = WSClean(
#         imsize=500,
#         cellsize="0.66asec",
#         fits_mask=image_path,
#     )

#     # Abstract version of MS object, only needs access to path attribute
#     mock_ms = mocker.Mock(path=ms_path)

#     with pytest.raises(ValueError):
#         wsclean.run(mock_ms, name="test")


def test_wsclean_fits_mask(mocker, ms_path, im_paths):
    image_path = im_paths["mask"]
    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        fits_mask=image_path,
    )

    # Abstract version of MS object, only needs access to path attribute
    mock_ms = mocker.Mock(path=ms_path)

    cmd = wsclean.run(mock_ms, name="test")

    assert f"-fits-mask {image_path.absolute()}" in cmd


def test_wsclean_run_command(mocker, ms_path):
    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
    )

    # Abstract version of MS object, only needs access to path attribute
    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-name test" in cmd


def test_get_pb_correction():
    pass
