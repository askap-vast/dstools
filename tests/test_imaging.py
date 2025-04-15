import astropy.units as u
import numpy as np
import pytest

from dstools.imaging import (
    Image,
    WSClean,
    WSCleanModel,
    get_pb_correction,
    make_pb_image,
)
from dstools.ms import MeasurementSet


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


def test_wsclean_fits_mask(mocker, ms_path, im_paths):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

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
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
    )

    # Abstract version of MS object, only needs access to path attribute
    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-name test" in cmd


def test_wsclean_run_command_multi_threads_mem_limit(mocker, ms_path):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        threads=2,
        abs_mem=50,
    )

    # Abstract version of MS object, only needs access to path attribute
    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-j 2" in cmd
    assert "-parallel-reordering 2" in cmd
    assert "-parallel-gridding 2" in cmd
    assert "-abs-mem 50" in cmd


def test_wsclean_run_command_multiscale(mocker, ms_path):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        multiscale=True,
    )

    # Abstract version of MS object, only needs access to path attribute
    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-multiscale" in cmd
    assert "-multiscale-scale-bias 0.7" in cmd
    assert "-multiscale-max-scales 8" in cmd


def test_get_pb_correction(ms):
    scale = get_pb_correction(
        ms,
        ms.phasecentre,
        pb_image=ms.path.with_suffix(".pb.fits"),
    )

    assert round(scale, 3) == 1


def test_get_pb_correction_existing_image(temp_environment):
    ms = MeasurementSet(temp_environment["onespw"])
    pb_image = temp_environment["pb"]

    scale = get_pb_correction(
        ms,
        ms.phasecentre,
        pb_image=pb_image,
    )

    assert round(scale, 3) == 1


def test_get_pb_correction_existing_image_offset_location(temp_environment):
    ms = MeasurementSet(temp_environment["onespw"])
    pb_image = temp_environment["pb"]
    position = ms.phasecentre.directional_offset_by(0 * u.deg, 30 * u.arcsec)

    scale = get_pb_correction(
        ms,
        position,
        pb_image=pb_image,
    )

    assert round(scale, 3) == 0.998


def test_get_pb_correction_provided_outside_image(temp_environment):
    ms = MeasurementSet(temp_environment["onespw"])
    pb_image = temp_environment["pb"]
    position = ms.phasecentre.directional_offset_by(0 * u.deg, 3 * u.arcmin)

    scale = get_pb_correction(
        ms,
        position,
        pb_image=pb_image,
    )

    assert scale == 1


def test_get_pb_correction_no_existing_image(mocker, temp_environment):
    ms = MeasurementSet(temp_environment["onespw"])
    pb_image = temp_environment["pb"]
    mocker.patch("dstools.imaging.make_pb_image", return_value=pb_image)

    scale = get_pb_correction(
        ms,
        ms.phasecentre,
        pb_image=pb_image.with_suffix(".noexist.fits"),
    )

    assert round(scale, 3) == 1


def test_make_pb_image(temp_environment):
    ms = MeasurementSet(temp_environment["onespw"])
    pb_image = temp_environment["pb"]

    pb_image_out = make_pb_image(
        pb_image,
        ms,
        ms.phasecentre,
    )

    assert pb_image == pb_image_out
