import astropy.units as u
import matplotlib
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

pytestmark = pytest.mark.fullstack

imaging = pytest.importorskip("dstools.imaging")
ms_module = pytest.importorskip("dstools.ms")

get_pb_correction = imaging.get_pb_correction
make_pb_image = imaging.make_pb_image
Image = imaging.Image
MeasurementSet = ms_module.MeasurementSet
WSClean = imaging.WSClean
WSCleanModel = imaging.WSCleanModel


matplotlib.use("Agg")


def test_wsclean_fits_mask(mocker, ms_path, imaging_workspace):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    image_path = imaging_workspace["clean_mask"]
    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        fits_mask=image_path,
    )

    mock_ms = mocker.Mock(path=ms_path)

    cmd = wsclean.run(mock_ms, name="test")

    assert f"-fits-mask {image_path.absolute()}" in cmd


def test_wsclean_run_command_args(mocker, ms_path):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        threads=4,
        abs_mem=50,
        parallel_deconvolution=100,
        parallel_gridding=2,
        parallel_reordering=3,
    )

    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-name test" in cmd
    assert "-j 4" in cmd
    assert "-abs-mem 50" in cmd
    assert "-parallel-reordering 3" in cmd
    assert "-parallel-gridding 2" in cmd
    assert "-parallel-deconvolution 100" in cmd


@pytest.mark.parametrize(
    "phasecentre",
    [
        (281.2719, -63.9632),
        ("281.2719", "-63.9632"),
        ("18:45:05.256", "-63:57:47.520"),
        ("18h45m05.256s", "-63d57m47.520s"),
    ],
)
def test_wsclean_coordinate_parsing(phasecentre, mocker, ms_path):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        phasecentre=phasecentre,
    )

    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-shift 18h45m05.256s -63d57m47.52s" in cmd


def test_wsclean_run_command_multiscale(mocker, ms_path):
    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        multiscale=True,
    )

    mock_ms = mocker.Mock(path=ms_path)
    cmd = wsclean.run(mock_ms, name="test")

    assert "-multiscale" in cmd
    assert "-multiscale-scale-bias 0.7" in cmd
    assert "-multiscale-max-scales 8" in cmd


def test_get_pb_correction(mocked_ms):
    scale = get_pb_correction(
        mocked_ms,
        mocked_ms.phasecentre,
        pb_image=mocked_ms.path.with_suffix(".pb.fits"),
    )

    assert round(scale, 3) == 1


def test_get_pb_correction_no_existing_image(
    mocker,
    workspace_onespw_ms,
    imaging_workspace,
):
    ms = MeasurementSet(workspace_onespw_ms)
    pb_image = imaging_workspace["pb"]
    mocker.patch("dstools.imaging.make_pb_image", return_value=pb_image)

    scale = get_pb_correction(
        ms,
        ms.phasecentre,
        pb_image=pb_image.with_suffix(".noexist.fits"),
    )

    assert round(scale, 3) == 1


def test_make_pb_image(workspace_onespw_ms, imaging_workspace, casa_task_mocks):
    ms = MeasurementSet(workspace_onespw_ms)
    pb_image = imaging_workspace["pb"]

    pb_image_out = make_pb_image(
        pb_image,
        ms,
        ms.phasecentre,
    )

    assert pb_image == pb_image_out


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


@pytest.mark.parametrize("radius", [None, 5 * u.arcsec])
def test_wsclean_model_get_circular_mask(radius, im_paths):
    model_path = im_paths["model"]
    model = WSCleanModel(model_dir=model_path.parent)

    position = SkyCoord(model.phasecentre, unit="hourangle,deg")

    mask = model.get_circular_mask(
        position=position,
        radius=radius,
    )

    assert mask.shape == Image(model.model).data.shape
    assert mask.sum() > 0


def test_wsclean_model_get_interactive_mask(im_paths, mocker):
    mocker.patch("matplotlib.pyplot.show")

    model_path = im_paths["model"]
    model = WSCleanModel(model_dir=model_path.parent)

    mask = model.get_interactive_mask()

    assert mask.shape == Image(model.model).data.shape
    assert (~mask).sum() == 0


def test_wsclean_model_applymask_alltrue(imaging_workspace):
    model_path = imaging_workspace["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((500, 500), True)
    model.apply_mask(mask)

    model_image = Image(model_path / "test-MFS-I-model.fits", name="model")

    assert np.any(model_image.data > 0)


def test_wsclean_model_applymask_allfalse(imaging_workspace):
    model_path = imaging_workspace["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((500, 500), False)
    model.apply_mask(mask)

    model_image = Image(model_path / "test-MFS-I-model.fits", name="model")

    assert np.all(model_image.data == 0)


def test_wsclean_model_applymask_invalid_shape_raises_error(imaging_workspace):
    model_path = imaging_workspace["model"]
    model = WSCleanModel(model_dir=model_path)

    mask = np.full((499, 499), True)

    with pytest.raises(ValueError):
        model.apply_mask(mask)


def test_get_pb_correction_existing_image(workspace_onespw_ms, imaging_workspace):
    ms = MeasurementSet(workspace_onespw_ms)
    pb_image = imaging_workspace["pb"]

    scale = get_pb_correction(
        ms,
        ms.phasecentre,
        pb_image=pb_image,
    )

    assert round(scale, 3) == 1


def test_get_pb_correction_existing_image_offset_location(
    workspace_onespw_ms,
    imaging_workspace,
):
    ms = MeasurementSet(workspace_onespw_ms)
    pb_image = imaging_workspace["pb"]
    position = ms.phasecentre.directional_offset_by(0 * u.deg, 30 * u.arcsec)

    scale = get_pb_correction(
        ms,
        position,
        pb_image=pb_image,
    )

    assert round(scale, 3) == 0.998


def test_get_pb_correction_provided_outside_image(
    workspace_onespw_ms,
    imaging_workspace,
):
    ms = MeasurementSet(workspace_onespw_ms)
    pb_image = imaging_workspace["pb"]
    position = ms.phasecentre.directional_offset_by(0 * u.deg, 3 * u.arcmin)

    scale = get_pb_correction(
        ms,
        position,
        pb_image=pb_image,
    )

    assert scale == 1
