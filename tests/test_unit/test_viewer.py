import matplotlib
import numpy as np
import pytest
from matplotlib.backend_bases import KeyEvent

pytestmark = pytest.mark.fullstack

imaging = pytest.importorskip("dstools.imaging")
viewer_module = pytest.importorskip("dstools.viewer")

Image = imaging.Image
Viewer = viewer_module.Viewer


matplotlib.use("Agg")


@pytest.fixture
def images(im_paths):
    im = Image(im_paths["image"], name="image")
    mod = Image(im_paths["model"], name="model")
    res = Image(im_paths["residual"], name="residual")

    im.data = im.data * 0 + 10
    mod.data = im.data * 0 + 5
    res.data = im.data * 0 + 1

    return [im, mod, res]


@pytest.mark.parametrize(
    ("switch_image", "val"), [("image", 10), ("model", 5), ("residual", 1)]
)
def test_viewer_button_switch_image(switch_image, val, images, mocker):
    mocker.patch("matplotlib.pyplot.show")

    viewer = Viewer(images=images)
    image = next(img for img in images if img.name == switch_image)
    viewer._switch_image(image)(event=mocker.Mock())

    assert np.all(viewer.image_object.get_array().data == val)


def test_keypress_draw_mask(images, mocker):
    mocker.patch("matplotlib.pyplot.show")

    viewer = Viewer(images=images)
    viewer.clicker = mocker.Mock()
    viewer.clicker.get_positions.return_value = {"mask": [(10, 10), (20, 10), (15, 20)]}
    viewer.coords = np.array([[x, y] for x in range(100) for y in range(100)]).reshape(
        (100 * 100, 2)
    )
    viewer.mask = np.ones((100, 100), dtype=bool)

    viewer._on_press(
        KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="x")
    )
    assert np.any(~viewer.mask)

    viewer.clicker.get_positions.return_value = {"mask": [(10, 10), (20, 10), (15, 20)]}
    viewer._on_press(
        KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="c")
    )
    assert np.all(viewer.mask)

    viewer._on_press(
        KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="v")
    )
    assert np.all(viewer.mask)


def test_update_colorscale(images, mocker):
    mocker.patch("matplotlib.pyplot.show")
    images[0].data = np.arange(16, dtype=float).reshape(4, 4)

    viewer = Viewer(images=images)
    expected_vmin = np.nanpercentile(images[0].data, 1)
    expected_vmax = np.nanpercentile(images[0].data, 99)
    draw_idle = mocker.patch.object(viewer.fig.canvas, "draw_idle")

    updater = viewer._update_colorscale(99)
    updater(None)

    assert viewer.image_object.norm.vmin == expected_vmin
    assert viewer.image_object.norm.vmax == expected_vmax
    draw_idle.assert_called_once_with()
