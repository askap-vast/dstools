import numpy as np
import pytest
from matplotlib.backend_bases import KeyEvent

from dstools.imaging import Image
from dstools.viewer import Viewer


@pytest.fixture
def images(im_paths):
    im = Image(im_paths["image"], name="image")
    mod = Image(im_paths["model"], name="model")
    res = Image(im_paths["residual"], name="residual")

    # Replace data with known constant values
    im.data = im.data * 0 + 10
    mod.data = im.data * 0 + 5
    res.data = im.data * 0 + 1

    return [im, mod, res]


@pytest.mark.parametrize(
    "switch_image, val",
    [
        ("image", 10),
        ("model", 5),
        ("residual", 1),
    ],
)
def test_viewer_button_switch_image(switch_image, val, images, mocker):
    mocker.patch("matplotlib.pyplot.show")

    viewer = Viewer(images=images)

    # Simulate button click to switch image
    image = next(img for img in images if img.name == switch_image)
    viewer._switch_image(image)(event=mocker.Mock())

    assert np.all(viewer.image_object.get_array().data == val)


def test_keypress_draw_mask(images, mocker):
    mocker.patch("matplotlib.pyplot.show")

    viewer = Viewer(images=images)

    # Simulate some click positions
    viewer.clicker = mocker.Mock()
    viewer.clicker.get_positions.return_value = {"mask": [(10, 10), (20, 10), (15, 20)]}

    coords = np.array([[x, y] for x in range(100) for y in range(100)])
    viewer.coords = coords.reshape((100 * 100, 2))
    viewer.mask = np.ones((100, 100), dtype=bool)

    # Simulate pressing the 'x' key
    event = KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="x")
    viewer._on_press(event)

    # Assert some pixels were masked
    assert np.any(~viewer.mask)

    # Simulate pressing 'c' to clear the same area
    viewer.clicker.get_positions.return_value = {"mask": [(10, 10), (20, 10), (15, 20)]}

    event = KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="c")
    viewer._on_press(event)

    # Assert those pixels are now unmasked
    assert np.all(viewer.mask)

    # Simulate pressing 'v' as an unsupported key does nothing
    event = KeyEvent(name="key_press_event", canvas=viewer.fig.canvas, key="v")
    viewer._on_press(event)

    # Assert pixels are still unmasked
    assert np.all(viewer.mask)


def test_update_colorscale(images, mocker):
    mocker.patch("matplotlib.pyplot.show")

    viewer = Viewer(images=images)

    # Grab updater and call it
    updater = viewer._update_colorscale(99)
    updater(None)

    # Ensure norm values changed
    assert viewer.image_object.norm.vmin is not None
    assert viewer.image_object.norm.vmax is not None
