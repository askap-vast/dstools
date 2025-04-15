import multiprocessing

import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from dstools.utils import (
    Array,
    get_available_cpus,
    parse_coordinates,
    prompt,
    rebin,
    rebin2D,
    slice_array,
)


def test_get_available_cpus_non_slurm():
    cpus = get_available_cpus()

    assert cpus == multiprocessing.cpu_count()


def test_get_available_cpus_slurm_cpus_per_task(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")

    cpus = get_available_cpus()

    assert cpus == 8


def test_get_available_cpus_slurm_cpus_on_node(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "64")

    cpus = get_available_cpus()

    assert cpus == 64


def test_get_available_cpus_slurm_cpus_on_node_and_cpus_per_task(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "64")

    cpus = get_available_cpus()

    assert cpus == 8


@pytest.mark.parametrize(
    "coord",
    [
        ("18:45:05.250", "-63:57:47.450"),
        ("18h45m05.250s", "-63d57m47.450s"),
    ],
)
def test_parse_coordinates(coord):
    parsed = parse_coordinates(coord)

    assert parsed == ("18h45m05.250s", "-63d57m47.450s")


def test_prompt_no_bypass_yes(mocker):
    mocker.patch("builtins.input", side_effect=["y"])

    assert prompt("continue")


def test_prompt_no_bypass_no(mocker):
    mocker.patch("builtins.input", side_effect=["n"])

    assert not prompt("continue")


def test_prompt_no_bypass_wrong_input(mocker):
    mocker.patch("builtins.input", side_effect=["d", "n"])

    assert not prompt("continue")


def test_prompt_bypass_yes():
    assert prompt("continue", bypass=True, default_response=True)


def test_prompt_bypass_no():
    assert not prompt("continue", bypass=True, default_response=False)


def test_prompt_bypass_message_warns():
    assert not prompt(
        "continue", bypass=True, bypass_msg="hello", default_response=False
    )


def test_rebin_row():
    bin_array = rebin(o=5, n=3, axis=0)
    expected = np.array(
        [
            [0.6, 0.4, 0.0, 0.0, 0.0],
            [0.0, 0.2, 0.6, 0.2, 0.0],
            [0.0, 0.0, 0.0, 0.4, 0.6],
        ]
    )

    assert np.allclose(bin_array, expected)


def test_rebin_column():
    bin_array = rebin(o=5, n=3, axis=1)
    expected = np.array(
        [
            [0.6, 0.4, 0.0, 0.0, 0.0],
            [0.0, 0.2, 0.6, 0.2, 0.0],
            [0.0, 0.0, 0.0, 0.4, 0.6],
        ]
    ).T

    assert np.allclose(bin_array, expected)


def test_rebin_empty():
    bin_array = rebin(o=0, n=0, axis=0)
    expected = np.zeros((0, 0))

    assert np.allclose(bin_array, expected)


@pytest.mark.parametrize(
    "orig_shape, new_shape",
    [
        ((100, 100), (5, 8)),
        ((50, 100), (5, 8)),
        ((80, 91), (8, 14)),
    ],
)
def test_rebin2D_conserves_total_flux(orig_shape, new_shape):
    array = np.random.random(orig_shape) + 1j * np.random.random(orig_shape)
    binned = rebin2D(array, new_shape)

    tolerance = 1e-12 + 1e-12j
    assert np.mean(array) - np.mean(binned) < tolerance


def test_rebin2D_masked_array():
    array = np.outer([1, 2, 3], [1, 2, 3]) * 1 + 0j
    masked_array = np.ma.array(array, mask=array > 5)
    binned = rebin2D(masked_array, (3, 3))

    assert np.isnan(binned).sum() == 3


@pytest.fixture
def array():
    return np.random.random((5, 3))


def test_rebin2D_same_shape(array):
    binned = rebin2D(array, (5, 3))

    assert np.allclose(array, binned)


def test_rebin2D_invalid_shape(array):
    with pytest.raises(ValueError):
        rebin2D(array, (3, 5))


def test_slice_array_1d(array):
    sliced = slice_array(array, ax1_min=1, ax1_max=4)

    assert sliced.shape == (3, 3)


def test_slice_array_2d(array):
    sliced = slice_array(array, ax1_min=1, ax1_max=4, ax2_min=1, ax2_max=2)

    assert sliced.shape == (3, 1)


def test_array_atca():
    array = Array(band="AT_L", config="6km")
    assert array.cell == 0.66
    assert array.imradius == 0.4125


def test_array_askap():
    array = Array(band="AK_low")
    assert array.cell == 2.5
    assert array.config == "AK"
    assert round(array.imradius, 2) == 2.13
