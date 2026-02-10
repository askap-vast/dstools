import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import SkyCoord

from dstools.ms import (
    CalTable,
    MeasurementSet,
    combine_spws,
    extract_baseline,
    extract_baselines,
    get_polslice,
    rotate_circular_feeds,
    rotate_linear_feeds,
    run_selfcal,
    swap_xy_feeds,
)
from dstools.utils import DataError

ms_properties = [
    ("nspws", 1),
    ("nbaselines", 3),
    ("integrations", 4),
    ("nchannels", 21),
    ("npols", 1),
    ("dimensions", (3, 4, 21, 1)),
    ("telescope", "ATCA"),
    ("feedtype", "linear"),
]


@pytest.mark.parametrize("prop, val", ms_properties)
def test_ms_basic_properties(prop, val, ms):
    assert getattr(ms, prop) == val


def test_ms_opens_with_str_path(ms_path):
    path = str(ms_path)
    ms = MeasurementSet(path=path)

    assert isinstance(path, str)
    assert isinstance(ms.path, Path)


def test_ms_not_exists_raises_error():
    with pytest.raises(FileNotFoundError):
        MeasurementSet("faketable.ms")


def test_ms_phasecentre(ms):
    assert ms.phasecentre.to_string("hmsdms") == "18h45m13.62000576s -63d57m34.3899846s"


columns = [
    "DATA",
    "CORRECTED_DATA",
    "MODEL_DATA",
    "SIGMA_SPECTRUM",
    "WEIGHT_SPECTRUM",
]


@pytest.mark.parametrize("column", columns)
def test_ms_column_exists(ms, column):
    assert ms.column_exists(column)


def test_combine_multi_spw_fails_if_not_split(ms):
    with pytest.raises(DataError):
        ms._combine_multi_spw()


def test_get_reference_antenna_auto(ms):
    refant = ms.get_reference_antenna(interactive=False)

    assert refant == "1"


def test_get_reference_antenna_interactive(ms, mocker):
    mocker.patch("builtins.input", side_effect="6")
    refant = ms.get_reference_antenna(interactive=True)

    assert refant == "6"


def test_get_reference_antenna_invalid_selection_prompts_user(
    ms,
    capfd,
    mocker,
):
    mocker.patch("builtins.input", side_effect=["41", "6"])
    refant = ms.get_reference_antenna(interactive=True)
    out, _ = capfd.readouterr()

    assert "Reference antenna must be in:" in out
    assert refant == "6"


def test_nspw_conversion_1to1(ms):
    assert ms.nspws == 1

    ms.to_nspws(1)

    assert ms.nspws == 1


def test_nspw_conversion_1to2(ms):
    assert ms.nspws == 1

    ms.to_nspws(2)

    assert ms.nspws == 2


def test_nspw_conversion_2to1(temp_environment):
    ms_path = temp_environment["twospw"]
    ms = MeasurementSet(path=ms_path)
    ms.original_path = temp_environment["onespw"]

    assert ms.nspws == 2

    ms.to_nspws(1)

    assert ms.nspws == 1


def test_nspw_conversion_sequence(ms):
    assert ms.nspws == 1

    ms.to_nspws(2)

    assert ms.nspws == 2

    ms.to_nspws(1)

    assert ms.nspws == 1


def test_invalid_nspw_conversion_raises_error(ms):
    with pytest.raises(ValueError):
        ms.to_nspws(-1)


def test_ms_average_baselines(ms):
    baseavg = ms.average_baselines()

    assert "baseavg" in str(baseavg.path)


def test_ms_rotate_phasecentre_inplace(ms):
    with pytest.raises(NotImplementedError):
        ms.rotate_phasecentre(ms.phasecentre, inplace=True)


def test_ms_rotate_phasecentre_new_location(ms):
    assert ms.phasecentre.to_string("hmsdms") == "18h45m13.62000576s -63d57m34.3899846s"

    phasecentre = SkyCoord(ra="18h45m00s", dec="-63d57m00s", unit="hourangle,deg")
    baseavg = ms.rotate_phasecentre(phasecentre)

    assert baseavg.phasecentre.to_string("hmsdms") == phasecentre.to_string("hmsdms")


def test_ms_rotate_phasecentre_same_location_does_nothing(ms):
    assert ms.phasecentre.to_string("hmsdms") == "18h45m13.62000576s -63d57m34.3899846s"

    baseavg = ms.rotate_phasecentre(ms.phasecentre)

    assert ms == baseavg


def test_subtract_model_no_model_column_raises_error(temp_environment):
    ms_path = temp_environment["rotated"]
    ms = MeasurementSet(path=ms_path)

    with pytest.raises(DataError):
        ms.subtract_model()


def test_subtract_model_no_split(ms):
    subtracted = ms.subtract_model()

    corrected = subtracted.getcolumn("CORRECTED_DATA")
    model = subtracted.getcolumn("MODEL_DATA")
    data = subtracted.getcolumn("DATA")

    assert subtracted.path == ms.path
    assert np.allclose(corrected, data - model)


def test_subtract_model_split(ms):
    no_split = ms.subtract_model()
    split = ms.subtract_model(split_ms=True)

    no_split_corrected = no_split.getcolumn("CORRECTED_DATA")
    split_corrected = split.getcolumn("DATA")

    assert np.allclose(no_split_corrected, split_corrected)


def test_ms_increment_selfcal_round_initial(ms):
    assert ms.increment_selfcal_round() == ms.path.with_suffix(".selfcal1.ms")


def test_ms_increment_selfcal_round_multiple(ms_path):
    ms = MeasurementSet(path=ms_path)
    selfcal_ms_path = ms.increment_selfcal_round()

    ms.path = selfcal_ms_path

    assert ms.increment_selfcal_round() == ms_path.with_suffix(".selfcal2.ms")


def test_ms_calc_flag_statistics(ms):
    df = ms.calc_flag_statistics()

    assert isinstance(df, pd.DataFrame)
    assert all(c in df.columns for c in ["antenna", "flagged", "total", "percentage"])


def test_ms_solve_gains(ms):
    ms.solve_gains(interval="10s", calmode="ap", gaintype="T")

    assert ms.caltable.path == ms.path.with_suffix(".cal")


ct_properties = [
    ("nspws", 16),
    ("npols", 1),
]


@pytest.mark.parametrize("prop, val", ct_properties)
def test_caltable_dimensions(prop, val, temp_environment):
    caltable_path = temp_environment["cal"]
    caltable = CalTable(path=caltable_path)

    assert getattr(caltable, prop) == val


@pytest.mark.parametrize("calmode", ["a", "p"])
def test_caltable_plot_solutions(calmode, temp_environment):
    caltable_path = temp_environment["cal"]
    caltable = CalTable(path=caltable_path)

    caltable.plot_solutions(calmode=calmode)


def test_ms_split_selfcal_round(temp_environment):
    ms_path = temp_environment["onespw"]
    ms = MeasurementSet(path=ms_path)
    ms2 = ms.split_selfcal_round()

    assert ms2.path == ms_path.with_suffix(".selfcal1.ms")


def test_run_selfcal_bad_interval_raises_error(mocker, temp_environment):
    ms_path = temp_environment["onespw"]
    ms = MeasurementSet(path=ms_path)

    with pytest.raises(ValueError):
        run_selfcal(
            ms,
            calmode="ap",
            gaintype="T",
            interval="10",
            split_data=False,
            interactive=False,
        )


def test_run_selfcal_no_model_raises_error(mocker, temp_environment):
    ms_path = temp_environment["rotated"]
    ms = MeasurementSet(path=ms_path)

    with pytest.raises(DataError):
        run_selfcal(
            ms,
            calmode="ap",
            gaintype="T",
            interval="10s",
            split_data=False,
            interactive=False,
        )


def test_run_selfcal_integration_cal_good(mocker, temp_environment):
    ms_path = temp_environment["onespw"]
    ms = MeasurementSet(path=ms_path)

    selfcal_ms = run_selfcal(
        ms,
        calmode="ap",
        gaintype="T",
        interval="10s",
        split_data=True,
        interactive=False,
        refant="6",
    )

    assert selfcal_ms.path == ms.path.with_suffix(".selfcal1.ms")
    assert selfcal_ms.nspws == 1


def test_run_selfcal_integration_cal_bad(mocker, temp_environment):
    ms_path = temp_environment["onespw"]

    mocker.patch.object(MeasurementSet, "to_nspws", return_value=None)
    mocker.patch.object(CalTable, "plot_solutions", return_value=None)
    mocker.patch("dstools.ms.prompt", return_value=False)
    mocker.patch("matplotlib.pyplot.show")

    ms = MeasurementSet(path=ms_path)

    selfcal_ms = run_selfcal(
        ms,
        calmode="ap",
        gaintype="T",
        interval="10s",
        split_data=True,
        interactive=True,
        refant="6",
    )

    assert selfcal_ms.path == ms.path
    assert selfcal_ms.nspws == 1


def test_combine_spws(temp_environment):
    ms_path = temp_environment["onespw"]
    ms = MeasurementSet(ms_path)

    combined = combine_spws(ms)

    assert combined.path == ms_path.with_suffix(".dstools-temp.comb.ms")


def test_extract_baseline(temp_environment):
    ms = MeasurementSet(temp_environment["minimal"])
    baseline = 0, (0, 5)
    data = extract_baseline(
        ms,
        baseline=baseline,
        datacolumn="DATA",
    )

    assert data["baseline"] == 0
    assert np.allclose(data["data_idx"], np.array([0, 1]))
    assert np.all(data["data"] == 0 + 0j)
    assert np.all(data["flags"])


@pytest.mark.parametrize("ncpus", [1, 2])
def test_extract_baselines_1baseline(ncpus, mocker, temp_environment):
    mocker.patch("dstools.ms.get_available_cpus", return_value=ncpus)

    ms = MeasurementSet(temp_environment["minimal"])
    vis, flags, uvws = extract_baselines(
        ms,
        datacolumn="DATA",
    )

    assert vis.shape == (1, 2, 3, 4)
    assert np.all(vis[:, :, :, 0] == 0 + 0j)
    assert np.all(np.isnan(vis[:, :, :, 1:]))
    assert np.all(flags)
    assert np.all(uvws > 0)


def test_swap_xy_feeds_array():
    nvis = 100
    nchan = 5

    data = np.random.random((nvis, nchan, 4)) + 1j * np.random.random((nvis, nchan, 4))
    swapped = swap_xy_feeds(data)

    assert np.all(data[:, :, 0] == swapped[:, :, 3])
    assert np.all(data[:, :, 1] == swapped[:, :, 2])
    assert np.all(data[:, :, 2] == swapped[:, :, 1])
    assert np.all(data[:, :, 3] == swapped[:, :, 0])


def test_swap_xy_feeds_ms(temp_environment):
    ms = MeasurementSet(temp_environment["mkt_3c286"])
    swapped_ms = MeasurementSet(temp_environment["mkt_3c286_fix"])

    swapped_ms.swap_xy_feeds(datacolumn="DATA")

    with ms.open_table() as t:
        data = t.getcol("DATA")

    with swapped_ms.open_table() as t:
        swapped = t.getcol("DATA")

    # Check that data axes have been swapped
    assert np.all(data[:, :, 0] == swapped[:, :, 3])
    assert np.all(data[:, :, 1] == swapped[:, :, 2])
    assert np.all(data[:, :, 2] == swapped[:, :, 1])
    assert np.all(data[:, :, 3] == swapped[:, :, 0])

    # Check that new feed angle is correct
    # This test portion of the band should align to +37 degrees
    # so that the median across the band is +33 degrees
    q = (swapped[:, :, 0] - swapped[:, :, 3]) / 2
    u = (swapped[:, :, 1] + swapped[:, :, 2]) / 2
    pa = 0.5 * np.atan2(u.real, q.real)

    assert np.round(np.rad2deg(np.nanmedian(pa)), 1) == 37.0


def test_rotate_linear_feeds():
    nvis = 100
    nchan = 1

    # Create pure Stokes Q signal in linear basis
    data = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    data[:, :, 0] = 2

    # Reshape visibility matrices to (nvis, 2, 2)
    V = np.zeros((nvis, 2, 2), dtype=np.complex128)
    V[:, 0, 0] = data[:, 0, 0]
    V[:, 0, 1] = data[:, 0, 1]
    V[:, 1, 0] = data[:, 0, 2]
    V[:, 1, 1] = data[:, 0, 3]

    # Construct rotation matrix
    chi = np.linspace(0, np.pi / 2, nvis).astype(np.float64)
    R = np.empty((nvis, 2, 2), dtype=np.complex128)
    cos_chi = np.cos(chi)
    sin_chi = np.sin(chi)
    R[:, 0, 0] = cos_chi
    R[:, 0, 1] = sin_chi
    R[:, 1, 0] = -sin_chi
    R[:, 1, 1] = cos_chi

    # Rotate to simulate observation
    V_rot = R @ V @ np.transpose(R.conj(), axes=(0, 2, 1))

    # Reshape back to 4 x nchan x nvis
    rotated = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    rotated[:, 0, 0] = V_rot[:, 0, 0]
    rotated[:, 0, 1] = V_rot[:, 0, 1]
    rotated[:, 0, 2] = V_rot[:, 1, 0]
    rotated[:, 0, 3] = V_rot[:, 1, 1]

    # Test de-rotation
    corrected = rotate_linear_feeds(rotated, chi)

    # Check recovery of original pure Q signal
    assert np.allclose(corrected, data, atol=1e-8)


def test_rotate_circular_feeds():
    nvis = 100
    nchan = 1

    # Create pure Stokes Q signal in circular basis
    # RR = 1, LL = 1, RL = 0, LR = 0
    data = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    data[:, 0, 0] = 1  # RR
    data[:, 0, 3] = 1  # LL

    # Rotate to simulate observation with circular feeds
    rotated = data.copy()

    chi = np.linspace(0, np.pi / 2, nvis).astype(np.float64)
    rotated[:, 0, 1] = data[:, 0, 1] * np.exp(2j * chi)  # RL
    rotated[:, 0, 2] = data[:, 0, 2] * np.exp(-2j * chi)  # LR

    # Test de-rotation
    corrected = rotate_circular_feeds(rotated, chi)

    # Check recovery of original pure Q signal
    assert np.allclose(corrected, data, atol=1e-8)


def test_correct_feed_rotation_askap_does_not_rotate(temp_environment, caplog):
    path = temp_environment["askap"]
    ms = MeasurementSet(path)

    ms.correct_feed_rotation(datacolumn="DATA")

    assert "Will not apply" in caplog.text


def test_correct_feed_rotation_linear(temp_environment, caplog):
    path = temp_environment["minimal"]
    ms = MeasurementSet(path)

    with caplog.at_level(logging.INFO):
        ms.correct_feed_rotation(datacolumn="DATA")

    assert "Correcting" in caplog.text


def test_correct_feed_rotation_circular(temp_environment, caplog):
    path = temp_environment["vla"]
    ms = MeasurementSet(path)

    with caplog.at_level(logging.INFO):
        ms.correct_feed_rotation(datacolumn="DATA")

    assert "Correcting" in caplog.text


@pytest.mark.parametrize(
    "pol_indices",
    [
        [5, 6, 7, 8],
        [9, 10, 11, 12],
    ],
)
def test_get_polslice_all_pols(pol_indices):
    polslice = get_polslice(pol_indices)

    assert polslice == slice(0, 4)


@pytest.mark.parametrize(
    "pol_indices",
    [
        [5, 8],
        [9, 12],
    ],
)
def test_get_polslice_parallel_hand_pols(pol_indices):
    polslice = get_polslice(pol_indices)

    assert polslice == slice(0, 4, 3)


@pytest.mark.parametrize(
    "pol_indices",
    [
        [5],
        [9],
    ],
)
def test_get_polslice_single_pol(pol_indices):
    polslice = get_polslice(pol_indices)

    assert polslice == slice(0, 1)


@pytest.mark.parametrize(
    "pol_indices",
    [
        [1, 2, 3, 4],  # Stokes IQUV
        [9, 12, 10, 11],  # XX, YY, XY, YX order
        [10],  # just YY
        [0],  # unregistered corr_type
    ],
)
def test_get_polslice_single_pol(pol_indices):
    with pytest.raises(DataError):
        get_polslice(pol_indices)
