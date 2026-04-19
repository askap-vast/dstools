import logging

import numpy as np
import pytest
from astropy.coordinates import SkyCoord

pytest.importorskip("casacore.tables")
pytestmark = pytest.mark.fullstack

from dstools.ms import (  # noqa: E402
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
from dstools.utils import DataError  # noqa: E402

MS_PROPERTIES = [
    ("nspws", 1),
    ("nbaselines", 3),
    ("integrations", 4),
    ("nchannels", 21),
    ("npols", 1),
    ("dimensions", (3, 4, 21, 1)),
    ("telescope", "ATCA"),
    ("feedtype", "linear"),
]
COLUMNS = ["DATA", "CORRECTED_DATA", "MODEL_DATA", "SIGMA_SPECTRUM", "WEIGHT_SPECTRUM"]


@pytest.mark.parametrize("prop, val", MS_PROPERTIES)
def test_ms_basic_properties(prop, val, ms):
    assert getattr(ms, prop) == val


def test_ms_opens_with_str_path(ms_path):
    ms = MeasurementSet(path=str(ms_path))

    assert ms.path == ms_path
    assert ms.column_exists("DATA")


def test_ms_phasecentre(ms):
    assert ms.phasecentre.to_string("hmsdms") == "18h45m13.62000576s -63d57m34.3899846s"


@pytest.mark.parametrize("column", COLUMNS)
def test_ms_column_exists(ms, column):
    assert ms.column_exists(column)


def test_combine_multi_spw_fails_if_not_split(ms):
    with pytest.raises(DataError):
        ms._combine_multi_spw()


def test_nspw_conversion_1to1(ms):
    assert ms.nspws == 1
    ms.to_nspws(1)
    assert ms.nspws == 1


def test_nspw_conversion_1to2(mocked_ms):
    assert mocked_ms.nspws == 1

    mocked_ms.to_nspws(2)
    assert mocked_ms.nspws == 2


def test_nspw_conversion_2to1(copy_paths_into_workspace, ms_sources, casa_task_mocks):
    copied = copy_paths_into_workspace(ms_sources, ("onespw", "twospw"))
    ms = MeasurementSet(path=copied["twospw"])
    ms.original_path = copied["onespw"]

    assert ms.nspws == 2

    ms.to_nspws(1)
    assert ms.nspws == 1


def test_nspw_conversion_sequence(mocked_ms):
    assert mocked_ms.nspws == 1

    mocked_ms.to_nspws(2)
    assert mocked_ms.nspws == 2

    mocked_ms.to_nspws(1)
    assert mocked_ms.nspws == 1


def test_invalid_nspw_conversion_raises_error(ms):
    with pytest.raises(ValueError):
        ms.to_nspws(-1)


def test_ms_rotate_phasecentre_inplace(ms):
    with pytest.raises(NotImplementedError):
        ms.rotate_phasecentre(ms.phasecentre, inplace=True)


def test_ms_rotate_phasecentre_same_location_does_nothing(ms):
    assert ms.phasecentre.to_string("hmsdms") == "18h45m13.62000576s -63d57m34.3899846s"
    assert ms.rotate_phasecentre(ms.phasecentre) == ms


def test_ms_increment_selfcal_round_initial(ms):
    assert ms.increment_selfcal_round() == ms.path.with_suffix(".selfcal1.ms")


def test_ms_increment_selfcal_round_multiple(ms_path):
    ms = MeasurementSet(path=ms_path)
    ms.path = ms.increment_selfcal_round()

    assert ms.increment_selfcal_round() == ms_path.with_suffix(".selfcal2.ms")


def test_combine_spws(workspace_onespw_ms, casa_task_mocks):
    ms = MeasurementSet(workspace_onespw_ms)
    combined = combine_spws(ms)

    assert combined.path == workspace_onespw_ms.with_suffix(".dstools-temp.comb.ms")


def test_ms_not_exists_raises_error():
    with pytest.raises(FileNotFoundError):
        MeasurementSet("faketable.ms")


def test_swap_xy_feeds_array():
    data = np.random.random((100, 5, 4)) + 1j * np.random.random((100, 5, 4))
    swapped = swap_xy_feeds(data)

    assert np.all(data[:, :, 0] == swapped[:, :, 3])
    assert np.all(data[:, :, 1] == swapped[:, :, 2])
    assert np.all(data[:, :, 2] == swapped[:, :, 1])
    assert np.all(data[:, :, 3] == swapped[:, :, 0])


def test_rotate_linear_feeds():
    nvis, nchan = 100, 1
    data = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    data[:, :, 0] = 2
    V = np.zeros((nvis, 2, 2), dtype=np.complex128)
    V[:, 0, 0], V[:, 0, 1], V[:, 1, 0], V[:, 1, 1] = (
        data[:, 0, 0],
        data[:, 0, 1],
        data[:, 0, 2],
        data[:, 0, 3],
    )
    chi = np.linspace(0, np.pi / 2, nvis).astype(np.float64)
    R = np.empty((nvis, 2, 2), dtype=np.complex128)
    cos_chi, sin_chi = np.cos(chi), np.sin(chi)
    R[:, 0, 0], R[:, 0, 1], R[:, 1, 0], R[:, 1, 1] = (
        cos_chi,
        sin_chi,
        -sin_chi,
        cos_chi,
    )
    V_rot = R @ V @ np.transpose(R.conj(), axes=(0, 2, 1))
    rotated = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    rotated[:, 0, 0], rotated[:, 0, 1], rotated[:, 0, 2], rotated[:, 0, 3] = (
        V_rot[:, 0, 0],
        V_rot[:, 0, 1],
        V_rot[:, 1, 0],
        V_rot[:, 1, 1],
    )
    corrected = rotate_linear_feeds(rotated, chi)

    assert np.allclose(corrected, data, atol=1e-8)


def test_rotate_circular_feeds():
    nvis, nchan = 100, 1
    data = np.zeros((nvis, nchan, 4), dtype=np.complex128)
    data[:, 0, 0] = 1
    data[:, 0, 3] = 1
    rotated = data.copy()
    chi = np.linspace(0, np.pi / 2, nvis).astype(np.float64)
    rotated[:, 0, 1] = data[:, 0, 1] * np.exp(2j * chi)
    rotated[:, 0, 2] = data[:, 0, 2] * np.exp(-2j * chi)
    corrected = rotate_circular_feeds(rotated, chi)

    assert np.allclose(corrected, data, atol=1e-8)


@pytest.mark.parametrize("pol_indices", [[5, 6, 7, 8], [9, 10, 11, 12]])
def test_get_polslice_all_pols(pol_indices):
    assert get_polslice(pol_indices) == slice(0, 4)


@pytest.mark.parametrize("pol_indices", [[5, 8], [9, 12]])
def test_get_polslice_parallel_hand_pols(pol_indices):
    assert get_polslice(pol_indices) == slice(0, 4, 3)


@pytest.mark.parametrize("pol_indices", [[5], [9]])
def test_get_polslice_single_pol(pol_indices):
    assert get_polslice(pol_indices) == slice(0, 1)


@pytest.mark.parametrize("pol_indices", [[1, 2, 3, 4], [9, 12, 10, 11], [10], [0]])
def test_get_polslice_unsupported_raises_error(pol_indices):
    with pytest.raises(DataError):
        get_polslice(pol_indices)


def test_get_reference_antenna_auto(mocked_ms):
    assert mocked_ms.get_reference_antenna(interactive=False) == "1"


def test_get_reference_antenna_interactive(mocked_ms, mocker):
    mocker.patch("builtins.input", side_effect="6")
    assert mocked_ms.get_reference_antenna(interactive=True) == "6"


def test_get_reference_antenna_invalid_selection_prompts_user(mocked_ms, capfd, mocker):
    mocker.patch("builtins.input", side_effect=["41", "6"])
    refant = mocked_ms.get_reference_antenna(interactive=True)
    out, _ = capfd.readouterr()

    assert "Reference antenna must be in:" in out
    assert refant == "6"


def test_ms_average_baselines(mocked_ms):
    assert "baseavg" in str(mocked_ms.average_baselines().path)


def test_ms_rotate_phasecentre_new_location(mocked_ms):
    original_centre = "18h45m13.62000576s -63d57m34.3899846s"

    assert mocked_ms.phasecentre.to_string("hmsdms") == original_centre

    new_centre = SkyCoord(ra="18h45m00s", dec="-63d57m00s", unit="hourangle,deg")
    baseavg = mocked_ms.rotate_phasecentre(new_centre)

    assert baseavg.phasecentre.to_string("hmsdms") == new_centre.to_string("hmsdms")


def test_subtract_model_no_model_column_raises_error(
    copy_paths_into_workspace,
    ms_sources,
    casa_task_mocks,
):
    ms = MeasurementSet(
        path=copy_paths_into_workspace(ms_sources, ("rotated",))["rotated"]
    )

    with pytest.raises(DataError):
        ms.subtract_model()


def test_subtract_model_no_split(mocked_ms):
    subtracted = mocked_ms.subtract_model()
    corrected = subtracted.getcolumn("CORRECTED_DATA")
    model = subtracted.getcolumn("MODEL_DATA")
    data = subtracted.getcolumn("DATA")

    assert subtracted.path == mocked_ms.path
    assert np.allclose(corrected, data - model)


def test_subtract_model_split(mocked_ms):
    no_split = mocked_ms.subtract_model()
    split = mocked_ms.subtract_model(split_ms=True)

    assert np.allclose(no_split.getcolumn("CORRECTED_DATA"), split.getcolumn("DATA"))


def test_ms_calc_flag_statistics(mocked_ms):
    df = mocked_ms.calc_flag_statistics()

    assert "antenna" in df.columns
    assert "flagged" in df.columns
    assert "total" in df.columns
    assert "percentage" in df.columns


def test_ms_solve_gains(mocked_ms):
    mocked_ms.solve_gains(interval="10s", calmode="ap", gaintype="T")
    assert mocked_ms.caltable.path == mocked_ms.path.with_suffix(".cal")


def test_ms_split_selfcal_round(workspace_onespw_ms, casa_task_mocks):
    ms = MeasurementSet(path=workspace_onespw_ms)

    assert ms.split_selfcal_round().path == workspace_onespw_ms.with_suffix(
        ".selfcal1.ms"
    )


def test_run_selfcal_bad_interval_raises_error(workspace_onespw_ms, casa_task_mocks):
    ms = MeasurementSet(path=workspace_onespw_ms)
    with pytest.raises(ValueError):
        run_selfcal(
            ms,
            calmode="ap",
            gaintype="T",
            interval="10",
            split_data=False,
            interactive=False,
        )


def test_run_selfcal_no_model_raises_error(
    copy_paths_into_workspace, ms_sources, casa_task_mocks
):
    ms = MeasurementSet(
        path=copy_paths_into_workspace(ms_sources, ("rotated",))["rotated"]
    )
    with pytest.raises(DataError):
        run_selfcal(
            ms,
            calmode="ap",
            gaintype="T",
            interval="10s",
            split_data=False,
            interactive=False,
        )


def test_extract_baseline(copy_paths_into_workspace, ms_sources):
    ms = MeasurementSet(copy_paths_into_workspace(ms_sources, ("minimal",))["minimal"])
    data = extract_baseline(ms, baseline=(0, (0, 5)), datacolumn="DATA")

    assert data["baseline"] == 0
    assert np.allclose(data["data_idx"], np.array([0, 1]))
    assert np.all(data["data"] == 0 + 0j)
    assert np.all(data["flags"])


@pytest.mark.parametrize("ncpus", [1, 2])
def test_extract_baselines_1baseline(
    ncpus,
    mocker,
    copy_paths_into_workspace,
    ms_sources,
):
    mocker.patch("dstools.ms.get_available_cpus", return_value=ncpus)
    ms = MeasurementSet(copy_paths_into_workspace(ms_sources, ("minimal",))["minimal"])
    vis, flags, uvws = extract_baselines(ms, datacolumn="DATA")

    assert vis.shape == (1, 2, 3, 4)
    assert np.all(vis[:, :, :, 0] == 0 + 0j)
    assert np.all(np.isnan(vis[:, :, :, 1:]))
    assert np.all(flags)
    assert np.all(uvws > 0)


def test_run_selfcal_cal_good(workspace_onespw_ms, casa_task_mocks):
    ms = MeasurementSet(path=workspace_onespw_ms)
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


def test_run_selfcal_cal_bad(mocker, workspace_onespw_ms, casa_task_mocks):
    mocker.patch.object(MeasurementSet, "to_nspws", return_value=None)
    mocker.patch.object(CalTable, "plot_solutions", return_value=None)
    mocker.patch("dstools.ms.prompt", return_value=False)
    mocker.patch("matplotlib.pyplot.show")

    ms = MeasurementSet(path=workspace_onespw_ms)
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


def test_swap_xy_feeds_ms(copy_paths_into_workspace, ms_sources):
    copied = copy_paths_into_workspace(
        {"original": ms_sources["mkt_3c286"], "fixed": ms_sources["mkt_3c286"]},
        ("original", "fixed"),
        destination_names={"fixed": "3C286.MKT_UHF.xyswapped.fix.ms"},
    )
    ms = MeasurementSet(copied["original"])
    swapped_ms = MeasurementSet(copied["fixed"])
    swapped_ms.swap_xy_feeds(datacolumn="DATA")

    with ms.open_table() as t:
        data = t.getcol("DATA")
    with swapped_ms.open_table() as t:
        swapped = t.getcol("DATA")

    assert np.all(data[:, :, 0] == swapped[:, :, 3])
    assert np.all(data[:, :, 1] == swapped[:, :, 2])
    assert np.all(data[:, :, 2] == swapped[:, :, 1])
    assert np.all(data[:, :, 3] == swapped[:, :, 0])

    q = (swapped[:, :, 0] - swapped[:, :, 3]) / 2
    u = (swapped[:, :, 1] + swapped[:, :, 2]) / 2
    pa = 0.5 * np.atan2(u.real, q.real)
    assert np.round(np.rad2deg(np.nanmedian(pa)), 1) == 37.0


def test_correct_feed_rotation_askap_does_not_rotate(
    copy_paths_into_workspace,
    ms_sources,
    caplog,
):
    ms = MeasurementSet(copy_paths_into_workspace(ms_sources, ("askap",))["askap"])
    ms.correct_feed_rotation(datacolumn="DATA")

    assert "Will not apply" in caplog.text


def test_correct_feed_rotation_linear(
    copy_paths_into_workspace,
    ms_sources,
    caplog,
):
    ms = MeasurementSet(copy_paths_into_workspace(ms_sources, ("rotated",))["rotated"])
    with caplog.at_level(logging.INFO):
        ms.correct_feed_rotation(datacolumn="DATA")

    assert "Correcting" in caplog.text


def test_correct_feed_rotation_circular(
    copy_paths_into_workspace,
    ms_sources,
    caplog,
):
    ms = MeasurementSet(copy_paths_into_workspace(ms_sources, ("vla",))["vla"])
    with caplog.at_level(logging.INFO):
        ms.correct_feed_rotation(datacolumn="DATA")

    assert "Correcting" in caplog.text
