from pathlib import Path
from unittest.mock import Mock

import astropy.units as u
import h5py
import numpy as np
import pytest
from click.testing import CliRunner

pytest.importorskip("casacore.tables")
pytestmark = pytest.mark.fullstack

from dstools.cli import (  # noqa: E402
    askap_preprocess,
    derotate_feeds,
    extract_ds,
    insert_model,
    selfcal,
    subtract_model,
)
from dstools.ms import MeasurementSet  # noqa: E402
from dstools.utils import DataError  # noqa: E402


def test_askap_preprocess_aborts_when_already_fixed(monkeypatch):
    runner = CliRunner()
    monkeypatch.setattr(askap_preprocess, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        askap_preprocess, "tableexists", Mock(return_value=True), raising=False
    )
    filtered_fixms = Mock()
    monkeypatch.setattr(
        askap_preprocess, "filtered_fixms", filtered_fixms, raising=False
    )

    result = runner.invoke(askap_preprocess.main, ["dummy.ms"])

    assert result.exit_code == 1
    filtered_fixms.assert_not_called()


def test_askap_preprocess_runs_filtered_fixms(monkeypatch):
    runner = CliRunner()
    monkeypatch.setattr(askap_preprocess, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        askap_preprocess, "tableexists", Mock(return_value=False), raising=False
    )
    filtered_fixms = Mock()
    monkeypatch.setattr(
        askap_preprocess, "filtered_fixms", filtered_fixms, raising=False
    )

    result = runner.invoke(askap_preprocess.main, ["dummy.ms"])

    assert result.exit_code == 0
    filtered_fixms.assert_called_once_with(Path("dummy.ms"))


def test_extract_ds_reports_missing_column(monkeypatch, tmp_path):
    runner = CliRunner()
    fake_ms = Mock()
    fake_ms.column_exists.return_value = False

    monkeypatch.setattr(extract_ds, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        extract_ds, "MeasurementSet", Mock(return_value=fake_ms), raising=False
    )

    result = runner.invoke(extract_ds.main, ["dummy.ms", str(tmp_path / "out.ds")])

    assert result.exit_code == 1
    fake_ms.column_exists.assert_called_once_with("DATA")


def test_selfcal_maps_combine_pols_to_gaintype(monkeypatch):
    runner = CliRunner()
    fake_ms = object()
    run_selfcal = Mock()

    monkeypatch.setattr(selfcal, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        selfcal, "MeasurementSet", Mock(return_value=fake_ms), raising=False
    )
    monkeypatch.setattr(selfcal, "run_selfcal", run_selfcal, raising=False)

    result = runner.invoke(selfcal.main, ["--no-combine-pols", "dummy.ms"])

    assert result.exit_code == 0
    assert run_selfcal.call_args.kwargs["gaintype"] == "G"


def test_subtract_model_forwards_split_ms(monkeypatch):
    runner = CliRunner()
    fake_ms = Mock()

    monkeypatch.setattr(subtract_model, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        subtract_model,
        "MeasurementSet",
        Mock(return_value=fake_ms),
        raising=False,
    )

    result = runner.invoke(subtract_model.main, ["--split-ms", "dummy.ms"])

    assert result.exit_code == 0
    fake_ms.subtract_model.assert_called_once_with(split_ms=True)


def test_subtract_model_reports_data_error_when_model_missing(monkeypatch):
    runner = CliRunner()
    fake_ms = Mock()
    fake_ms.subtract_model.side_effect = DataError("MODEL_DATA missing")

    monkeypatch.setattr(subtract_model, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        subtract_model,
        "MeasurementSet",
        Mock(return_value=fake_ms),
        raising=False,
    )

    result = runner.invoke(subtract_model.main, ["dummy.ms"])

    assert result.exit_code == 1
    fake_ms.subtract_model.assert_called_once_with(split_ms=False)


def test_derotate_feeds_cli_updates_rotated_ms(copy_paths_into_workspace, ms_sources):
    ms_path = copy_paths_into_workspace(ms_sources, ("rotated",))["rotated"]
    before = MeasurementSet(ms_path).getcolumn("DATA").copy()

    runner = CliRunner()
    result = runner.invoke(derotate_feeds.main, [str(ms_path)])

    after = MeasurementSet(ms_path).getcolumn("DATA")

    assert result.exit_code == 0
    assert not np.allclose(before, after)


def test_subtract_model_cli_success(workspace_onespw_ms, casa_task_mocks):
    ms = MeasurementSet(workspace_onespw_ms)
    runner = CliRunner()

    result = runner.invoke(subtract_model.main, [str(workspace_onespw_ms)])

    corrected = ms.getcolumn("CORRECTED_DATA")
    model = ms.getcolumn("MODEL_DATA")
    data = ms.getcolumn("DATA")

    assert result.exit_code == 0
    assert np.allclose(corrected, data - model)


def test_selfcal_cli_creates_new_selfcal_round(workspace_onespw_ms, casa_task_mocks):
    runner = CliRunner()
    expected = workspace_onespw_ms.with_suffix(".selfcal1.ms")

    result = runner.invoke(
        selfcal.main, ["--no-interactive", "--split-data", str(workspace_onespw_ms)]
    )

    assert result.exit_code == 0
    assert expected.exists()


def test_insert_model_requires_existing_model_dir():
    runner = CliRunner()

    result = runner.invoke(insert_model.main, ["missing-model", "dummy.ms"])

    assert result.exit_code == 1


def test_extract_ds_cli_writes_output(workspace_onespw_ms, tmp_path, casa_task_mocks):
    runner = CliRunner()
    output = tmp_path / "extracted.ds"

    result = runner.invoke(extract_ds.main, [str(workspace_onespw_ms), str(output)])

    assert result.exit_code == 0
    with h5py.File(output, "r") as handle:
        assert "flux" in handle
        assert "time" in handle
        assert "frequency" in handle
        assert "uvdist" in handle


def test_insert_model_reports_missing_model_dir():
    runner = CliRunner()

    result = runner.invoke(insert_model.main, ["missing-model", "dummy.ms"])

    assert result.exit_code == 1


def test_insert_model_automatic_mask_path(monkeypatch, tmp_path):
    runner = CliRunner()
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    fake_ms = Mock(phasecentre="PHASECENTRE")
    fake_model = Mock()
    fake_model.get_circular_mask.return_value = "MASK"
    parsed = object()

    monkeypatch.setattr(insert_model, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        insert_model, "MeasurementSet", Mock(return_value=fake_ms), raising=False
    )
    monkeypatch.setattr(
        insert_model, "WSCleanModel", Mock(return_value=fake_model), raising=False
    )
    monkeypatch.setattr(
        insert_model, "parse_coordinates", Mock(return_value=parsed), raising=False
    )

    result = runner.invoke(
        insert_model.main,
        [
            "--mask-pos",
            "00:00:00",
            "-30:00:00",
            "--mask-radius",
            "3",
            str(model_dir),
            "dummy.ms",
        ],
    )

    assert result.exit_code == 0
    insert_model.parse_coordinates.assert_called_once_with(("00:00:00", "-30:00:00"))
    fake_model.get_circular_mask.assert_called_once_with(parsed, 3.0 * u.arcsec)
    fake_model.apply_mask.assert_called_once_with("MASK")
    fake_model.insert_into.assert_called_once_with(fake_ms)


def test_insert_model_interactive_mask(monkeypatch, tmp_path):
    runner = CliRunner()
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    fake_ms = Mock()
    fake_model = Mock()
    fake_model.get_interactive_mask.return_value = "INTERACTIVE_MASK"

    monkeypatch.setattr(insert_model, "HAS_CASA_SUPPORT", True)
    monkeypatch.setattr(
        insert_model, "MeasurementSet", Mock(return_value=fake_ms), raising=False
    )
    monkeypatch.setattr(
        insert_model, "WSCleanModel", Mock(return_value=fake_model), raising=False
    )

    result = runner.invoke(
        insert_model.main, ["--interactive", str(model_dir), "dummy.ms"]
    )

    assert result.exit_code == 0
    fake_model.get_interactive_mask.assert_called_once_with()
    fake_model.get_circular_mask.assert_not_called()
    fake_model.apply_mask.assert_called_once_with("INTERACTIVE_MASK")
    fake_model.insert_into.assert_called_once_with(fake_ms)
