import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord

import dstools

PACKAGE_ROOT = Path(dstools.__path__[0]).parent
TEST_DATA_ROOT = PACKAGE_ROOT / "tests" / "data"
TEST_MPLCONFIGDIR = Path(tempfile.gettempdir()) / "dstools-mpl"

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", str(TEST_MPLCONFIGDIR))
TEST_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)


def _copy_path(src: Path, dst: Path) -> Path:
    if src.is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)

    return dst


def _copy_named_paths(
    sources: dict[str, Path],
    workspace: Path,
    names: tuple[str, ...],
    destination_names: dict[str, str] | None = None,
) -> dict[str, Path]:
    copied = {}
    destination_names = destination_names or {}

    for name in names:
        dst_name = destination_names.get(name, sources[name].name)
        copied[name] = _copy_path(sources[name], workspace / dst_name)

    return copied


def _replace_path(src: Path, dst: Path) -> Path:
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()

    return _copy_path(src, dst)


def _set_ms_phasecentre(ms_path: Path, phasecenter: str) -> None:
    casacore = pytest.importorskip("casacore.tables")
    position = SkyCoord(phasecenter.replace("J2000 ", ""), unit=("hourangle", "deg"))

    with casacore.table(
        (ms_path / "FIELD").as_posix(),
        readonly=False,
        ack=False,
    ) as table:
        for column in ("PHASE_DIR", "REFERENCE_DIR", "DELAY_DIR"):
            if column not in table.colnames():
                continue

            data = table.getcol(column)
            if data.shape[0] == 2:
                data[0, ...] = position.ra.rad
                data[1, ...] = position.dec.rad
            elif data.shape[-1] == 2:
                data[..., 0] = position.ra.rad
                data[..., 1] = position.dec.rad
            else:
                raise ValueError(f"Unexpected FIELD/{column} shape: {data.shape}")

            table.putcol(column, data)


def _set_data_from_corrected(ms_path: Path) -> None:
    casacore = pytest.importorskip("casacore.tables")

    with casacore.table(ms_path.as_posix(), readonly=False, ack=False) as table:
        corrected = table.getcol("CORRECTED_DATA")
        table.putcol("DATA", corrected)


@pytest.fixture
def disable_jit(monkeypatch):
    monkeypatch.setenv("NUMBA_DISABLE_JIT", "1")

    return


@pytest.fixture
def dispersed_pulse():
    return np.load(
        TEST_DATA_ROOT / "ds" / "dispersed_pulse_dm6000.npy", allow_pickle=True
    )


@pytest.fixture
def ms_path():
    return TEST_DATA_ROOT / "msets" / "fred.atca.ms"


@pytest.fixture
def ms_sources(ms_path):
    return {
        "onespw": ms_path,
        "twospw": ms_path.with_suffix(".2spw.ms"),
        "minimal": ms_path.with_suffix(".minimal.ms"),
        "rotated": ms_path.with_suffix(".dstools-temp.rotated.ms"),
        "averaged": ms_path.with_suffix(".dstools-temp.baseavg.ms"),
        "subtracted": ms_path.with_suffix(".subtracted.ms"),
        "vla": TEST_DATA_ROOT / "msets" / "gpm.vla.ms",
        "mkt_3c286": TEST_DATA_ROOT / "msets" / "3C286.MKT_UHF.xyswapped.ms",
        "askap": TEST_DATA_ROOT / "msets" / "j1755.askap.ms",
    }


@pytest.fixture
def caltable_sources(ms_path):
    return {"fred": ms_path.with_suffix(".cal")}


@pytest.fixture
def im_paths():
    return {
        "image": TEST_DATA_ROOT / "images" / "test-MFS-I-image.fits",
        "model": TEST_DATA_ROOT / "images" / "test-MFS-I-model.fits",
        "residual": TEST_DATA_ROOT / "images" / "test-MFS-I-residual.fits",
    }


@pytest.fixture
def image_sources():
    return {
        "model_dir": TEST_DATA_ROOT / "images",
        "pb": TEST_DATA_ROOT / "images" / "fred.atca.pb.fits",
        "target_mask": TEST_DATA_ROOT / "images" / "target_mask.fits",
        "clean_mask": TEST_DATA_ROOT / "images" / "clean_mask.fits",
        "final_mask": TEST_DATA_ROOT / "images" / "final_mask.fits",
    }


@pytest.fixture
def temp_workspace(tmp_path):
    return tmp_path


@pytest.fixture
def workspace_onespw_ms(copy_paths_into_workspace, ms_sources):
    return copy_paths_into_workspace(
        ms_sources,
        ("onespw",),
        destination_names={"onespw": "test.ms"},
    )["onespw"]


@pytest.fixture
def imaging_workspace(copy_paths_into_workspace, image_sources):
    copied = copy_paths_into_workspace(
        image_sources,
        ("model_dir", "pb", "target_mask", "clean_mask", "final_mask"),
        destination_names={
            "model_dir": "model",
            "pb": "test.pb.fits",
            "target_mask": "target_mask.fits",
            "clean_mask": "clean_mask.fits",
            "final_mask": "final_mask.fits",
        },
    )

    return {
        "model": copied["model_dir"],
        "pb": copied["pb"],
        "target_mask": copied["target_mask"],
        "clean_mask": copied["clean_mask"],
        "final_mask": copied["final_mask"],
    }


@pytest.fixture
def copy_paths_into_workspace(temp_workspace):
    def _copy(
        sources: dict[str, Path],
        names: tuple[str, ...],
        destination_names: dict[str, str] | None = None,
    ) -> dict[str, Path]:
        return _copy_named_paths(
            sources=sources,
            workspace=temp_workspace,
            names=names,
            destination_names=destination_names,
        )

    return _copy


@pytest.fixture
def casa_task_mocks(
    mocker,
    ms_sources,
    caltable_sources,
    image_sources,
):
    flagstats = {
        "antenna": {
            "1": {"flagged": 240.0, "total": 1620.0},
            "6": {"flagged": 244.0, "total": 1620.0},
        },
    }

    def _copy_ms_output(src: str, dst: str, fallback: Path | None = None) -> None:
        src_path = Path(src)
        dst_path = Path(dst)
        _replace_path(fallback or src_path, dst_path)

    def _mstransform_side_effect(vis, outputvis, **kwargs):
        if kwargs.get("nspw") == 2:
            _copy_ms_output(vis, outputvis, fallback=ms_sources["twospw"])
        elif str(outputvis).endswith(".dstools-temp.baseavg.ms"):
            _copy_ms_output(vis, outputvis, fallback=ms_sources["averaged"])
        else:
            _copy_ms_output(vis, outputvis)

    def _cvel_side_effect(vis, outputvis, **kwargs):
        _copy_ms_output(vis, outputvis, fallback=ms_sources["onespw"])

    def _phaseshift_side_effect(vis, outputvis, phasecenter, **kwargs):
        _copy_ms_output(vis, outputvis)
        _set_ms_phasecentre(Path(outputvis), phasecenter)

    def _split_side_effect(vis, outputvis, datacolumn, **kwargs):
        _copy_ms_output(vis, outputvis)
        if datacolumn == "corrected":
            _set_data_from_corrected(Path(outputvis))

    def _gaincal_side_effect(vis, caltable, **kwargs):
        _replace_path(caltable_sources["fred"], Path(caltable))

    def _tablecopy_side_effect(tablename, newtablename):
        _replace_path(Path(tablename), Path(newtablename))

    def _exportuvfits_side_effect(imagename, fitsimage, **kwargs):
        _replace_path(image_sources["pb"], Path(fitsimage))

    mocker.patch("dstools.ms.mstransform", side_effect=_mstransform_side_effect)
    mocker.patch("dstools.ms.cvel", side_effect=_cvel_side_effect)
    mocker.patch("dstools.ms.phaseshift", side_effect=_phaseshift_side_effect)
    mocker.patch("dstools.ms.uvsub")
    mocker.patch("dstools.ms.flagdata", return_value=flagstats)
    mocker.patch("dstools.ms.split", side_effect=_split_side_effect)
    mocker.patch("dstools.ms.gaincal", side_effect=_gaincal_side_effect)
    mocker.patch("dstools.ms.applycal")
    mocker.patch("dstools.imaging.tclean")
    mocker.patch("dstools.imaging.exportuvfits", side_effect=_exportuvfits_side_effect)
    mocker.patch("dstools.imaging.parse_stdout_stderr")
    mocker.patch("dstools.ms.tablecopy", side_effect=_tablecopy_side_effect)
    mocker.patch("os.system")
    mocker.patch("os.chdir")

    return flagstats


@pytest.fixture
def ms(workspace_onespw_ms):
    MeasurementSet = pytest.importorskip("dstools.ms").MeasurementSet

    return MeasurementSet(workspace_onespw_ms)


@pytest.fixture
def mocked_ms(workspace_onespw_ms, casa_task_mocks):
    MeasurementSet = pytest.importorskip("dstools.ms").MeasurementSet

    return MeasurementSet(workspace_onespw_ms)


@pytest.fixture
def ds_paths():
    return {
        "atca_pulse": str(TEST_DATA_ROOT / "ds" / "fred.atca.pulse.ds"),
        "atca_calscan": str(TEST_DATA_ROOT / "ds" / "fred.atca.calscan.ds"),
        "vla_pulse": str(TEST_DATA_ROOT / "ds" / "gpm1839.vla.pulse.ds"),
        "askap_pulse": str(TEST_DATA_ROOT / "ds" / "j1755.askap.pulse.ds"),
    }
