import os
from pathlib import Path

import numpy as np
import pytest

import dstools
from dstools.ms import MeasurementSet

package_root = Path(dstools.__path__[0]).parent


@pytest.fixture
def dispersed_pulse():
    return np.load(
        f"{package_root}/tests/data/ds/dispersed_pulse_dm3000.npy", allow_pickle=True
    )


@pytest.fixture
def ms_path():
    return package_root / "tests/data/msets/fred.atca.ms"


@pytest.fixture
def im_paths():
    images = {
        "mask": package_root / "tests/data/images/mask.fits",
        "image": package_root / "tests/data/images/test-MFS-I-image.fits",
        "model": package_root / "tests/data/images/test-MFS-I-model.fits",
        "residual": package_root / "tests/data/images/test-MFS-I-residual.fits",
    }

    return images


@pytest.fixture
def temp_environment(tmp_path_factory, mocker, ms_path):
    """Set up temporary filesystem environment to test I/O operations with CASA tasks mocked."""

    # Set up temporary directory for test
    tmp_path = tmp_path_factory.mktemp("temp")

    # Copy MSets over
    nspw_ms_path = ms_path.with_suffix(".2spw.ms")
    caltable_path = ms_path.with_suffix(".cal")
    rotated_ms_path = ms_path.with_suffix(".dstools-temp.rotated.ms")
    averaged_ms_path = ms_path.with_suffix(".dstools-temp.baseavg.ms")
    subbed_ms_path = ms_path.with_suffix(".subtracted.ms")

    model_path = package_root / "tests/data/images"

    tmp_ms_path = tmp_path / "test.ms"
    tmp_caltable_path = tmp_path / "test.cal"
    tmp_selfcal_ms_path = tmp_path / "test.selfcal1.ms"
    tmp_rotated_ms_path = tmp_path / "test.dstools-temp.rotated.ms"
    tmp_averaged_ms_path = tmp_path / "test.dstools-temp.baseavg.ms"
    tmp_combined_ms_path = tmp_path / "test.dstools-temp.comb.ms"
    tmp_onespw_ms_path = tmp_path / "test.1spw.ms"
    tmp_twospw_ms_path = tmp_path / "test.2spw.ms"
    tmp_subbed_ms_path = tmp_path / "test.subtracted.ms"
    tmp_model_path = tmp_path / "model"

    os.system(f"cp -r {ms_path} {tmp_ms_path}")
    os.system(f"cp -r {ms_path} {tmp_selfcal_ms_path}")
    os.system(f"cp -r {caltable_path} {tmp_caltable_path}")
    os.system(f"cp -r {ms_path} {tmp_onespw_ms_path}")
    os.system(f"cp -r {nspw_ms_path} {tmp_twospw_ms_path}")
    os.system(f"cp -r {rotated_ms_path} {tmp_rotated_ms_path}")
    os.system(f"cp -r {averaged_ms_path} {tmp_averaged_ms_path}")
    os.system(f"cp -r {ms_path} {tmp_combined_ms_path}")
    os.system(f"cp -r {subbed_ms_path} {tmp_subbed_ms_path}")
    os.system(f"cp -r {model_path} {tmp_model_path}")

    # Mock CASA tasks and file-system operations, as we will directly compare
    # the MS state to the temporary MS files
    flagstats = {
        "antenna": {
            "1": {"flagged": 240.0, "total": 1620.0},
            "6": {"flagged": 244.0, "total": 1620.0},
        },
    }
    mocker.patch("dstools.casa.casatasks.mstransform")
    mocker.patch("dstools.casa.casatasks.cvel")
    mocker.patch("dstools.casa.casatasks.phaseshift")
    mocker.patch("dstools.casa.casatasks.uvsub")
    mocker.patch("dstools.casa.casatasks.flagdata", return_value=flagstats)
    mocker.patch("dstools.casa.casatasks.split")
    mocker.patch("dstools.casa.casatasks.gaincal")
    mocker.patch("dstools.casa.casatasks.applycal")
    mocker.patch("dstools.imaging.parse_stdout_stderr")
    mocker.patch("casatools.table.table.copy")
    mocker.patch("subprocess.Popen")
    mocker.patch("os.system")
    mocker.patch("os.chdir")

    # Pass temporary paths to test
    tmp_ms_paths = {
        "onespw": tmp_ms_path,
        "twospw": tmp_twospw_ms_path,
        "rotated": tmp_rotated_ms_path,
        "averaged": tmp_averaged_ms_path,
        "cal": tmp_caltable_path,
        "model": tmp_model_path,
    }
    yield tmp_ms_paths

    # Clean up temp path
    os.system(f"rm -r {tmp_ms_path}")

    return


@pytest.fixture
def ms(temp_environment):
    return MeasurementSet(temp_environment["onespw"])


@pytest.fixture
def ds_paths():
    paths = {
        "atca_pulse": f"{package_root}/tests/data/ds/fred.atca.pulse.ds",
        "atca_calscan": f"{package_root}/tests/data/ds/fred.atca.calscan.ds",
        "vla_pulse": f"{package_root}/tests/data/ds/gpm1839.vla.pulse.ds",
        "askap_pulse": f"{package_root}/tests/data/ds/j1755.askap.pulse.ds",
    }

    return paths
