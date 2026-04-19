import matplotlib
import numpy as np
import pytest
from astropy.io import fits

pytestmark = pytest.mark.fullstack

imaging = pytest.importorskip("dstools.imaging")
ms_module = pytest.importorskip("dstools.ms")
WSClean = imaging.WSClean


matplotlib.use("Agg")


def test_wsclean_fits_mask_without_target(mocker, ms_path, imaging_workspace):
    """Check that the masked pixels in this path agrees with expectations.

    We have initialised wsclean with fits_mask containing non-zero pixels on
    both the field source to be cleaned and the target, and target_mask
    with non-zero pixels only covering the target. We expect that after running
    _get_fits_mask, the final fits_mask should only include the field source pixels.
    """

    mocker.patch("subprocess.Popen")
    mocker.patch("dstools.imaging.parse_stdout_stderr")

    clean_mask_path = imaging_workspace["clean_mask"]
    target_mask_path = imaging_workspace["target_mask"]

    wsclean = WSClean(
        imsize=500,
        cellsize="0.66asec",
        fits_mask=clean_mask_path,
        target_mask=target_mask_path,
    )
    wsclean._get_fits_mask()

    final_mask_path = imaging_workspace["final_mask"]
    with fits.open(clean_mask_path) as hdul:
        clean_mask = hdul[0].data
    with fits.open(final_mask_path) as hdul:
        fits_mask = hdul[0].data

    assert np.allclose(clean_mask, fits_mask)
