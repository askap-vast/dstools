import numpy as np
import pytest

from dstools.dynamic_spectrum import DynamicSpectrum
from dstools.polarisation import (
    RMTimeSeries,
    derotate_dynamic_spectrum,
    dynamic_rmts_from_fdf,
    peak_rm_from_fdf,
    rm_error,
    rm_synthesis,
    smooth_rm_timeseries,
)


def test_rm_synthesis(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    fdf = rm_synthesis(
        ds.data["L"],
        ds.freq,
    )

    assert fdf.fdf.shape == (len(ds.time), 40001)


def test_peak_rm_from_fdf(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    fdf = rm_synthesis(
        ds.data["L"],
        ds.freq,
    )

    rm = peak_rm_from_fdf(fdf, ds.data["I"])

    assert round(rm, 1) == -830.5


def test_peak_rm_from_fdf_handles_calibrator_scan(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    # Test that RM extraction runs on the pulse at later time index
    fdf = rm_synthesis(ds.data["L"], ds.freq)

    # Simulate early cal scan
    ds.data["L"][5, :] = np.nan

    rm = peak_rm_from_fdf(fdf, ds.data["I"])

    assert round(rm, 1) == -830.5


def test_dynamic_rmts_from_fdf(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j
    Li = ds.data["Q"].imag + ds.data["U"].imag * 1j

    fdf = rm_synthesis(
        L,
        ds.freq,
    )
    fdf_im = rm_synthesis(
        Li,
        ds.freq,
    )

    rmts = dynamic_rmts_from_fdf(fdf, fdf_im, ds.data["I"])

    assert isinstance(rmts, RMTimeSeries)
    assert len(rmts.data) == len(ds.time)


def test_dynamic_rmts_from_fdf_all_masked(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j
    Li = ds.data["Q"].imag + ds.data["U"].imag * 1j

    fdf = rm_synthesis(
        L,
        ds.freq,
    )
    fdf_im = rm_synthesis(
        Li,
        ds.freq,
    )

    rmts = dynamic_rmts_from_fdf(fdf, fdf_im, ds.data["I"], snr_min=1000)
    assert np.all(np.isnan(rmts.model))
    assert np.all(rmts.mask)

    rmts = smooth_rm_timeseries(rmts)

    assert np.all(np.isnan(rmts.model))
    assert np.all(rmts.mask)


@pytest.mark.parametrize(
    "mode",
    ["constant", "polynomial", "spline", "periodic"],
)
def test_smooth_rm_timeseries(ds_paths, mode):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j
    Li = ds.data["Q"].imag + ds.data["U"].imag * 1j

    fdf = rm_synthesis(
        L,
        ds.freq,
    )
    fdf_im = rm_synthesis(
        Li,
        ds.freq,
    )

    rmts = dynamic_rmts_from_fdf(fdf, fdf_im, ds.data["I"])
    rmts = smooth_rm_timeseries(rmts, mode=mode)

    assert isinstance(rmts, RMTimeSeries)
    assert len(rmts.data) == len(ds.time)
    assert len(rmts.error) == len(ds.time)
    assert len(rmts.model) == len(ds.time)


def test_smooth_rm_timeseries_unknown_mode_raises_error(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j
    Li = ds.data["Q"].imag + ds.data["U"].imag * 1j

    fdf = rm_synthesis(
        L,
        ds.freq,
    )
    fdf_im = rm_synthesis(
        Li,
        ds.freq,
    )

    rmts = dynamic_rmts_from_fdf(fdf, fdf_im, ds.data["I"])

    with pytest.raises(ValueError):
        smooth_rm_timeseries(rmts, mode="poly")


def test_rm_error(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j
    Li = ds.data["Q"].imag + ds.data["U"].imag * 1j

    fdf = rm_synthesis(
        L,
        ds.freq,
    )
    fdf_im = rm_synthesis(
        Li,
        ds.freq,
    )

    rm_err = rm_error(fdf, fdf_im)

    assert len(rm_err) == len(ds.time)


def test_derotate_l_dynamic_spectrum(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    L = ds.data["Q"].real + ds.data["U"].real * 1j

    derot_L = derotate_dynamic_spectrum(
        L,
        ds.freq,
        rm=-800,
    )

    assert L.shape == derot_L.shape
