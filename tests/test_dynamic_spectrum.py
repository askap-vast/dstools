import os

import astropy.units as u
import numpy as np
import pytest

import dstools
from dstools.dynamic_spectrum import DynamicSpectrum, LightCurve, Spectrum


def test_testdata_updated(ds_paths):
    for _, ds_path in ds_paths.items():
        ds = DynamicSpectrum(ds_path)

        assert f"dstools_version: {dstools.__version__}" in str(ds)


def test_ds_trim(ds_paths):
    ds_path = ds_paths.get("atca_calscan")
    ds_notrim = DynamicSpectrum(ds_path, trim=False)
    ds_trim = DynamicSpectrum(ds_path, trim=True)

    assert ds_notrim.header["channels"] == 101
    assert ds_trim.header["channels"] == 59


def test_ds_minuv_when_already_baseline_averaged(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    DynamicSpectrum(ds_path, minuvdist=500)


def test_ds_fold_without_period_raises_error(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    with pytest.raises(ValueError):
        DynamicSpectrum(ds_path, fold=True)


def make_pulse_array(period, tres):
    # Initialise array
    tsamples = 1000
    shape = tsamples, 50
    array = np.zeros(shape, dtype=np.complex128)

    # Make a pulse profile
    sigma = 90
    pulse_halfwidth = sigma // tres
    x = np.tile(np.arange(-pulse_halfwidth, pulse_halfwidth).reshape(-1, 1), 50)
    pulse = 100 * np.exp(-(x**2) / (2 * pulse_halfwidth**2)) + 1j * 0 * x

    # Inject pulses
    for pulse_index in range(10, tsamples - 10, period // tres):
        array[pulse_index - pulse_halfwidth : pulse_index + pulse_halfwidth, :] += pulse

    time = np.arange(0, tres * tsamples, tres)

    return array, time


@pytest.mark.parametrize("period", [1700, 3224])
def test_ds_fold_conserves_flux(period, ds_paths):
    ds_path = ds_paths.get("atca_pulse")

    tres = 10
    phase_bins = round(period / tres)

    ds = DynamicSpectrum(
        ds_path,
        period=period,
        tunit=u.s,
        fold_periods=1,
        phase_bins=phase_bins,
    )

    array, time = make_pulse_array(period, tres)
    ds.time = time

    folded = ds._fold(array)

    # Check that folded pulse is within 1% of original
    percentage_error = abs(folded.max() - array.max()) / array.max()
    assert percentage_error < 0.01


def test_ds_crop(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(
        ds_path,
        mintime=1,
        maxtime=3,
        minfreq=2800,
        maxfreq=2900,
        tunit=u.min,
    )

    assert ds.data["I"].shape == (14, 101)


def test_flag_channels(ds_paths):
    ds_path = ds_paths.get("askap_pulse")

    flag_chans = [(950, 1000)]
    ds = DynamicSpectrum(ds_path, flag_channels=flag_chans)
    I_lc = np.nanmean(ds.data["I"].real, axis=0)

    assert np.isnan(I_lc).sum() == 50


@pytest.mark.parametrize(
    "flag_times",
    [
        [(1, 2)],
        [("00:59:00", "01:00:00")],
    ],
)
def test_flag_times(ds_paths, flag_times):
    ds_path = ds_paths.get("askap_pulse")

    ds = DynamicSpectrum(ds_path, flag_times=flag_times, tunit=u.min)
    I_lc = np.nanmean(ds.data["I"].real, axis=1)

    assert np.isnan(I_lc).sum() == 6


def test_dedispersion(ds_paths, dispersed_pulse):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, DM=6000, dedisperse=True, tunit=u.s)

    # Insert fake dispersed data
    tbins, fbins = (1000, 100)
    ds.time = np.linspace(0, 100, tbins)
    ds.freq = np.linspace(500, 1000, fbins)

    # Build noise + pulse array
    dedispersed = ds._dedisperse(dispersed_pulse)

    # De-dispersed array should have same dimensions as input
    # and recover a boosted flux in the channel-averaged data
    assert dispersed_pulse.shape == dedispersed.shape
    assert round(np.max(dispersed_pulse.real.mean(axis=1)), 2) == 0.62
    assert round(np.max(dedispersed.real.mean(axis=1)), 1) == 5.2


def test_barycentric_correction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, barycentre=True)

    assert ds.header["time_start"] == "2022-09-22 06:14:33.495"
    assert ds.header["time_scale"] == "tdb"


def test_rm_synthesis(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, derotate=True)

    assert round(ds.RM, 1) == -830.5


def test_rm_synthesis_handles_calibrator_scan(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    # Simulate early cal scan
    ds.data["I"][5, :] = np.nan
    ds.data["L"][5, :] = np.nan

    # Test that RM extraction runs on the pulse at later time index
    I = ds.data["I"]
    L = ds.data["L"].T
    RM = ds.rm_synthesis(I, L)

    assert round(RM, 1) == -830.5


def test_ds_acf(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    acf2d = ds.acf(stokes="I")

    assert acf2d.shape == (301, 47)


def test_lc_construction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    lc = LightCurve(ds)

    assert lc.column == "time"


def test_lc_folded_construction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, fold=True, period=90, tunit=u.s)
    lc = LightCurve(ds)

    assert lc.column == "time"


def test_lc_imag_construction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    lc = LightCurve(ds, imag=True)

    assert lc.column == "time"


def test_sp_construction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    sp = Spectrum(ds)

    assert sp.column == "frequency"


def test_save_spectrum(ds_paths, tmp_path_factory):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    sp = Spectrum(ds)
    savepath = tmp_path_factory.mktemp("temp") / "test_spec.csv"

    sp.save(savepath, include_pols=True)

    assert os.path.exists(savepath)


def test_save_lightcurve(ds_paths, tmp_path_factory):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    sp = LightCurve(ds)
    savepath = tmp_path_factory.mktemp("temp") / "test_lc.csv"

    sp.save(savepath, include_pols=True)

    assert os.path.exists(savepath)
