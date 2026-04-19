import astropy.units as u
import numpy as np
import pandas as pd
import pytest

import dstools
from dstools.dynamic_spectrum import DynamicSpectrum, LightCurve, Spectrum
from tests.builders import (
    make_pulse_array,
    make_stokes_data,
    make_synthetic_dynamic_spectrum,
    make_synthetic_dynamic_spectrum_with_stokes_data,
)


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
    flagged_channels = np.all(np.isnan(ds.data["I"].real), axis=0)

    assert flagged_channels.sum() == 50


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
    flagged_times = np.all(np.isnan(ds.data["I"].real), axis=1)

    assert flagged_times.sum() == 6


def test_barycentric_correction(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, barycentre=True)

    assert ds.header["time_start"] == "2022-09-22 06:14:33.495"
    assert ds.header["time_scale"] == "tdb"


def test_derotate_faraday(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)

    rm = ds.derotate_faraday()

    assert round(rm, 1) == -830.5


def test_ds_acf(ds_paths):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path)
    acf2d = ds.acf(stokes="I")

    assert acf2d.shape == (301, 47)


@pytest.mark.parametrize("period", [1700, 3224])
def test_ds_fold_conserves_flux(period):
    tres = 10
    phase_bins = round(period / tres)

    ds = make_synthetic_dynamic_spectrum(
        tunit=u.s,
        period=period,
        fold_periods=1,
        phase_bins=phase_bins,
    )

    array, time = make_pulse_array(period, tres)
    ds.time = time

    folded = ds._fold(array)

    percentage_error = abs(folded.max() - array.max()) / array.max()
    assert percentage_error < 0.01


def test_dedispersion_recovers_channel_averaged_peak(dispersed_pulse):
    # Create synthetic DS to generate required properties
    ds = make_synthetic_dynamic_spectrum(
        tunit=u.s,
        DM=6000,
        time=np.linspace(0, 100, 1000),
        freq=np.linspace(500, 1000, 100),
    )

    # Desdisperse synthetic DM=6000 dispersed pulse
    dedispersed = ds._dedisperse(dispersed_pulse)

    assert dispersed_pulse.shape == dedispersed.shape
    assert round(np.max(dispersed_pulse.real.mean(axis=1)), 2) == 0.62
    assert round(np.max(dedispersed.real.mean(axis=1)), 1) == 5.2


def test_rebin_updates_time_frequency_axes_and_data_shape():
    ds = make_synthetic_dynamic_spectrum(
        time=np.linspace(0, 7, 8),
        freq=np.linspace(900, 1000, 6),
        tavg=2,
        favg=3,
    )

    shape = (8, 6)
    XX = np.ones(shape, dtype=np.complex128)
    XY = np.ones(shape, dtype=np.complex128)
    YX = np.ones(shape, dtype=np.complex128)
    YY = np.ones(shape, dtype=np.complex128)

    XX, XY, YX, YY = ds._rebin(XX, XY, YX, YY)

    assert XX.shape == (4, 2)
    assert XY.shape == (4, 2)
    assert YX.shape == (4, 2)
    assert YY.shape == (4, 2)
    assert ds.time.shape == (4,)
    assert ds.freq.shape == (2,)


def test_flag_channels_masks_selected_frequency_range():
    freq = np.linspace(900, 1000, 11)
    ds = make_synthetic_dynamic_spectrum(
        freq=freq,
        flag_channels=[(940, 970)],
    )

    shape = (4, freq.size)
    XX = np.ones(shape, dtype=np.complex128)
    XY = np.ones(shape, dtype=np.complex128)
    YX = np.ones(shape, dtype=np.complex128)
    YY = np.ones(shape, dtype=np.complex128)

    XX, XY, YX, YY = ds._flag_channels(XX, XY, YX, YY)

    masked = (freq > 940) & (freq <= 970)
    assert np.all(np.isnan(XX[:, masked]))
    assert np.all(np.isnan(XY[:, masked]))
    assert np.all(np.isnan(YX[:, masked]))
    assert np.all(np.isnan(YY[:, masked]))
    assert np.all(np.isfinite(XX[:, ~masked]))


@pytest.mark.parametrize(
    ("flag_times", "expected_start", "expected_end"),
    [
        ([(1, 2)], 2, 3),
        ([("00:01:00", "00:02:00")], 1, 2),
    ],
)
def test_flag_times_masks_selected_range(flag_times, expected_start, expected_end):
    time = np.arange(0, 6)
    ds = make_synthetic_dynamic_spectrum(
        time=time,
        tunit=u.min,
        flag_times=flag_times,
        header={
            "time_start": "2020-01-01 00:00:00",
            "time_scale": "utc",
            "phasecentre": "00:00:00 +00:00:00",
            "telescope": "ATCA",
            "feeds": "linear",
            "channels": 4,
        },
    )

    shape = (time.size, 4)
    XX = np.ones(shape, dtype=np.complex128)
    XY = np.ones(shape, dtype=np.complex128)
    YX = np.ones(shape, dtype=np.complex128)
    YY = np.ones(shape, dtype=np.complex128)

    XX, XY, YX, YY = ds._flag_times(XX, XY, YX, YY)

    assert np.all(np.isnan(XX[expected_start:expected_end, :]))
    assert np.all(np.isnan(XY[expected_start:expected_end, :]))
    assert np.all(np.isnan(YX[expected_start:expected_end, :]))
    assert np.all(np.isnan(YY[expected_start:expected_end, :]))
    assert np.all(np.isfinite(XX[:expected_start, :]))
    assert np.all(np.isfinite(XX[expected_end:, :]))


def test_acf_shape_matches_input_quadrant():
    data = make_stokes_data(shape=(6, 4))
    ds = make_synthetic_dynamic_spectrum(data=data)

    acf2d = ds.acf(stokes="I")

    assert acf2d.shape == (4, 6)
    assert np.nanmax(acf2d) == 1


def test_lightcurve_construction_from_synthetic_ds():
    ds = make_synthetic_dynamic_spectrum()
    ds.data = make_stokes_data(shape=(5, 10))

    lc = LightCurve(ds)

    assert lc.column == "time"
    assert lc.x.shape == (5,)
    assert lc.flux["I"].shape == (5,)


def test_lightcurve_folded_construction_from_synthetic_ds():
    ds = make_synthetic_dynamic_spectrum(fold=True, absolute_times=False)
    ds.data = make_stokes_data(shape=(5, 10))

    lc = LightCurve(ds)

    assert lc.column == "time"
    assert lc.x.shape == (5,)


def test_lightcurve_imag_construction_from_synthetic_ds():
    ds = make_synthetic_dynamic_spectrum()
    ds.data = make_stokes_data(shape=(5, 10))
    ds.data["I"] = ds.data["I"] + 1j * np.ones((5, 10))

    lc = LightCurve(ds, imag=True)

    assert lc.column == "time"
    assert lc.flux["I"].shape == (5,)


def test_spectrum_construction_from_synthetic_ds():
    ds = make_synthetic_dynamic_spectrum()
    ds.data = make_stokes_data(shape=(5, 10))

    sp = Spectrum(ds)

    assert sp.column == "frequency"
    assert sp.x.shape == (10,)
    assert sp.flux["I"].shape == (10,)


def test_save_spectrum_writes_csv(tmp_path):
    ds = make_synthetic_dynamic_spectrum_with_stokes_data()

    sp = Spectrum(ds)
    savepath = tmp_path / "test_spec.csv"

    sp.save(savepath, include_pols=True)

    saved = pd.read_csv(savepath)

    assert list(saved.columns[:5]) == [
        "frequency",
        "flux_density_I",
        "flux_density_I_err",
        "flux_density_Q",
        "flux_density_Q_err",
    ]
    assert len(saved) == sp.x.size
    assert np.allclose(saved["flux_density_I"], 10.0)
    assert "polarisation_angle" in saved
    assert saved["polarisation_angle"].notna().all()


def test_save_lightcurve_writes_csv(tmp_path):
    ds = make_synthetic_dynamic_spectrum_with_stokes_data()

    lc = LightCurve(ds)
    savepath = tmp_path / "test_lc.csv"

    lc.save(savepath, include_pols=True)

    saved = pd.read_csv(savepath)

    assert list(saved.columns[:5]) == [
        "time",
        "flux_density_I",
        "flux_density_I_err",
        "flux_density_Q",
        "flux_density_Q_err",
    ]
    assert len(saved) == lc.x.size
    assert np.allclose(saved["flux_density_I"], 10.0)
    assert "linear_fraction" in saved
    assert saved["linear_fraction"].notna().all()
