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


@pytest.mark.parametrize(
    "period, time_res",
    [
        (60, "39.149 s"),
        (90, "25.532 s"),
        (100, "22.979 s"),
    ],
)
def test_ds_fold(ds_paths, period, time_res):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, fold=True, period=period, tunit=u.s)

    assert ds.header["time_resolution"] == time_res


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


def test_dedispersion(ds_paths, dispersed_pulse):
    ds_path = ds_paths.get("atca_pulse")
    ds = DynamicSpectrum(ds_path, DM=3000, dedisperse=True, tunit=u.s)

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
