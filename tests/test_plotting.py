from types import SimpleNamespace

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest

from dstools.dynamic_spectrum import DynamicSpectrum
from dstools.plotting import (
    _plot_polarisations,
    _plot_timefreqseries,
    format_timeaxis,
    plot_acf,
    plot_ds,
    plot_lightcurve,
    plot_polarisation_lightcurve,
    plot_polarisation_spectrum,
    plot_spectrum,
    plot_summary,
)


@pytest.fixture(autouse=True)
def cleanup_pyplot():
    yield

    plt.close("all")


def create_fake_ds(
    absolute_times: bool = True,
    fold: bool = False,
    tmax: int = 6,
):
    # Simulate 6 hours of 10s samples and 5 frequency channels
    tunit = u.hour
    fmin, fmax = 100, 105
    fold_periods = 1

    header = {
        "time_start": "2020-01-01 20:00:00",
        "time_scale": "utc",
    }

    rng = np.random.default_rng()
    shape = 5, 10

    data = {
        "I": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "Q": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "U": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "V": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "L": np.abs(rng.normal(0, 1, shape)),
    }

    label = "Time (s)" if not fold else "Phase (deg)"

    ds = SimpleNamespace(
        data=data,
        tunit=tunit,
        fmin=fmin,
        fmax=fmax,
        tmin=0,
        tmax=tmax,
        fold=fold,
        fold_periods=fold_periods,
        absolute_times=absolute_times,
        header=header,
        _timelabel=label,
    )

    return ds


def create_fake_tfseries(
    absolute_times: bool,
    column: str,
    fold: bool,
):
    ds = create_fake_ds(absolute_times=absolute_times, fold=fold)

    tf = SimpleNamespace(
        x=np.array([1, 2, 3]),
        flux={
            "I": np.array([1.0, 1.0, 3.0]),
            "Q": np.array([1.0, 2.0, 3.0]),
        },
        flux_err={
            "I": np.array([0.1, 0.1, 0.1]),
            "Q": np.array([0.1, 0.4, 0.1]),
        },
        column=column,
        polangle=np.array([20, 50, 30]),
        polangle_err=np.array([20, 50, 30]),
        ellipticity=np.array([20, 50, 30]),
        ellipticity_err=np.array([20, 50, 30]),
        circular_fraction=np.array([0.4, 0.5, 0.3]),
        circular_fraction_err=np.array([0.1, 0.1, 0.1]),
        linear_fraction=np.array([0.5, 0.3, 0.4]),
        linear_fraction_err=np.array([0.1, 0.1, 0.1]),
        ds=ds,
    )

    return tf


def test_format_timeaxis_absolute_one_date():
    fake_ds = create_fake_ds(tmax=2)

    _, ax = plt.subplots()

    x = np.arange(fake_ds.tmax)
    y = [np.random.random() for _ in x]
    ax.plot(x, y)

    format_timeaxis(fake_ds, ax)

    label = ax.get_xlabel()

    # Single timestamp is set as xlabel
    assert "UTC" in label


def test_format_timeaxis_absolute_two_dates():
    fake_ds = create_fake_ds(tmax=6)

    _, ax = plt.subplots()

    x = np.arange(fake_ds.tmax)
    y = [np.random.random() for _ in x]
    ax.plot(x, y)

    format_timeaxis(fake_ds, ax)

    label = ax.get_xlabel()

    # Timestamps are added as text and no label is necessary
    assert "" in label


@pytest.mark.parametrize("stokes", ["I", "L"])
def test_plot_ds_without_fig_ax(stokes):
    fake_ds = create_fake_ds()

    fig, ax = plot_ds(fake_ds, stokes=stokes)

    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)


@pytest.mark.parametrize("stokes", ["I", "L"])
def test_plot_ds_with_existing_fig(stokes):
    fake_ds = create_fake_ds()

    fig, ax = plt.subplots()
    fig_out, ax_out = plot_ds(fake_ds, stokes=stokes, fig=fig, ax=ax)

    assert fig == fig_out
    assert ax == ax_out


def test_plot_ds_with_all_nan():
    fake_ds = create_fake_ds()
    fake_ds.data["I"] = np.full(fake_ds.data["I"].shape, fill_value=np.nan)

    fig, ax = plot_ds(fake_ds, stokes="I")

    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)


def test_plot_ds_with_folding_no_absolute_times():
    fake_ds = create_fake_ds(fold=True, absolute_times=False)

    _, ax = plot_ds(fake_ds, stokes="I")

    assert "Phase (deg)" in ax.get_xlabel()


def test_plot_ds_with_folding_and_absolute_times():
    fake_ds = create_fake_ds(fold=True, absolute_times=True)

    _, ax = plot_ds(fake_ds, stokes="I")

    assert "Phase (deg)" in ax.get_xlabel()


@pytest.mark.parametrize("column", ["time", "frequency"])
def test_plot_timefreqseries(column):
    tf = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column=column,
    )

    fig, ax = _plot_timefreqseries(tf, stokes="I")

    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)


def test_plot_timefreqseries_multiple_stokes():
    tf = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="time",
    )
    _, ax = _plot_timefreqseries(tf, stokes="IQ")

    lines = ax.get_lines()

    assert len(lines) == 2


def test_plot_timefreqseries_spectrum():
    tf = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="frequency",
    )
    _, ax = _plot_timefreqseries(tf, stokes="IQ")

    assert ax.get_xlabel() == "Frequency (MHz)"


def test_plot_polarisations():
    lc = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="time",
    )

    _, (data_ax, pa_ax, ell_ax, pol_ax) = _plot_polarisations(
        lc,
        stokes="IQ",
        error_alpha=0.4,
    )

    assert len(data_ax.get_lines()) == 2
    assert "P.A. (deg)" == pa_ax.get_ylabel()
    assert "Ellipticity (deg)" == ell_ax.get_ylabel()
    assert "Fractional Polarisation" == pol_ax.get_ylabel()


def test_plot_lightcurve():
    lc = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="time",
    )

    _, ax = plot_lightcurve(lc, stokes="IQ")

    assert "UTC" in ax.get_xlabel()


def test_plot_spectrum():
    sp = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="frequency",
    )

    _, ax = plot_spectrum(sp, stokes="IQ")

    assert "Frequency (MHz)" in ax.get_xlabel()


def test_plot_polarisation_lightcurve():
    lc = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="time",
    )

    _, (data_ax, _, _, _) = plot_polarisation_lightcurve(lc, stokes="IQ")

    assert "UTC" in data_ax.get_xlabel()


def test_plot_polarisation_spectrum():
    sp = create_fake_tfseries(
        absolute_times=True,
        fold=False,
        column="frequency",
    )

    _, (data_ax, _, _, _) = plot_polarisation_spectrum(sp, stokes="IQ")

    assert "Frequency (MHz)" in data_ax.get_xlabel()


def test_plot_summary_integration_default_cmax(ds_paths):
    ds = DynamicSpectrum(ds_paths["atca_pulse"])

    _, axes = plot_summary(ds, stokes="IQUVL")

    assert len(axes) == 6


def test_plot_summary_integration(ds_paths):
    ds = DynamicSpectrum(ds_paths["atca_pulse"])
    cmax = {stokes: 30 for stokes in "IQUV"}

    _, axes = plot_summary(ds, stokes="IQUVL", cmax=cmax)

    assert len(axes) == 6


def test_plot_fdf(ds_paths):
    pass


def test_plot_acf_integration(ds_paths):
    ds = DynamicSpectrum(ds_paths["atca_pulse"])

    acf_fig, acf_ax, acfz_fig, acfz_ax = plot_acf(ds, stokes="I")

    assert isinstance(acf_fig, plt.Figure)
    assert isinstance(acfz_fig, plt.Figure)
    assert isinstance(acf_ax, plt.Axes)
    assert isinstance(acfz_ax, plt.Axes)
