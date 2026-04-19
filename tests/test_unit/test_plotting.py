import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from dstools.plotting.plotting import (
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
from tests.builders import (
    make_synthetic_dynamic_spectrum,
    make_synthetic_dynamic_spectrum_with_stokes_data,
    make_synthetic_lightcurve,
    make_synthetic_spectrum,
)

matplotlib.use("Agg")


def assert_axis_label_contains(ax, text, axis="x"):
    label = ax.get_xlabel() if axis == "x" else ax.get_ylabel()

    assert text in label


def assert_colorbar_label(fig, expected_label_fragment):
    assert len(fig.axes) >= 2
    assert expected_label_fragment in fig.axes[-1].get_ylabel()


def assert_dynamic_spectrum_annotation(ax, expected_text):
    labels = [text.get_text() for text in ax.texts]

    assert expected_text in labels


@pytest.fixture(autouse=True)
def cleanup_pyplot():
    yield

    plt.close("all")


def test_format_timeaxis_absolute_one_date():
    fake_ds = make_synthetic_dynamic_spectrum(tmax=2)

    _, ax = plt.subplots()

    x = np.arange(fake_ds.tmax)
    y = [np.random.random() for _ in x]
    ax.plot(x, y)

    format_timeaxis(fake_ds, ax)

    label = ax.get_xlabel()

    assert "UTC" in label


def test_format_timeaxis_absolute_two_dates():
    fake_ds = make_synthetic_dynamic_spectrum(tmax=6)

    _, ax = plt.subplots()

    x = np.arange(fake_ds.tmax)
    y = [np.random.random() for _ in x]
    ax.plot(x, y)

    format_timeaxis(fake_ds, ax)

    assert ax.get_xlabel() == ""


@pytest.mark.parametrize("stokes", ["I", "L"])
def test_plot_ds_without_fig_ax(stokes):
    fake_ds = make_synthetic_dynamic_spectrum()

    fig, ax = plot_ds(fake_ds, stokes=stokes)

    assert len(ax.images) == 1
    assert_axis_label_contains(ax, "Frequency", axis="y")
    assert_colorbar_label(fig, "Flux Density")
    expected = "Stokes I" if stokes == "I" else r"L = $\sqrt{Q^2 + U^2}$"
    assert_dynamic_spectrum_annotation(ax, expected)


@pytest.mark.parametrize("stokes", ["I", "L"])
def test_plot_ds_with_existing_fig(stokes):
    fake_ds = make_synthetic_dynamic_spectrum()

    fig, ax = plt.subplots()
    fig_out, ax_out = plot_ds(fake_ds, stokes=stokes, fig=fig, ax=ax)

    assert fig == fig_out
    assert ax == ax_out
    assert len(ax.images) == 1
    assert_colorbar_label(fig, "Flux Density")


def test_plot_ds_with_all_nan():
    fake_ds = make_synthetic_dynamic_spectrum()
    fake_ds.data["I"] = np.full(fake_ds.data["I"].shape, fill_value=np.nan)

    fig, ax = plot_ds(fake_ds, stokes="I")

    assert len(ax.images) == 1
    assert np.isnan(np.asarray(ax.images[0].get_array())).all()
    assert_axis_label_contains(ax, "Frequency", axis="y")
    assert_colorbar_label(fig, "Flux Density")


def test_plot_ds_with_folding_no_absolute_times():
    fake_ds = make_synthetic_dynamic_spectrum(fold=True, absolute_times=False)

    _, ax = plot_ds(fake_ds, stokes="I")

    assert "Phase (deg)" in ax.get_xlabel()


def test_plot_ds_with_folding_and_absolute_times():
    fake_ds = make_synthetic_dynamic_spectrum(fold=True, absolute_times=True)

    _, ax = plot_ds(fake_ds, stokes="I")

    assert "Phase (deg)" in ax.get_xlabel()


@pytest.mark.parametrize("column", ["time", "frequency"])
def test_plot_timefreqseries(column):
    tf = (
        make_synthetic_lightcurve(absolute_times=True, fold=False, shape=(3, 10))
        if column == "time"
        else make_synthetic_spectrum(absolute_times=True, fold=False)
    )

    _, ax = _plot_timefreqseries(tf, stokes="I")

    assert len(ax.lines) == 1
    assert_axis_label_contains(ax, "Flux Density", axis="y")
    expected = "UTC" if column == "time" else "Frequency"
    assert_axis_label_contains(ax, expected)


def test_plot_timefreqseries_multiple_stokes():
    tf = make_synthetic_lightcurve(absolute_times=True, fold=False)
    _, ax = _plot_timefreqseries(tf, stokes="IQ")

    assert len(ax.get_lines()) == 2
    assert {text.get_text() for text in ax.get_legend().get_texts()} == {"I", "Q"}


def test_plot_timefreqseries_spectrum():
    tf = make_synthetic_spectrum(absolute_times=True, fold=False)
    _, ax = _plot_timefreqseries(tf, stokes="IQ")

    assert_axis_label_contains(ax, "Frequency")


def test_plot_timefreqseries_relative_times_no_folding():
    tf = make_synthetic_lightcurve(absolute_times=False, fold=False)
    _, ax = _plot_timefreqseries(tf, stokes="I")

    assert_axis_label_contains(ax, "Time")


def test_plot_polarisations():
    lc = make_synthetic_lightcurve(absolute_times=True, fold=False)

    _, (data_ax, pa_ax, ell_ax, pol_ax) = _plot_polarisations(
        lc,
        stokes="IQ",
        error_alpha=0.4,
    )

    assert len(data_ax.get_lines()) == 2
    assert_axis_label_contains(pa_ax, "P.A.", axis="y")
    assert_axis_label_contains(ell_ax, "Ellipticity", axis="y")
    assert_axis_label_contains(pol_ax, "Polarisation", axis="y")


def test_plot_lightcurve():
    lc = make_synthetic_lightcurve(absolute_times=True, fold=False, shape=(3, 10))

    _, ax = plot_lightcurve(lc, stokes="IQ")

    assert_axis_label_contains(ax, "UTC")
    assert_axis_label_contains(ax, "Flux Density", axis="y")


def test_plot_spectrum():
    sp = make_synthetic_spectrum(absolute_times=True, fold=False)

    _, ax = plot_spectrum(sp, stokes="IQ")

    assert_axis_label_contains(ax, "Frequency")
    assert_axis_label_contains(ax, "Flux Density", axis="y")


def test_plot_polarisation_lightcurve():
    lc = make_synthetic_lightcurve(absolute_times=True, fold=False, shape=(3, 10))

    _, (data_ax, _, _, _) = plot_polarisation_lightcurve(lc, stokes="IQ")

    assert_axis_label_contains(data_ax, "UTC")
    assert_axis_label_contains(data_ax, "Flux Density", axis="y")


def test_plot_polarisation_spectrum():
    sp = make_synthetic_spectrum(absolute_times=True, fold=False)

    _, (data_ax, _, _, _) = plot_polarisation_spectrum(sp, stokes="IQ")

    assert_axis_label_contains(data_ax, "Frequency")
    assert_axis_label_contains(data_ax, "Flux Density", axis="y")


def test_plot_acf():
    ds = make_synthetic_dynamic_spectrum_with_stokes_data(shape=(20, 8))
    periodic = np.sin(2 * np.pi * np.arange(20) / 5).reshape(-1, 1)
    ds.data["I"] = np.tile(periodic, (1, 8)) + 0j

    acf_fig, acf_ax, _, acfz_ax = plot_acf(ds, stokes="I")

    assert len(acf_ax.images) == 1
    assert_axis_label_contains(acf_ax, "Time Lag")
    assert_axis_label_contains(acf_ax, "Frequency Lag", axis="y")
    assert_colorbar_label(acf_fig, "ACF")
    assert len(acfz_ax.lines) == 2
    assert acfz_ax.lines[1].get_linestyle() == "--"
    assert_axis_label_contains(acfz_ax, "Time Lag")
    assert_axis_label_contains(acfz_ax, "ACF", axis="y")
    assert hasattr(ds, "peak_lags")
    assert len(ds.peak_lags) > 0


def test_plot_summary(ds_paths):
    ds = make_synthetic_dynamic_spectrum()

    _, axes = plot_summary(ds, stokes="IQUVL")

    assert len(axes) == 6
