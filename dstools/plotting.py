import logging
from typing import Optional

import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.visualization import ImageNormalize, ZScaleInterval
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks

from dstools.dynamic_spectrum import (
    DynamicSpectrum,
    LightCurve,
    Spectrum,
    TimeFreqSeries,
)

logger = logging.getLogger(__name__)

COLORS = {
    "I": "firebrick",
    "Q": "lightgreen",
    "U": "darkorchid",
    "V": "darkorange",
    "L": "cornflowerblue",
    "P": "",
    "PA": "",
}

DS_LABELS = {
    "I": "Stokes I",
    "Q": "Stokes Q",
    "U": "Stokes U",
    "V": "Stokes V",
    "L": r"L = $\sqrt{Q^2 + U^2}$",
    "P": r"P = $\sqrt{Q^2 + U^2 + V^2}$",
    "PA": r"P.A. = $\atan{U/Q}$",
}


def format_timeaxis(ds: DynamicSpectrum, ax):
    """Apply custom timestamp formatter to time axis."""

    timescale = ds.header["time_scale"]
    timestart = ds.header["time_start"]
    t0 = Time(timestart, format="iso", scale=timescale)

    # Set tick locations
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    ax.xaxis.set_major_locator(locator)

    # Determine first occurrences of unique dates
    tick_dates = [Time(t0 + t * ds.tunit).strftime("%Y-%m-%d") for t in ax.get_xticks()]
    unique_dates, unique_date_idx = list(np.unique(tick_dates, return_index=True))
    single_date = len(unique_dates) == 1

    # Set x-axis label to UTC date if data confined to a single day
    # otherwise we add the date as text overlays at the first xtick occurence in a day
    xlabel = f"{timescale.upper()} {unique_dates[0]}" if single_date else None
    ax.set_xlabel(xlabel)

    # Add custom timestamp formatter
    def timestamp_formatter(x, pos):
        t = t0 + x * ds.tunit
        mpl_date = mdates.num2date(t.plot_date)

        if not single_date and pos in unique_date_idx:
            fmt = f"%H:%M\n{timescale.upper()} %Y-%m-%d"
        else:
            fmt = "%H:%M"

        return mpl_date.strftime(fmt)

    ax.xaxis.set_major_formatter(timestamp_formatter)

    return


def plot_ds(
    ds: DynamicSpectrum,
    stokes,
    cmax=20,
    imag=False,
    fig=None,
    ax=None,
):
    """Plot dynamic spectrum for single Stokes parameter."""

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Select polarisation
    if stokes == "L":
        data = np.abs(ds.data[stokes])
    else:
        data = ds.data[stokes].imag if imag else ds.data[stokes].real

    # Produce normalisation for products with valid data
    if not np.isnan(data).all():
        norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
    else:
        norm = None

    # Configure plot parameters
    cmap = "plasma" if stokes in ["I", "L"] else "coolwarm"
    cmin = -2 if stokes in ["I", "L"] else -cmax

    # Set time axis to phase if folding
    phasemax = 0.5 * ds.fold_periods
    tmin, tmax = (-phasemax, phasemax) if ds.fold else (ds.tmin, ds.tmax)

    im = ax.imshow(
        data.T,
        extent=[tmin, tmax, ds.fmin, ds.fmax],
        aspect="auto",
        origin="lower",
        norm=norm,
        clim=(cmin, cmax),
        cmap=cmap,
    )

    ax.set_xlabel(ds._timelabel)
    ax.set_ylabel("Frequency (MHz)")

    if ds.absolute_times and not ds.fold:
        format_timeaxis(ds, ax)

    ax.text(
        0.05,
        0.95,
        DS_LABELS[stokes],
        color="white",
        weight="heavy",
        path_effects=[pe.withStroke(linewidth=2, foreground="black")],
        transform=ax.transAxes,
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Flux Density (mJy)")

    fig.tight_layout()

    return fig, ax


def plot_lightcurve(
    lc: LightCurve,
    stokes: str,
    fig: Optional = None,
    ax: Optional = None,
):
    return _plot_timefreqseries(lc, stokes, fig, ax)


def plot_spectrum(
    sp: Spectrum,
    stokes: str,
    fig: Optional = None,
    ax: Optional = None,
):
    return _plot_timefreqseries(sp, stokes, fig, ax)


def _plot_timefreqseries(
    tf: TimeFreqSeries,
    stokes: str,
    fig: Optional = None,
    ax: Optional = None,
):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_xlabel(tf.ds._timelabel)

    # Overplot each specified polarisation
    for s in stokes:
        ax.errorbar(
            tf.x,
            y=tf.flux[s],
            yerr=tf.flux_err[s],
            lw=1,
            color=COLORS[s],
            marker="o",
            markersize=1,
            label=s,
        )

    ax.set_ylabel("Flux Density (mJy)")
    ax.legend()

    pad = (tf.x.max() - tf.x.min()) * 0.05
    ax.set_xlim([tf.x.min() - pad, tf.x.max() + pad])

    if tf.column == "time" and tf.ds.absolute_times:
        format_timeaxis(tf.ds, ax)
    else:
        ax.set_xlabel("Frequency (MHz)")

    fig.tight_layout()

    return fig, ax


def plot_summary(
    ds: DynamicSpectrum,
    stokes: str,
    cmax: Optional[dict[float]] = None,
    imag: bool = False,
):
    """Plot all Stokes dynamic spectra and averaged light-curve / spectrum."""

    if cmax is None:
        cmax = {stokes: 30 for stokes in "IQUV"}

    fig = plt.figure(figsize=(14, 15))
    gs = GridSpec(3, 2, figure=fig)

    I_ax = fig.add_subplot(gs[0, 0])
    Q_ax = fig.add_subplot(gs[0, 1])
    U_ax = fig.add_subplot(gs[1, 0])
    V_ax = fig.add_subplot(gs[1, 1])
    lc_ax = fig.add_subplot(gs[2, 0])
    sp_ax = fig.add_subplot(gs[2, 1])

    fig, I_ax = plot_ds(ds, stokes="I", cmax=cmax["I"], fig=fig, ax=I_ax, imag=imag)
    fig, Q_ax = plot_ds(ds, stokes="Q", cmax=cmax["Q"], fig=fig, ax=Q_ax, imag=imag)
    fig, U_ax = plot_ds(ds, stokes="U", cmax=cmax["U"], fig=fig, ax=U_ax, imag=imag)
    fig, V_ax = plot_ds(ds, stokes="V", cmax=cmax["V"], fig=fig, ax=V_ax, imag=imag)

    lc = LightCurve(ds, imag=imag)
    sp = Spectrum(ds, imag=imag)

    fig, sp_ax = plot_spectrum(
        sp,
        stokes=stokes,
        fig=fig,
        ax=sp_ax,
    )
    fig, lc_ax = plot_lightcurve(
        lc,
        stokes=stokes,
        fig=fig,
        ax=lc_ax,
    )

    fig.subplots_adjust(
        left=0.06,
        top=0.98,
        right=0.96,
        bottom=0.05,
        hspace=0.18,
        wspace=0.24,
    )

    axes = [I_ax, Q_ax, U_ax, V_ax, lc_ax, sp_ax]

    return fig, axes

def plot_fdf(ds: DynamicSpectrum, fig=None, ax=None):
    """Plot Faraday dispersion function."""

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
def plot_polarisation_lightcurve(lc: LightCurve, stokes: str, error_alpha: float = 0.4):
    return _plot_polarisations(lc, stokes=stokes, error_alpha=error_alpha)


def plot_polarisation_spectrum(lc: LightCurve, stokes: str, error_alpha: float = 0.4):
    return _plot_polarisations(lc, stokes=stokes, error_alpha=error_alpha)


    if not ds.polobs:
        I = ds.data["I"]
        Q = ds.data["Q"]
        U = ds.data["U"]
        _ = ds.rm_synthesis(I, Q, U)
def _plot_polarisations(tf: TimeFreqSeries, stokes: str, error_alpha: float):
    """Plot lightcurve / spectrum with polarisation parameters."""

    fig = plt.figure(figsize=(10, 10))
    gs = GridSpec(5, 1, figure=fig)

    data_ax = fig.add_subplot(gs[3:, 0])
    pa_ax = fig.add_subplot(gs[2, 0])
    ell_ax = fig.add_subplot(gs[1, 0])
    pol_ax = fig.add_subplot(gs[0, 0])

    ax.plot(
        ds.polobs.rmsf_phi,
        np.abs(ds.polobs.rmsf),
    # Plot lightcurve / spectrum
    fig, data_ax = _plot_timefreqseries(
        tf,
        stokes=stokes,
        fig=fig,
        ax=data_ax,
    )

    # Plot polarisation angle
    pa_ax.errorbar(
        tf.x,
        y=tf.polangle,
        yerr=tf.polangle_err,
        alpha=error_alpha,
        color="k",
        label="RMSF",
        marker="o",
        markersize=1,
        ls="none",
    )
    ax.plot(
        ds.polobs.phi,
        np.abs(ds.polobs.fdf),
        color="r",
        label="FDF",
    pa_ax.axhline(
        0,
        ls=":",
        color="k",
        alpha=0.5,
    )
    ax.plot(
        ds.polobs.phi,
        np.abs(ds.polobs.rm_cleaned),
        color="b",
        label="Clean",

    ell_ax.errorbar(
        tf.x,
        y=tf.ellipticity,
        yerr=tf.ellipticity_err,
        alpha=error_alpha,
        color="k",
        marker="o",
        markersize=1,
        ls="none",
    )
    ax.plot(
        ds.polobs.phi,
        np.abs(ds.polobs.rm_comps),
        color="g",
        label="Model",
    ell_ax.axhline(
        0,
        ls=":",
        color="k",
        alpha=0.5,
    )

    ax.set_xlabel(r"RM ($rad/m^2$)")
    ax.set_ylabel("Amplitude")
    pol_ax.errorbar(
        tf.x,
        y=tf.linear_fraction,
        yerr=tf.linear_fraction_err,
        color="dodgerblue",
        alpha=error_alpha,
        marker="o",
        markersize=1,
        label="$|L/I|$",
        ls="none",
    )
    pol_ax.errorbar(
        tf.x,
        y=tf.circular_fraction,
        yerr=tf.circular_fraction_err,
        color="darkorange",
        alpha=error_alpha,
        marker="o",
        markersize=1,
        label="$|V/I|$",
        ls="none",
    )

    ax.legend()
    # Apply formatting
    pa_ax.set_xticklabels([])
    ell_ax.set_xticklabels([])
    pol_ax.set_xticklabels([])

    return fig, ax
    pad = (tf.x.max() - tf.x.min()) * 0.05
    pa_ax.set_xlim([tf.x.min() - pad, tf.x.max() + pad])
    ell_ax.set_xlim([tf.x.min() - pad, tf.x.max() + pad])
    pol_ax.set_xlim([tf.x.min() - pad, tf.x.max() + pad])

    pa_ax.set_ylim(-100, 100)
    ell_ax.set_ylim(-100, 100)
    pol_ax.set_ylim(-0.1, 1.1)
    pa_ax.set_yticks([-90, 0, 90])
    ell_ax.set_yticks([-90, 0, 90])
    pol_ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    pa_ax.set_ylabel("P.A. (deg)")
    ell_ax.set_ylabel("Ellipticity (deg)")
    pol_ax.set_ylabel("Fractional Polarisation")

    pol_ax.legend()

    fig.tight_layout()

    axes = [data_ax, pa_ax, ell_ax, pol_ax]

    return fig, axes


def plot_acf(ds, stokes="I", contrast=0.4):
    """Plot 2D auto-correlation function of dynamic spectrum."""

    acf2d = ds.acf(stokes)

    # Plot 2D ACF
    acf_fig, acf_ax = plt.subplots(figsize=(7, 5))

    norm = ImageNormalize(acf2d, interval=ZScaleInterval(contrast=contrast))
    im = acf_ax.imshow(
        acf2d,
        extent=[0, ds.tmax - ds.tmin, 0, ds.fmax - ds.fmin],
        aspect="auto",
        norm=norm,
        cmap="plasma",
    )
    cb = acf_fig.colorbar(
        im,
        ax=acf_ax,
        fraction=0.05,
        pad=0.02,
    )
    cb.formatter.set_powerlimits((0, 0))
    cb.set_label("ACF")

    acf_ax.set_xlabel(f"Time Lag ({ds.tunit})")
    acf_ax.set_ylabel("Frequency Lag (MHz)")

    # Plot zero frequency lag trace
    acfz_fig, acfz_ax = plt.subplots(figsize=(7, 5))

    zero_trace_acf = acf2d[-1, 1:]
    time_lag = np.linspace(0, ds.tmax - ds.tmin, len(zero_trace_acf))
    acfz_ax.plot(
        time_lag,
        zero_trace_acf,
        color="k",
    )

    acfz_ax.set_xlabel(f"Time Lag ({ds.tunit})")
    acfz_ax.set_ylabel("ACF")

    acf_peaks, props = find_peaks(zero_trace_acf, prominence=(None, None))

    max_prom = np.argsort(props["prominences"])[::-1]
    ds.peak_lags = time_lag[acf_peaks[max_prom]]

    acfz_ax.axvline(
        ds.peak_lags[0],
        color="darkorange",
        ls="--",
    )

    peak_lag = ds.peak_lags[0] * ds.tunit
    logger.debug(f"Stokes {stokes} ACF peak at {peak_lag:.3f}")

    return acf_fig, acf_ax, acfz_fig, acfz_ax
