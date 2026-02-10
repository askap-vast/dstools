import logging
from dataclasses import dataclass, field
from typing import Literal, Optional

import astropy.units as u
import numpy as np
from astropy import constants as c
from astropy.timeseries import LombScargle
from rm_lite.utils.synthesis import freq_to_lambda2, make_phi_arr, rmsynth_nufft
from scipy.interpolate import UnivariateSpline

logger = logging.getLogger(__name__)


@dataclass
class FDF2D:
    """Result of RM synthesis on a dynamic spectrum.

    phi      : [n_phi] Faraday depths (rad/m^2)
    fdf      : [n_time, n_phi] complex Faraday dispersion function
    lam_sq_0 : reference lambda^2 used in RM-synthesis (m^2)
    meta     : arbitrary extras (RMSF, weights, etc.)
    """

    phi: np.ndarray
    fdf: np.ndarray
    lam_sq: np.ndarray
    lam_sq_0: float
    meta: dict = field(default_factory=dict)


@dataclass
class RMTimeSeries:
    """Time-resolved RM measurement and model."""

    rm_data: np.ndarray
    rm_error: Optional[np.ndarray]
    rm_model: np.ndarray
    mask: np.ndarray
    snr_min: float = 5.0

    @property
    def data(self):
        a = self.rm_data.copy()
        a[self.mask] = np.nan

        return a

    @property
    def error(self):
        a = self.rm_error.copy()
        a[self.mask] = np.nan

        return a

    @property
    def model(self):
        a = self.rm_model.copy()
        a[self.mask] = np.nan

        return a


def rm_synthesis(
    L: np.ndarray,
    freq: np.ndarray,
    n_phi: int = 2000,
    dphi: float = 0.1,
    lam_sq_0: float | None = None,
) -> FDF2D:
    """Run RM-synthesis on a 2D IQU dynamic spectrum and return the FDF."""

    # Prepare data for RM synthesis
    phi_arr = make_phi_arr(n_phi, dphi)
    freq_hz = freq * 1e6

    lam_sq = freq_to_lambda2(freq_hz)
    if lam_sq_0 is None:
        lam_sq_0_m2 = np.mean(lam_sq)

    weights_arr = np.ones_like(freq_hz)

    # Perform RM synthesis on the 2D complex polarisation dynamic spectrum
    fdf_timeseries = rmsynth_nufft(
        complex_pol_arr=L.T,
        lambda_sq_arr_m2=freq_to_lambda2(freq_hz),
        phi_arr_radm2=phi_arr,
        weight_arr=weights_arr,
        lam_sq_0_m2=lam_sq_0_m2,
    )

    fdf = FDF2D(
        phi=phi_arr,
        fdf=fdf_timeseries.T,
        lam_sq=lam_sq,
        lam_sq_0=lam_sq_0,
        meta={},
    )

    return fdf


def peak_rm_from_fdf(
    fdf: FDF2D,
    ds_array: np.ndarray,
) -> float:
    """Measure a single RM from a 2D FDF at lightcurve peak."""

    # Get time of where ds_array (typically Stokes I or L) peaks
    lc = np.nanmean(ds_array.real, axis=1)
    tslice = np.nanargmax(lc)

    # Get max Faraday depth at lightcurve peak
    fdf_lc_peak = np.abs(fdf.fdf[tslice])
    phi_index = np.argmax(fdf_lc_peak)

    return fdf.phi[phi_index]


def fit_rm_peak(rm_raw, n_time):
    # rm_raw_fit = rm_raw.copy()
    # for t in range(n_time):
    #     j = int(phi_idx[t])
    #     if j == 0 or j == n_phi - 1:
    #         continue

    #     x = phi[j - 1 : j + 2]
    #     y = abs_fdf[t, j - 1 : j + 2]
    #     # Quadratic fit y = a x^2 + b x + c
    #     coeff = np.polyfit(x, y, deg=2)
    #     a = coeff[0]
    #     b = coeff[1]
    #     if a == 0:
    #         continue
    #     x_peak = -b / (2 * a)
    #     rm_raw_fit[t] = x_peak

    # rm_raw = rm_raw_fit
    # method = "fdf_peak_parabola"

    return


def fwhm_rmsf_from_lam_sq(lam_sq: np.ndarray) -> float:
    """Approximate FWHM of the RMSF from lambda^2 coverage (Brentjens & de Bruyn 2005)."""

    delta_lam_sq = np.nanmax(lam_sq) - np.nanmin(lam_sq)

    return 2 * np.sqrt(3) / delta_lam_sq


def rm_error(fdf: FDF2D, fdf_im: FDF2D) -> np.ndarray:
    # Get peak RM from real FDF
    abs_F_sig = np.abs(fdf.fdf)
    peak_amp = np.nanmax(abs_F_sig, axis=1)

    # Get noise from imaginary FDF
    sigma_re = fdf_im.fdf.real.std(axis=1, ddof=1)
    sigma_amp = np.sqrt(2) * sigma_re
    snr = peak_amp / sigma_amp

    # Compute RM error from RMSF and SNR
    fwhm_rmsf = fwhm_rmsf_from_lam_sq(fdf.lam_sq)
    rm_err = fwhm_rmsf / (2 * snr)

    return rm_err


def dynamic_rmts_from_fdf(
    fdf: FDF2D,
    fdf_im: FDF2D,
    mask_array: np.ndarray,
    snr_min: float = 5.0,
    use_parabolic_fit: bool = True,
) -> RMTimeSeries:
    """Estimate raw per-time-bin RM from FDF(t, phi)."""

    abs_fdf = np.abs(fdf.fdf)
    phi_idx = np.argmax(abs_fdf, axis=1)
    rm_raw = fdf.phi[phi_idx]

    # if use_parabolic_fit:
    #     rm_raw = fit_rm_peak(fdf)

    sqrtn = np.sqrt(np.isfinite(mask_array).sum(axis=1))
    noise = np.nanstd(mask_array, axis=1) / sqrtn
    snr_mask = np.nanmean(mask_array, axis=1) < snr_min * noise

    rm_model = rm_raw.copy()
    rm_err = rm_error(fdf, fdf_im)

    return RMTimeSeries(
        rm_data=rm_raw,
        rm_error=rm_err,
        rm_model=rm_model,
        mask=snr_mask,
    )


def _sin_model(t, A, omega, phi, C):
    return A * np.sin(omega * t + phi) + C


def fit_periodic_rm(
    t_fit: np.ndarray,
    rm_fit: np.ndarray,
    min_period: float | None = None,
    max_period: float | None = None,
) -> tuple[float, float, float, float]:
    """Fit sinusoidal RM using Lomb-Scargle for frequency
    estimation and linear least-squares for amplitude/phase/offset.
    """

    # Detrend only for the period search; final fit uses original rm_fit
    rm0 = rm_fit - rm_fit.mean()

    # Convert period bounds to frequency bounds if given
    min_freq = 1.0 / max_period if max_period is not None else None
    max_freq = 1.0 / min_period if min_period is not None else None

    ls = LombScargle(t_fit, rm0, nterms=1)
    freq, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
    )

    best_freq = freq[np.argmax(power)]  # cycles / time-unit
    omega = 2.0 * np.pi * best_freq  # rad / time-unit

    # Linear LS at fixed ω: RM(t) = a*sin(ωt) + b*cos(ωt) + C
    s = np.sin(omega * t_fit)
    c = np.cos(omega * t_fit)
    X = np.column_stack((s, c, np.ones_like(t_fit)))

    a, b, C = np.linalg.lstsq(X, rm_fit, rcond=None)[0]

    # Map (a, b) -> (A, φ) so we can write A*sin(ω t + φ) + C
    A = np.hypot(a, b)
    phi = np.arctan2(b, a)  # because a = A*cos(phi), b = A*sin(phi)

    return A, omega, phi, C


def smooth_rm_timeseries(
    rm_ts: RMTimeSeries,
    mode: Literal["constant", "polynomial", "spline", "periodic"] = "constant",
    poly_deg: int = 1,
) -> RMTimeSeries:
    """Replace rm_ts.rm_model with a smoothed / modeled RM(t)."""

    rm = rm_ts.rm_data
    mask = rm_ts.mask

    n_time = rm.shape[0]
    t = np.arange(n_time, dtype=float)

    if np.all(mask):
        rm_ts.rm_model = np.full_like(rm, np.nan)
        return rm_ts

    t_fit = t[~mask]
    rm_fit = rm[~mask]

    match mode:
        case "constant":
            val = np.nanmedian(rm_fit)
            rm_model = np.full_like(rm, val)

        case "polynomial":
            coeff = np.polyfit(t_fit, rm_fit, deg=poly_deg)
            rm_model = np.polyval(coeff, t)

        case "spline":
            spline = UnivariateSpline(t_fit, rm_fit)
            rm_model = spline(t)

        case "periodic":
            params = fit_periodic_rm(t_fit, rm_fit)
            rm_model = _sin_model(t, *params)
        case _:
            raise ValueError(f"Unknown mode: {mode}")

    rm_model[mask] = np.nan

    # Return a new RMTimeseries to avoid mutating original
    rmts = RMTimeSeries(
        rm_data=rm_ts.data,
        rm_error=rm_ts.error,
        rm_model=rm_model,
        mask=rm_ts.mask,
    )

    return rmts


def derotate_dynamic_spectrum(
    L: np.ndarray,
    freq_mhz: np.ndarray,
    rm: float | np.ndarray,
) -> np.ndarray:
    """Derotate a complex L(t, nu) dynamic spectrum for Faraday rotation.

    Parameters
    ----------
    L
        Complex linear polarisation dynamic spectrum, shape [n_time, n_chan].
    freq_mhz
        Channel centre frequencies in MHz, shape [n_chan].
    rm
        Rotation measure in rad/m^2. Can be:
        - scalar (constant RM for all times),
        - [n_time] array (RM(t)),
        - or something already broadcastable to [n_time, n_chan].

    Returns
    -------
    np.ndarray
        Derotated complex L' with same shape as L.
    """

    rm_arr = np.asarray(rm, dtype=np.float32)

    # Replace NaN with 0 for derotation
    # to keep polarised / low RM time bins
    rm_arr[np.isnan(rm_arr)] = 0

    lam = (c.c / (freq_mhz * u.MHz)).to(u.m).value
    lam_sq = lam * lam

    # Broadcast phase to full DS shape
    phase = -2 * rm_arr[..., None] * lam_sq[None, :]

    return L * np.exp(1j * phase)
