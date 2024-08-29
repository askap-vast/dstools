import logging
import warnings
from abc import ABC
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import astropy.constants as c
import astropy.units as u
import h5py
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.visualization import ImageNormalize, ZScaleInterval
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import correlate, find_peaks

from dstools.rm import PolObservation

logger = logging.getLogger(__name__)

COLORS = {
    "I": "firebrick",
    "Q": "lightgreen",
    "U": "darkorchid",
    "V": "darkorange",
    "L": "cornflowerblue",
}


def snr_mask(data, noise, n_sigma):
    """Mask data array below n_sigma based on imaginary component of stokes array."""

    mask = np.abs(noise.real) < n_sigma * np.nanstd(noise.imag)
    data[mask] = np.nan

    return data


def rebin(o, n, axis):
    """Create unitary array compression matrix from o -> n length.

    if rebinning along row axis we want:
        - (o // n) + 1 entries in each row that sum to unity,
        - each column to sum to the compression ratio o / n
        - values distributed along the row in units of o / n until expired

        >>> rebin(5, 3)
        array([[0.6, 0.4, 0. , 0. , 0. ],
               [0. , 0.2, 0.6, 0.2, 0. ],
               [0. , 0. , 0. , 0.4, 0.6]])

        - transpose of this for column rebinning

    The inner product of this compressor with an array will rebin
    the array conserving the total intensity along the given axis.
    """

    compressor = np.zeros((n, o))

    # Exit early with empty array if chunk is empty
    if compressor.size == 0:
        return compressor

    comp_ratio = n / o

    nrow = 0
    ncol = 0

    budget = 1
    overflow = 0

    # While loop to avoid visiting n^2 zero-value cells
    while nrow < n and ncol < o:
        # Use overflow if just spilled over from last row
        if overflow > 0:
            value = overflow
            overflow = 0
            budget -= value
            row_shift = 0
            col_shift = 1
        # Use remaining budget if at end of current row
        elif budget < comp_ratio:
            value = budget
            overflow = comp_ratio - budget
            budget = 1
            row_shift = 1
            col_shift = 0
        # Otherwise spend n / o and move to next column
        else:
            value = comp_ratio
            budget -= value
            row_shift = 0
            col_shift = 1

        compressor[nrow, ncol] = value
        nrow += row_shift
        ncol += col_shift

    return compressor if axis == 0 else compressor.T


def rebin2D(array, new_shape):
    """Re-bin along time / frequency axes conserving flux."""

    if new_shape == array.shape:
        array[array == 0 + 0j] = np.nan
        return array

    if new_shape[0] > array.shape[0] or new_shape[1] > array.shape[1]:
        raise ValueError(
            "New shape should not be greater than old shape in either dimension"
        )

    time_comp = rebin(array.shape[0], new_shape[0], axis=0)
    freq_comp = rebin(array.shape[1], new_shape[1], axis=1)
    array[np.isnan(array)] = 0 + 0j
    result = time_comp @ np.array(array) @ freq_comp
    result[result == 0 + 0j] = np.nan

    return result


def slice_array(a, ax1_min, ax1_max, ax2_min=None, ax2_max=None):
    """Slice 1D or 2D array with variable lower and upper boundaries."""

    if ax2_min is None and ax2_max is None:
        a = a[ax1_min:] if ax1_max == 0 else a[ax1_min:ax1_max]
    else:
        a = a[ax1_min:, :] if ax1_max == 0 else a[ax1_min:ax1_max, :]
        a = a[:, ax2_min:] if ax2_max == 0 else a[:, ax2_min:ax2_max]

    return a


def make_summary_plot(ds, stokes, cmax, imag):
    """Plot all-stokes dynamic spectra and averaged light-curve / spectrum."""

    fig = plt.figure(figsize=(14, 15))
    gs = GridSpec(3, 2, figure=fig)

    I_ax = fig.add_subplot(gs[0, 0])
    Q_ax = fig.add_subplot(gs[0, 1])
    U_ax = fig.add_subplot(gs[1, 0])
    V_ax = fig.add_subplot(gs[1, 1])
    lc_ax = fig.add_subplot(gs[2, 0])
    sp_ax = fig.add_subplot(gs[2, 1])

    fig, I_ax = ds.plot_ds(stokes="I", cmax=cmax["I"], fig=fig, ax=I_ax, imag=imag)
    fig, Q_ax = ds.plot_ds(stokes="Q", cmax=cmax["Q"], fig=fig, ax=Q_ax, imag=imag)
    fig, U_ax = ds.plot_ds(stokes="U", cmax=cmax["U"], fig=fig, ax=U_ax, imag=imag)
    fig, V_ax = ds.plot_ds(stokes="V", cmax=cmax["V"], fig=fig, ax=V_ax, imag=imag)

    fig, sp_ax = ds.plot_spectrum(
        stokes=stokes,
        fig=fig,
        ax=sp_ax,
        imag=imag,
    )
    fig, lc_ax = ds.plot_lightcurve(
        stokes=stokes,
        fig=fig,
        ax=lc_ax,
        polangle=False,
        imag=imag,
    )

    fig.subplots_adjust(
        left=0.06,
        top=0.94,
        right=0.96,
        bottom=0.05,
    )

    return fig


@dataclass
class DynamicSpectrum:
    ds_path: str

    favg: int = 1
    tavg: int = 1

    minfreq: Optional[float] = None
    maxfreq: Optional[float] = None
    mintime: Optional[float] = None
    maxtime: Optional[float] = None
    minuvdist: float = 0
    maxuvdist: float = np.inf
    minuvwave: float = 0
    maxuvwave: float = np.inf

    tunit: u.Quantity = u.hour
    corr_dumptime: float = 10.1

    derotate: bool = False

    fold: bool = False
    period: Optional[float] = None
    period_offset: float = 0.0
    fold_periods: int = 2

    calscans: bool = True
    trim: bool = True

    def __post_init__(self):

        # Load instrumental polarisation time/frequency/uvdist arrays
        XX, XY, YX, YY = self._load_data()

        # Insert calibrator scan breaks
        XX, XY, YX, YY = self._stack_cal_scans(XX, XY, YX, YY)

        # Store time and frequency resolution
        timebins = len(self.time) / self.tavg
        freqbins = len(self.freq) / self.favg
        self.time_res = (self.tmax - self.tmin) * self.tunit / timebins
        self.freq_res = (self.fmax - self.fmin) * u.MHz / freqbins
        self.header.update(
            {
                "time_resolution": f"{self.time_res.to(u.s):.1f}",
                "freq_resolution": f"{self.freq_res.to(u.MHz):.1f}",
            }
        )

        # Fold data to selected period
        if self.fold:
            if not self.period:
                raise ValueError("Must pass period argument when folding.")

            XX = self._fold(XX)
            XY = self._fold(XY)
            YX = self._fold(YX)
            YY = self._fold(YY)

            self.time = rebin(len(self.time), len(XX), axis=0) @ self.time

        # Average data in time and frequency
        XX, XY, YX, YY = self._rebin(XX, XY, YX, YY)

        # Compute Stokes products and store in data attribute
        self._make_stokes(XX, XY, YX, YY)

    def __str__(self):
        str_rep = ""
        for attr in self.header:
            str_rep += f"{attr}: {self.header[attr]}\n"
        return str_rep

    def _fold(self, data):
        """Average chunks of data folding at specified period."""

        # Calculate number of pixels in each chunk
        pixel_duration = self.time[1] - self.time[0]
        chunk_length = min(int(self.period // pixel_duration), len(data))

        # Create left-padded nans, derived from period phase offset
        offset = (0.5 + self.period_offset) * self.period
        leftpad_length = int(offset // pixel_duration)
        leftpad_chunk = np.full((leftpad_length, data.shape[1]), np.nan)

        # Create right-padded nans
        rightpad_length = chunk_length - (leftpad_length + len(data)) % chunk_length
        rightpad_chunk = np.full((rightpad_length, data.shape[1]), np.nan)

        # Stack and split data
        data = np.vstack((leftpad_chunk, data, rightpad_chunk))
        numsplits = int(data.shape[0] // chunk_length)
        arrays = np.split(data, numsplits)

        # Compute average along stack axis
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            data = np.nanmean(arrays, axis=0)

        return np.tile(data, (self.fold_periods, 1))

    def _get_scan_intervals(self):
        """Find indices of start/end of each calibrator scan cycle."""

        dts = [0]
        dts.extend([self.time[i] - self.time[i - 1] for i in range(1, len(self.time))])
        dts = np.array(dts)

        # Locate indices signaling beginning of cal-scan
        # (scan intervals longer than correlator dump time)
        scan_start_idx = np.where(np.abs(dts) > self.corr_dumptime)[0]

        # End indices are just prior to the next start index, then
        scan_end_idx = scan_start_idx - 1

        # Insert first scan start index and last scan end index
        scan_start_idx = np.insert(scan_start_idx, 0, 0)
        scan_end_idx = np.append(scan_end_idx, len(self.time) - 1)

        return scan_start_idx, scan_end_idx

    def _validate(self, datafile):

        # Check if baselines have been pre-averaged and disable uvdist selection if so
        default_uv_params = [
            self.minuvdist == 0,
            self.maxuvdist == np.inf,
            self.minuvwave == 0,
            self.maxuvwave == np.inf,
        ]
        made_uvdist_selection = not all(default_uv_params)
        baseline_averaged = len(datafile["uvdist"][:]) == 1

        if made_uvdist_selection and baseline_averaged:
            logger.warning(
                f"DS is already baseline averaged, disabling uvdist selection."
            )
            self.minuvdist = 0
            self.maxuvdist = np.inf
            self.minuvwave = 0
            self.maxuvwave = np.inf

        return

    def _load_data(self):
        """Load instrumental pols and uvdist/time/freq data, converting to MHz, s, and mJy."""

        # Import instrumental polarisations and time/frequency/uvdist arrays
        with h5py.File(self.ds_path, "r") as f:

            self._validate(f)

            # Read header
            self.header = dict(f.attrs)

            # Read uvdist, time, frequency, and flux arrays
            uvdist = f["uvdist"][:]
            time = f["time"][:]
            freq = f["frequency"][:] / 1e6
            flux = f["flux"][:] * 1e3

            # Make baseline selection using UV distance
            blmask = (uvdist >= self.minuvdist) & (uvdist <= self.maxuvdist)
            uvdist = uvdist[blmask]
            flux = flux[blmask, :, :, :]

            # Construct array of UV distance in units of wavelength
            wavelength = (freq * u.MHz).to(u.m, equivalencies=u.spectral()).value
            uvdist_expanded = uvdist[:, np.newaxis, np.newaxis, np.newaxis]
            wavelength_expanded = wavelength[np.newaxis, np.newaxis, :, np.newaxis]
            uvwave = np.tile(
                uvdist_expanded / wavelength_expanded,
                (1, len(time), 1, 4),
            )

            uvwave_mask = (uvwave <= self.minuvwave) | (uvwave >= self.maxuvwave)

            # Apply uvwave limit mask
            flux[uvwave_mask] = np.nan
            uvwave[uvwave_mask] = np.nan

            # Average over baseline axis
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                flux = np.nanmean(flux, axis=0)

            # Read out instrumental polarisations
            XX = flux[:, :, 0]
            XY = flux[:, :, 1]
            YX = flux[:, :, 2]
            YY = flux[:, :, 3]

        # Set timescale
        time_scale_factor = self.tunit.to(u.s)
        time /= time_scale_factor
        self.corr_dumptime /= time_scale_factor
        self._timelabel = "Phase" if self.fold else f"Time ({self.tunit})"

        # Flip ATCA L-band frequency axis to intuitive order
        if freq[0] > freq[-1]:
            XX = np.flip(XX, axis=1)
            XY = np.flip(XY, axis=1)
            YX = np.flip(YX, axis=1)
            YY = np.flip(YY, axis=1)

            freq = np.flip(freq)

        # Optionally remove flagged channels at top/bottom of band
        if self.trim:
            # Create binary mask identifying non-nan values across all polarisations
            full = np.nansum((XX + XY + YX + YY), axis=0)
            full[full == 0.0 + 0.0j] = np.nan
            allpols = np.isfinite(full)

            # Set minimum and maximum non-nan channel indices
            minchan = np.argmax(allpols)
            if np.isnan(allpols[-1]):
                maxchan = -np.argmax(allpols[::-1]) + 1
            else:
                maxchan = 0
        else:
            minchan = 0
            maxchan = 0

        # Select channel range
        if self.minfreq:
            minchan = -np.argmax((freq < self.minfreq)[::-1]) - 1
        if self.maxfreq:
            maxchan = np.argmax(freq > self.maxfreq)

        # Select time range
        if self.mintime:
            mintime = np.argmax(time - time[0] > self.mintime)
        else:
            mintime = 0

        if self.maxtime:
            maxtime = -np.argmax((time - time[0] < self.maxtime)[::-1]) + 1
        else:
            maxtime = 0

        # Identify start time and set observation start to t=0
        time_start = Time(
            time[0] * time_scale_factor / 3600 / 24,
            format="mjd",
            scale="utc",
        )
        time_start.format = "iso"
        self.header["time_start"] = time_start
        time -= time[0]

        # Make data selection
        XX = slice_array(XX, mintime, maxtime, minchan, maxchan)
        XY = slice_array(XY, mintime, maxtime, minchan, maxchan)
        YX = slice_array(YX, mintime, maxtime, minchan, maxchan)
        YY = slice_array(YY, mintime, maxtime, minchan, maxchan)

        self.uvdist = uvdist
        self.freq = slice_array(freq, minchan, maxchan)
        self.time = slice_array(time, mintime, maxtime)

        self.tmin = self.time[0]
        self.tmax = self.time[-1]
        self.fmin = self.freq[0]
        self.fmax = self.freq[-1]

        self.header.update(
            {
                "integrations": len(self.time),
                "channels": len(self.freq),
            }
        )

        return XX, XY, YX, YY

    def _stack_cal_scans(self, XX, XY, YX, YY):
        """Insert null data representing off-source time."""

        scan_start_idx, scan_end_idx = self._get_scan_intervals()

        # Calculate number of cycles in each calibrator/stow break
        time_end_break = self.time[scan_start_idx[1:]]
        time_start_break = self.time[scan_end_idx[:-1]]

        dt = self.time[1] - self.time[0]
        num_break_cycles = np.append((time_end_break - time_start_break), 0) / dt
        num_channels = self.header["channels"]

        # Create initial time-slice to start stacking target and calibrator scans together
        new_data_XX = new_data_XY = new_data_YX = new_data_YY = np.zeros(
            (1, num_channels),
            dtype=complex,
        )
        new_time = np.zeros(1)

        for start_index, end_index, num_scans in zip(
            scan_start_idx, scan_end_idx, num_break_cycles
        ):
            # Select each contiguous on-target chunk of data
            XX_chunk = XX[start_index : end_index + 1, :]
            XY_chunk = XY[start_index : end_index + 1, :]
            YX_chunk = YX[start_index : end_index + 1, :]
            YY_chunk = YY[start_index : end_index + 1, :]
            time_chunk = self.time[start_index : end_index + 1]

            # Make array of complex NaN's for subsequent calibrator / stow gaps
            # and append to each on-target chunk of data
            if self.calscans and num_scans > 0:
                num_nans = (int(round(num_scans) - 1), num_channels)
                nan_chunk = np.full(num_nans, np.nan + np.nan * 1j)

                XX_chunk = np.ma.vstack([XX_chunk, nan_chunk])
                XY_chunk = np.ma.vstack([XY_chunk, nan_chunk])
                YX_chunk = np.ma.vstack([YX_chunk, nan_chunk])
                YY_chunk = np.ma.vstack([YY_chunk, nan_chunk])

                time_break_start = self.time[end_index] + dt
                time_break_scans = time_break_start + np.arange(num_nans[0]) * dt
                time_chunk = np.append(time_chunk, time_break_scans)

            new_data_XX = np.ma.vstack([new_data_XX, XX_chunk])
            new_data_XY = np.ma.vstack([new_data_XY, XY_chunk])
            new_data_YX = np.ma.vstack([new_data_YX, YX_chunk])
            new_data_YY = np.ma.vstack([new_data_YY, YY_chunk])
            new_time = np.append(new_time, time_chunk)

        new_data_XX = new_data_XX[1:]
        new_data_XY = new_data_XY[1:]
        new_data_YX = new_data_YX[1:]
        new_data_YY = new_data_YY[1:]

        self.time = new_time[1:]

        return new_data_XX, new_data_XY, new_data_YX, new_data_YY

    def _rebin(self, XX, XY, YX, YY):
        num_tsamples, num_channels = XX.shape
        tbins = num_tsamples // self.tavg
        fbins = num_channels // self.favg

        XX = rebin2D(XX, (tbins, fbins))
        XY = rebin2D(XY, (tbins, fbins))
        YX = rebin2D(YX, (tbins, fbins))
        YY = rebin2D(YY, (tbins, fbins))

        self.time = rebin(num_tsamples, tbins, axis=0) @ self.time
        self.freq = self.freq @ rebin(num_channels, fbins, axis=1)

        return XX, XY, YX, YY

    def _make_stokes(self, XX, XY, YX, YY):
        """Convert instrumental polarisations to Stokes products."""

        # Compute Stokes products from instrumental pols
        feedtype = self.header["feeds"]
        if feedtype == "linear":
            I = (XX + YY) / 2
            Q = (XX - YY) / 2
            U = (XY + YX) / 2
            V = 1j * (YX - XY) / 2
        elif feedtype == "circular":
            I = (XX + YY) / 2
            Q = (XY + YX) / 2
            U = 1j * (XY - YX) / 2
            V = (XX - YY) / 2
        else:
            raise ValueError(
                f"Feed type {feedtype} not recognised, should be either 'linear' or 'circular'."
            )

        L = Q.real + 1j * U.real

        if self.derotate:

            # Compute RM
            RM = self.rm_synthesis(I, Q, U)

            # Build L from imaginary components
            Li = Q.imag + 1j * U.imag

            # Derotate real and imaginary L
            L = self.derotate_faraday(L, RM)
            Li = self.derotate_faraday(Li, RM)

            # Compute complex Q and U from L
            Q = L.real + 1j * Li.real
            U = L.imag + 1j * Li.imag
        else:
            self.polobs = None

        P = np.sqrt(Q.real**2 + U.real**2 + V.real**2) / I.real

        PA = 0.5 * np.arctan2(U.real, Q.real) * u.rad.to(u.deg)

        self.data = {
            "XX": XX,
            "XY": XY,
            "YX": YX,
            "YY": YY,
            "I": I,
            "Q": Q,
            "U": U,
            "V": V,
            "L": L,
            "P": P,
            "PA": PA,
        }

        return

    def acf(self, stokes):
        """Generate a 2D auto-correlation of the dynamic spectrum."""

        # Replace NaN with zeros to calculate auto-correlation
        data = self.data[stokes].real.copy()
        data[np.isnan(data)] = 0.0

        # Compute auto-correlation and select upper-right quadrant
        acf2d = correlate(data, data)
        acf2d = acf2d[acf2d.shape[0] // 2 :, acf2d.shape[1] // 2 :]

        # Reorder time-frequency axes and normalise
        acf2d = np.flip(acf2d, axis=1).T
        acf2d /= np.nanmax(acf2d)

        return acf2d

    def rm_synthesis(self, I, Q, U):

        # Zero out null values in Stokes arrays
        I[np.isnan(I)] = 0 + 0j
        Q[np.isnan(Q)] = 0 + 0j
        U[np.isnan(U)] = 0 + 0j

        # Compute RM along brightest time-sample
        tslice = np.argmax(np.nanmean(I.real, axis=1))

        # RM synthesis
        stokes_arrays = (I[tslice, :].real, Q[tslice, :].real, U[tslice, :].real)
        self.polobs = PolObservation(
            self.freq * 1e6,
            stokes_arrays,
            verbose=False,
        )

        phi_axis = np.arange(-2000.0, 2000.1, 0.1)
        self.polobs.rmsynthesis(
            phi_axis,
            verbose=False,
        )
        self.polobs.rmclean(
            cutoff=1.0,
            verbose=False,
        )

        RM = self.polobs.phi[np.argmax(abs(self.polobs.fdf))]
        logger.debug(f"Peak RM of {RM:.1f} rad/m2")

        return RM

    def plot_fdf(self, fig=None, ax=None):
        """Plot Faraday dispersion function derived with RMclean (Heald 2009)."""

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 5))

        if not self.polobs:
            I = self.data["I"]
            Q = self.data["Q"]
            U = self.data["U"]
            _ = self.rm_synthesis(I, Q, U)

        ax.plot(
            self.polobs.rmsf_phi,
            np.abs(self.polobs.rmsf),
            color="k",
            label="RMSF",
        )
        ax.plot(
            self.polobs.phi,
            np.abs(self.polobs.fdf),
            color="r",
            label="FDF",
        )
        ax.plot(
            self.polobs.phi,
            np.abs(self.polobs.rm_cleaned),
            color="b",
            label="Clean",
        )
        ax.plot(
            self.polobs.phi,
            np.abs(self.polobs.rm_comps),
            color="g",
            label="Model",
        )

        ax.set_xlabel(r"RM ($rad/m^2$)")
        ax.set_ylabel("Amplitude")

        ax.legend()

        return fig, ax

    def derotate_faraday(self, L, RM):
        """Correct linear polarisation DS for Faraday rotation."""

        lam = (c.c / (self.freq * u.MHz)).to(u.m).value
        L = L * np.exp(-2j * RM * lam**2)

        return L

    def plot_lightcurve(self, stokes, polangle, fig=None, ax=None, imag=False):
        """Plot channel-averaged lightcurve."""

        lc = LightCurve(self, stokes, imag=imag)

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 5))
        fig, ax = lc.plot(fig, ax, polangle=polangle)

        return fig, ax

    def plot_spectrum(self, stokes, fig=None, ax=None, imag=False):
        """Plot time-averaged spectrum."""

        sp = Spectrum(self, stokes, imag=imag)

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 5))

        fig, ax = sp.plot(fig, ax)

        return fig, ax

    def _plot_ds(self, data, cmin, cmax, cmap, fig, ax):
        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))

        phasemax = 0.5 * self.fold_periods
        tmin, tmax = (-phasemax, phasemax) if self.fold else (self.tmin, self.tmax)
        norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))

        im = ax.imshow(
            np.transpose(data),
            extent=[tmin, tmax, self.fmin, self.fmax],
            aspect="auto",
            origin="lower",
            norm=norm,
            clim=(cmin, cmax),
            cmap=cmap,
        )

        ax.set_xlabel(self._timelabel)
        ax.set_ylabel("Frequency (MHz)")

        return fig, ax, im

    def plot_ds(self, stokes, cmax=20, imag=False, fig=None, ax=None):
        """Plot dynamic spectrum."""

        if stokes == "L":
            data = np.abs(self.data[stokes])
        else:
            data = self.data[stokes].imag if imag else self.data[stokes].real

        cmap = "plasma" if stokes in ["I", "L"] else "coolwarm"
        cmin = -2 if stokes in ["I", "L"] else -cmax

        fig, ax, im = self._plot_ds(data, cmin, cmax, cmap, fig, ax)

        ax.text(
            0.05,
            0.95,
            f"Stokes {stokes}",
            color="white",
            weight="heavy",
            path_effects=[pe.withStroke(linewidth=2, foreground="black")],
            transform=ax.transAxes,
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label("Flux Density (mJy)")

        return fig, ax

    def plot_pol_ds(self, fig=None, ax=None, mask_sigma=1):
        """Plot dynamic spectrum of polarisation fraction."""

        # Mask based on Stokes I RMS
        P = snr_mask(
            data=self.data["P"],
            noise=self.data["I"],
            n_sigma=mask_sigma,
        )

        fig, ax, im = self._plot_ds(P, 0, 100, "plasma", fig, ax)

        ax.text(
            0.05,
            0.95,
            r"$\sqrt{Q^2 + U^2 + V^2}/I$",
            color="white",
            weight="heavy",
            path_effects=[pe.withStroke(linewidth=2, foreground="black")],
            transform=ax.transAxes,
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label("Fractional Polarisation")

        return fig, ax

    def plot_polangle_ds(
        self,
        fig=None,
        ax=None,
        cmin=-180,
        cmax=180,
        mask_sigma=1,
    ):
        """Plot dynamic spectrum of polarisation angle."""

        # Mask based on Stokes I RMS
        PA = snr_mask(
            data=self.data["PA"],
            noise=self.data["L"],
            n_sigma=mask_sigma,
        )

        fig, ax, im = self._plot_ds(PA, cmin, cmax, "coolwarm", fig, ax)

        ax.text(
            0.05,
            0.95,
            "PA",
            color="white",
            weight="heavy",
            path_effects=[pe.withStroke(linewidth=2, foreground="black")],
            transform=ax.transAxes,
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label("PA (deg)")

        return fig, ax

    def plot_acf(self, stokes="I", contrast=0.4):
        """Plot 2D auto-correlation function of the dynamic spectrum."""

        acf2d = self.acf(stokes)

        # Plot 2D ACF
        acf_fig, acf_ax = plt.subplots(figsize=(7, 5))

        norm = ImageNormalize(acf2d, interval=ZScaleInterval(contrast=contrast))
        im = acf_ax.imshow(
            acf2d,
            extent=[0, self.tmax - self.tmin, 0, self.fmax - self.fmin],
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

        acf_ax.set_xlabel(f"Time Lag ({self.tunit})")
        acf_ax.set_ylabel("Frequency Lag (MHz)")

        # Plot zero frequency lag trace
        acfz_fig, acfz_ax = plt.subplots(figsize=(7, 5))

        zero_trace_acf = acf2d[-1, 1:]
        time_lag = np.linspace(0, self.tmax - self.tmin, len(zero_trace_acf))
        acfz_ax.plot(
            time_lag,
            zero_trace_acf,
            color="k",
        )

        acfz_ax.set_xlabel(f"Time Lag ({self.tunit})")
        acfz_ax.set_ylabel("ACF")

        acf_peaks, props = find_peaks(zero_trace_acf, prominence=(None, None))

        max_prom = np.argsort(props["prominences"])[::-1]
        self.peak_lags = time_lag[acf_peaks[max_prom]]

        acfz_ax.axvline(
            self.peak_lags[0],
            color="darkorange",
            ls="--",
        )

        logger.debug(f"Stokes {stokes} ACF peak at {self.peak_lags[0]*self.tunit:.3f}")

        return acf_fig, acf_ax, acfz_fig, acfz_ax


class TimeFreqSeries(ABC):

    def _construct_yaxis(self, avg_axis):
        """Construct y-axis averaging flux over x-axis."""

        # Catch RuntimeWarning that occurs when averaging empty time/freq slices
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            y = defaultdict()
            yerr = defaultdict()
            sqrtn = np.sqrt(self.ds.data["I"].shape[avg_axis])

            for stokes in self.stokes:
                data = self.ds.data[stokes]
                if stokes == "L":
                    # Compute L as complex magnitude of averaged DS
                    y[stokes] = np.abs(np.nanmean(data, axis=avg_axis))
                    yerr[stokes] = np.abs(np.nanstd(data, axis=avg_axis)) / sqrtn
                else:
                    ydata = data.imag if self.imag else data.real

                    y[stokes] = np.nanmean(ydata, axis=avg_axis)
                    yerr[stokes] = np.nanstd(data.imag, axis=avg_axis) / sqrtn

        return y, yerr

    def plot(self):

        # Overplot each specified polarisation
        for stokes in self.stokes:

            # Plot time/frequency series
            self.ax.errorbar(
                self.x,
                y=self.y[stokes],
                yerr=self.yerr[stokes],
                lw=1,
                color=COLORS[stokes],
                marker="o",
                markersize=1,
                label=stokes,
            )

        self.ax.set_ylabel("Flux Density (mJy)")
        self.ax.legend()

        return

    def save(self, savepath):
        values = self.valstart + self.x * self.unit

        df = pd.DataFrame({self.column: values})
        for stokes in self.stokes:
            df[f"flux_density_{stokes}"] = self.y[stokes].real.reshape(1, -1)[0]
            df[f"flux_density_{stokes}_err"] = self.yerr[stokes]

        df.dropna().to_csv(savepath, index=False)

        return


@dataclass
class LightCurve(TimeFreqSeries):
    ds: DynamicSpectrum
    stokes: str
    pa_sigma: int = 2
    imag: bool = False

    def __post_init__(self):
        self.column = "time"
        self.unit = self.ds.tunit
        self.valstart = self.ds.header["time_start"]

        # Set time/phase axis limits if folded
        phasemax = 0.5 * self.ds.fold_periods
        valmin, valmax = (
            (-phasemax, phasemax) if self.ds.fold else (self.ds.tmin, self.ds.tmax)
        )

        # Construct time and flux axes
        bins = self.ds.data["I"].shape[0]
        interval = (valmax - valmin) / bins
        self.x = np.array([valmin + i * interval for i in range(bins)])
        self.y, self.yerr = self._construct_yaxis(avg_axis=1)

        Q = self.ds.data["Q"]
        U = self.ds.data["U"]
        L = self.ds.data["L"]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            signal_L = np.abs(np.nanmean(L, axis=1))
            noise_L = np.abs(np.nanstd(L, axis=1)) / np.sqrt(L.shape[1])
            Q = np.nanmean(Q.real, axis=1)
            U = np.nanmean(U.real, axis=1)

        # Mask low signifigance PA
        self.polangle = 0.5 * np.arctan2(U, Q) * u.rad.to(u.deg)
        mask = signal_L < self.pa_sigma * noise_L
        self.polangle[mask] = np.nan

    def plot(self, fig, ax, polangle=False):
        self.fig = fig
        self.ax = ax

        # Set time axis label
        self.ax.set_xlabel(self.ds._timelabel)

        # Plot with lightcurve/spectrum independent parameters
        super().plot()

        pad = (self.x.max() - self.x.min()) * 0.05
        self.ax.set_xlim([self.x.min() - pad, self.x.max() + pad])

        if polangle:
            divider = make_axes_locatable(self.ax)
            ax2 = divider.append_axes("top", size="25%", pad=0.1)

            ax2.scatter(
                self.x,
                self.polangle,
                color="k",
                marker="o",
                s=1,
            )
            ax2.axhline(
                0,
                ls=":",
                color="k",
                alpha=0.5,
            )
            ax2.set_xticklabels([])
            ax2.set_xlim([self.x.min() - pad, self.x.max() + pad])
            ax2.set_ylim(-120, 120)
            ax2.set_yticks([-90, 0, 90])
            ax2.set_ylabel(r"$\chi$ (deg)")

        return self.fig, self.ax


@dataclass
class Spectrum(TimeFreqSeries):
    ds: DynamicSpectrum
    stokes: str
    imag: bool = False

    def __post_init__(self):
        self.column = "frequency"
        self.unit = u.MHz
        self.valstart = 0

        # Construct frequency axis
        bins = self.ds.data["I"].shape[1]
        interval = (self.ds.fmax - self.ds.fmin) / bins
        self.x = np.array([self.ds.fmin + i * interval for i in range(bins)])
        self.y, self.yerr = self._construct_yaxis(avg_axis=0)

    def plot(self, fig, ax):
        self.fig = fig
        self.ax = ax

        # Set frequency axis label
        ax.set_xlabel("Frequency (MHz)")

        # Plot with lightcurve/spectrum independent parameters
        super().plot()

        return self.fig, self.ax
