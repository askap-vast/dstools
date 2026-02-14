import logging
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from importlib.metadata import version
from typing import Optional, Sequence

import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from scipy.signal import correlate

from dstools.polarisation import (
    RMTimeSeries,
    derotate_dynamic_spectrum,
    dynamic_rmts_from_fdf,
    peak_rm_from_fdf,
    rm_synthesis,
    smooth_rm_timeseries,
)
from dstools.utils import LOCATIONS, parse_time, rebin, rebin2D, slice_array

logger = logging.getLogger(__name__)

FlagRanges = Sequence[tuple[float, float]]


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
    flag_channels: Optional[FlagRanges] = None
    flag_times: Optional[FlagRanges] = None
    flag_imag_snr: Optional[float] = 10000

    tunit: u.Quantity = u.hour
    corr_dumptime: float = 10.1

    barycentre: bool = False
    derotate: bool = False
    rm_smooth_method: Optional[str] = None
    dedisperse: bool = False
    RM: Optional[float] = None
    RM_reffreq: Optional[u.Quantity] = None
    DM: Optional[float] = None
    DM_reffreq: Optional[u.Quantity] = None

    fold: bool = False
    period: Optional[float] = None
    period_offset: float = 0.0
    fold_periods: int = 2
    phase_bins: Optional[int] = None

    absolute_times: bool = True
    calscans: bool = True
    trim: bool = True

    def __post_init__(self):
        # Load instrumental polarisation time/frequency/uvdist arrays
        XX, XY, YX, YY = self._load_data()

        # Flag specified channels / times
        if self.flag_channels is not None:
            XX, XY, YX, YY = self._flag_channels(XX, XY, YX, YY)
        if self.flag_times is not None:
            XX, XY, YX, YY = self._flag_times(XX, XY, YX, YY)
        if self.flag_imag_snr is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                XX, XY, YX, YY = self._flag_imag_snr(XX, XY, YX, YY)

        # Insert calibrator scan breaks
        XX, XY, YX, YY = self._stack_cal_scans(XX, XY, YX, YY)

        # Incoherently dedisperse at DM
        if self.dedisperse and self.DM is not None:
            XX = self._dedisperse(XX)
            XY = self._dedisperse(XY)
            YX = self._dedisperse(YX)
            YY = self._dedisperse(YY)

        # Average data in time and frequency
        XX, XY, YX, YY = self._rebin(XX, XY, YX, YY)

        # Fold data to selected period
        if self.fold:
            if not self.period:
                raise ValueError("Must pass period argument when folding.")

            # Compute folding stats
            dt = self.time[1] - self.time[0]

            if self.phase_bins is None:
                self.phase_bins = round(self.period / dt)

            samples_per_period = self.period / dt
            oversample = self.phase_bins / samples_per_period
            num_folds = (self.tmax - self.tmin) / self.period

            logger.info(f"Folding at {self.period * self.tunit} period.")
            logger.info(
                f"Oversampling by factor of {oversample:.1f} "
                f"with {num_folds:.1f} folds and {self.phase_bins:.1f} phase bins"
            )

            # Phase fold instrumental pols
            XX = self._fold(XX)
            XY = self._fold(XY)
            YX = self._fold(YX)
            YY = self._fold(YY)

            # Disable plotting with absolute times as we will plot phase instead
            self.absolute_times = False
            self.time = rebin(len(self.time), len(XX), axis=0) @ self.time

        # Store time and frequency resolution
        if len(self.time) > 1:
            self.time_res = (self.time[1] - self.time[0]) * self.tunit
        if len(self.freq) > 1:
            self.freq_res = (self.freq[1] - self.freq[0]) * u.MHz

        self.header.update(
            {
                "time_resolution": f"{self.time_res.to(u.s):.3f}",
                "freq_resolution": f"{self.freq_res.to(u.MHz):.2f}",
            }
        )

        # Compute Stokes products and store in data attribute
        self._make_stokes(XX, XY, YX, YY)

        return

    def __str__(self):
        str_rep = ""
        for attr in self.header:
            str_rep += f"{attr}: {self.header[attr]}\n"
        return str_rep

    def _fold(self, data):
        """Fold data at specified period with linearly interpolated phase binning."""

        # Compute fractional phases of each sampled timestep
        phases = ((self.time - self.time[0]) / self.period + self.period_offset) % 1.0

        # Create bins spanning the phase range
        bin_floats = phases * self.phase_bins
        bin_lower = np.floor(bin_floats).astype(int)
        bin_upper = (bin_lower + 1) % self.phase_bins

        # Compute weights for neighbouring bins linearly interpolate sampled point to bin centres
        weights_upper = bin_floats - bin_lower
        weights_lower = 1.0 - weights_upper

        folded = np.zeros((self.phase_bins, data.shape[1]), dtype=np.complex128)
        counts = np.zeros((self.phase_bins, data.shape[1]), dtype=np.float64)

        # Iterate through each timestep assigning weighted contribution to neighbouring phase bins
        for time in range(data.shape[0]):
            valid = ~np.isnan(data[time])

            folded[bin_lower[time], valid] += data[time, valid] * weights_lower[time]
            folded[bin_upper[time], valid] += data[time, valid] * weights_upper[time]
            counts[bin_lower[time], valid] += weights_lower[time]
            counts[bin_upper[time], valid] += weights_upper[time]

        # Normalise binned counts
        with np.errstate(invalid="ignore", divide="ignore"):
            folded = folded / counts

        # Tile number of fold_periods together for display
        folded = np.tile(folded, (self.fold_periods, 1))

        return folded

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
        """Validate the HDF5 DS file."""

        # Check if baselines have been pre-averaged and disable uvdist selection if so.
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
                "DS is already baseline averaged, disabling uvdist selection."
            )
            self.minuvdist = 0
            self.maxuvdist = np.inf
            self.minuvwave = 0
            self.maxuvwave = np.inf

        # Check extraction and library versions of DStools
        ds_version = datafile.attrs.get("dstools_version", "1.0.0")
        dstools_version = version("radio-dstools")

        if ds_version != dstools_version:
            logger.warning(
                f"Using DStools v{dstools_version} to open DS extracted with DStools v{ds_version}."
            )

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

        # Set an initial time / freq resolution.
        # This will be updated after folding / data selection,
        # but we set a default to fall back on in case further
        # processing restricts to single channel / integration
        self.time_res = (time[1] - time[0]) * self.tunit
        self.freq_res = (freq[1] - freq[0]) * u.MHz

        # Flip ATCA L-band frequency axis to intuitive order
        if freq[0] > freq[-1]:
            XX = np.flip(XX, axis=1)
            XY = np.flip(XY, axis=1)
            YX = np.flip(YX, axis=1)
            YY = np.flip(YY, axis=1)

            freq = np.flip(freq)

        # Optionally remove flagged channels at top/bottom of band
        if self.trim:
            # Create binary mask identifying non-nan values
            # after summing across time and polarisation axes
            full = np.nansum([XX, XY, YX, YY], axis=(0, 1))
            full[full == 0.0 + 0.0j] = np.nan
            allpols = np.isfinite(full)

            # Set minimum and maximum non-nan channel indices
            minchan = np.argmax(allpols)
            maxchan = 0 if allpols[-1] else -np.argmax(allpols[::-1]) + 1
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
            maxtime = -(np.argmax((time - time[0] < self.maxtime)[::-1])) + 1
        else:
            maxtime = 0

        # Convert times to sensible format
        t = Time(
            (time * self.tunit).to(u.day),
            format="mjd",
            scale="utc",
        )

        # Correct to barycentric dynamic time
        if self.barycentre:
            t = self._barycentre_times(t)

        # Set start time of observation in appopriate timescale (UTC or TDB)
        time_start = getattr(t[0], t[0].scale)
        self.header["time_start"] = time_start.iso
        self.header["time_scale"] = t[0].scale

        # Set time array relative to time_start
        time = (t - t[0]).value * u.day.to(self.tunit)

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

    def _flag_imag_snr(self, XX, XY, YX, YY):
        I = (XX + YY) / 2
        Q = (XX - YY) / 2
        U = (XY + YX) / 2
        V = 1j * (YX - XY) / 2

        chan_mask = np.full(XX.shape[1], False)
        time_mask = np.full(XX.shape[0], False)

        for vis in [I.imag]:  # , Q.imag, U.imag, V.imag]:
            sp = np.nanmean(vis, axis=0)
            lc = np.nanmean(vis, axis=1)
            sqrt_f = np.sqrt(vis.shape[0])
            sqrt_t = np.sqrt(vis.shape[0])

            sp_snr = np.abs(self.flag_imag_snr * np.nanstd(vis, axis=0) / sqrt_f)
            lc_snr = np.abs(self.flag_imag_snr * np.nanstd(vis, axis=1) / sqrt_t)

            chan_mask = chan_mask | (sp > sp_snr) | (sp < -sp_snr)
            time_mask = time_mask | (lc > lc_snr) | (lc < -lc_snr)

        XX[:, chan_mask] = np.nan
        XY[:, chan_mask] = np.nan
        YX[:, chan_mask] = np.nan
        YY[:, chan_mask] = np.nan

        XX[time_mask, :] = np.nan
        XY[time_mask, :] = np.nan
        YX[time_mask, :] = np.nan
        YY[time_mask, :] = np.nan

        # rfi_mask = self.flag_imag_snr
        # print(sp_snr, lc_snr)

        return XX, XY, YX, YY

    def _flag_channels(self, XX, XY, YX, YY):
        for flagrange in self.flag_channels:
            minfreq, maxfreq = flagrange

            logger.debug(f"Flagging channels from {minfreq}-{maxfreq} MHz")
            minfreq = np.argmax(self.freq > minfreq)
            maxfreq = np.argmax(self.freq > maxfreq)

            XX[:, minfreq:maxfreq] = np.nan
            XY[:, minfreq:maxfreq] = np.nan
            YX[:, minfreq:maxfreq] = np.nan
            YY[:, minfreq:maxfreq] = np.nan

        return XX, XY, YX, YY

    def _flag_times(self, XX, XY, YX, YY):
        for flagrange in self.flag_times:
            tstart = self.header["time_start"]

            mintime = parse_time(flagrange[0], self.tunit, tstart)
            maxtime = parse_time(flagrange[1], self.tunit, tstart)

            logger.debug(
                f"Flagging time range from {mintime:.2f}-{maxtime:.2f} {self.tunit}"
            )
            mintime = np.argmax(self.time > mintime)
            maxtime = np.argmax(self.time > maxtime)

            XX[mintime:maxtime, :] = np.nan
            XY[mintime:maxtime, :] = np.nan
            YX[mintime:maxtime, :] = np.nan
            YY[mintime:maxtime, :] = np.nan

        return XX, XY, YX, YY

    def _barycentre_times(self, time: Time):
        """Apply corrections to Barycentric Dynamical Timescale."""

        ra, dec = self.header.get("phasecentre").split()
        location = LOCATIONS.get(self.header.get("telescope"))

        target_coord = SkyCoord(
            ra=ra,
            dec=dec,
            unit="hourangle,deg",
            frame="icrs",
        )
        time = Time(time, format="mjd", scale="utc", location=location)
        ltt_bary = time.light_travel_time(target_coord)
        time = time.tdb + ltt_bary

        self._timescale = "TDB"

        return time

    def _dedisperse(self, array):
        """Incoherently dedisperse using Fourier shift."""

        # Compute time-domain delays
        a = (c.e.si**2 / (8 * np.pi**2 * c.eps0 * c.m_e * c.c)).to(
            u.GHz**2 * u.cm**3 * u.pc**-1 * u.ms
        )
        DM = self.DM * u.pc / u.cm**3
        reffreq = self.freq[-1] if self.DM_reffreq is None else self.DM_reffreq
        tau = (self.freq * u.MHz) ** -2 - (reffreq * u.MHz) ** -2

        time_resolution = (self.time[1] - self.time[0]) * self.tunit
        dt = (a * DM * tau / time_resolution).to(1)

        # FFT
        array = np.fft.fft(array, axis=0)
        fsamp = np.fft.fftfreq(array.shape[0])

        # Phase shift
        phasor = np.exp(2j * np.pi * np.outer(fsamp, dt))
        array = array * phasor.value

        # Invert back to time/freq space
        array = np.fft.ifft(array, axis=0)

        return array

    def _stack_cal_scans(self, XX, XY, YX, YY):
        """Insert null data representing off-source time."""

        scan_start_idx, scan_end_idx = self._get_scan_intervals()
        dt = self.time_res.value

        # Calculate number of cycles in each calibrator/stow break
        time_end_break = self.time[scan_start_idx[1:]]
        time_start_break = self.time[scan_end_idx[:-1]]

        # Count number of samples within each calibrator / stow break
        num_break_cycles = np.append((time_end_break - time_start_break), 0) / dt
        num_channels = self.header["channels"]

        # Create initial time-slice to start stacking target and calibrator scans together
        stacked_XX = stacked_XY = stacked_YX = stacked_YY = np.zeros(
            (1, num_channels),
            dtype=complex,
        )
        stacked_time = np.zeros(1)

        for start_index, end_index, num_scans in zip(
            scan_start_idx,
            scan_end_idx,
            num_break_cycles,
        ):
            # Select each contiguous on-target chunk of data
            XX_chunk = XX[start_index : end_index + 1, :]
            XY_chunk = XY[start_index : end_index + 1, :]
            YX_chunk = YX[start_index : end_index + 1, :]
            YY_chunk = YY[start_index : end_index + 1, :]
            time_chunk = self.time[start_index : end_index + 1]

            # Make array of complex NaN's for subsequent calibrator / stow gaps
            # and append to each on-target chunk of data.
            if self.calscans and num_scans > 0:
                # We round down to the nearest integer number of correlator cycles
                # to populate the nan-break. The final cycle will therefore be slightly
                # longer than the rest, but this only affects the visual presentation
                # of the lightcurve / dynamic spectrum, not the timestamps.
                num_timesteps = int(round(num_scans) - 1)
                num_nans = (num_timesteps, num_channels)
                nan_chunk = np.full(num_nans, np.nan + np.nan * 1j)

                XX_chunk = np.ma.vstack([XX_chunk, nan_chunk])
                XY_chunk = np.ma.vstack([XY_chunk, nan_chunk])
                YX_chunk = np.ma.vstack([YX_chunk, nan_chunk])
                YY_chunk = np.ma.vstack([YY_chunk, nan_chunk])

                time_break_start = self.time[end_index] + dt
                time_break_scans = time_break_start + np.arange(num_nans[0]) * dt
                time_chunk = np.append(time_chunk, time_break_scans)

            stacked_XX = np.ma.vstack([stacked_XX, XX_chunk])
            stacked_XY = np.ma.vstack([stacked_XY, XY_chunk])
            stacked_YX = np.ma.vstack([stacked_YX, YX_chunk])
            stacked_YY = np.ma.vstack([stacked_YY, YY_chunk])
            stacked_time = np.append(stacked_time, time_chunk)

        stacked_XX = stacked_XX[1:]
        stacked_XY = stacked_XY[1:]
        stacked_YX = stacked_YX[1:]
        stacked_YY = stacked_YY[1:]

        self.time = stacked_time[1:]
        self.dts = self.time[1:] - self.time[:-1]

        return stacked_XX, stacked_XY, stacked_YX, stacked_YY

    def _rebin(self, XX, XY, YX, YY):
        """Bin data in time and frequency."""

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
            I = 0.5 * np.nansum([XX, +YY], axis=0)
            Q = 0.5 * np.nansum([XX, -YY], axis=0)
            U = 0.5 * np.nansum([XY, +YX], axis=0)
            V = 0.5 * np.nansum([YX, -XY], axis=0) * 1j
        elif feedtype == "circular":
            I = 0.5 * np.nansum([XX, +YY], axis=0)
            Q = 0.5 * np.nansum([XY, +YX], axis=0)
            U = 0.5 * np.nansum([XY, -YX], axis=0) * 1j
            V = 0.5 * np.nansum([XX, -YY], axis=0)
        else:
            raise ValueError(
                f"Feed type {feedtype} not recognised, should be either 'linear' or 'circular'."
            )

        # Identically zero values should only arise from
        # np.nansum converting an all-NaN range to 0 + 0j
        # Set these back to NaN for clarity
        I[I == 0 + 0j] = np.nan
        Q[Q == 0 + 0j] = np.nan
        U[U == 0 + 0j] = np.nan
        V[V == 0 + 0j] = np.nan

        L = Q.real + 1j * U.real

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
        }

        return

    def _get_rmts(
        self,
        L: np.ndarray,
        Li: np.ndarray,
        model: str,
        mask_array: np.ndarray,
        snr_min: float,
        poly_deg: int,
    ):
        # Compute 2D FDF with RM synthesis
        fdf = rm_synthesis(L, self.freq)
        fdf_im = rm_synthesis(Li, self.freq)

        # The 'peak' model returns the Faraday depth at the timestep
        # of the maximum value in mask_array.
        if model == "peak":
            rm = peak_rm_from_fdf(fdf, mask_array)
            return rm

        # Otherwise we get a per-timestep RM time series
        rmts = dynamic_rmts_from_fdf(fdf, fdf_im, mask_array, snr_min=snr_min)

        # Smooth the RM time series with a model, default to constant median RM
        rmts = smooth_rm_timeseries(rmts, mode=model, poly_deg=poly_deg).model

        return rmts

    def derotate_faraday(
        self,
        RM: Optional[float | RMTimeSeries] = None,
        model: str = "peak",
        mask_pol: str = "I",
        snr_min: float = 10.0,
        poly_deg: int = 1,
    ) -> float | RMTimeSeries:
        # Default to masking with Stokes I
        mask_array = self.data.get(mask_pol, "I")

        # Compute complex L and
        Q = self.data["Q"]
        U = self.data["U"]
        L = Q.real + 1j * U.real
        Li = Q.imag + 1j * U.imag

        if RM is None:
            RM = self._get_rmts(
                L,
                Li,
                model=model,
                mask_array=mask_array,
                snr_min=snr_min,
                poly_deg=poly_deg,
            )

        # Derotate complex L with RM time series
        L = derotate_dynamic_spectrum(L, self.freq, RM)
        Li = derotate_dynamic_spectrum(Li, self.freq, RM)

        # Compute complex Q and U from L
        Q = L.real + 1j * Li.real
        U = L.imag + 1j * Li.imag

        # Update data arrays
        self.data["Q"] = Q
        self.data["U"] = U
        self.data["L"] = L

        return RM

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


class TimeFreqSeries(ABC):
    """Abstract base class for construction of common elements of lightcurves / 1D spectra."""

    @abstractmethod
    def x(self):
        """An array representing the x-axis of the averaged data (time or frequency)."""

    def _construct_yaxis(self, avg_axis):
        # Catch RuntimeWarning that occurs when averaging empty time/freq slices
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            self.flux = defaultdict()
            self.flux_err = defaultdict()

            # Compute flux / errors in each Stokes parameter averaging over the DS
            for stokes in "IQUV":
                data = self.ds.data[stokes]
                sqrtn = np.sqrt(np.isfinite(data).sum(axis=avg_axis))

                ydata = data.imag if self.imag else data.real

                self.flux[stokes] = np.nanmean(ydata, axis=avg_axis)
                self.flux_err[stokes] = np.nanstd(data.imag, axis=avg_axis) / sqrtn

        with np.errstate(invalid="ignore", divide="ignore"):
            self._calc_polarisation_params()

        return

    def _calc_polarisation_params(self):
        """Compute polarisation angle, ellipticity, and fractional polarisation."""

        I = self.flux["I"]
        Ierr = self.flux_err["I"]
        Q = self.flux["Q"]
        Qerr = self.flux_err["Q"]
        U = self.flux["U"]
        Uerr = self.flux_err["U"]
        V = self.flux["V"]
        Verr = self.flux_err["V"]
        L = self.flux["L"] = np.sqrt(Q**2 + U**2)
        Lerr = self.flux_err["L"] = 1 / L * np.sqrt((Q * Qerr) ** 2 + (U * Uerr) ** 2)

        P = np.sqrt(L**2 + V**2)
        Perr = (1 / P) * np.sqrt((L * Lerr) ** 2 + (V * Verr) ** 2)

        self.pol_fraction = P / I

        self.polangle = 0.5 * np.arctan2(U, Q) * u.rad.to(u.deg)
        self.ellipticity = 0.5 * np.arctan2(V, L) * u.rad.to(u.deg)
        self.linear_fraction = np.abs(L / I)
        self.circular_fraction = np.abs(V / I)

        # Propagate errors
        qu_err = (Q * Uerr) ** 2 + (U * Qerr) ** 2
        vl_err = (V * Lerr) ** 2 + (L * Verr) ** 2
        pi_err = (Perr / P) ** 2 + (Ierr / I) ** 2
        li_err = (Lerr / L) ** 2 + (Ierr / I) ** 2
        vi_err = (Verr / V) ** 2 + (Ierr / I) ** 2
        self.polangle_err = (0.5 * np.sqrt(qu_err) / L**2) * u.rad.to(u.deg)
        self.ellipticity_err = (0.5 * np.sqrt(vl_err) / P**2) * u.rad.to(u.deg)
        self.pol_fraction_err = np.abs(P / I) * np.sqrt(pi_err)
        self.linear_fraction_err = np.abs(L / I) * np.sqrt(li_err)
        self.circular_fraction_err = np.abs(V / I) * np.sqrt(vi_err)

        # Mask low signifiance points
        L_mask = L < self.pol_sigma * Lerr
        V_mask = np.abs(V) < self.pol_sigma * Verr
        P_mask = P < self.pol_sigma * Perr

        # Remove any isolated unmasked values (likely noise)
        isolated = L_mask[:-2] & L_mask[2:]
        L_mask[1:-1][isolated] = True
        isolated = V_mask[:-2] & V_mask[2:]
        V_mask[1:-1][isolated] = True
        isolated = P_mask[:-2] & P_mask[2:]
        P_mask[1:-1][isolated] = True

        self.polangle[L_mask] = np.nan
        self.polangle_err[L_mask] = np.nan
        self.ellipticity[P_mask] = np.nan
        self.ellipticity_err[P_mask] = np.nan
        self.pol_fraction[P_mask] = np.nan
        self.pol_fraction_err[P_mask] = np.nan
        self.linear_fraction[L_mask] = np.nan
        self.linear_fraction_err[L_mask] = np.nan
        self.circular_fraction[V_mask] = np.nan
        self.circular_fraction_err[V_mask] = np.nan

        return

    def save(self, savepath, stokes="IQUVL", include_pols: bool = False):
        values = self.valstart + self.x * self.unit

        df = pd.DataFrame({self.column: values})
        for s in stokes:
            df[f"flux_density_{s}"] = self.flux[s].real.reshape(1, -1)[0]
            df[f"flux_density_{s}_err"] = self.flux_err[s]

        if include_pols:
            df["polarisation_angle"] = self.polangle
            df["polarisation_angle_err"] = self.polangle_err
            df["ellipticity"] = self.ellipticity
            df["ellipticity_err"] = self.ellipticity_err
            df["pol_fraction"] = self.pol_fraction
            df["pol_fraction_err"] = self.pol_fraction_err
            df["linear_fraction"] = self.linear_fraction
            df["linear_fraction_err"] = self.linear_fraction_err
            df["circular_fraction"] = self.circular_fraction
            df["circular_fraction_err"] = self.circular_fraction_err

        df.dropna().to_csv(savepath, index=False)

        return


@dataclass
class LightCurve(TimeFreqSeries):
    ds: DynamicSpectrum
    imag: bool = False
    pol_sigma: float = 4

    def __post_init__(self):
        self.column = "time"
        self.unit = self.ds.tunit
        self.valstart = Time(self.ds.header["time_start"])

        # Construct time and flux axes, using phase if folding enabled
        if self.ds.fold:
            phasemax = 0.5 * self.ds.fold_periods
            phasebins = self.ds.data["I"].shape[0]
            self.time = np.linspace(-phasemax, phasemax, phasebins)
        else:
            self.time = self.ds.time

        self._construct_yaxis(avg_axis=1)

        return

    @property
    def x(self):
        return self.time


@dataclass
class Spectrum(TimeFreqSeries):
    ds: DynamicSpectrum
    imag: bool = False
    pol_sigma: float = 4

    def __post_init__(self):
        self.column = "frequency"
        self.unit = u.MHz
        self.valstart = 0

        # Construct frequency axis
        bins = self.ds.data["I"].shape[1]
        interval = (self.ds.fmax - self.ds.fmin) / bins
        self.frequency = np.array([self.ds.fmin + i * interval for i in range(bins)])
        self._construct_yaxis(avg_axis=0)

    @property
    def x(self):
        return self.frequency
