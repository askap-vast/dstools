import logging
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from importlib.metadata import version
from typing import Optional

import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import pandas as pd
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
from scipy.signal import correlate

from dstools.rm import PolObservation
from dstools.utils import rebin, rebin2D, slice_array

logger = logging.getLogger(__name__)

LOCATIONS = {
    "ATCA": EarthLocation(
        lat=Angle("-30:18:46.385", unit=u.degree),
        lon=Angle(149.5501388, unit=u.degree),
        height=236.87 * u.m,
    ),
    "GMRT": EarthLocation.of_site("GMRT"),
    "VLA": EarthLocation.of_site("vla"),
    "MeerKAT": EarthLocation.of_site("MeerKAT"),
    "ASKAP": EarthLocation.of_site("ASKAP"),
}


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

    barycentre: bool = False
    derotate: bool = False
    dedisperse: bool = False
    RM: float = None
    DM: float = None
    DM_reffreq: u.Quantity = None

    fold: bool = False
    period: Optional[float] = None
    period_offset: float = 0.0
    fold_periods: int = 2

    absolute_times: bool = True
    calscans: bool = True
    trim: bool = True

    def __post_init__(self):
        # Load instrumental polarisation time/frequency/uvdist arrays
        XX, XY, YX, YY = self._load_data()

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

            XX = self._fold(XX)
            XY = self._fold(XY)
            YX = self._fold(YX)
            YY = self._fold(YY)

            # Disable plotting with absolute times as we will plot phase instead
            self.absolute_times = False
            self.time = rebin(len(self.time), len(XX), axis=0) @ self.time

        # Store time and frequency resolution
        self.time_res = (self.time[1] - self.time[0]) * self.tunit
        self.freq_res = (self.freq[1] - self.freq[0]) * u.MHz
        self.header.update(
            {
                "time_resolution": f"{self.time_res.to(u.s):.3f}",
                "freq_resolution": f"{self.freq_res.to(u.MHz):.2f}",
            }
        )

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
                "DS is already baseline averaged, disabling uvdist selection."
            )
            self.minuvdist = 0
            self.maxuvdist = np.inf
            self.minuvwave = 0
            self.maxuvwave = np.inf

        ds_version = datafile.attrs.get("dstools_version", "1.0.0")
        dstools_version = version("radio-dstools")

        if ds_version != dstools_version:
            logger.warning(
                f"Attempting to open v{ds_version} DS with DStools v{dstools_version}."
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

    def _barycentre_times(self, time: Time):
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
        a = (c.e.si**2 / (4 * np.pi**2 * c.eps0 * c.m_e * c.c)).to(
            u.GHz**2 * u.cm**3 * u.pc**-1 * u.ms
        )
        DM = self.DM * u.pc / u.cm**3
        reffreq = self.freq[-1] if self.DM_reffreq is None else self.DM_reffreq
        tau = (self.freq * u.MHz) ** -2 - (reffreq * u.MHz) ** -2
        dt = (a * DM * tau / self.time_res).to(1)

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

        # Calculate number of cycles in each calibrator/stow break
        time_end_break = self.time[scan_start_idx[1:]]
        time_start_break = self.time[scan_end_idx[:-1]]

        # Count number of samples within each calibrator / stow break
        dt = self.time[1] - self.time[0]
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
            if self.RM is None:
                self.RM = self.rm_synthesis(I, Q, U)

            # Build L from imaginary components
            Li = Q.imag + 1j * U.imag

            # Derotate real and imaginary L
            L = self.derotate_faraday(L)
            Li = self.derotate_faraday(Li)

            # Compute complex Q and U from L
            Q = L.real + 1j * Li.real
            U = L.imag + 1j * Li.imag
        else:
            self.polobs = None

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


    def derotate_faraday(self, L):
        """Correct linear polarisation DS for Faraday rotation."""

        lam = (c.c / (self.freq * u.MHz)).to(u.m).value
        L = L * np.exp(-2j * self.RM * lam**2)

        return L


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
            sqrtn = np.sqrt(self.ds.data["I"].shape[avg_axis])

            # Compute flux / errors in each Stokes parameter averaging over the DS
            for stokes in "IQUV":
                data = self.ds.data[stokes]
                ydata = data.imag if self.imag else data.real

                self.flux[stokes] = np.nanmean(ydata, axis=avg_axis)
                self.flux_err[stokes] = np.nanstd(data.imag, axis=avg_axis) / sqrtn

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

        P = self.pol_fraction = np.sqrt(L**2 + V**2)
        self.pol_fraction_err = 1 / P * np.sqrt((L * Lerr) ** 2 + (V * Verr) ** 2)

        self.polangle = 0.5 * np.arctan2(U, Q) * u.rad.to(u.deg)
        self.ellipticity = 0.5 * np.arctan2(V, L) * u.rad.to(u.deg)
        self.linear_fraction = np.abs(L / I)
        self.circular_fraction = np.abs(V / I)

        # Propagate errors
        qu_err = (Q * Qerr) ** 2 + (U * Uerr) ** 2
        vl_err = (V * Verr) ** 2 + (L * Lerr) ** 2
        li_err = (Lerr / L) ** 2 + (Ierr / I) ** 2
        vi_err = (Verr / V) ** 2 + (Ierr / I) ** 2
        self.polangle_err = (0.5 * np.sqrt(qu_err) / L**2) * u.rad.to(u.deg)
        self.ellipticity_err = (0.5 * np.sqrt(vl_err) / P**2) * u.rad.to(u.deg)
        self.linear_fraction_err = np.abs(L / I) * np.sqrt(li_err)
        self.circular_fraction_err = np.abs(V / I) * np.sqrt(vi_err)

        # Mask low signifiance points
        mask = I < self.pa_sigma * Ierr

        # Remove any isolated unmasked values (likely noise)
        isolated = mask[:-2] & mask[2:]
        mask[1:-1][isolated] = True

        self.polangle[mask] = np.nan
        self.polangle_err[mask] = np.nan
        self.ellipticity[mask] = np.nan
        self.ellipticity_err[mask] = np.nan
        self.linear_fraction[mask] = np.nan
        self.linear_fraction_err[mask] = np.nan
        self.circular_fraction[mask] = np.nan
        self.circular_fraction_err[mask] = np.nan

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
    pa_sigma: int = 5

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
    pa_sigma: int = 5

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
