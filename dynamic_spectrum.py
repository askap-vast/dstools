import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.visualization import ZScaleInterval, ImageNormalize


class DynamicSpectrum:
    def __init__(self, project, band="L", calscans=True, prefix=""):
        self.src = project.split("/")[-1]
        self.band = band

        self.ds_path = f"{project}/dynamic_spectra/{band}"
        self.prefix = prefix

        self._load_data()

        if calscans:
            self._stack_cal_scans()
        self._make_stokes()

        self.tmin = self.time[0]
        self.tmax = self.time[-1]
        self.fmin = self.freq[0]
        self.fmax = self.freq[-1]

    def _get_scan_intervals(self):
        """Find indices of start/end of each calibrator scan cycle"""

        t0 = Time(self.time[0] / 24.0, format="mjd", scale="utc")
        tN = Time(self.time[-1] / 24.0, format="mjd", scale="utc")
        t0.format = "iso"
        tN.format = "iso"

        dts = [0]
        dts.extend([self.time[i] - self.time[i - 1] for i in range(1, len(self.time))])
        self.dts = np.array(dts)

        # Locate indices signaling beginning of cal-scan (scan intervals longer than ~10 seconds)
        scan_start_indices = np.where(np.abs(self.dts) > 10.1 / 3600)[0]

        # End indices are just prior to the next start index, then insert first scan start index
        scan_end_indices = scan_start_indices - 1
        self.scan_start_indices = np.insert(scan_start_indices, 0, 0)
        self.scan_end_indices = np.append(scan_end_indices, len(self.time) - 1)

    def _load_data(self):
        """Load instrumental pols and time/freq data, converting to MHz, s, and mJy"""

        file_prefix = f"{self.ds_path}/{self.prefix}"

        self.freq = np.load(f"{file_prefix}freq.npy", allow_pickle=True) / 1e6
        self.time = np.load(f"{file_prefix}time.npy", allow_pickle=True) / 3600
        self.time_start = Time(self.time[0] / 24.0, format="mjd", scale="utc")
        self.time_start.format = "iso"

        self.time -= self.time[0]

        self.XX = np.load(f"{file_prefix}dynamic_spectra_XX.npy", allow_pickle=True) * 1e3
        self.XY = np.load(f"{file_prefix}dynamic_spectra_XY.npy", allow_pickle=True) * 1e3
        self.YX = np.load(f"{file_prefix}dynamic_spectra_YX.npy", allow_pickle=True) * 1e3
        self.YY = np.load(f"{file_prefix}dynamic_spectra_YY.npy", allow_pickle=True) * 1e3

    def _stack_cal_scans(self):
        """Insert null data representing off-source time"""

        self._get_scan_intervals()
        avg_scan_dt = np.median(self.dts[: self.scan_end_indices[0]])
        num_scan_intervals = (
            self.time[self.scan_start_indices[1:]] - self.time[self.scan_end_indices[:-1]]
        ) / avg_scan_dt
        num_channels = self.XX.shape[1]

        # Create initial time-slice to start stacking target and calibrator scans together
        new_data_XX = new_data_XY = new_data_YX = new_data_YY = np.zeros(
            (1, num_channels), dtype=complex
        )

        for start_index, end_index, num_scans in zip(
            self.scan_start_indices[:-1], self.scan_end_indices[:-1], num_scan_intervals
        ):

            # Select each contiguous on-target chunk of data
            XX_chunk = self.XX[start_index:end_index, :]
            XY_chunk = self.XY[start_index:end_index, :]
            YX_chunk = self.YX[start_index:end_index, :]
            YY_chunk = self.YY[start_index:end_index, :]

            # Make array of complex NaN's for subsequent calibrator / stow gaps
            # and append to each on-target chunk of data
            nan_chunk = np.full((int(num_scans), num_channels), np.nan + np.nan * 1j)

            new_data_XX = np.vstack([new_data_XX, np.vstack([XX_chunk, nan_chunk])])
            new_data_XY = np.vstack([new_data_XY, np.vstack([XY_chunk, nan_chunk])])
            new_data_YX = np.vstack([new_data_YX, np.vstack([YX_chunk, nan_chunk])])
            new_data_YY = np.vstack([new_data_YY, np.vstack([YY_chunk, nan_chunk])])

        # Mask out NaN values
        XX = np.ma.masked_invalid(new_data_XX[1:])
        XY = np.ma.masked_invalid(new_data_XY[1:])
        YX = np.ma.masked_invalid(new_data_YX[1:])
        YY = np.ma.masked_invalid(new_data_YY[1:])

        self.XX = np.ma.masked_where(XX == 0, XX)
        self.XY = np.ma.masked_where(XY == 0, XY)
        self.YX = np.ma.masked_where(YX == 0, YX)
        self.YY = np.ma.masked_where(YY == 0, YY)

    def _make_stokes(self):
        """Convert instrumental polarisations to Stokes products.

        NOTE:
        For some inexplicable reason the I, Q, and U products result in
        an array full of NaNs if both calibrator scans are added and
        the array averaged/reshaped. This -1j*1j hack seems to fix it.
        Object and data types are <class 'numpy.ma.core.MaskedArray'>
        and complex128 in both cases
        """

        I = -1j * 1j * (self.XX + self.YY) / 2
        Q = 1j * 1j * (self.XY + self.YX) / 2
        U = 1j * 1j * (self.XX - self.YY) / 2
        V = 1j * (self.XY - self.YX) / 2

        # Flip frequency axis of data to intuitive order
        I = np.flip(I, axis=1)
        Q = np.flip(Q, axis=1)
        U = np.flip(U, axis=1)
        V = np.flip(V, axis=1)
        self.data = {"I": I, "Q": Q, "U": U, "V": V}

    def rebin(self, o, n, axis):
        """Create unitary array compression matrix from o -> n length.

        if rebinning along row axis we want:
         - (o // n) + 1 entries in each row that sum to unity,
         - each column to sum to the compression ratio o / n
         - values distributed along the row in units of o / n
           until expired

           >>> rebin(5, 3)
           array([[0.6, 0.4, 0. , 0. , 0. ],
                  [0. , 0.2, 0.6, 0.2, 0. ],
                  [0. , 0. , 0. , 0.4, 0.6]])

         - transpose of this for column rebinning

        The inner product of this compressor with an array will rebin
        the array conserving the total intensity along the given axis.
        """

        compressor = np.zeros((n, o))
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
                dnrow = 0
                dncol = 1
            # Use remaining budget if at end of current row
            elif budget < comp_ratio:
                value = budget
                overflow = comp_ratio - budget
                budget = 1
                dnrow = 1
                dncol = 0
            # Otherwise spend n / o and move to next column
            else:
                value = comp_ratio
                budget -= value
                dnrow = 0
                dncol = 1

            compressor[nrow, ncol] = value
            nrow += dnrow
            ncol += dncol

        return compressor if axis == 0 else compressor.T

    def rebin2D(self, array, new_shape):
        """Re-bin along time / frequency axes conserving flux."""

        if new_shape[0] > array.shape[0] or new_shape[1] > array.shape[1]:
            raise ValueError("New shape should not be greater than old shape in either dimension")

        time_comp = self.rebin(array.shape[0], new_shape[0], axis=0)
        freq_comp = self.rebin(array.shape[1], new_shape[1], axis=1)
        result = time_comp @ np.array(array) @ freq_comp

        return result

    def plot_ds(self, favg, tavg, stokes, cmax=20, save=False, imag=False):

        fbins = self.data["I"].shape[1] // favg
        tbins = self.data["I"].shape[0] // tavg

        data = self.data[stokes].imag if imag else self.data[stokes].real
        data = self.rebin2D(data, (tbins, fbins))

        norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
        cmin = -2 if stokes == "I" else -cmax
        cmax = cmax

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(
            np.transpose(data),
            extent=[self.tmin, self.tmax, self.fmax, self.fmin],
            aspect="auto",
            origin="lower",
            norm=norm,
            clim=(cmin, cmax),
            cmap="plasma",
        )
        ax.text(
            0.05, 0.95, f"Stokes {stokes}", color="white", weight="bold", transform=ax.transAxes
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label("Flux Density (mJy)")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Frequency (MHz)")

        if save:
            path_template = "{}/{}_stokes{}_subbed_ds_fbins{}-tbins{}.png"
            fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), stokes.lower(), fbins, tbins),
                bbox_inches="tight",
                format="png",
                dpi=300,
            )

        return fig, ax

    def plot_spectrum(self, favg, save=False):

        sp_fig, sp_ax = plt.subplots(figsize=(7, 5))

        fbins = self.data["I"].shape[1] // favg
        df = (self.fmax - self.fmin) / fbins
        x = [(i * df + self.fmin) for i in range(fbins)]

        for pol in ["I", "Q", "U", "V"]:
            y = self.rebin2D(self.data[pol], (1, fbins))
            y = np.flip(np.transpose(y))

            sp_ax.plot(x, y.real, label=pol)

            if pol == "I":
                rms = np.sqrt(np.mean(np.square(y.imag)))

        sp_ax.set_ylabel("Flux Density (mJy)")
        sp_ax.set_xlabel("Frequency (MHz)")
        sp_ax.axhline(rms, ls=":", color="r", label=f"rms={rms:.1f} mJy")
        sp_ax.axhline(-rms, ls=":", color="r")
        sp_ax.legend()

        if save:
            path_template = "{}/{}_spectrum_fbins{}.png"
            sp_fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), fbins),
                bbox_inches="tight",
                format="png",
                dpi=300,
            )

        return sp_fig, sp_ax

    def plot_lightcurve(self, tavg, save=False):
        lc_fig, lc_ax = plt.subplots(figsize=(7, 5))

        tbins = self.data["I"].shape[0] // tavg
        dt = (self.tmax - self.tmin) / tbins
        x = np.array([i * dt for i in range(tbins)])

        for pol in ["I", "Q", "U", "V"]:
            y = self.rebin2D(self.data[pol], (tbins, 1))

            lc_ax.plot(x, y.real, label=pol)

            df = pd.DataFrame({"time": x, "flux_density": y.real.reshape(1, -1)[0]})

            df["time"] = Time("2021-09-16T08:31:34.927736") + x * u.hour

            if pol == "I":
                rms = np.sqrt(np.mean(np.square(y.imag)))

        lc_ax.axhline(rms, ls=":", color="r", label=f"rms={rms:.1f} mJy")
        lc_ax.axhline(-rms, ls=":", color="r")
        lc_ax.legend()
        lc_ax.set_ylabel("Flux Density (mJy)")
        lc_ax.set_xlabel("Time (hours)")

        if save:
            path_template = "{}/{}_lc_tbins{}.png"
            lc_fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), tbins),
                bbox_inches="tight",
                format="png",
                dpi=300,
            )

        return lc_fig, lc_ax
