import logging
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.visualization import ZScaleInterval, ImageNormalize
from scipy.signal import correlate

logger = logging.getLogger(__name__)

class DynamicSpectrum:
    def __init__(self, project, band='L', tavg=1, favg=1, calscans=True, prefix=''):
        self.src = project.split('/')[-1]
        self.band = band


        self.ds_path = f'{project}/dynamic_spectra/{band}'
        self.prefix = prefix

        self._load_data()
        self._stack_cal_scans(tavg, favg)

        self._make_stokes()

        self.tmin = self.time[0]
        self.tmax = self.time[-1]
        self.fmin = self.freq[0]
        self.fmax = self.freq[-1]

    def _get_scan_intervals(self):
        '''Find indices of start/end of each calibrator scan cycle'''

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
        '''Load instrumental pols and time/freq data, converting to MHz, s, and mJy'''

        file_prefix = f'{self.ds_path}/{self.prefix}'

        self.freq = np.load(f'{file_prefix}freq.npy', allow_pickle=True) / 1e6
        self.time = np.load(f'{file_prefix}time.npy', allow_pickle=True) / 3600

        self.time_start = Time(self.time[0] / 24.0, format='mjd', scale='utc')
        self.time_start.format = 'iso'

        self.time -= self.time[0]

        self.XX = np.load(f'{file_prefix}dynamic_spectra_XX.npy', allow_pickle=True) * 1e3
        self.XY = np.load(f'{file_prefix}dynamic_spectra_XY.npy', allow_pickle=True) * 1e3
        self.YX = np.load(f'{file_prefix}dynamic_spectra_YX.npy', allow_pickle=True) * 1e3
        self.YY = np.load(f'{file_prefix}dynamic_spectra_YY.npy', allow_pickle=True) * 1e3

    def _stack_cal_scans(self, tavg, favg):
        '''Insert null data representing off-source time'''

        self._get_scan_intervals()

        # Calculate average correlator cycle time from first scan
        avg_scan_dt = np.median(self.dts[: self.scan_end_indices[0]])

        # Calculate number of cycles in each calibrator/stow break
        time_end_break = self.time[self.scan_start_indices[1:]]
        time_start_break = self.time[self.scan_end_indices[:-1]]
        num_break_cycles = np.append((time_end_break - time_start_break), 0) / avg_scan_dt
        num_channels = self.XX.shape[1]

        fbins = num_channels // favg

        # Create initial time-slice to start stacking target and calibrator scans together
        new_data_XX = new_data_XY = new_data_YX = new_data_YY = np.zeros(
            (1, fbins), dtype=complex
        )

        # Check that data is not being truncated by stacking insufficient time chunks together
        # This is a diagnostic of errors in data combination, where duplicated times are
        # dropped in dstools.make_dspec if data is combined incorrectly.
        if self.scan_end_indices[-1] + 1 < self.XX.shape[0]:
            logger.warning("More time samples in data array than time array")
            logger.warning("Check results with -C to see full data")

        for start_index, end_index, num_scans in zip(
            self.scan_start_indices, self.scan_end_indices, num_break_cycles
        ):

            tbins = (end_index - start_index) // tavg

            # Select each contiguous on-target chunk of data
            XX_chunk = self.XX[start_index:end_index, :]
            XY_chunk = self.XY[start_index:end_index, :]
            YX_chunk = self.YX[start_index:end_index, :]
            YY_chunk = self.YY[start_index:end_index, :]
            
            XX_chunk = self.rebin2D(XX_chunk, (tbins, fbins))
            XY_chunk = self.rebin2D(XY_chunk, (tbins, fbins))
            YX_chunk = self.rebin2D(YX_chunk, (tbins, fbins))
            YY_chunk = self.rebin2D(YY_chunk, (tbins, fbins))

            # Make array of complex NaN's for subsequent calibrator / stow gaps
            # and append to each on-target chunk of data
            nan_chunk = np.full((int(num_scans) // tavg, fbins), 0 + 0 * 1j)

            new_data_XX = np.vstack([new_data_XX, np.vstack([XX_chunk, nan_chunk])])
            new_data_XY = np.vstack([new_data_XY, np.vstack([XY_chunk, nan_chunk])])
            new_data_YX = np.vstack([new_data_YX, np.vstack([YX_chunk, nan_chunk])])
            new_data_YY = np.vstack([new_data_YY, np.vstack([YY_chunk, nan_chunk])])

        # Mask out NaN values
        self.XX = np.ma.masked_invalid(new_data_XX[1:])
        self.XY = np.ma.masked_invalid(new_data_XY[1:])
        self.YX = np.ma.masked_invalid(new_data_YX[1:])
        self.YY = np.ma.masked_invalid(new_data_YY[1:])

    def _make_stokes(self):
        '''Convert instrumental polarisations to Stokes products.'''

        I = (self.XX + self.YY) / 2
        Q = - (self.XY + self.YX) / 2
        U = - (self.XX - self.YY) / 2
        V = 1j * (self.XY - self.YX) / 2

        # Flip frequency axis of data to intuitive order
        I = np.flip(I, axis=1)
        Q = np.flip(Q, axis=1)
        U = np.flip(U, axis=1)
        V = np.flip(V, axis=1)

        self.data = {'I': I, 'Q': Q, 'U': U, 'V': V}

    def rebin(self, o, n, axis):
        '''Create unitary array compression matrix from o -> n length.

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
        '''

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
        '''Re-bin along time / frequency axes conserving flux.'''

        if new_shape[0] > array.shape[0] or new_shape[1] > array.shape[1]:
            raise ValueError('New shape should not be greater than old shape in either dimension')

        time_comp = self.rebin(array.shape[0], new_shape[0], axis=0)
        freq_comp = self.rebin(array.shape[1], new_shape[1], axis=1)
        result = time_comp @ np.array(array) @ freq_comp

        return result

    def plot_crosspol_ds(self, cmax=20, save=False):

        data = np.sqrt(np.abs(self.data['U']**2 + self.data['V']**2))

        norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
        cmax = cmax
        cmin = 0
        
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(
            np.transpose(data),
            extent=[self.tmin, self.tmax, self.fmax, self.fmin],
            aspect='auto',
            origin='lower',
            norm=norm,
            clim=(cmin, cmax),
            cmap='plasma',
        )
        ax.text(
            0.05, 0.95, r'sqrt(U^2 + V^2)', color='white', weight='bold', transform=ax.transAxes
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label('Flux Density (mJy)')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Frequency (MHz)')

        if save:
            data.dump(f'{self.ds_path}/{self.src.lower()}_crosspol_ds.npy')

            path_template = '{}/{}_crosspol_subbed_ds_fbins{}-tbins{}.png'
            fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), 'REPL', 'REPL'),
                bbox_inches='tight',
                format='png',
                dpi=300,
            )

        return fig, ax

    def plot_ds(self, stokes, cmax=20, save=False, imag=False):

        data = self.data[stokes].imag if imag else self.data[stokes].real

        norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
        cmin = -2 if stokes == 'I' else -cmax
        cmax = cmax

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(
            np.transpose(data),
            extent=[self.tmin, self.tmax, self.fmax, self.fmin],
            aspect='auto',
            origin='lower',
            norm=norm,
            clim=(cmin, cmax),
            cmap='plasma',
        )
        ax.text(
            0.05, 0.95, f'Stokes {stokes}', color='white', weight='bold', transform=ax.transAxes
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label('Flux Density (mJy)')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Frequency (MHz)')

        if save:
            data.dump(f'{self.ds_path}/{self.src.lower()}_{stokes.lower()}_ds.npy')

            path_template = '{}/{}_stokes{}_subbed_ds_fbins{}-tbins{}.png'
            fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), stokes.lower(), 'REPL', 'REPL'),
                bbox_inches='tight',
                format='png',
                dpi=300,
            )

        return fig, ax

    def plot_spectrum(self, save=False):

        sp_fig, sp_ax = plt.subplots(figsize=(7, 5))

        favg = 1
        fbins = self.data['I'].shape[1] // favg
        df = (self.fmax - self.fmin) / fbins
        x = [(i * df + self.fmin) for i in range(fbins)]

        for pol in ['I', 'Q', 'U', 'V']:
            # y = self.rebin2D(self.data[pol], (1, fbins))
            y = np.flip(np.transpose(y))

            sp_ax.plot(x, y.real, label=pol)

            if pol == 'I':
                rms = np.sqrt(np.mean(np.square(y.imag)))

        sp_ax.set_ylabel('Flux Density (mJy)')
        sp_ax.set_xlabel('Frequency (MHz)')
        sp_ax.axhline(rms, ls=':', color='r', label=f'rms={rms:.1f} mJy')
        sp_ax.axhline(-rms, ls=':', color='r')
        sp_ax.legend()

        if save:
            path_template = '{}/{}_spectrum_fbins{}.png'
            sp_fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), fbins),
                bbox_inches='tight',
                format='png',
                dpi=300,
            )

        return sp_fig, sp_ax

    def plot_acf(self, save=False):

        acf_fig, acf_ax = plt.subplots(figsize=(7, 5))

        acf2d = correlate(self.data['V'].real, self.data['V'].real)
        # acf2d = acf2d[acf2d.shape[0]//2:, acf2d.shape[1]//2:]
        # acf2d = np.flip(acf2d, axis=1).T

        norm = ImageNormalize(acf2d, interval=ZScaleInterval(contrast=0.2))

        acf_ax.imshow(
            acf2d,
            extent=[0, self.tmax-self.tmin, 0, self.fmax-self.fmin],
            aspect='auto',
            norm=norm,
            cmap='plasma',
        )
        acf_ax.set_xlabel('Time Lag (h)')
        acf_ax.set_ylabel('Frequency Lag (MHz)')

        if save:
            acf_fig.savefig(f'{self.ds_path}/{self.src.lower()}_acf.png')

        return acf_fig, acf_ax
        
    def plot_lightcurve(self, save=False):
        lc_fig, lc_ax = plt.subplots(figsize=(7, 5))

        tavg = 1
        tbins = self.data['I'].shape[0] // tavg
        dt = (self.tmax - self.tmin) / tbins
        x = np.array([i * dt for i in range(tbins)])

        for pol in ['I', 'Q', 'U', 'V']:
            y = self.data[pol].mean(axis=1)

            lc_ax.plot(x, y.real, label=pol)

            df = pd.DataFrame({
                'time': self.time_start + x * u.hour,
                'flux_density': y.real.reshape(1, -1)[0]
            })
            df.loc[df.flux_density.isin([-1.000000, 1.000000, 0.000000]), 'flux_density'] = np.nan
            df.dropna().to_csv(f'{self.ds_path}/lc_stokes{pol.lower()}.csv')

            if pol == 'I':
                rms = np.sqrt(np.mean(np.square(y.imag)))

        lc_ax.axhline(rms, ls=':', color='r', label=f'rms={rms:.1f} mJy')
        lc_ax.axhline(-rms, ls=':', color='r')
        lc_ax.legend()
        lc_ax.set_ylabel('Flux Density (mJy)')
        lc_ax.set_xlabel('Time (hours)')

        if save:
            path_template = '{}/{}_lc_tbins{}.png'
            lc_fig.savefig(
                path_template.format(self.ds_path, self.src.lower(), tbins),
                bbox_inches='tight',
                format='png',
                dpi=300,
            )

        return lc_fig, lc_ax
