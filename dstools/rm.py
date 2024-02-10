"""RM module that implements 1-D RM Synthesis/RMCLEAN."""

"""
MIT License

Copyright (c) 2017 George Heald

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# Written by George Heald
# v1.0, 14 November 2017

from numpy import *
from pylab import *
from scipy.optimize import curve_fit


class PolObservation:
    """Class to describe an observation & perform polarimetry operations"""

    def __init__(self, freq, IQU, IQUerr=None, verbose=True):
        """
        Initialise a PolObservation object.

        The object needs to know the frequency setup and
        obtain Stokes IQU data.
        If you don't have I, but only QU, just provide ones.
        Optionally, errors can be provided for I, Q, and U.

        Parameters
        ----------
        freq : array
            Frequencies in Hz
        IQU : array
            Stokes IQU values in a 3xN array where
            N is the length of freq.
        IQUerr : array, optional (default None)
            Stokes IQU uncertainties in a 3xN array
        verbose : boolean, optional (default True)
            Print some output?

        When initialised, a powerlaw will be fit to Stokes I
        for possible later use.
        """

        self.freq = freq
        self.i = IQU[0]
        self.q = IQU[1]
        self.u = IQU[2]
        if IQUerr is None:
            self.ierr = None
            self.qerr = None
            self.uerr = None
            p, pcov = curve_fit(
                lambda f, s, a: s * (f / mean(f)) ** a,
                freq,
                self.i,
                p0=(mean(self.i), -0.7),
            )
        else:
            self.ierr = IQUerr[0]
            self.qerr = IQUerr[1]
            self.uerr = IQUerr[2]
            p, pcov = curve_fit(
                lambda f, s, a: s * (f / mean(f)) ** a,
                freq,
                self.i,
                p0=(mean(self.i), -0.7),
                sigma=self.ierr,
            )
        self.i_model = p[0] * (freq / mean(freq)) ** p[1]
        if verbose:
            print("Fitted spectral index is %f" % p[1])
        self.model_s = p[0]
        self.model_a = p[1]
        self.rmsynth_done = False
        self.rmclean_done = False

    def rmsynthesis(
        self,
        phi,
        norm_mod=False,
        norm_vals=False,
        double=True,
        clip=None,
        pclip=None,
        weightmode="none",
        verbose=True,
    ):
        """
        Perform RM Synthesis.

        This function performs RM Synthesis on the provided
        IQU data and computes the RMSF.
        It can adopt weights based on the IQU errors if they
        are provided.

        Parameters
        ----------
        phi : array
           RM values to use (should be contiguous and regularly spaced)
        norm_mod : boolean, optional (default False)
           Normalise QU values with the Stokes I fitted power-law model?
        norm_vals : boolean, optional (default False)
           Normalise QU values with the Stokes I data points
           per channel?
        double : boolean, optional (default True)
           Create the RMSF double the length along the RM axis?
           This should be kept as True unless you are not planning
           to RMCLEAN.
        clip : float, optional (default -inf)
           Ignore channels with Stokes I S/N ratio < clip
        pclip : float, optional (default -inf)
           Ignore channels with sqrt(Q**2+U**2) S/N ratio < pclip
        weightmode : str, optional (default 'none')
           How to do weighting in the Fourier transform. Options are:
             'none'	= no weighting (uniform weights)
             'varwt'	= inverse variance weights based on
                          Stokes I noise
           Further options may be added later.
        verbose : boolean, optional (default True)
            Print some output?

        """

        # Constants
        c = 2.998e8
        ci = complex(0.0, 1.0)

        # Compute lamba-squared values
        l2 = (c / self.freq) ** 2

        # Clip values if needed
        if clip is None:
            clip = -inf
        if pclip is None:
            pclip = -inf
        if self.ierr is None:
            gp = logical_and(self.i > clip, sqrt(self.q**2 + self.u**2) > pclip)
        else:
            gp = logical_and(
                self.i / self.ierr > clip,
                sqrt((self.q / self.qerr) ** 2 + (self.u / self.uerr) ** 2) > pclip,
            )
        if verbose:
            print("Using %d of %d channels" % (sum(gp), len(l2)))

        # Work out weights if needed and overall normalisation
        if weightmode == "varwt":
            if self.ierr is None:
                print(
                    "Error: you need to provide Stokes I errors to use inverse variance weighting"
                )
                return
            if verbose:
                print("Inverse variance weighting mode activated")
            weights = 1.0 / self.ierr**2
        else:
            if verbose:
                print("Using uniform weights")
            weights = ones(self.freq.shape)
        K = 1.0 / sum(weights[gp])  # Normalisation
        if verbose:
            print("Normalisation value = %f" % (K))

        # Calculate and report FWHM and other metrics
        min_l2 = min(l2[gp])
        min_dl2 = min(abs(diff(l2[gp])))
        fwhm = 2.0 * sqrt(3.0) / (max(l2[gp]) - min(l2[gp]))
        maxscale = pi / min_l2
        rm_max = sqrt(3.0) / min_dl2
        if verbose:
            print("FWHM of RMSF is %f rad/m2" % fwhm)
            print("Max RM scale is %f rad/m2" % maxscale)
            print("Max RM value is %f rad/m2" % rm_max)
        l20 = mean(l2[gp])

        # QU values to operate on, normalized if needed
        quvec = self.q + ci * self.u
        if norm_mod:
            if verbose:
                print("Normalising QU by Stokes I model")
            quvec /= self.i_model
        if norm_vals:
            if verbose:
                print("Normalising QU by Stokes I values")
            quvec /= self.i

        # Work out dimensions of RMSF
        dphi = min(diff(phi))
        if double:
            nrmsf = len(phi) * 2 + 1
        else:
            nrmsf = len(phi)
        rmsf_phi = linspace(-float(nrmsf / 2) * dphi, float(nrmsf / 2) * dphi, nrmsf)
        rmsf = zeros(rmsf_phi.shape, dtype=complex)
        fdf = zeros(phi.shape, dtype=complex)

        # Now do the work
        if verbose:
            print("Calculating RMSF...")
        for j, phival in enumerate(rmsf_phi):
            rmsf[j] = K * sum(weights[gp] * exp(-2.0 * ci * phival * (l2[gp] - l20)))
        if verbose:
            print("Calculating FDF...")
        for j, phival in enumerate(phi):
            fdf[j] = K * sum(
                quvec[gp] * weights[gp] * exp(-2.0 * ci * phival * (l2[gp] - l20))
            )

        # For later use
        self.l20 = l20
        self.weights = weights
        self.K = K
        self.rmsf = rmsf
        self.rmsf_phi = rmsf_phi
        self.fdf = fdf
        self.rmsf_fwhm = fwhm
        self.maxscale = maxscale
        self.rm_max = rm_max
        self.dphi = dphi
        self.rmsf_fwhm_pix = int(fwhm / dphi + 0.5)
        self.phi = phi
        self.nphi = len(phi)
        self.rmsynth_done = True
        self.rmclean_done = False

    def plot_fdf(self, display=True, save=None, rescale=False, plot_rmsf=True):
        """
        Plot the FDF and RMSF

        The plot will show the RMSF in black, and the FDF in red.
        If RMCLEAN has already been performed then the cleaned spectrum
        will be shown in blue, and the clean model in green.

        Parameters
        ----------
        display : boolean, optional (default True)
               Show plot on screen?
        save : str, optional (default None)
               Save figure to disk? Provide filename if desired.
        rescale : boolean, optional (default False)
               Rescale RMSF peak to match that of FDF?
        plot_rmsf : boolean, optional (default True)
               Plot RMSF?

        """

        if self.rmsynth_done:
            figure()
            if rescale:
                scfac = max(abs(self.fdf))
            else:
                scfac = 1.0
            if plot_rmsf:
                plot(self.rmsf_phi, scfac * abs(self.rmsf), "k-")
                plot(self.phi, abs(self.fdf), "r-")
                if self.rmclean_done:
                    plot(self.phi, abs(self.rm_cleaned), "b-")
                    plot(self.phi, abs(self.rm_comps), "g-")
                    legend(("RMSF", "FDF", "Clean", "Model"), loc="best")
                else:
                    legend(("RMSF", "FDF"), loc="best")
            else:
                plot(self.phi, abs(self.fdf), "r-")
                if self.rmclean_done:
                    plot(self.phi, abs(self.rm_cleaned), "b-")
                    plot(self.phi, abs(self.rm_comps), "g-")
                    legend(("FDF", "Clean", "Model"), loc="best")
                else:
                    legend(("FDF"), loc="best")
            xlabel("RM (rad/m2)")
            ylabel("Amplitude")
            if save is not None:
                savefig(save, bbox_inches="tight")
            if display:
                show()
        # close('all')

    def plot_stokesi(self, display=True, save=None):
        """
        Plot Stokes I

        The plot will show the Stokes I values (and errors if available).

        Parameters
        ----------
        display : boolean, optional (default True)
               Show plot on screen?
        save : str, optional (default None)
               Save figure to disk? Provide filename if desired.

        """

        figure()
        if self.ierr is None:
            plot(self.freq, self.i, marker="o", ls="none")
        else:
            errorbar(self.freq, self.i, yerr=self.ierr, marker="o", ls="none")
        plot(self.freq, self.i_model, "k--")
        xlabel("Frequency (Hz)")
        ylabel("Stokes I")
        if save is not None:
            savefig(save, bbox_inches="tight")
        if display:
            show()
        # close('all')

    def get_fdf_peak(self, verbose=True):
        """
        Obtain the peak of the FDF and its Faraday depth

        This will fit the peak of the FSF and report results
        if requested.
        Values that are calculated and reported are:
                Absolute value at peak, PA at peak, peak RM value
        This is done for the dirty FDF.
        If RMCLEAN has been performed then this is also done for the
        deconvolved Faraday spectrum.

        Parameters
        ----------
        verbose : boolean, optional (default True)
            Print some output?

        """
        if self.rmsynth_done:
            x0 = where(abs(self.fdf) == max(abs(self.fdf)))[0][0]
            if x0 == 0 or x0 == len(self.phi) - 1:
                self.fdf_peak_rm = self.phi[x0]
                self.fdf_peak = max(abs(self.fdf))
                self.fdf_peak_err = 0.0
            else:
                y = abs(self.fdf[x0 - 1 : x0 + 2])
                x = self.phi[x0 - 1 : x0 + 2]
                p = polyfit(x, y, 2)
                self.fdf_peak_rm = -p[1] / (2.0 * p[0])
                self.fdf_peak = polyval(p, self.fdf_peak_rm)
                self.fdf_peak_rm_err = self.rmsf_fwhm / (
                    2.355 * self.fdf_peak / std(real(self.fdf))
                )
                rotpeak = self.fdf[x0] * exp(
                    -2.0 * complex(0.0, 1.0) * self.fdf_peak_rm * self.l20
                )
                pa = 0.5 * arctan2(imag(rotpeak), real(rotpeak))
                self.pa = (pa * 180.0 / pi) % 360.0
            if verbose:
                print(
                    "FDF peaks at an amplitude %f, at RM=%f +/- %f rad/m2"
                    % (self.fdf_peak, self.fdf_peak_rm, self.fdf_peak_rm_err)
                )
            if verbose:
                print("RM-corrected PA of FDF peak is %f degrees" % (self.pa))
        if self.rmclean_done:
            x0 = where(abs(self.rm_cleaned) == max(abs(self.rm_cleaned)))[0][0]
            if x0 == 0 or x0 == len(self.phi) - 1:
                self.cln_fdf_peak_rm = self.phi[x0]
                self.cln_fdf_peak = max(abs(self.fdf))
            else:
                y = abs(self.rm_cleaned[x0 - 1 : x0 + 2])
                x = self.phi[x0 - 1 : x0 + 2]
                p = polyfit(x, y, 2)
                self.cln_fdf_peak_rm = -p[1] / (2.0 * p[0])
                self.cln_fdf_peak = polyval(p, self.fdf_peak_rm)
                self.cln_fdf_peak_rm_err = self.rmsf_fwhm / (
                    2.355 * self.fdf_peak / std(real(self.rm_resid))
                )
                cln_rotpeak = self.rm_cleaned[x0] * exp(
                    -2.0 * complex(0.0, 1.0) * self.fdf_peak_rm * self.l20
                )
                cln_pa = 0.5 * arctan2(imag(cln_rotpeak), real(cln_rotpeak))
                self.cln_pa = (cln_pa * 180.0 / pi) % 360.0
            if verbose:
                print(
                    "Cleaned FDF peaks at an amplitude %f, at RM=%f +/- %f rad/m2"
                    % (
                        self.cln_fdf_peak,
                        self.cln_fdf_peak_rm,
                        self.cln_fdf_peak_rm_err,
                    )
                )
            if verbose:
                print(
                    "RM-corrected PA of cleaned FDF peak is %f degrees" % (self.cln_pa)
                )

    def rmclean(self, niter=1000, gain=0.1, cutoff=2.0, mask=False, verbose=True):
        """
        Perform RMCLEAN

        The FDF will be deconvolved using the RMSF.

        Parameters
        ----------
        niter : int, optional (default 1000)
             Maximum number of clean iterations
        gain : float, optional (default 0.1)
             Clean gain
        cutoff : float, optional (default 2)
             Clean cutoff, in units of S/N
             The default stops at 2*sigma above the mean
        mask : boolean (default False)
             If True, all clean components must be within an
             RMSF FWHM of the first peak
        verbose : boolean, optional (default True)
            Print some output?

        """

        noise = std(real(self.fdf))
        zerolev = median(abs(self.fdf))
        cleanlim = cutoff * noise + zerolev
        if verbose:
            print("CLEAN will proceed down to %f" % (cleanlim))
        num = 0
        res = self.fdf.copy()
        modcomp = zeros(res.shape, dtype=complex)
        resp = abs(res)
        mr = range(len(resp))
        while max(resp[mr]) > cleanlim and num < niter:
            maxloc = where(resp[mr] == max(resp[mr]))[0] + mr[0]
            if num == 0 and verbose:
                print("First component found at %f rad/m2" % (self.phi[maxloc]))
            if num == 0 and mask:
                mr = range(
                    maxloc - self.rmsf_fwhm_pix / 2, maxloc + self.rmsf_fwhm_pix / 2 + 1
                )
                if verbose:
                    print(
                        "Masking: Clean components must fall within mask of %d/%d pixels"
                        % (len(mr), len(resp))
                    )
                    print(
                        "(i.e. within RM range %f - %f rad/m2)"
                        % (self.phi[mr[0]], self.phi[mr[-1]])
                    )
            num += 1
            if num % 10 ** int(log10(num)) == 0 and verbose:
                print("Iteration %d: max residual = %f" % (num, max(resp)))
            srmsf = roll(self.rmsf, maxloc - self.nphi)
            modcomp[maxloc] += res[maxloc] * gain
            subtr = res[maxloc] * gain * srmsf[: self.nphi]
            res -= subtr
            resp = abs(res)
        if verbose:
            print("Convolving clean components...")
        if 10 * self.rmsf_fwhm_pix > len(self.phi):
            kernel = exp(
                -((self.phi - mean(self.phi)) ** 2)
                / (2.0 * (self.rmsf_fwhm / 2.355) ** 2)
            )
        else:
            kernel = exp(
                -arange(-self.rmsf_fwhm * 5.0, self.rmsf_fwhm * 5.0, self.dphi) ** 2
                / (2.0 * (self.rmsf_fwhm / 2.355) ** 2)
            )
        self.rm_model = convolve(modcomp, kernel, mode="same")
        if verbose:
            print("Restoring convolved clean components...")
        self.rm_cleaned = self.rm_model + res
        self.rm_comps = modcomp
        self.niters = num
        self.rm_resid = res
        self.rmclean_done = True

    def print_rmstats(self):
        """
        Print some stats about the results

        Some basic statistics will be reported to the terminal:
                mean of RM clean components,
                dispersion of RM clean components

        """

        if self.rmclean_done:
            cp = where(abs(self.rm_comps) > 0.0)
            # cabs = abs(self.rm_comps[cp])
            # crm = self.phi[cp]
            cabs = abs(self.rm_comps)
            crm = self.phi.copy()
            mcrm = sum(cabs * crm) / sum(cabs)
            crmdisp = sqrt(sum(cabs * (crm - mcrm) ** 2) / sum(cabs))
            print("Weighted mean RM of clean components is %f" % mcrm)
            print("Weighted RM dispersion of clean components is %f" % crmdisp)
