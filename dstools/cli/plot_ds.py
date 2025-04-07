import logging
import warnings
from itertools import chain, combinations

import astropy.units as u
import click
import matplotlib.pyplot as plt
import numpy as np
from erfa import ErfaWarning

from dstools.dynamic_spectrum import DynamicSpectrum, LightCurve, Spectrum
from dstools.logger import setupLogger
from dstools.plotting import (
    _plot_polarisations,
    plot_acf,
    plot_ds,
    plot_lightcurve,
    plot_polarisation_lightcurve,
    plot_polarisation_spectrum,
    plot_spectrum,
    plot_summary,
)

warnings.filterwarnings("ignore", category=ErfaWarning, append=True)
warnings.filterwarnings("ignore", category=UserWarning, append=True)

logger = logging.getLogger(__name__)

stokes_choices = [
    "".join(stokes)
    for stokes in chain(
        *(list(combinations(["I", "Q", "U", "V", "L", "P"], i)) for i in range(1, 6))
    )
]


@click.command(context_settings={"show_default": True})
@click.option(
    "-f",
    "--favg",
    default=1,
    type=int,
    help="Averaging factor across frequency axis.",
)
@click.option(
    "-t",
    "--tavg",
    default=1,
    type=int,
    help="Averaging factor across time axis.",
)
@click.option(
    "--fmin",
    default=None,
    type=float,
    help="Selection of minimum frequency in MHz.",
)
@click.option(
    "--fmax",
    default=None,
    type=float,
    help="Selection of maximum frequency in MHz.",
)
@click.option(
    "--uvmin",
    default=0,
    type=float,
    help="Selection of minimum projected baseline distance in m.",
)
@click.option(
    "--uvmax",
    default=np.inf,
    type=float,
    help="Selection of maximum projected baseline distance in m.",
)
@click.option(
    "--uvwavemin",
    default=0,
    type=float,
    help="Selection of minimum frequency-dependent projected baseline distance in units of wavelengths.",
)
@click.option(
    "--uvwavemax",
    default=np.inf,
    type=float,
    help="Selection of maximum frequency-dependent projected baseline distance in units of wavelengths.",
)
@click.option(
    "--tmin",
    default=None,
    type=float,
    help="Selection of minimum time in hours.",
)
@click.option(
    "--tmax",
    default=None,
    type=float,
    help="Selection of maximum time in hours.",
)
@click.option(
    "-u",
    "--tunit",
    default="hour",
    type=click.Choice(["h", "hour", "min", "minute", "s", "second"]),
    help="Selection of time axis unit.",
)
@click.option(
    "-I",
    "--cmax_i",
    default=15,
    type=float,
    help="Maximum colormap normalisation in Stokes I.",
)
@click.option(
    "-L",
    "--cmax_l",
    default=15,
    type=float,
    help="Maximum colormap normalisation in linear polarisations.",
)
@click.option(
    "-V",
    "--cmax_v",
    default=15,
    type=float,
    help="Maximum colormap normalisation in Stokes V.",
)
@click.option(
    "-i",
    "--imag",
    is_flag=True,
    default=False,
    help="Toggle plotting of imaginary component of visibilities in dynamic spectra.",
)
@click.option(
    "-s",
    "--stokes",
    default="IQUV",
    type=click.Choice(stokes_choices),
    help="Stokes parameters that will be included in each plot.",
)
@click.option(
    "-d",
    "--dspec",
    is_flag=True,
    default=False,
    help="Plot dynamic spectrum.",
)
@click.option(
    "-l",
    "--lightcurve",
    is_flag=True,
    default=False,
    help="Plot channel-averaged lightcurve.",
)
@click.option(
    "-p",
    "--spectrum",
    is_flag=True,
    default=False,
    help="Plot time-averaged spectrum.",
)
@click.option(
    "-P",
    "--polarisations",
    is_flag=True,
    default=False,
    help="Include polarisation fraction / angle / ellipticity in lightcurve / spectrum plot.",
)
@click.option(
    "--fdf",
    is_flag=True,
    default=False,
    help="Plot RM synthesis Faraday dispersion function.",
)
@click.option(
    "-a",
    "--acf",
    is_flag=True,
    default=False,
    help="Plot 2D auto-correlation function.",
)
@click.option(
    "-k",
    "--trim",
    is_flag=True,
    default=True,
    help="Remove flagged channels at top/bottom of band.",
)
@click.option(
    "--fold",
    is_flag=True,
    default=False,
    help="Toggle to enable folding of data.",
)
@click.option(
    "-R",
    "--derotate",
    is_flag=True,
    default=False,
    help="Toggle RM de-rotation of Stokes Q/U.",
)
@click.option(
    "--RM",
    type=float,
    default=None,
    help="Rotation measure in units of rad/m^2.",
)
@click.option(
    "-D",
    "--dedisperse",
    is_flag=True,
    default=False,
    help="Toggle de-dispersion.",
)
@click.option(
    "--DM",
    type=float,
    default=None,
    help="Dispersion measure in units of pc/cm^3.",
)
@click.option(
    "-T",
    "--period",
    default=None,
    type=float,
    help="Period to use when folding.",
)
@click.option(
    "-o",
    "--period_offset",
    default=0,
    type=float,
    help="Period phase offset to use when folding.",
)
@click.option(
    "-B",
    "--barycentre",
    is_flag=True,
    default=False,
    help="Toggle Barycentric correction.",
)
@click.option(
    "--absolute-times",
    is_flag=True,
    default=True,
    help="Toggle plotting time axes with absolute vs relative times.",
)
@click.option(
    "-C",
    "--calscans",
    is_flag=True,
    default=True,
    help="Toggle inclusion of null-valued time chunks while off-source.",
)
@click.option(
    "-Y",
    "--summary",
    is_flag=True,
    default=False,
    help="Plot all Stokes dynamic spectrum, lightcurve, and time-averaged spectrum.",
)
@click.option("-v", "--verbose", is_flag=True, default=False)
@click.argument("ds_path")
def main(
    favg,
    tavg,
    fmin,
    fmax,
    uvmin,
    uvmax,
    uvwavemin,
    uvwavemax,
    tmin,
    tmax,
    tunit,
    cmax_i,
    cmax_l,
    cmax_v,
    imag,
    stokes,
    dspec,
    lightcurve,
    spectrum,
    polarisations,
    fdf,
    rm,
    acf,
    fold,
    barycentre,
    dedisperse,
    dm,
    derotate,
    trim,
    period,
    period_offset,
    absolute_times,
    calscans,
    summary,
    verbose,
    ds_path,
):
    setupLogger(verbose)

    tunit = u.Unit(tunit)

    cmax = {
        "I": cmax_i,
        "Q": cmax_l,
        "U": cmax_l,
        "V": cmax_v,
        "L": cmax_l,
    }

    ds = DynamicSpectrum(
        ds_path=ds_path,
        tavg=tavg,
        favg=favg,
        minfreq=fmin,
        maxfreq=fmax,
        minuvdist=uvmin,
        maxuvdist=uvmax,
        minuvwave=uvwavemin,
        maxuvwave=uvwavemax,
        mintime=tmin,
        maxtime=tmax,
        tunit=tunit,
        trim=trim,
        absolute_times=absolute_times,
        calscans=calscans,
        barycentre=barycentre,
        derotate=derotate,
        dedisperse=dedisperse,
        DM=dm,
        RM=rm,
        fold=fold,
        period=period,
        period_offset=period_offset,
    )

    if verbose:
        logger.debug(f"Dynamic spectrum attributes:\n{ds}")

    # Dynamic Spectrum
    # --------------------------------------
    if dspec:
        for s in stokes:
            plot_ds(ds, stokes=s, cmax=cmax[s], imag=imag)

    # Spectrum
    # --------------------------------------
    if spectrum:
        sp = Spectrum(ds)
        if polarisations:
            plot_polarisation_spectrum(sp, stokes=stokes, error_alpha=0.4)
        else:
            plot_spectrum(sp, stokes=stokes)

    # Light Curve
    # --------------------------------------
    if lightcurve:
        lc = LightCurve(ds)
        if polarisations:
            plot_polarisation_lightcurve(lc, stokes=stokes, error_alpha=0.4)
        else:
            plot_lightcurve(lc, stokes=stokes)

    # Summary plot
    # --------------------------------------
    if summary:
        plot_summary(
            ds,
            stokes,
            cmax,
            imag,
        )

    # Dynamic Spectrum 2D Auto-correlation Function
    # --------------------------------------
    if acf:
        for s in stokes:
            plot_acf(ds, stokes=s, contrast=0.2)

    # RM FDF and Lightcurve/Dynamic Spectrum of Polarisation Angle
    # --------------------------------------
    if fdf:
        plot_fdf(ds)

    with warnings.catch_warnings():
        plt.show()


if __name__ == "__main__":
    main()
