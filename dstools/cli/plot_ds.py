import logging
import sys
import warnings
from itertools import chain, combinations

import astropy.units as u
import click
import matplotlib.pyplot as plt
import numpy as np
from astroutils.logger import setupLogger
from erfa import ErfaWarning
from matplotlib.gridspec import GridSpec

from dstools.dynamic_spectrum import DynamicSpectrum
from dstools.utils import BANDS

warnings.filterwarnings("ignore", category=ErfaWarning, append=True)

logger = logging.getLogger(__name__)

stokes_choices = [
    "".join(stokes)
    for stokes in chain(
        *(list(combinations(["I", "Q", "U", "V", "L"], i + 1)) for i in range(5))
    )
]


def make_summary_plot(ds, stokes, cmax, imag):

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

    fig, sp_ax = ds.plot_spectrum(stokes=stokes, fig=fig, ax=sp_ax)
    fig, lc_ax = ds.plot_lightcurve(stokes=stokes, fig=fig, ax=lc_ax, polangle=False)

    fig.subplots_adjust(
        left=0.06,
        top=0.94,
        right=0.96,
        bottom=0.05,
    )

    return fig


@click.command()
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
    default=30,
    type=float,
    help="Maximum colormap normalisation in Stokes I.",
)
@click.option(
    "-L",
    "--cmax_l",
    default=50,
    type=float,
    help="Maximum colormap normalisation in linear polarisations.",
)
@click.option(
    "-V",
    "--cmax_v",
    default=10,
    type=float,
    help="Maximum colormap normalisation in Stokes V.",
)
@click.option(
    "-r",
    "--real",
    is_flag=True,
    default=True,
    help="Toggle plotting of real component of visibilities in dynamic spectra.",
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
    "-P",
    "--linpols",
    is_flag=True,
    default=False,
    help="Plot linear polarisation fraction and angle dynamic spectra.",
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
    "-x",
    "--polangle",
    is_flag=True,
    default=False,
    help="Include polarisation angle in lightcurve plot.",
)
@click.option(
    "-R",
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
    "-F",
    "--fold",
    is_flag=True,
    default=False,
    help="Toggle to enable folding of data. Must also provide period with -T.",
)
@click.option(
    "-E",
    "--derotate",
    is_flag=True,
    default=False,
    help="Toggle RM de-rotation of Stokes Q/U.",
)
@click.option(
    "-T",
    "--period",
    default=None,
    type=float,
    help="Period to use when folding data.",
)
@click.option(
    "-o",
    "--period_offset",
    default=0,
    type=float,
    help="Period phase offset to use when folding data.",
)
@click.option(
    "-C",
    "--calscans",
    is_flag=True,
    default=True,
    help="Toggle inclusion of null-valued time chunks while off-source (e.g. calibrator scans, wind stows)",
)
@click.option(
    "-B",
    "--band",
    default="AT_L",
    type=click.Choice(BANDS),
    help="Frequency band. Must correspond to a sub-directory of <project>/dynamic_spectra/",
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
    real,
    imag,
    linpols,
    stokes,
    dspec,
    lightcurve,
    spectrum,
    polangle,
    fdf,
    acf,
    fold,
    derotate,
    trim,
    period,
    period_offset,
    calscans,
    band,
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
        band=band,
        calscans=calscans,
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
        fold=fold,
        trim=trim,
        derotate=derotate,
        period=period,
        period_offset=period_offset,
    )

    # Dynamic Spectrum
    # --------------------------------------
    if dspec:
        for s in stokes:
            if real:
                ds.plot_ds(stokes=s, cmax=cmax[s])
            if imag:
                ds.plot_ds(stokes=s, cmax=cmax[s], imag=True)

    if linpols:
        ds.plot_pol_ds()
        ds.plot_polangle_ds()

    # Spectrum
    # --------------------------------------
    if spectrum:
        ds.plot_spectrum(stokes=stokes)

    # Light Curve
    # --------------------------------------
    if lightcurve:
        ds.plot_lightcurve(stokes=stokes, polangle=polangle)

    # Summary plot
    # --------------------------------------
    if summary:
        make_summary_plot(
            ds,
            stokes,
            cmax,
            imag,
        )

    # Dynamic Spectrum 2D Auto-correlation Function
    # --------------------------------------
    if acf:
        for s in stokes:
            ds.plot_acf(stokes=s, contrast=0.2)

    # RM FDF and Lightcurve/Dynamic Spectrum of Polarisation Angle
    # --------------------------------------
    if fdf:
        ds.plot_fdf()

    plt.show()


if __name__ == "__main__":
    main()