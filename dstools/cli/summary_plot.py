import logging
import os
import warnings

import astropy.units as u
import click
import matplotlib.pyplot as plt
from astropy.utils.exceptions import ErfaWarning
from astroutils.logger import setupLogger
from matplotlib.gridspec import GridSpec

from dstools.dynamic_spectrum import DynamicSpectrum
from dstools.utils import BANDS

warnings.filterwarnings("ignore", category=ErfaWarning, append=True)

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "-f", "--favg", default=1, type=int, help="Averaging factor across frequency axis.."
)
@click.option(
    "-t", "--tavg", default=1, type=int, help="Averaging factor across time axis."
)
@click.option(
    "--fmin", default=None, type=float, help="Selection of minimum frequency in MHz."
)
@click.option(
    "--fmax", default=None, type=float, help="Selection of maximum frequency in MHz."
)
@click.option(
    "--tmin", default=None, type=float, help="Selection of minimum time in hours."
)
@click.option(
    "--tmax", default=None, type=float, help="Selection of maximum time in hours."
)
@click.option(
    "-u",
    "--tunit",
    default="hour",
    type=click.Choice(["h", "hour", "min", "minute", "s", "second"]),
    help="Selection of time axis unit.",
)
@click.option(
    "-s",
    "--stokes",
    default="IQUV",
    help="Stokes parameters that will be included in each plot.",
)
@click.option(
    "-I",
    "--cmax_i",
    default=50,
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
    default=50,
    type=float,
    help="Maximum colormap normalisation in Stokes V.",
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
    "-N", "--versionname", default=None, help="Prefix for different processing versions"
)
@click.option("-D", "--title", default="", type=str)
@click.option("-P", "--show_plot", default=True, is_flag=True)
@click.option("-v", "--verbose", is_flag=True, default=False)
@click.option("-S", "--savepath", type=click.Path(), default=None)
@click.argument("project")
def main(
    favg,
    tavg,
    fmin,
    fmax,
    tmin,
    tmax,
    tunit,
    stokes,
    cmax_i,
    cmax_l,
    cmax_v,
    fold,
    derotate,
    trim,
    period,
    period_offset,
    calscans,
    band,
    versionname,
    title,
    show_plot,
    verbose,
    savepath,
    project,
):

    setupLogger(verbose)

    tunit = u.Unit(tunit)
    prefix = f"{versionname}_" if versionname else ""

    ds = DynamicSpectrum(
        project=project,
        band=band,
        calscans=calscans,
        tavg=tavg,
        favg=favg,
        minfreq=fmin,
        maxfreq=fmax,
        mintime=tmin,
        maxtime=tmax,
        tunit=tunit,
        prefix=prefix,
        fold=fold,
        derotate=derotate,
        trim=trim,
        period=period,
        period_offset=period_offset,
    )

    fig = plt.figure(figsize=(14, 15))
    gs = GridSpec(3, 2, figure=fig)

    I_ax = fig.add_subplot(gs[0, 0])
    Q_ax = fig.add_subplot(gs[0, 1])
    U_ax = fig.add_subplot(gs[1, 0])
    V_ax = fig.add_subplot(gs[1, 1])
    lc_ax = fig.add_subplot(gs[2, 0])
    sp_ax = fig.add_subplot(gs[2, 1])

    fig, I_ax = ds.plot_ds(stokes="I", cmax=cmax_i, fig=fig, ax=I_ax)
    fig, Q_ax = ds.plot_ds(stokes="Q", cmax=cmax_l, fig=fig, ax=Q_ax)
    fig, U_ax = ds.plot_ds(stokes="U", cmax=cmax_l, fig=fig, ax=U_ax)
    fig, V_ax = ds.plot_ds(stokes="V", cmax=cmax_v, fig=fig, ax=V_ax)

    fig, sp_ax = ds.plot_spectrum(stokes=stokes, fig=fig, ax=sp_ax)
    fig, lc_ax = ds.plot_lightcurve(stokes=stokes, fig=fig, ax=lc_ax, polangle=False)

    fig.suptitle(title, size=16, fontweight="bold")

    fig.subplots_adjust(
        left=0.06,
        top=0.94,
        right=0.96,
        bottom=0.05,
    )

    if savepath:
        fig.savefig(savepath)

    if show_plot:
        plt.show()


if __name__ == "__main__":
    main()
