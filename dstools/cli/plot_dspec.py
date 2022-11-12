import click
import logging
import matplotlib.pyplot as plt
from astroutils.logger import setupLogger
from dstools.dynamic_spectrum import DynamicSpectrum

logger = logging.getLogger(__name__)

@click.command()
@click.option('-f', '--favg', default=1, type=int)
@click.option('-t', '--tavg', default=1, type=int)
@click.option('-I', '--cmax_i', default=50, type=float)
@click.option('-Q', '--cmax_qu', default=10, type=float)
@click.option('-V', '--cmax_v', default=50, type=float)
@click.option('-r', '--real', is_flag=True, default=True)
@click.option('-i', '--imag', is_flag=True, default=False)
@click.option('-s', '--stokes', default='IQUV')
@click.option('-x', '--crosspols', is_flag=True, default=False)
@click.option('-S', '--save', is_flag=True, default=False)
@click.option('-a', '--acf', is_flag=True, default=False)
@click.option('-l', '--lightcurve', is_flag=True, default=False)
@click.option('-p', '--spectrum', is_flag=True, default=False)
@click.option('-C', '--calscans', is_flag=True, default=True)
@click.option('-B', '--band', default='L', type=click.Choice(['low', 'mid', 'L', 'C', 'X']))
@click.option('-N', '--versionname', default=None, help='Prefix for different processing versions')
@click.option('-v', '--verbose', is_flag=True, default=False)
@click.argument('project')
def main(
    favg,
    tavg,
    cmax_i,
    cmax_qu,
    cmax_v,
    real,
    imag,
    stokes,
    crosspols,
    save,
    acf,
    lightcurve,
    spectrum,
    calscans,
    band,
    versionname,
    verbose,
    project,
):

    setupLogger(verbose)

    prefix = f'{versionname}_' if versionname else ''

    cmax = {
        'I': cmax_i,
        'Q': cmax_qu,
        'U': cmax_qu,
        'V': cmax_v,
    }

    ds = DynamicSpectrum(project, band=band, calscans=calscans, tavg=tavg, favg=favg, prefix=prefix)

    # Dynamic Spectrum
    # --------------------------------------
    for s in stokes:
        if real:
            fig, ax = ds.plot_ds(stokes=s, cmax=cmax[s], save=save)
        if imag:
            fig, ax = ds.plot_ds(stokes=s, cmax=cmax[s], save=save, imag=True)

    if crosspols:
        fig, ax = ds.plot_crosspol_ds(cmax=cmax['I'], save=save)

    # Spectrum
    # --------------------------------------
    if spectrum:
        fig, ax = ds.plot_spectrum(save=save)

    # Light Curve
    # --------------------------------------
    if lightcurve:
        fig, ax = ds.plot_lightcurve(save=save)

    # Dynamic Spectrum 2D Auto-correlation Function
    # --------------------------------------
    if acf:
        fig, ax = ds.plot_acf(save=save)

    plt.show()


if __name__ == '__main__':
    main()
