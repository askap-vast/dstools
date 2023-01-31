import click
import logging
import matplotlib.pyplot as plt
from dstools.utils import setupLogger
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
@click.option('-d', '--dspec', is_flag=True, default=False)
@click.option('-l', '--lightcurve', is_flag=True, default=False)
@click.option('-p', '--spectrum', is_flag=True, default=False)
@click.option('-x', '--crosspols', is_flag=True, default=False)
@click.option('-a', '--acf', is_flag=True, default=False)
@click.option('-F', '--fold', is_flag=True, default=False)
@click.option('-T', '--period', default=None, type=float)
@click.option('-o', '--period_offset', default=0, type=float)
@click.option('-C', '--calscans', is_flag=True, default=True)
@click.option('-B', '--band', default='L', type=click.Choice(['low', 'mid', 'L', 'C', 'X']))
@click.option('-S', '--save', is_flag=True, default=False)
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
    dspec,
    lightcurve,
    spectrum,
    crosspols,
    acf,
    fold,
    period,
    period_offset,
    calscans,
    band,
    save,
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

    ds = DynamicSpectrum(
        project=project,
        band=band,
        calscans=calscans,
        tavg=tavg,
        favg=favg,
        prefix=prefix,
        fold=fold,
        period=period,
        period_offset=period_offset,
        save=save,
    )

    # Dynamic Spectrum
    # --------------------------------------
    if dspec:
        for s in stokes:
            if real:
                ds.plot_ds(stokes=s, cmax=cmax[s])
            if imag:
                ds.plot_ds(stokes=s, cmax=cmax[s], imag=True)

    # Cross-polarisations sqrt(U^2 + V^2)
    if crosspols:
        ds.plot_crosspol_ds(cmax=cmax['I'])

    # Spectrum
    # --------------------------------------
    if spectrum:
        ds.plot_spectrum(stokes=stokes)

    # Light Curve
    # --------------------------------------
    if lightcurve:
        ds.plot_lightcurve(stokes=stokes)

    # Dynamic Spectrum 2D Auto-correlation Function
    # --------------------------------------
    if acf:
        for s in stokes:
            ds.plot_acf(stokes=s)

    plt.show()


if __name__ == '__main__':
    main()
