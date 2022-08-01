import click
import subprocess

from dstools.utils import parse_casa_args

@click.command()
@click.option('-C', '--config', type=click.Choice(['6km', '750_no6', '750_6', 'H168']), default='6km',
              help='Array configuration, used to calculate image and pixel sizes. ASKAP is equivalent to 6km.')
@click.option('-B', '--band', type=click.Choice(['low', 'mid', 'L', 'C', 'X']), default='L',
              help='Observing band, used to calculate image and pixel sizes.')
@click.option('-N', '--iterations', default=2000,
              help='Maximum number of clean iterations.')
@click.option('-t', '--threshold', default=8e-5,
              help='Clean threshold in Jy.')
@click.option('-n', '--nterms', default=3,
              help='Number of Taylor terms to use in imaging.')
@click.option('-S', '--scale', default=[0, 4, 8], type=int, multiple=True,
              help='Multi-scale clean pixel smoothing scales.')
@click.option('-w', '--widefield/--no-widefield', default=False,
              help='Enable w-projection to correct for widefield artefacts in off-axis sources.')
@click.option('-r', '--robust', default=0.5,
              help='Briggs weighting robust parameter.')
@click.option('-p', '--phasecenter', default='',
              help='Imaging phasecentre (in format J2000 hms dms).')
@click.option('-l', '--pblim', default=-0.1,
              help='Primary beam fractional power cutoff, default is no cutoff.')
@click.option('-P', '--proj_dir', type=click.Path(), default=None,
              help='Project namespace / sub-directory, default is to infer from data path')
@click.option('-I', '--interactive/--no-interactive', default=True,
              help='Run tclean in interactive mode')
@click.option('--log2term/--no-log2term', default=True,
              help='Disable CASA logger GUI and log straight to terminal.')
@click.argument('data')
def main(**kwargs):

    logconfig = ' --log2term --nologger --logfile ./casa.log' if kwargs.pop('log2term') else ''

    # Read off multiple scale arguments from tuple and rewrite as individual flags
    scales = kwargs.pop('scale')
    scaleargs = [f' -S {clean_scale}' for clean_scale in scales]

    # Read off phasecentre separately so that multi-word string can be parsed 
    phasecenter = kwargs.pop('phasecenter')
    phasecenter = ['-p', phasecenter] if phasecenter is not None else []

    # Construct string call signature to pass on to CASA
    path, argstr, kwargstr = parse_casa_args(main, 'model_field_casa.py', kwargs, args=['data'])
    kwargstr += ''.join(scaleargs)

    call = f'casa{logconfig} -c {path} {argstr} {kwargstr}'.split(' ') + phasecenter
    subprocess.run(call)

if __name__ == '__main__':
    main()
