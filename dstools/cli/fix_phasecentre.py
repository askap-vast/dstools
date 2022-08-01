import click
import subprocess

import dstools


@click.command()
@click.argument('msname')
@click.argument('phasecenter')
def main(msname, phasecenter):

    path = dstools.__path__[0]
    path = f'{path}/cli/fix_phasecentre_casa.py'
    
    call = f'casa -c {path} {msname}'.split(' ') + [phasecenter]
    subprocess.run(call)

if __name__ == '__main__':
    main()
