import click
import os
import dstools

@click.command()
@click.option('-m', '--mfinterval', default=1.0,
              help='Time interval to solve for antenna gains in bandpass calibration')
@click.option('-b', '--bpinterval', default=1.0,
              help='Time interval to solve for bandpass in bandpass calibration')
@click.option('-g', '--gpinterval', default=0.1,
              help='Time interval to solve for antenna gains in gain calibration')
@click.option('-r', '--refant', type=click.Choice(['1', '2', '3', '4', '5', '6']), default='1',
              help='Reference antenna.')
@click.option('-d' ,'--data_dir', default='data',
              help='Path to directory containing raw miriad RPFITS visibilities.')
@click.argument('project_dir')
@click.argument('project_code')
def main(project_dir, project_code, mfinterval, bpinterval, gpinterval, refant, data_dir):

    path = dstools.__path__[0]

    call = '{}/atca_cal.sh {} {} {} {} {} {} {} {}'.format(
        path,
        path,
        project_dir,
        data_dir,
        project_code,
        refant,
        mfinterval,
        bpinterval,
        gpinterval,
    )

    os.system(call)

if __name__ == '__main__':
    main()
