import click


@click.command()
@click.argument('msname')
@click.argument('phasecenter')
def main(msname, phasecenter):

    rotated_ms = msname.replace('.ms', '.rotated.target.ms')
    phaseshift(vis=msname, outputvis=rotated_ms, phasecenter=phasecenter)

if __name__ == '__main__':
    main()
