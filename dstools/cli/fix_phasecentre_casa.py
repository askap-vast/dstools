import click

from dstools.utils import parse_coordinates


@click.command()
@click.option("-p", "--phasecenter", nargs=2)
@click.argument("msname")
def main(msname, phasecenter):
    rotated_ms = msname.replace(".ms", ".rotated.target.ms")
    ra, dec = parse_coordinates(phasecenter)
    phaseshift(vis=msname, outputvis=rotated_ms, phasecenter=f"J2000 {ra} {dec}")


if __name__ == "__main__":
    main()
