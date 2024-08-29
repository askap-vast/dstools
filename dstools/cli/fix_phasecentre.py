import shutil
import subprocess

import click

import dstools


@click.command()
@click.argument("msname")
@click.argument("phasecenter", nargs=2)
def main(msname, phasecenter):
    path = dstools.__path__[0]
    path = f"{path}/cli/fix_phasecentre_casa.py"

    ra, dec = phasecenter

    casa_bin = shutil.which("casa")
    call = f"{casa_bin} --nologger --nologfile -c {path} {msname} -p {ra} {dec}".split(
        " "
    )
    subprocess.run(call)


if __name__ == "__main__":
    main()
