import shutil
import click
import subprocess
from astropy.coordinates import SkyCoord

import dstools


@click.command()
@click.argument("msname")
@click.argument("ra")
@click.argument("dec")
def main(msname, ra, dec):
    path = dstools.__path__[0]
    path = f"{path}/cli/fix_phasecentre_casa.py"

    ra_unit = "hourangle" if ":" in ra or "h" in ra else "deg"
    c = SkyCoord(ra=ra, dec=dec, unit=(ra_unit, "deg"))
    c = c.to_string(style="hmsdms", precision=3)

    phasecenter = f"J2000 {c}"

    casa_bin = shutil.which("casa")
    call = f"{casa_bin} --nologger -c {path} {msname}".split(" ") + [phasecenter]
    subprocess.run(call)


if __name__ == "__main__":
    main()
