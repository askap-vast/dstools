import shutil
import subprocess

import click

import dstools


@click.command()
@click.argument("ms")
@click.argument("outputvis")
def main(ms, outputvis):
    path = dstools.__path__[0]
    path = f"{path}/cli/combine_spws_casa.py"

    casa_bin = shutil.which("casa")

    call = f"{casa_bin} --nologger --nologfile -c {path} {ms} {outputvis}".split(" ")
    subprocess.run(call)


if __name__ == "__main__":
    main()
