import subprocess

import click

from dstools.utils import parse_casa_args


@click.command()
@click.option(
    "-u",
    "--minuvdist",
    type=float,
    default=0,
    help="Minimum UV distance in meters to retain if averaging over baseline axis.",
)
@click.argument("msname")
def main(**kwargs):
    # Construct string call signature to pass on to CASA
    path, argstr, kwargstr = parse_casa_args(
        main,
        "avg_baselines_casa.py",
        kwargs,
        args=["msname"],
    )

    call = f"casa --nologger --nologfile -c {path} {argstr} {kwargstr}".split(" ")
    subprocess.run(call)


if __name__ == "__main__":
    main()
