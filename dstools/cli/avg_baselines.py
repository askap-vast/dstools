import subprocess

import click

from dstools.utils import parse_casa_args


@click.command()
@click.option(
    "-u",
    "--uvrange",
    type=str,
    default="0m",
    help="Minimum uv distance to select from data column (e.g. -u 200m).",
)
@click.option(
    "-c",
    "--datacolumn",
    type=click.Choice(["data", "corrected"]),
    default="data",
    help="Selection of DATA or CORRECTED_DATA column.",
)
@click.argument("msname")
def main(**kwargs):
    # Construct string call signature to pass on to CASA
    path, argstr, kwargstr = parse_casa_args(
        main, "avg_baselines_casa.py", kwargs, args=["msname"]
    )

    call = f"casa -c {path} {argstr} {kwargstr}".split(" ")
    subprocess.run(call)


if __name__ == "__main__":
    main()
