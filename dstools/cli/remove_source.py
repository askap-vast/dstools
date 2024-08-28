import shutil
import click
import subprocess

from dstools.utils import BANDS, CONFIGS, parse_casa_args


@click.command()
@click.option(
    "-C",
    "--config",
    type=click.Choice(CONFIGS),
    default="6km",
    help="Array configuration, used to calculate image and pixel sizes. ASKAP is equivalent to 6km.",
)
@click.option(
    "-B",
    "--band",
    default="AT_L",
    type=click.Choice(BANDS),
    help="Observing band, used to calculate image and pixel sizes.",
)
@click.option(
    "-N",
    "--iterations",
    default=10000,
    help="Maximum number of clean iterations.",
)
@click.option(
    "-t",
    "--threshold",
    default=8e-5,
    help="Clean threshold in Jy.",
)
@click.option(
    "-n",
    "--nterms",
    default=3,
    help="Number of Taylor terms to use in imaging.",
)
@click.option(
    "-S",
    "--scale",
    default=[0, 4, 8],
    type=int,
    multiple=True,
    help="Multi-scale clean pixel smoothing scales.",
)
@click.option(
    "-w",
    "--widefield/--no-widefield",
    default=False,
    help="Enable w-projection to correct for widefield artefacts in off-axis sources.",
)
@click.option(
    "-r",
    "--robust",
    default=0.5,
    help="Briggs weighting robust parameter.",
)
@click.option(
    "-l",
    "--pblim",
    default=-0.1,
    help="Primary beam fractional power cutoff, default is no cutoff.",
)
@click.option(
    "-I",
    "--interactive/--no-interactive",
    default=True,
    help="Run tclean in interactive mode",
)
@click.option(
    "-i",
    "--imsize",
    default=400,
    help="Image size in pixels.",
)
@click.option(
    "-m",
    "--mpinodes",
    default=1,
    help="Set greater than 1 to run casa in mpi mode with mpinodes available nodes.",
)
@click.argument("data")
@click.argument("killcoord", nargs=2)
def main(**kwargs):

    # Read off multiple scale arguments from tuple and rewrite as individual flags
    scales = kwargs.pop("scale")
    scaleargs = [f" -S {clean_scale}" for clean_scale in scales]

    mpinodes = kwargs.pop("mpinodes")
    killcoord = kwargs.pop("killcoord")

    # Construct string call signature to pass on to CASA
    path, argstr, kwargstr = parse_casa_args(
        main,
        "remove_source_casa.py",
        kwargs,
        args=["data"],
    )

    kwargstr += "".join(scaleargs)

    casa_bin = shutil.which("casa")
    if mpinodes > 1:
        casa_cmd = f"mpicasa -n {mpinodes} {casa_bin}"
        argstr += " --mpi"
    else:
        casa_cmd = casa_bin

    call = f"{casa_cmd} --nologger -c {path} {argstr} {kwargstr}".split(" ") + [
        "-k",
        *killcoord,
    ]
    subprocess.run(call)


if __name__ == "__main__":
    main()
