import shutil
import subprocess

import click

from dstools.utils import parse_casa_args


@click.command()
@click.option(
    "-p",
    "--phasecenter",
    type=str,
    nargs=2,
    default=None,
    help="Imaging phasecenter",
)
@click.option(
    "-r",
    "--robust",
    default=0.5,
    help="Briggs weighting robust parameter.",
)
@click.option(
    "-t",
    "--target_position",
    type=str,
    nargs=2,
    default=None,
    help="Coordinates around which to mask model images before subtraction.",
)
@click.option(
    "-I",
    "--interactive/--no-interactive",
    is_flag=True,
    default=True,
    help="Mask out target source in interactive tclean.",
)
@click.option(
    "-w",
    "--widefield/--no-widefield",
    default=False,
    help="Enable w-projection to correct for widefield artefacts in off-axis sources.",
)
@click.option(
    "-m",
    "--mpinodes",
    default=1,
    help="Set greater than 1 to run casa in mpi mode with mpinodes available nodes.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms")
@click.argument("model", nargs=-1)
def main(**kwargs):

    # Read off phasecentre separately so that multi-word string can be parsed
    phasecenter = kwargs.pop("phasecenter")
    phasecenter = ["-p", *phasecenter] if phasecenter is not None else []

    target_position = kwargs.pop("target_position")
    target_position = ["-t", *target_position] if target_position is not None else []

    mpinodes = kwargs.pop("mpinodes")
    model = kwargs.pop("model")

    # Construct string call signature to pass on to CASA
    path, argstr, kwargstr = parse_casa_args(
        main,
        "subtract_model_casa.py",
        kwargs,
        args=["ms"],
    )
    argstr += " " + " ".join(list(model))

    casa_bin = shutil.which("casa")
    if mpinodes > 1:
        casa_bin = f"mpicasa -n {mpinodes} {casa_bin}"
        argstr += " --mpi"

    call = (
        f"{casa_bin} --nologger -c {path} {argstr} {kwargstr}".split(" ")
        + phasecenter
        + target_position
    )
    subprocess.run(call)


if __name__ == "__main__":
    main()
