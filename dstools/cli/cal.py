import subprocess

import click

import dstools


@click.command()
@click.option(
    "-m",
    "--mfinterval",
    default="1.0",
    type=str,
    help="Time interval to solve for antenna gains in bandpass calibration.",
)
@click.option(
    "-b",
    "--bpinterval",
    default="1.0",
    type=str,
    help="Time interval to solve for bandpass in bandpass calibration.",
)
@click.option(
    "-g",
    "--gpinterval",
    default="0.1",
    type=str,
    help="Time interval to solve for antenna gains in gain calibration.",
)
@click.option(
    "-F",
    "--noflag",
    is_flag=True,
    default=False,
    help="Disable birdie and rfiflag options in atlod.",
)
@click.option(
    "-r",
    "--refant",
    type=click.Choice(["1", "2", "3", "4", "5", "6"]),
    default="3",
    help="Reference antenna.",
)
@click.option(
    "-d",
    "--data_dir",
    default="data",
    help="Path to directory containing raw miriad RPFITS visibilities.",
)
@click.option(
    "-A",
    "--auto",
    is_flag=True,
    default=False,
)
@click.argument("project_dir")
@click.argument("project_code")
def main(
    project_dir,
    project_code,
    mfinterval,
    bpinterval,
    gpinterval,
    noflag,
    refant,
    data_dir,
    auto,
):
    path = dstools.__path__[0]

    script = "atca_quickcal.sh" if auto else "atca_cal.sh"
    call = [
        f"{path}/{script}",
        path,
        project_dir,
        data_dir,
        project_code,
        refant,
        str(mfinterval),
        str(bpinterval),
        str(gpinterval),
        str(noflag).lower(),
    ]

    subprocess.run(call)

    return


if __name__ == "__main__":
    main()
