import logging

import click

from dstools.logger import setupLogger
from dstools.ms import MeasurementSet, run_selfcal
from dstools.utils import DataError

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-m",
    "--calmode",
    type=click.Choice(["p", "ap"]),
    default="p",
    help="Solve for phase only (p) or amplitude + phase (ap) gain corrections.",
)
@click.option(
    "-i",
    "--interval",
    type=str,
    default="10s",
    help="Time interval over which to solve for gains in format <int>[min/s].",
)
@click.option(
    "-r",
    "--refant",
    type=str,
    default=None,
    help="Name of reference antenna.",
)
@click.option(
    "-P",
    "--combine-pols/--no-combine-pols",
    is_flag=True,
    default=False,
    help="Solve for calibration solutions on combined orthogonal polarisations.",
)
@click.option(
    "-S",
    "--split-data/--no-split-data",
    is_flag=True,
    default=False,
    help="Split self-calibrated data into DATA column of output MS with .selfcal.ms suffix.",
)
@click.option(
    "--interactive/--no-interactive",
    is_flag=True,
    default=True,
    help="Run interactively to select solution in",
)
@click.option(
    "-n",
    "--nspws",
    default=1,
    type=int,
    help="Number of spectral windows to derive indepdendent gain solutions within.",
)
@click.argument("ms", type=MeasurementSet)
def main(ms, calmode, interval, refant, combine_pols, split_data, interactive, nspws):
    setupLogger(verbose=False)

    gaintype = "T" if combine_pols else "G"
    try:
        run_selfcal(
            ms,
            calmode=calmode,
            gaintype=gaintype,
            interval=interval,
            refant=refant,
            split_data=split_data,
            interactive=interactive,
            nspws=nspws,
        )
    except (ValueError, DataError) as exc:
        logger.error(exc)
        exit(1)

    return


if __name__ == "__main__":
    main()
