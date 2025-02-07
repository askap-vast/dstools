import logging

import click

from dstools.logger import setupLogger
from dstools.ms import MeasurementSet
from dstools.utils import DataError

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-S",
    "--split-ms/--no-split-ms",
    default=False,
    help="Split subtracted data into DATA column of output MS with .subbed.ms suffix.",
)
@click.argument("ms", type=MeasurementSet)
def main(ms, split_ms):

    setupLogger(verbose=False)

    # Perform field model subtraction
    # ------------------------------
    try:
        ms.subtract_model(split_ms=split_ms)
    except DataError as e:
        logger.error(e)
        exit(1)

    return


if __name__ == "__main__":
    main()
