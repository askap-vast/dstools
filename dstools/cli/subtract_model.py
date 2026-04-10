import logging
from pathlib import Path

import click

from dstools.logger import setupLogger
from dstools.utils import DataError

try:
    from dstools.ms import MeasurementSet

    HAS_CASA_SUPPORT = True
except ImportError:
    HAS_CASA_SUPPORT = False

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-S",
    "--split-ms/--no-split-ms",
    default=False,
    help="Split subtracted data into DATA column of output MS with .subtracted.ms suffix.",
)
@click.argument("ms", type=Path)
def main(ms, split_ms):
    setupLogger(verbose=False)

    if not HAS_CASA_SUPPORT:
        logger.error("Model subtraction is not supported on this system.")
        raise SystemExit(1)

    ms = MeasurementSet(ms)

    # Perform field model subtraction
    # ------------------------------
    try:
        ms.subtract_model(split_ms=split_ms)
    except DataError as e:
        logger.error(e)
        raise SystemExit(1)

    return


if __name__ == "__main__":
    main()
