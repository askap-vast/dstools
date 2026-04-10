import logging
from pathlib import Path

import click

from dstools.logger import setupLogger

try:
    from dstools.ms import MeasurementSet

    HAS_CASA_SUPPORT = True
except ImportError:
    HAS_CASA_SUPPORT = False

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-d",
    "--datacolumn",
    type=click.Choice(["data", "corrected", "model"]),
    default="data",
    help="Selection of DATA, CORRECTED_DATA, or MODEL column.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms", type=Path)
def main(ms, datacolumn, verbose):
    setupLogger(verbose=verbose)

    if not HAS_CASA_SUPPORT:
        logger.error("Feed derotation is not supported on this system.")
        raise SystemExit(1)

    ms = MeasurementSet(ms)

    columns = {
        "data": "DATA",
        "corrected": "CORRECTED_DATA",
        "model": "MODEL_DATA",
    }
    datacolumn = columns[datacolumn]

    # Swap incorrectly labelled X and Y feeds for MeerKAT
    if ms.telescope == "MeerKAT":
        ms.swap_xy_feeds(datacolumn=datacolumn)

    ms.correct_feed_rotation(datacolumn=datacolumn)

    return


if __name__ == "__main__":
    main()
