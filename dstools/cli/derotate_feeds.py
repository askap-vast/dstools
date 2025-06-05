import logging

import click

from dstools.logger import setupLogger
from dstools.ms import MeasurementSet

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
@click.argument("ms", type=MeasurementSet)
def main(ms, datacolumn, verbose):
    setupLogger(verbose=verbose)

    columns = {
        "data": "DATA",
        "corrected": "CORRECTED_DATA",
        "model": "MODEL_DATA",
    }
    datacolumn = columns[datacolumn]

    ms.correct_feed_rotation(datacolumn=datacolumn)

    return


if __name__ == "__main__":
    main()
