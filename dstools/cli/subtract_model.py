import logging
import os
from pathlib import Path

import astropy.units as u
import click
import numpy as np
from dstools.casa import split, uvsub
from dstools.imaging import CASAModel, WSCleanModel
from dstools.logger import setupLogger
from dstools.utils import column_exists, parse_coordinates

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-S",
    "--split-data/--no-split-data",
    default=False,
    help="Split subtracted data into DATA column of output MS with .subbed.ms suffix.",
)
@click.option(
    "-v",
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
@click.argument("ms")
def main(
    ms,
    split_data,
    verbose,
):

    setupLogger(verbose=verbose)

    if not column_exists(ms, column="MODEL_DATA"):
        logger.error(
            f"{ms} does not contain a MODEL_DATA column. Create or insert a model first!"
        )
        exit(1)

    # Perform field model subtraction
    # ------------------------------

    uvsub(vis=ms)

    # Split subtracted visibilities into new file
    # ------------------------------------------

    if split_data:
        split(
            vis=ms,
            outputvis=ms.replace(".ms", ".subbed.ms"),
            datacolumn="corrected",
        )

    return


if __name__ == "__main__":
    main()
