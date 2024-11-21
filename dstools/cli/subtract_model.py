import logging
import os
from pathlib import Path

import astropy.units as u
import click
import numpy as np
from casaconfig import config

config.logfile = "/dev/null"
from casatasks import split, uvsub
from dstools.imaging import CASAModel, WSCleanModel
from dstools.logger import setupLogger
from dstools.utils import parse_coordinates

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
    split_data,
    verbose,
    ms,
):

    setupLogger(verbose=verbose)

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
