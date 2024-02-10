import logging
import os

import click
from astroutils.logger import setupLogger

import dstools
from dstools.utils import BANDS

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "-d",
    "--datacolumn",
    default="data",
    type=click.Choice(["data", "corrected"]),
    help="Use DATA or CORRECTED_DATA column.",
)
@click.option("-B", "--band", default="AT_L", type=click.Choice(BANDS))
@click.option(
    "-a",
    "--askap",
    is_flag=True,
    default=False,
    help="Apply beam direction and flux scale corrections to ASKAPsoft visibilities.",
)
@click.option("-v", "--verbose", is_flag=True, default=False)
@click.argument("ms")
@click.argument("target_coord")
def main(
    ms,
    target_coord,
    datacolumn,
    band,
    askap,
    verbose,
):
    setupLogger(False)

    root = dstools.__path__[0]

    logger.info(f"Processing {ms}")

    if askap:
        # Fix beam pointing
        os.system(f"python {root}/fix_dir.py {ms}")

        # Convert instrumental pol visibilities from average to total flux
        os.system(f"python {root}/rescale_askapsoft.py {ms}")

    # Rotate to target
    os.system(f'dstools-rotate {ms} "J2000 {target_coord}"')
    ms = ms.replace(".ms", ".rotated.target.ms")

    # Average over baselines
    os.system(f"dstools-avg-baselines -u 500m -c {datacolumn} {ms}")
    ms = ms.replace(".ms", ".baseavg.ms")

    # Create DS arrays
    os.system(f"dstools-make-dspec -c data {ms} -B {band}")


if __name__ == "__main__":
    main()
