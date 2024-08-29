import logging
import os

import click
from casacore.tables import tableexists

import dstools
from dstools.logger import setupLogger

logger = logging.getLogger(__name__)
setupLogger(verbose=False)


@click.command()
@click.argument("ms")
def main(ms):

    root = dstools.__path__[0]

    if tableexists(f"{ms}/FIELD_OLD"):
        logger.error(f"ASKAP beam pointing and flux re-scaling already applied.")
        exit(1)

    # Fix beam pointing
    os.system(f"python {root}/fix_dir.py {ms}")

    # Convert instrumental pol visibilities from average to total flux
    os.system(f"python {root}/rescale_askapsoft.py {ms}")

    return


if __name__ == "__main__":
    main()
