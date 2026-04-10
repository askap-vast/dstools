import logging
from pathlib import Path

import click

from dstools.logger import filter_stdout, setupLogger

try:
    from casacore.tables import tableexists
    from fixms.fix_ms_corrs import fix_ms_corrs
    from fixms.fix_ms_dir import fix_ms_dir

    HAS_CASA_SUPPORT = True
except ImportError:
    HAS_CASA_SUPPORT = False

logger = logging.getLogger(__name__)


@filter_stdout("Successful read/write open of default-locked table")
def filtered_fixms(ms: Path):
    # Fix beam pointing
    fix_ms_dir(str(ms))

    # Convert instrumental pol visibilities from average to total flux
    fix_ms_corrs(ms, chunksize=250)

    return


@click.command(context_settings={"show_default": True})
@click.argument("ms", type=Path)
def main(ms):
    setupLogger(verbose=False)

    if not HAS_CASA_SUPPORT:
        logger.error("ASKAP preprocessing is not supported on this system.")
        raise SystemExit(1)

    if tableexists(f"{ms}/FIELD_OLD"):
        logger.error("ASKAP beam pointing and flux re-scaling already applied.")
        raise SystemExit(1)

    filtered_fixms(ms)

    return


if __name__ == "__main__":
    main()
