import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt

from dstools.ms import MeasurementSet
from dstools.utils import DataError, prompt

logger = logging.getLogger(__name__)


def run_selfcal(
    ms: MeasurementSet,
    calmode: str,
    gaintype: str,
    interval: str,
    split_data: bool,
    interactive: bool,
    refant: str = None,
    nspws: int = 1,
) -> Path:
    """Perform self-calibration on MS with field model in the MODEL_DATA column."""

    try:
        unit = "min" if "min" in interval else "s" if "s" in interval else ""
        int(interval.replace(unit, ""))
    except ValueError:
        raise ValueError(
            "Argument 'interval' must have format <int>[min/s] (e.g. 10s, 1min, 5min)."
        )

    if not ms.column_exists("MODEL_DATA"):
        raise DataError(f"{ms} does not contain a MODEL_DATA column.")

    # Produce MS with multiple spectral windows
    ms.to_nspws(nspws)

    # Solve for self calibration solutions
    ms.solve_gains(
        interval,
        calmode,
        gaintype,
        refant=refant,
        interactive=interactive,
    )

    # Generate phase and amplitude calibration plots
    for mode in calmode:
        ms.caltable.plot_solutions(mode)

    if interactive:
        plt.show(block=False)

    # Confirm solution is good before applying
    cal_good = prompt(
        msg="Apply gain solutions?",
        bypass=not interactive,
        default_response=True,
    )

    # If unacceptable, remove calibration tables, plots, and multi-spw MS and return
    if not cal_good:
        os.system(f"rm -r {ms.caltable.path}")
        os.system(f"rm {ms.path.stem}*.png")

        if nspws > 1:
            os.system(f"rm -r {ms.path} ")
            ms.path = ms.original_path
            ms.original_path = None

        return ms

    # Otherwise proceed with applying calibration solutions
    ms.applycal()

    # Transform back to single SPW MS
    ms.to_nspws(1)

    # Split out calibrated MS
    if split_data:
        ms = ms.split_selfcal_round()

    plt.close("all")

    return ms
