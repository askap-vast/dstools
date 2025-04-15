import logging
import multiprocessing
import os
from dataclasses import dataclass
from typing import Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from numpy.typing import ArrayLike

CONFIGS = ["6km", "750_no6", "750_6", "H168"]
BANDS = [
    "AK_low",
    "AK_mid",
    "AK_high",
    "AT_L",
    "AT_C",
    "AT_X",
    "AT_K",
    "MKT_UHF",
    "MKT_L",
]


logger = logging.getLogger(__name__)


class DataError(Exception):
    pass


def parse_coordinates(coord: tuple[str, str]) -> tuple[str, str]:
    """Convert decimal degrees or hexagesimal coordinates to hms dms format."""

    ra, dec = coord
    raunit = "hourangle" if ":" in ra or "h" in ra else "deg"
    pos = SkyCoord(ra=ra, dec=dec, unit=(raunit, "deg"))
    ra, dec = pos.to_string(style="hmsdms", precision=3).split()

    return ra, dec
def get_available_cpus():
    """Returns the number of CPUs allocated by SLURM or falls back to system count."""

    if "SLURM_CPUS_PER_TASK" in os.environ:
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    elif "SLURM_CPUS_ON_NODE" in os.environ:
        return int(os.environ["SLURM_CPUS_ON_NODE"])
    else:
        return multiprocessing.cpu_count()


def prompt(msg, bypass=False, bypass_msg=None, default_response=True):
    if bypass:
        if bypass_msg is not None:
            logger.warning(bypass_msg)
        return default_response

    msg = f"{msg} (y/n)\n"

    resp = input(msg)
    if resp not in ["y", "n"]:
        resp = input(msg)

    return True if resp == "y" else False


def rebin(o: int, n: int, axis: int) -> ArrayLike:
    """Create l1-norm preserving array compression matrix from o -> n length.

    if rebinning along row axis we want:
        - (o // n) + 1 entries in each row that sum to unity and preserve l1-norm,
        - each column to sum to the compression ratio o / n
        - values distributed along the row in units of o / n until expired

        >>> rebin(5, 3, axis=0)
        array([[0.6, 0.4, 0. , 0. , 0. ],
               [0. , 0.2, 0.6, 0.2, 0. ],
               [0. , 0. , 0. , 0.4, 0.6]])

        - transpose of this for column rebinning

    The inner product of this compressor with an array will rebin
    the array conserving the total intensity along the given axis.
    """

    compressor = np.zeros((n, o))

    # Exit early with empty array if chunk is empty
    if compressor.size == 0:
        return compressor

    comp_ratio = n / o

    nrow = 0
    ncol = 0

    budget = 1
    overflow = 0

    # While loop to avoid visiting n^2 zero-value cells
    while nrow < n and ncol < o:
        # Use overflow if just spilled over from last row
        if overflow > 0:
            value = overflow
            overflow = 0
            budget -= value
            row_shift = 0
            col_shift = 1
        # Use remaining budget if at end of current row
        elif budget < comp_ratio:
            value = budget
            overflow = comp_ratio - budget
            budget = 1
            row_shift = 1
            col_shift = 0
        # Otherwise spend n / o and move to next column
        else:
            value = comp_ratio
            budget -= value
            row_shift = 0
            col_shift = 1

        compressor[nrow, ncol] = value
        nrow += row_shift
        ncol += col_shift

    return compressor if axis == 0 else compressor.T


def rebin2D(array: ArrayLike, new_shape: tuple[int, int]) -> ArrayLike:
    """Re-bin along time / frequency axes conserving flux."""

    # Convert from masked array to pure numpy array
    if isinstance(array, np.ma.MaskedArray):
        array[array.mask] = np.nan
        array = array.data

    if new_shape == array.shape:
        array[array == 0 + 0j] = np.nan
        return array

    if new_shape[0] > array.shape[0] or new_shape[1] > array.shape[1]:
        raise ValueError(
            "New shape should not be greater than old shape in either dimension"
        )

    time_comp = rebin(array.shape[0], new_shape[0], axis=0)
    freq_comp = rebin(array.shape[1], new_shape[1], axis=1)
    array[np.isnan(array)] = 0 + 0j
    result = time_comp @ np.array(array) @ freq_comp
    result[result == 0 + 0j] = np.nan

    return result


def slice_array(
    a: ArrayLike,
    ax1_min: int,
    ax1_max: int,
    ax2_min: Optional[int] = None,
    ax2_max: Optional[int] = None,
) -> ArrayLike:
    """Slice 1D or 2D array with variable lower and upper boundaries."""

    if ax2_min is None and ax2_max is None:
        a = a[ax1_min:] if ax1_max == 0 else a[ax1_min:ax1_max]
    else:
        a = a[ax1_min:, :] if ax1_max == 0 else a[ax1_min:ax1_max, :]
        a = a[:, ax2_min:] if ax2_max == 0 else a[:, ax2_min:ax2_max]

    return a


@dataclass
class Array:
    band: str = "AT_L"
    config: str = "6km"

    def __post_init__(self):
        telescope = self.band.split("_")[0]
        self.config = self.config if telescope == "AT" else telescope

        # Frequencies are taken at centre of band for Taylor expansion
        frequencies = {
            "AK_low": "888.49",
            "AK_mid": "1367.49",
            "AK_high": "1655.49",
            "AT_L": "2100",
            "AT_C": "5500",
            "AT_X": "9000",
            "AT_K": "17000",
            "MKT_UHF": "797.5",
            "MKT_L": "1285",
        }

        imsize = {
            "AK_low": 6144,
            "AK_mid": 4500,
            "AK_high": 4500,
            "AT_L": 4500,
            "AT_C": 2048,
            "AT_X": 2048,
            "AT_K": 1200,
            "MKT_UHF": 4500,
            "MKT_L": 6144,
        }

        cellsize = {
            "AK_low": 2.5,
            "AK_mid": 1.5,
            "AK_high": 1.25,
            "AT_L": 0.66,
            "AT_C": 0.32,
            "AT_X": 0.21,
            "AT_K": 0.10,
            "MKT_UHF": 1.5,
            "MKT_L": 1,
        }

        self.frequency = frequencies[self.band]
        self.cell = cellsize[self.band]
        self.imsize = imsize[self.band]

    @property
    def imradius(self):
        return self.imsize * self.cell * u.arcsec.to(u.deg) / 2

    def __str__(self):
        return str(
            {
                "band": self.band,
                "config": self.config,
                "frequency": self.frequency,
                "cell": self.cell,
                "imradius": self.imradius,
                "imsize": self.imsize,
            }
        )
