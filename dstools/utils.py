import logging
from dataclasses import dataclass

import astropy.units as u
from astropy.coordinates import SkyCoord

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
