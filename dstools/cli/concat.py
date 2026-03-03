import logging
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import astropy.units as u
import click
import h5py
import numpy as np
from astropy.time import Time

from dstools.logger import setupLogger
from dstools.utils import DataError

logger = logging.getLogger(__name__)

DSTOOLS_VERSION = version("radio-dstools")


@dataclass
class HDF5DS:
    header: dict
    flux: np.ndarray
    time: np.ndarray
    freq: np.ndarray
    uvdist: np.ndarray


def validate(f: h5py.File):
    # Do not allow combination of dynamic spectra extracted with old versions
    ds_version = f.attrs.get("dstools_version", "1.0.0")
    if not ds_version == DSTOOLS_VERSION:
        msg = f"Can not concat v{ds_version} DS {f.filename}. Re-extract with v{DSTOOLS_VERSION}"
        raise DataError(msg)


def read_ds(ds: Path) -> HDF5DS:
    with h5py.File(ds, "r") as f:
        validate(f)

        header = dict(f.attrs)
        flux = f["flux"][:]
        time = f["time"][:]
        freq = f["frequency"][:]
        uvdist = f["uvdist"][:]

    hdf5_ds = HDF5DS(
        header=header,
        flux=flux,
        time=time,
        freq=freq,
        uvdist=uvdist,
    )

    return hdf5_ds


def valid_concat_list(ds_list: list[HDF5DS]) -> bool:
    # Must have equal number of channels
    # nchan = np.array([len(ds.freq) for ds in ds_list])
    # all_chan_eq = np.all(nchan == nchan[0])
    # chan0 = np.array([ds.freq[0] for ds in ds_list])
    # all_chan0_eq = np.all(chan0 == chan0[0])

    # if not (all_chan_eq and all_chan0_eq):
    #     logger.error("All datasets must share the same channel range.")
    #     return False

    # All must be baseline-averaged
    nbaselines = np.array([len(ds.uvdist) for ds in ds_list])
    if not np.all(nbaselines == 1):
        logger.error("All datasets must be baseline averaged")
        return False

    # Same time and frequency resolution
    chan_res = [abs(round(ds.freq[1] - ds.freq[0])) for ds in ds_list]
    time_res = [round(ds.time[1] - ds.time[0], 1) for ds in ds_list]

    all_chan_res_eq = np.allclose(chan_res, chan_res[0])
    all_time_res_eq = np.allclose(time_res, time_res[0])

    if not (all_chan_res_eq and all_time_res_eq):
        logger.error("All datasets must share the same time/frequency resolution.")
        return False

    # Must have no overlapping timestamps
    times = np.concat([np.array(ds.time) for ds in ds_list]).ravel()
    no_dupe_times = len(times) == len(np.unique(times))

    if not no_dupe_times:
        logger.error("Datasets must not share overlapping timestamps.")
        return False

    return True


@click.command()
@click.argument("dynamic-spectra", nargs=-1, type=Path)
@click.argument("out-ds", type=Path)
def main(dynamic_spectra, out_ds):
    setupLogger(verbose=True)

    dynamic_spectra_files = sorted(
        [read_ds(ds) for ds in dynamic_spectra],
        key=lambda ds: ds.time[0],
    )

    valid = valid_concat_list(dynamic_spectra_files)

    if not valid:
        exit(1)

    if len(dynamic_spectra_files) < 2:
        logger.error(f"Only {len(dynamic_spectra_files)} DS...")

    with h5py.File(out_ds, "w", track_order=True) as f:
        g = f.create_group("observations")
        for i, ds in enumerate(dynamic_spectra_files):
            telescope = ds.header["telescope"]
            t0 = Time((ds.time[0] * u.s).to(u.day), format="mjd").isot
            obs_id = f"{telescope}-{t0}"

            obs = g.create_group(obs_id)

            for attr in ds.header:
                obs.attrs[attr] = ds.header[attr]

            obs.create_dataset("flux", compression="gzip", data=ds.flux)
            obs.create_dataset("time", compression="gzip", data=ds.time)
            obs.create_dataset("frequency", compression="gzip", data=ds.freq)
            obs.create_dataset("uvdist", compression="gzip", data=ds.uvdist)


if __name__ == "__main__":
    main()
