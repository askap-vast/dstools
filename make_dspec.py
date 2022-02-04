#!/usr/bin/env python
import click
import os
import numpy as np
from casacore.tables import table, taql


@click.command()
@click.option("-F", "--noflag", is_flag=True, default=False, help="Remove flagging mask.")
@click.option("-B", "--band", default="L", type=click.Choice(["L", "C", "X", "low"]))
@click.argument("ms")
def main(noflag, band, ms):

    project = "/".join(ms.split("/")[:-1])

    ds_path = "{}/dynamic_spectra/{}".format(project, band)
    os.system("mkdir -p {}".format(ds_path))

    pols = ["XX", "XY", "YX", "YY"]

    for polidx, pol in enumerate(pols):

        outfile = "{}/dynamic_spectra_{}.npy".format(ds_path, pol)

        t = table(ms)
        tf = table("{}/SPECTRAL_WINDOW".format(ms))

        # Write time and frequency arrays, discarding duplicate timesamples
        times = np.unique(t.getcol("TIME"))
        times.dump("{}/time.npy".format(ds_path, project))
        tf[0]["CHAN_FREQ"].dump("{}/freq.npy".format(ds_path, project))

        val = len(times)

        try:
            waterfall = t.getcol("CORRECTED_DATA")[:, :, polidx]
        except RuntimeError as e:
            print("No corrected data, using uncorrected.")
            waterfall = t.getcol("DATA")[:, :, polidx]

        vis_flag = t.getcol("FLAG")[:, :, polidx]
        if not noflag:
            waterfall = np.ma.masked_where(vis_flag, waterfall)

        # Multiply any ASKAP visibilities by 2 to convert
        # from average to total flux density
        if band == "low":
            waterfall *= 2

        waterfall.dump(outfile)
        t.close()


if __name__ == "__main__":
    main()
