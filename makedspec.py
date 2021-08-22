#!/usr/bin/env python
import click
import os
import numpy as np
from casacore.tables import table, taql


@click.command()
@click.option('-F', '--noflag', is_flag=True, default=False, help='Remove flagging mask.')
@click.option('-N', '--versionname', default=None,
              help='Prefix for different processing versions')
@click.argument('ms')
def main(noflag, versionname, ms):

    src_name = '/'.join(ms.split('/')[:-1])

    if versionname:
        prefix = f'{versionname}_'
    else:
        prefix = ''

    ds_path = "{}/dynamic_spectra/".format(src_name)
    os.system("mkdir {}".format(ds_path))

    pols = ["XX", "XY", "YX", "YY"]

    for polidx, pol in enumerate(pols):

        outfile = "{}/{}dynamic_spectra_{}.npy".format(ds_path, prefix, pol)

        t = table(ms)
        tf = table("{}/SPECTRAL_WINDOW".format(ms))

        # Write time and frequency arrays, discarding duplicate timesamples
        times = np.unique(t.getcol('TIME'))
        times.dump("{}/{}time.npy".format(ds_path, prefix))
        tf[0]["CHAN_FREQ"].dump("{}/{}freq.npy".format(ds_path, prefix))

        val = len(times)
        
        try:
            waterfall = t.getcol("CORRECTED_DATA")[:, :, polidx]
        except RuntimeError:
            print("No corrected data, using uncorrected.")
            waterfall = t.getcol("DATA")[:, :, polidx]

        vis_flag = t.getcol("FLAG")[:, :, polidx]
        if not noflag:
            waterfall = np.ma.masked_where(vis_flag, waterfall)

        waterfall.dump(outfile)
        t.close()

if __name__ == '__main__':
    main()
