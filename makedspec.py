import numpy as np
import click
from casacore.tables import table, taql
import time
import os
import re


@click.command()
@click.option('-F', '--noflag', is_flag=True, default=False, help='Remove flagging mask.')
@click.argument('ms')
def main(noflag, ms):

    pattern = re.compile(r'\S*(SB\d{4,5}_POSSUM_)(\d{4}[-+]\d{2}.beam\d{2})\S*.ms')
    src_name = pattern.sub(r'\1\2', ms)
    sbid = ms.split('/')[0]

    ds_path = "{}/dynamic_spectra2/{}".format(sbid, src_name)
    os.system("mkdir {}/dynamic_spectra2".format(sbid))
    os.system("mkdir {}".format(ds_path))

    pols = ["XX", "XY", "YX", "YY"]

    for polidx, pol in enumerate(pols):

        outfile = "{}/dynamic_spectra_{}.npy".format(ds_path, pol)

        t = table(ms)
        tf = table("{}/SPECTRAL_WINDOW".format(ms))

        # Write time and frequency arrays, discarding duplicate timesamples
        times = np.unique(t.getcol('TIME'))
        times.dump("{}/time.npy".format(ds_path, src_name))
        tf[0]["CHAN_FREQ"].dump("{}/freq.npy".format(ds_path, src_name))

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
    t1 = time.time()
    main()
    t2 = time.time()

    print(t2-t1)
