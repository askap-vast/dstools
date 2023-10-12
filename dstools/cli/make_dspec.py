#!/usr/bin/env python
import click
import os
import numpy as np
from casacore.tables import table, taql

from dstools.utils import BANDS

@click.command()
@click.option('-F', '--noflag', is_flag=True, default=False, help='Remove flagging mask.')
@click.option('-B', '--band', default='AT_L', type=click.Choice(BANDS))
@click.option('-N', '--versionname', default=None, help='Prefix for different processing versions')
@click.option('-c', '--datacolumn', type=click.Choice(['data', 'corrected']), default='data',
              help='Selection of DATA or CORRECTED_DATA column.')
@click.argument('ms')
def main(noflag, band, versionname, datacolumn, ms):

    # Parse path into target directory to create dynamic spectra products
    if ms[-1] == '/':
        ms = ms[:-1]
   
    project = '/'.join(ms.split('/')[:-1])
    project = project if len(project) > 0 else '.'

    prefix = f'{versionname}_' if versionname else ''
    ds_path = f'{project}/dynamic_spectra/{band}'
    os.system(f'mkdir -p {ds_path}')

    pols = ['XX', 'XY', 'YX', 'YY']

    for polidx, pol in enumerate(pols):

        file_prefix = f'{ds_path}/{prefix}'
        outfile = f'{file_prefix}dynamic_spectra_{pol}.npy'

        t = table(ms)
        tf = table(f'{ms}/SPECTRAL_WINDOW')

        # Write time and frequency arrays, discarding duplicate timesamples
        times = np.unique(t.getcol('TIME'))
        times.dump(f'{file_prefix}time.npy')
        tf[0]['CHAN_FREQ'].dump(f'{file_prefix}freq.npy')

        val = len(times)

        datacolumn = 'CORRECTED_DATA' if datacolumn == 'corrected' else 'DATA'
        waterfall = t.getcol(datacolumn)[:, :, polidx]

        if not noflag:
            vis_flag = t.getcol('FLAG')[:, :, polidx]
            waterfall = np.ma.masked_where(vis_flag, waterfall)

        waterfall.dump(outfile)
        t.close()


if __name__ == '__main__':
    main()
