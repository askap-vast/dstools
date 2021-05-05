import numpy as np
import sys

msname = sys.argv[-2]
fixphase = sys.argv[-1]
savename = msname.replace('.ms', '.rotated.target.baseavg.ms')
rotated_ms = msname.replace('.ms', '.rotated.target.ms')

if fixphase == 'y':
    # phasecen = "J2000 21h16m06.23s +29d51m56.49s"
    phasecen = 'J2000 21h09m24.632s -56d18m43.851s'
    fixvis(vis=msname, outputvis=rotated_ms, phasecenter=phasecen)
    msname = rotated_ms

intab = tbtool()
intab.open(msname, nomodify=False)

ant1 = intab.getcol('ANTENNA1')
ant2 = intab.getcol('ANTENNA2')

# Set all antenna pairs equal for baseline averaging
nrows = intab.nrows()
intab.putcol('ANTENNA1', np.zeros(nrows))
intab.putcol('ANTENNA2', np.ones(nrows))

interval = intab.getcol('INTERVAL')
timebin = '{}s'.format(min(interval) * 1e-2)

# average over baselines greater than 200m
mstransform(vis=msname, outputvis=savename, datacolumn='data', uvrange='>200m',
            timeaverage=True, timebin=timebin, keepflags=False)

intab.putcol('ANTENNA1', ant1)
intab.putcol('ANTENNA2', ant2)

intab.unlock()
intab.close()
