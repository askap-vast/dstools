import numpy as np
import sys


msname = sys.argv[1]
if len(sys.argv) > 2:
    uvrange = ">{}".format(sys.argv[2])
else:
    uvrange = ">0m"

savename = msname.replace(".ms", ".baseavg.ms")

intab = tbtool()
intab.open(msname, nomodify=False)

ant1 = intab.getcol("ANTENNA1")
ant2 = intab.getcol("ANTENNA2")

# Set all antenna pairs equal for baseline averaging
nrows = intab.nrows()
intab.putcol("ANTENNA1", np.zeros(nrows))
intab.putcol("ANTENNA2", np.ones(nrows))

interval = intab.getcol("INTERVAL")
timebin = "{}s".format(min(interval) * 1e-2)

# average over baselines greater than 200m
mstransform(
    vis=msname,
    outputvis=savename,
    datacolumn="corrected",
    uvrange=uvrange,
    timeaverage=True,
    timebin=timebin,
    keepflags=False,
)

intab.putcol("ANTENNA1", ant1)
intab.putcol("ANTENNA2", ant2)

intab.unlock()
intab.close()
