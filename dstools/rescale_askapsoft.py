#!/usr/bin/env python

import os
import sys

import numpy as np
from casacore.tables import table

ms_file = sys.argv[-2]
ms_file_out = sys.argv[-1]
os.system("cp -R %s %s" % (ms_file, ms_file_out))
t = table(ms_file_out, readonly=False, ack=False)
nrows = t.nrows()
for row in range(nrows):
    if row % 1000 == 0:
        print("%d/%d" % (row, nrows))
    cdata = t.getcol("DATA", startrow=row, nrow=1)
    cdata *= 2.0
    t.putcol("DATA", cdata, startrow=row, nrow=1)
t.close()
