#!/usr/bin/env python

import sys

from casacore.tables import table

ms_file = sys.argv[-1]

tab = table(ms_file, readonly=False, ack=False)

data = tab.getcol("DATA")
data *= 2

tab.putcol("DATA", data)
tab.close()
