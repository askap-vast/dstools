#!/usr/bin/env python

def colored(msg):

    return '\033[91m{}\033[0m'.format(msg)

msname = sys.argv[-2]
phasecen = sys.argv[-1]

if msname.endswith('cal'):
    print(colored("Converting miriad UVFITS to CASA MeasurementSet"))
    mirfile = msname
    proj_dir = msname.split('/')[0]
    msname = '{}/{}.ms'.format(proj_dir, proj_dir.lower())
    importmiriad(mirfile=mirfile, vis=msname)

rotated_ms = msname.replace('.ms', '.rotated.target.ms')
fixvis(vis=msname, outputvis=rotated_ms, phasecenter=phasecen)
