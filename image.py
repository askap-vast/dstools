#!/usr/bin/env python

import os
import sys

proj_dir = sys.argv[-1]
if proj_dir[-1] != '/':
    proj_dir += '/'
source = proj_dir.split('/')[0].lower()

os.system("rm -r {}{}.ms*")
os.system("rm -r {}{}.1.im*".format(proj_dir, source))
os.system("rm -r {}{}.2.im*".format(proj_dir, source))
os.system("rm -r {}{}_phase_selfcal*".format(proj_dir, source))
os.system("mkdir {}diagnostic/".format(proj_dir))

    
for int_freq in [1,2]:
    mirfile = '{}.2100.cal.{}'.format(source, int_freq)
    msname = '{}.ms'.format(source)

    # DATA IMPORT / FLAGGING
    reimport = raw_input("Redo data import and flagging? (y/n)\n")
    if reimport == 'y':
        if os.path.exists(proj_dir + msname):
            os.system('rm -r {}{}'.format(proj_dir, msname))

        importmiriad(mirfile='{}{}'.format(proj_dir, mirfile),
                     vis='{}{}'.format(proj_dir, msname))

        # Run some additional flagging
        flagdata(vis=proj_dir+msname, mode='manual', antenna='3&4')
        flagdata(vis=proj_dir+msname, mode='rflag')
        flagdata(vis=proj_dir+msname, mode='tfcrop')
        flagdata(vis=proj_dir+msname, mode='tfcrop')


    # Run selfcal
    i = 0
    cont = 'y'
    while cont == 'y':
        tclean(vis=proj_dir+msname,
               field='0', 
               cell=['1.0arcsec', '1.0arcsec'], 
               imsize=[3024, 3024], 
               savemodel='modelcolumn', 
               threshold=5.5E-5, 
               niter=1000, 
               imagename='{}{}.{}.im_I'.format(proj_dir, source, int_freq),
               nterms=2, 
               deconvolver='mtmfs',
               reffreq='2099.999MHz', 
               mask=proj_dir+'{}.mask'.format(source),
               weighting='briggs', 
               robust=0.5, 
               interactive=True,
               pblimit=0.1)

        cal_table = '{}diagnostic/{}.{}.phase_selfcal.{}.cal'.format(proj_dir, source, int_freq, i)
        gaincal(vis=proj_dir+msname,
                caltable=cal_table,
                solint="10 min",
                calmode="p",
                gaintype="G")
        plotms(vis=cal_table,
               xaxis="time",
               yaxis="phase",
               plotrange=[0, 0, -30, 30],
               plotfile='{}{}.{}.phase_selfcal.{}.png'.format(proj_dir, msname, int_freq, i),
               showgui=True)
        applycal(vis=proj_dir+msname, gaintable=[cal_table], interp="linear")   

        # ans = raw_input("Proceed with more selfcal (y/n)? ")
        # while ans not in ['y', 'n']:
        #     if cont == 'n':
        #         break
        #     elif cont == 'y':
        #         i += 1
        #         break

        cont = 'n'

    # Copy mask to separate location for reuse and flush image data
    mask_loc = "{}/{}.mask".format(proj_dir, source)
    if not os.path.exists(mask_loc):
        os.system("cp -r {}_im_I.mask {}".format(proj_dir+source, mask_loc))

    # Clean with 20000 iterations
    tclean(vis=proj_dir+msname,
           field='0', 
           cell=['1.0arcsec', '1.0arcsec'], 
           imsize=[3024, 3024], 
           savemodel='modelcolumn', 
           threshold=5.5E-5, 
           niter=20000, 
           imagename='{}{}.{}.im_I'.format(proj_dir, source, int_freq),
           nterms=2, 
           deconvolver='mtmfs',
           reffreq='2099.999MHz', 
           mask='',
           weighting='briggs', 
           robust=0.5, 
           pblimit=0.1)

    for tt in ['tt0', 'tt1']:

        exportfits(imagename='{}{}.{}.im_I.image.{}'.format(proj_dir, source, int_freq, tt),
                   fitsimage='{}{}.{}.im_I.{}.fits'.format(proj_dir, source, int_freq, tt),
                   overwrite=True)



    tclean(vis=proj_dir+msname,
           field='0',
           cell=['1.0arcsec', '1.0arcsec'],
           imsize=[3024, 3024],
           savemodel='none',
           threshold=5.5E-5,
           niter=0,
           imagename='{}{}.{}.im_V'.format(proj_dir, source, int_freq),
           nterms=2,
           deconvolver='mtmfs',
           reffreq='2099.999MHz',
           stokes='V',
           weighting='briggs',
           robust=0.5,
           pblimit=0.1)

    for tt in ['tt0', 'tt1']:

        exportfits(imagename='{}{}.{}.im_V.image.{}'.format(proj_dir, source, int_freq, tt),
                   fitsimage='{}{}.{}.im_V.{}.fits'.format(proj_dir, source, int_freq, tt),
                   overwrite=True)
