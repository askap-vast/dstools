#!/usr/bin/env python

import os
import sys
import glob

def colored(msg):

    return '\033[91m{}\033[0m'.format(msg)

def prompt(msg):

    msg = msg[:-1] + ' (y/n)?\n'
    
    resp = raw_input(colored(msg))
    if resp not in ['y', 'n']:
        resp = raw_input(colored(msg))

    return resp
        

proj_dir = sys.argv[-1]
if proj_dir[-1] != '/':
    proj_dir += '/'

if 'obs' in proj_dir:
    source = proj_dir.split('_')[0].lower()
else:
    source = proj_dir.split('/')[0].lower()


imsize = 3584
masksize = 25
iterations = 2000
interactive = True
cellsize = ['1.0arcsec']
clean_scales = [0, 3, 8]
pblim = 0.01
threshold = 8E-5

clean_mask = '{}{}.mask'.format(proj_dir, source)

for int_freq in [1]:
    mirfile = '{}.2100.{}.cal'.format(source, int_freq)
    msname = '{}.ms'.format(source)

    # Check if version exists with field sources removed
    if os.path.exists(proj_dir + mirfile + '.zapped'):
        mirfile += '.zapped'

    # DATA IMPORT / FLAGGING
    if os.path.exists(proj_dir + msname):
        reimport = prompt('Redo data import?')

        if reimport == 'y':
            os.system('rm -r {}{}'.format(proj_dir, msname))
            os.system('rm -r {}{}.flagversions'.format(proj_dir, msname))
            os.system('rm -r {}selfcal/'.format(proj_dir))
            os.system('mkdir {}selfcal/'.format(proj_dir))

            importmiriad(mirfile='{}{}'.format(proj_dir, mirfile),
                         vis='{}{}'.format(proj_dir, msname))
    else:
        importmiriad(mirfile='{}{}'.format(proj_dir, mirfile),
                     vis='{}{}'.format(proj_dir, msname))


    reflag = prompt('Do automatic flagging?')
        
    if reflag == 'y':
        # Run some additional flagging
        # flagdata(vis=proj_dir + msname, mode='manual', antenna='3&4')
        flagdata(vis=proj_dir + msname, mode='rflag')
        flagdata(vis=proj_dir + msname, mode='tfcrop')
        flagdata(vis=proj_dir + msname, mode='tfcrop')
        

    test_params = prompt('Experiment with imsize / flagging?')

    # Check if default imsize is acceptable (e.g. no strong out of field sources)
    if test_params == 'y':
        os.system('mkdir {}test_params'.format(proj_dir))
        while True:
        
            tclean(vis=proj_dir + msname,
                   field='0',
                   cell=cellsize,
                   imsize=[imsize],
                   savemodel='modelcolumn',
                   threshold=threshold,
                   niter=2000,
                   imagename='{}test_params/test_im'.format(proj_dir),
                   nterms=2,
                   deconvolver='mtmfs',
                   scales=clean_scales,
                   reffreq='2099.999MHz',
                   weighting='briggs',
                   robust=0.5,
                   stokes='IQUV',
                   interactive=interactive,
                   pblimit=pblim)

            accept_imsize = prompt('Is current imsize of {} and pblim of {} acceptable?'.format(imsize, pblim))
            # flag_cmd = input('Enter flagging command (empty to skip): ')
            # if flag_cmd != '':
            #     try:
            #         eval(flag_cmd)
            #     except Exception:
            #         print('{ } \n ... is not a valid flagging command.'.format(flag_cmd))
                
            if accept_imsize == 'y':
                break
            else:
                imsize = int(raw_input("Enter new imsize: "))
                pblim = float(raw_input("Enter new pblim: "))
    os.system('rm -r {}test_params'.format(proj_dir))
    
    selfcal = prompt('Run selfcal?')

    if selfcal == 'y':
        
        # Run selfcal
        i = 0
        i = len(glob.glob('{}/{}.*.ms'.format(proj_dir, source)))
        cont = 'y'
        while cont == 'y':
            i += 1
            selfcal_template = '{}selfcal/{}.{}.selfcal_im{}_I'
            selfcal_im = selfcal_template.format(proj_dir, source, int_freq, i)

            # Use freshly split ms each iteration after first loop
            selfcal_ms = proj_dir + msname.replace('.ms', '.selfcal.{}.ms'.format(i-1))
            selfcal_ms = selfcal_ms if os.path.exists(selfcal_ms) else proj_dir + msname

            # Set mask from previous iteration as input to next tclean loop
            selfcal_mask = selfcal_template.format(proj_dir, source, int_freq, i-1) + '.mask'

            # If no selfcal mask exists, check for a backed up mask in main directory
            if os.path.exists(selfcal_mask):
                init_mask = selfcal_mask
            else:
                init_mask = clean_mask if os.path.exists(clean_mask) else ''

            tclean(vis=selfcal_ms,
                   field='0',
                   cell=cellsize,
                   imsize=[imsize],
                   savemodel='modelcolumn',
                   threshold=threshold,
                   niter=2000,
                   imagename=selfcal_im,
                   nterms=2,
                   deconvolver='mtmfs',
                   scales=clean_scales,
                   reffreq='2099.999MHz',
                   weighting='briggs',
                   robust=0.5,
                   stokes='IQUV',
                   interactive=interactive,
                   mask=init_mask,
                   pblimit=pblim)

            # Trial self-cal solution intervals
            while True:
                interval = raw_input("Select solution interval (in min/s): ")
                try:
                    unit = 'min' if 'min' in interval else 's' if 's' in interval else ''
                    num = int(interval.replace(unit, ''))
                except ValueError:
                    print("Invalid solution interval entered, must be format <int>[min/s].")
                    continue

                # Save self-cal plots
                cal_file = '{}selfcal/{}.{}.phase_selfcal_{}.{}'.format(proj_dir, source, int_freq, interval, i)
                cal_table = cal_file + '.cal'
                cal_plot = cal_file + '.png'
                gaincal(vis=selfcal_ms,
                        caltable=cal_table,
                        solint=interval,
                        calmode='p',
                        gaintype='G')
                plotms(vis=cal_table,
                       xaxis='time',
                       yaxis='phase',
                       plotrange=[0, 0, -30, 30],
                       plotfile=cal_plot,
                       showgui=True)

                # Confirm solution is good before applying, else trial another solution
                cal_good = prompt('Is self-cal a good solution?')

                if cal_good == 'y':
                    applycal(vis=selfcal_ms, gaintable=[cal_table], interp='linear')
                    split(vis=selfcal_ms,
                          outputvis=proj_dir + msname.replace('.ms', '.selfcal.{}.ms'.format(i)),
                          datacolumn="corrected")
                    break
                else:
                    os.system("rm -r {} {}".format(cal_table, cal_plot))

            # Break from loop once solution is sufficient
            cont = prompt('Proceed with more selfcal?')

        # Backup clean mask to main project directory
        backup_mask = selfcal_template.format(proj_dir, source, int_freq, i) + '.mask'
        if os.path.exists(clean_mask):
            replace_mask = prompt('Update the existing clean mask with current mask?')

            if replace_mask == 'y':
                os.system('cp -r {}.mask {}'.format(backup_mask, clean_mask))
        else:
            os.system('cp -r {}.mask {}'.format(backup_mask, clean_mask))

    elif selfcal == 'n':
        i = len(glob.glob("{}/*selfcal.*.ms".format(proj_dir)))
        selfcal_ms = proj_dir + msname.replace('.ms', '.selfcal.{}.ms'.format(i))
        
    field_model_path = '{}field_model'.format(proj_dir)
    if os.path.exists(field_model_path):
        os.system('rm -r {}'.format(field_model_path))
    os.system('mkdir {}'.format(field_model_path))

    # Deep clean to produce field model
    tclean(vis=selfcal_ms,
           field='0',
           cell=cellsize,
           imsize=[imsize],
           threshold=threshold,
           niter=iterations*5,
           imagename='{}/{}.{}.im_deep'.format(field_model_path, source, int_freq),
           nterms=2,
           deconvolver='mtmfs',
           scales=clean_scales,
           reffreq='2099.999MHz',
           weighting='briggs',
           stokes='IQUV',
           robust=0.5,
           mask=clean_mask,
           interactive=interactive,
           pblimit=pblim)

    # Mask out the source
    circle_cen = 'circle[[{}pix, {}pix], {}pix]'.format(imsize // 2, imsize // 2, masksize*3)
    model = '{}/{}.{}.im_deep.model'.format(field_model_path, source, int_freq)
    bgmodel = '{}/{}.{}.im_deep.bgmodel'.format(field_model_path, source, int_freq)

    os.system('rm -r {}*'.format(bgmodel))
    for tt in ['tt0', 'tt1']:
        mask = '{}/{}.mask.{}'.format(field_model_path, source, tt)
        modelfile = '{}.{}'.format(model, tt)
        makemask(mode='copy', inpimage=modelfile, output=mask, inpmask=circle_cen, overwrite=True)
        immath(imagename=[modelfile, mask], outfile=bgmodel + '.' + tt, expr='IM0*(1-IM1)')

    os.system('rm -r {}/{}.{}.im_presub*'.format(field_model_path, source, int_freq))
    # Insert masked background model into visibilities and subtract
    tclean(vis=selfcal_ms,
           field='0',
           cell=cellsize,
           imsize=[imsize],
           startmodel=[bgmodel + '.tt0', bgmodel + '.tt1'],
           savemodel='modelcolumn',
           niter=0,
           imagename='{}/{}.{}.im_presub'.format(field_model_path, source, int_freq),
           nterms=2,
           deconvolver='mtmfs',
           scales=clean_scales,
           reffreq='2099.999MHz',
           weighting='briggs',
           stokes='IQUV',
           robust=0.5,
           pblimit=pblim)

    uvsub(vis=selfcal_ms)

    os.system('rm -r {}/{}.{}.im_subbed*'.format(field_model_path, source, int_freq))
    # Reimage to confirm field subtraction
    tclean(vis=selfcal_ms,
           field='0',
           datacolumn='corrected',
           cell=cellsize,
           imsize=[imsize],
           threshold=threshold,
           niter=iterations // 4,
           imagename='{}/{}.{}.im_subbed'.format(field_model_path, source, int_freq),
           nterms=2,
           deconvolver='mtmfs',
           scales=clean_scales,
           reffreq='2099.999MHz',
           weighting='briggs',
           stokes='IQUV',
           mask=clean_mask,
           robust=0.5,
           pblimit=pblim)

    for tt in ['tt0', 'tt1']:
        for stokes in ['I', 'V']:
            for imtype in ['deep', 'subbed']:

                os.system('rm -r {}/{}.{}.im_{}.{}.image*'.format(field_model_path, source, int_freq, imtype, stokes))
                
                imsubimage(imagename='{}/{}.{}.im_{}.image.{}/'.format(field_model_path,
                                                                       source,
                                                                       int_freq,
                                                                       imtype,
                                                                       tt),
                           outfile='{}/{}.{}.im_{}.{}.image.{}'.format(field_model_path,
                                                                       source,
                                                                       int_freq,
                                                                       imtype,
                                                                       stokes,
                                                                       tt),
                           stokes=stokes)
                exportfits(imagename='{}/{}.{}.im_{}.{}.image.{}'.format(field_model_path,
                                                                         source,
                                                                         int_freq,
                                                                         imtype,
                                                                         stokes,
                                                                         tt),
                           fitsimage='{}/{}.{}.im_{}.{}.{}.fits'.format(field_model_path,
                                                                        source,
                                                                        int_freq,
                                                                        stokes,
                                                                        imtype,
                                                                        tt),
                           overwrite=True)
