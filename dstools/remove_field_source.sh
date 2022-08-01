#!/bin/bash

source functions.sh

export proj_dir=$1
export pcode=C3369
export ifsel=2
export freq=2100
export mfs=0

files=$(ls $proj_dir/ | grep $ifsel.cal | awk -F . '{print $1}')

echo -e "\e[38;5;166mChoose visibility to image\e[39m"
select file in $files
    do
	if [[ "$REPLY" == stop ]]; then break; fi

	if [[ "$REPLY" == "" ]]; then
	    echo -e "\e[38;5;196m'$REPLY' \e[38;5;166mis not a valid choice.\e[39m"
	    continue
	fi

	export src=$proj_dir/$file

	break
done

# Clear previous files
rm -r $src.sp.*
rm -r $src.imap $src.ibeam $src.irestor $src.imodel $src.iresidual
rm -r $src.$freq.$ifsel.cal.zapped $src.zapped.imap $src.zapped.ibeam $src.zapped.irestor $src.zapped.imodel $src.zapped.iresidual
rm -r $proj_dir/cgcurs*.region

# Initial clean to get image and cell sizes
invert vis=$src.$freq.$ifsel.cal map=$src.imap beam=$src.ibeam robust=0.5 stokes=i options=mfs,double imsize=3,3,beam
clean map=$src.imap beam=$src.ibeam out=$src.imodel options=negstop,positive cutoff=5e-5 niters=10000
restor model=$src.imodel beam=$src.ibeam map=$src.imap out=$src.irestor

if $mfs; then
    # Run statistics to get beam / pixel sizes
    imlist in=$src.irestor options=statistics
    impos in=$src.irestor coord=1,1 type=abspix

    print "Enter image size in pixels (width then height):"
    read imwidth
    read imheight
    print "Enter pix (1,1) pixel coordinate position (width then height)"
    read pixwidth
    read pixheight
    print "Enter pix (1,1) world coordinate position (RA then Dec)"
    read pixra
    read pixdec

    # Calculate imaging parameters
    export imwidth=`expr 3 \* $imwidth`
    export imheight=`expr 3 \* $imheight`
    export cellx=`calc $pixra/$pixwidth`
    export celly=`calc $pixdec/$pixheight`

    invert vis=$src.$freq.$ifsel.cal map=$src.sp.imap beam=$src.sp.ibeam robust=0.5 stokes=i options=mfs,sdb imsize=$imwidth,$imheight cell=$cellx,$celly
    mfclean map=$src.sp.imap beam=$src.sp.ibeam out=$src.sp.imodel cutoff=5e-5 niters=10000 region="relcenter,boxes(-$pixwidth,-$pixheight,$pixwidth,$pixheight)"
    restor model=$src.sp.imodel beam=$src.sp.ibeam map=$src.sp.imap out=$src.sp.irestor
    cgdisp in=$src.sp.irestor type=p device=/xs labtyp=hms,dms
fi

# Select bright off-axis sources to isolate in field model
cgcurs in=$src.imap type=p device=/xs labtyp=hms,dms options=wedge,region range=0,0,heq
immask in=$src.imodel region=@cgcurs.region flag=false logic=not
mv cgcurs.region $proj_dir/cgcurs_mask.$i.region

# Subtract off-axis source from visibilities
uvmodel vis=$src.$freq.$ifsel.cal model=$src.imodel options=subtract,mfs out=$src.$freq.$ifsel.cal.zapped

# Verify removal from residual map
invert vis=$src.$freq.$ifsel.cal.zapped map=$src.zapped.imap beam=$src.zapped.ibeam robust=0.5 stokes=i options=mfs,double imsize=3,3,beam
clean map=$src.zapped.imap beam=$src.zapped.ibeam out=$src.zapped.imodel options=negstop,positive cutoff=5e-5 niters=10000
restor model=$src.zapped.imodel beam=$src.zapped.ibeam map=$src.zapped.imap out=$src.zapped.irestor

cgdisp in=$src.zapped.irestor type=p device=/xs labtyp=hms,dms range=0,0,log

# mv $src.i* $proj_dir/miriad/.
# mv $src.zapped.i* $proj_dir/miriad/.
# mv $proj_dir/*cgcurs* $proj_dir/miriad/.
