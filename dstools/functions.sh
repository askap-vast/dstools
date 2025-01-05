function load_data {

    proj_dir=$1
    data_dir=$2
    noflag=$3
    pcode=$4
    shiftra=$5
    shiftdec=$6

    mkdir -p $proj_dir/miriad
    cd $proj_dir/miriad
    
    if $noflag; then
	atlod_options=noauto,xycorr,notsys
    else
	atlod_options=birdie,rfiflag,noauto,xycorr,notsys
    fi
    
    # Identify RPFITS files from top-level data directory so that backup scans (e.g. 1934)
    # can sit in subdirectories of the data directory without being auto-imported
    infiles=$(find -L $data_dir/* -maxdepth 1 -type f | grep $pcode | tr '\n' ',' | head -c -1)
    
    atlod in=$infiles out=$pcode.uv ifsel=$ifsel options=$atlod_options
    
    # Optionally shift phasecenter. This is to be used when you have offset the phasecenter
    # during an observation (e.g. by -120 arcsec in declination) to avoid DC correlator
    # errors. Correction would be to shift by +120 arcsec here.
    if [[ $shiftra -ne 0 ]] || [[ $shiftdec -ne 0 ]]; then
	echo "Shifting phasecentre by RA=$shiftra Dec=$shiftdec"
	uvedit vis=$pcode.uv ra=$shiftra dec=$shiftdec out=$pcode.fix.uv
	rm -r $pcode.uv
	mv $pcode.fix.uv $pcode.uv
    fi
    
    uvsplit vis=$pcode.uv
    
}

function uvtofits {

    uv=$1
    fits=$2

    fits in=$uv out=$fits op=uvout

}

function manflag {

    vis=$1
    x=$2
    y=$3
    options=$4

    blflag vis=$vis device=/xs stokes=xx,yy,xy,yx axis=$x,$y options=$options

}

function autoflag {

    vis=$1

    pgflag vis=$vis command="<b" device=/xs stokes=xx,yy,xy,yx options=nodisp
    pgflag vis=$vis command="<b" device=/xs stokes=xx,yy,yx,xy options=nodisp

}

function flag_timerange {

    vis=$1
    start_time=$2
    end_time=$3

    uvflag vis=$vis select="time($start_time,$end_time)" flagval=flag

}

function cal_bandpass {

    vis=$1
    interpolate=$2

    if $interpolate; then
	mfcal vis=$vis interval=$mfinterval,$mfinterval,$bpinterval options=interpolate refant=$refant;
    else
	mfcal vis=$vis interval=$mfinterval,$mfinterval,$bpinterval refant=$refant;
    fi

}

function cal_gains {

    vis=$1
    options=$2

    gpcal vis=$vis interval=$gpinterval options=$options minants=3 nfbin=$nfbin spec=$spec refant=$refant;
} 

function copy_cal {

    pcal=$1
    scal=$2

    gpcopy vis=$pcal out=$scal
}

function bootstrap {

    scal=$1
    pcal=$2

    gpboot vis=$scal cal=$pcal;
}

function average_gains {

    cal=$1

    gpaver vis=$cal interval=2;

}

function apply_gains {

    cal=$1

    uvaver vis=$cal out=$cal.cal

}

