#!/bin/bash

source $1/functions.sh

proj_dir=reduced/$2/miriad/
data_dir=$(pwd)/$3
pcode=$4
refant=$5
mfinterval=$6
bpinterval=$7
gpinterval=$8
noflag=$9

# Load data
prompt "Reload data?"
read reload
case $reload in

    [Yy]* )
	print "Flushing atlod files and reloading"

	mkdir -p $proj_dir/
	cd $proj_dir

	if $noflag; then
	    atlod_options=noauto,xycorr,notsys
	else
	    atlod_options=birdie,rfiflag,noauto,xycorr,notsys
	fi

	# Identify RPFITS files from top-level data directory so that backup scans (e.g. 1934)
	# can sit in subdirectories of the data directory without being auto-imported
	infiles=$(find -L $data_dir/* -maxdepth 1 -type f | grep $pcode | tr '\n' ',')

	atlod in=$infiles out=$pcode.uv options=$atlod_options

	# Optionally shift phasecenter. This is to be used when you have offset the phasecenter
	# during an observation (e.g. by -120 arcsec in declination) to avoid DC correlator
	# errors. Correction would be to shift by +120 arcsec here.
	print "Shift phasecenter?"
	read fix_phasecenter
	case $fix_phasecenter in

	    [Yy]* )
		print "Shift in right ascension (arcsec):"
		read shiftra
		print "Shift in declination (arcsec):"
		read shiftdec

		uvedit vis=$pcode.uv ra=$shiftra dec=$shiftdec out=$pcode.fix.uv
		rm -r $pcode.uv
		mv $pcode.fix.uv $pcode.uv

		;;

	    [Nn]* )
	        ;;

	esac

	uvsplit vis=$pcode.uv ;;

    [Nn]* )
	print "Skipping data reload"
	cd $proj_dir ;;
esac

# Choose frequency
freqs=$(ls | grep -E '2100|5500|9000|17000' | sed 's/^[^\.]*\.//g' | sort | uniq)
print "Choose frequency / IF"
select f in $freqs; do
    if [[ "$REPLY" == stop ]]; then break; fi

    if [[ "$REPLY" == "" ]]; then
	print "'$REPLY' is not a valid choice"
	continue
    fi

    # TODO: improve this logic
    if [ $(echo $f | grep 2100 | wc -l) -gt 0 ]; then
        freq=$f
        spec=2.1
	band=AT_L
    elif [ $(echo $f | grep 5500 | wc -l) -gt 0 ]; then
        freq=$f
        spec=5.5
	band=AT_C
    elif [ $(echo $f | grep 9000 | wc -l) -gt 0 ]; then
        freq=$f
        spec=9
	band=AT_X
    elif [ $(echo $f | grep 17000 | wc -l) -gt 0 ]; then
        freq=$f
        spec=17
	band=AT_K

    fi
    break
done

# Choose from files which sources to use as calibrators and target
files=$(ls | grep $freq | awk -F . '{print $1}' | uniq)

declare -a orders
orders+=('pcal')
orders+=('scal')
orders+=('target')

declare -A sources
sources['pcal']="primary calibrator"
sources['scal']="secondary calibrator"
sources['target']="science target"

for src in "${orders[@]}"; do
    print "Choose ${sources[$src]}"
    select file in $files
    do
	if [[ "$REPLY" == stop ]]; then break; fi

	if [[ "$REPLY" == "" ]]; then
	    print "'$REPLY' is not a valid choice"
	    continue
	fi

	export $src=$file

	break

    done
done

# ------------------------------
# PRIMARY / BANDPASS CALIBRATION
# ------------------------------

# Solve for initial bandpass on primary calibrator
cal_bandpass $pcal.$freq
autoflag $pcal.$freq 
cal_bandpass $pcal.$freq 

# Gain calibrations for primary
cal_flux $pcal.$freq

# ----------------------------
# SECONDARY / GAIN CALIBRATION
# ----------------------------

# Copy calibration to secondary calibrator
if [ $pcal != $scal ]; then
    gpcopy vis=$pcal.$freq out=$scal.$freq options=nocal;
    autoflag $scal.$freq
else
    print "Using primary as secondary, skipping calibration copy."
fi

# Gain calibrations for secondary calibrator
gpcal vis=$scal.$freq interval=$gpinterval options="xyvary,qusolve" minants=3 nfbin=4 refant=$refant;

# Transfer flux scale from primary to secondary
if [ $pcal != $scal ]; then
    print "Propagating flux scale to secondary calibrator."
    gpboot vis=$scal.$freq cal=$pcal.$freq;
fi

# Transfer gain calibrations to target
print "Transferring calibration tables to science target."
gpcopy vis=$scal.$freq out=$target.$freq;

# Average gain phase solutions over 2 minutes for better interpolation
print "Averaging 2 minute calibration samples."
gpaver vis=$target.$freq interval=2;

print "Applying calibration to science target."
uvaver vis=$scal.$freq out=$scal.$freq.cal
uvaver vis=$target.$freq out=$target.$freq.cal

# Clean up miriad files before moving to CASA (just keep calibrated target in MS format)
target_dir=$(echo $target | tr '[:lower:]' '[:upper:]')
mkdir ../$target_dir 2>/dev/null

fits in=$target.$freq.cal out=$target.$freq.cal.fits op=uvout
msfile=../$target_dir/$target.$band.ms
casa --nologger --nologfile -c "importuvfits(fitsfile='$target.$freq.cal.fits', vis='$msfile')" 1>/dev/null

print "DONE!"
