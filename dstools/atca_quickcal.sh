#!/bin/bash

source $1/functions.sh

export proj_dir=reduced/$2/miriad/
export data_dir=$(pwd)/$3
export pcode=$4
export refant=$5
export mfinterval=$6
export bpinterval=$7
export gpinterval=$8

# Load data
prompt "Reload data?"
read reload
case $reload in

    [Yy]* )
	print "Flushing atlod files and reloading"
	export reflag="y"

	# Move pre-flagged primary scan back in
	cp -r $proj_dir/*flagged_backup . 2>/dev/null
	rm -rf $proj_dir/ 2>/dev/null
	mkdir -p $proj_dir/
	cd $proj_dir

	# Identify RPFITS files from top-level data directory so that backup scans (e.g. 1934)
	# can sit in subdirectories of the data directory without being auto-imported
	export infiles=$(find -L $data_dir/* -maxdepth 1 -type f | grep $pcode | tr '\n' ',')

	atlod in=$infiles out=$pcode.uv options=birdie,rfiflag,noauto,xycorr,notsys
	uvflag vis=$pcode.uv edge=40 flagval=flag

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
	export reflag="n"
	print "Skipping data reload"
	cd $proj_dir ;;
esac

# Choose frequency
export freqs=$(ls | grep -E '2100|5500|9000|17000' | sed 's/^[^\.]*\.//g' | sort | uniq)
print "Choose frequency / IF"
select f in $freqs; do
    if [[ "$REPLY" == stop ]]; then break; fi

    if [[ "$REPLY" == "" ]]; then
	print "'$REPLY' is not a valid choice"
	continue
    fi

    # TODO: improve this logic
    if [ $(echo $f | grep 2100 | wc -l) -gt 0 ]; then
        export freq=$f
        export spec=2.1
    elif [ $(echo $f | grep 5500 | wc -l) -gt 0 ]; then
        export freq=$f
        export spec=5.5
    elif [ $(echo $f | grep 9000 | wc -l) -gt 0 ]; then
        export freq=$f
        export spec=9
    elif [ $(echo $f | grep 17000 | wc -l) -gt 0 ]; then
        export freq=$f
        export spec=17

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

# Flag known bad channels for this band
flag_channels $pcal.$freq $freq

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
else
    print "Using primary as secondary, skipping calibration copy."
fi

autoflag $scal.$freq

# Gain calibrations for secondary calibrator
gpcal vis=$scal.$freq interval=$gpinterval options="xyvary,qusolve" minants=3 nfbin=4 refant=$refant;

# Transfer flux scale from primary to secondary
if [ $pcal != $scal ]; then
    print "Propagating flux scale to secondary calibrator"
    gpboot vis=$scal.$freq cal=$pcal.$freq;
fi

# Transfer gain calibrations to target
print "Transferring calibration tables to science target"
gpcopy vis=$scal.$freq out=$target.$freq;

# Average gain phase solutions over 2 minutes for better interpolation
print "Averaging 2 minute calibration samples"
gpaver vis=$target.$freq interval=2;

print "Applying calibration to science target"
uvaver vis=$scal.$freq out=$scal.$freq.cal
uvaver vis=$target.$freq out=$target.$freq.cal

# Clean up miriad files before moving to CASA (just keep calibrated target)
export target_dir=$(echo $target | tr '[:lower:]' '[:upper:]')
mkdir ../$target_dir 2>/dev/null
mv $target.$freq.cal ../$target_dir/.

print "DONE!"
