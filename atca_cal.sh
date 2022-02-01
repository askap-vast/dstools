#!/bin/bash

source functions.sh

# Editable Script Parameters
export refant=1
export mfinterval=2.0
export bpinterval=2.0
# --------------------------

export data_dir=$(pwd)/data/
export proj_dir=reduced/$1/miriad/
export pcode=$2

# Load data
prompt "Reload data?"
read reload
case $reload in

    [Yy]* )
	print "Flushing atlod files and reloading"
	export reflag="y"

	# Move pre-flagged primary scan back in
	cp -r $proj_dir/*flagged_backup .
	rm -rf $proj_dir/
	mkdir -p $proj_dir/
	cd $proj_dir

	# Identify RPFITS files from top-level data directory so that backup scans (e.g. 1934)
	# can sit in subdirectories of the data directory without being auto-imported
	export infiles=$(find -L $data_dir/* -maxdepth 1 -type f | grep $pcode | tr '\n' ',')

	atlod in=$infiles out=$pcode.uv options=birdie,rfiflag,noauto,xycorr
	uvflag vis=$pcode.uv edge=40 flagval=flag
	uvsplit vis=$pcode.uv ;;

    [Nn]* )
	export reflag="n"
	print "Skipping data reload"
	cd $proj_dir ;;
esac

# Choose frequency
export freqs=$(ls | grep -E '2100|5500|9000' | sed 's/^[^\.]*\.//g' | sort | uniq)
print "Choose frequency / IF"
select f in $freqs; do
    if [[ "$REPLY" == stop ]]; then break; fi

    if [[ "$REPLY" == "" ]]; then
	print "'$REPLY' is not a valid choice"
	continue
    fi

    if [ $(echo $f | grep 2100 | wc -l) > 0 ]; then
        export freq=$f
        export spec=2.1
    else
        # TODO: add logic here to return 5500 and 9000 MHz files
        export freq=$f
        export spec=5.5
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

prompt "Proceed with primary flagging?"
read reflag

case $reflag in

    [Yy]* )

	# Flag initial extreme values
	# flag_extreme $pcal.$freq 200
	
	# Flag known bad channels for this band
	flag_channels $pcal.$freq $freq

	# Solve for initial bandpass on primary calibrator
	print "Running mfcal on $pcal.$freq"
	mfcal vis=$pcal.$freq interval=$mfinterval,$mfinterval,$bpinterval refant=$refant

	# Iterate between RFI flagging and bandpass solution until spectrum is clean
	while true; do

	    # Run blflag to clean up lower level RFI 
	    print "Running blflag on $pcal.$freq"
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav,nobase
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav
	    # blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav

	    # Solve for bandpass
	    print "Running mfcal on $pcal.$freq"
	    mfcal vis=$pcal.$freq interval=$mfinterval,$mfinterval,$bpinterval refant=$refant

	    # Automatic flagging
	    prompt "Run automatic flagging?"
	    read autoflag
	    case $autoflag in
		[Yy]* )
		    print "Running pgflag on $pcal.$freq"
		    pgflag vis=$pcal.$freq command="<b" device=/xs stokes=xx,yy,xy,yx
		    pgflag vis=$pcal.$freq command="<b" device=/xs stokes=xx,yy,yx,xy
		    print "Running mfcal on $pcal.$freq"
		    mfcal vis=$pcal.$freq interval=$mfinterval,$mfinterval,$bpinterval refant=$refant
		    ;;
		[Nn]* )
		    : ;;
	    esac

	    # Check bandpass calibration
	    # uvspec vis=$pcal.$freq interval=$mfinterval axis=channel,amp device=/xs nxy=1,1
	    # uvspec vis=$pcal.$freq interval=$mfinterval axis=channel,phase device=/xs nxy=1,1

	    prompt "Do more interactive flagging?"
	    read flagmore
	    case $flagmore in
		[Yy]* ) : ;;
		[Nn]* )
		    cp -r $pcal.$freq $pcal.$freq.flagged_backup
		    break ;;
	    esac

	done ;;

    [Nn]* )
	print "Skipping flagging" ;;

esac

prompt "Apply primary gain calibration?"
read primary_cal
case $primary_cal in

    [Yy]* ) 

	# Gain calibrations for primary
	gpcal vis=$pcal.$freq interval=0.1 options=xyvary minants=3 nfbin=16 spec=$spec refant=$refant;

	# Flag outliers in real/imaginary space
	blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=real,imag options=nofqav,nobase
	blflag vis=$pcal.$freq device=/xs stokes=xx,yy,yx,xy axis=real,imag options=nofqav,nobase

	# Check primary calibration
	uvplt vis=$pcal.$freq stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs
	;;

    [Nn]* )
	print "Skipping primary gain calibration" ;;

esac

# ----------------------------
# SECONDARY / GAIN CALIBRATION
# ----------------------------

# Copy calibration to secondary calibrator
if [ $pcal != $scal ]; then
gpcopy vis=$pcal.$freq out=$scal.$freq options=nocal;
    # flag_extreme $scal.$freq 200;
else
    print "Using primary as secondary, skipping calibration copy."
fi

prompt "Proceed with secondary flagging?"
read calibrate

case $calibrate in
	
	[Yy]* )

	    # Iterate between gain solutions and RFI flagging
	    print "Flagging secondary calibrator"
	    while true; do

		# Run blflag on secondary prior to calibrating
		print "Running blflag on $scal.$freq"
		blflag vis=$scal.$freq device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav,nobase
		blflag vis=$scal.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
		blflag vis=$scal.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav
		blflag vis=$scal.$freq device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav

		# Run automatic flagging on secondary
		prompt "Run automatic flagging?"
		read autoflag
		case $autoflag in
		    [Yy]* )
			print "Running pgflag on $scal.$freq"
			pgflag vis=$scal.$freq command="<b" device=/xs stokes=xx,yy,xy,yx
			pgflag vis=$scal.$freq command="<b" device=/xs stokes=xx,yy,yx,xy
			;;
		    [Nn]* ) : ;;
		esac
		
		# Gain calibrations for secondary calibrator
		gpcal vis=$scal.$freq interval=0.1 options="xyvary,qusolve" minants=3 nfbin=4 refant=$refant;

		# Check secondary calibration
		uvplt vis=$scal.$freq stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs;

		prompt "Do more interactive flagging?"
		read flagmore
		case $flagmore in
		    [Yy]* ) : ;;
		    [Nn]* )
			cp -r $scal.$freq $scal.$freq.flagged_backup
			break ;;
		esac
		
	    done ;;
	
	[Nn]* )
	    print "Skipping secondary flagging"

esac

# Transfer flux scale from primary to secondary
if [ $pcal != $scal ]; then
    print "Propagating flux scale to secondary calibrator"
    gpboot vis=$scal.$freq cal=$pcal.$freq;
fi

# Transfer gain calibrations to target
print "Transferring calibration tables to science target"
gpcopy vis=$scal.$freq out=$target.$freq;
# flag_extreme $target.$freq 20

# Average gain phase solutions over 2 minutes for better interpolation
print "Averaging 2 minute calibration samples"
gpaver vis=$target.$freq interval=2;

prompt "Proceed with target flagging?"
read calibrate
case $calibrate in
	
    [Yy]* )

	print "Flagging science target"
	while true; do

	    # Run blflag on science target prior to calibrating
	    print "Running blflag on $target.$freq"
	    # blflag vis=$target.$freq device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav
	    blflag vis=$target.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav

	    # Run automatic flagging on secondary
	    prompt "Run automatic flagging?"
	    read autoflag
	    case $autoflag in
		[Yy]* )
		    print "Running pgflag on $target.$freq"
		    pgflag vis=$target.$freq command="<b" device=/xs stokes=xx,yy,xy,yx
		    pgflag vis=$target.$freq command="<b" device=/xs stokes=xx,yy,xy,yx
		    ;;
		[Nn]* )
		    : ;;
	    esac

	    # Check target calibration
	    # uvplt vis=$target.$freq stokes=xx,yy axis=time,amp options=nofqav device=/xs
	    # uvplt vis=$target.$freq stokes=xx,yy axis=time,phase options=nofqav device=/xs
	    uvplt vis=$target.$freq stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs

	    prompt "Do more interactive flagging?"
	    read flagmore
	    case $flagmore in
		[Yy]* ) : ;;
		[Nn]* )
		    cp -r $target.$freq $target.$freq.flagged_backup
		    break ;;
	    esac
	    
	done
	;;
    
    [Nn]* ) :
	    ;;
    
esac

print "Applying calibration to science target"
uvaver vis=$scal.$freq out=$scal.$freq.cal
uvaver vis=$target.$freq out=$target.$freq.cal

# Clean up miriad files before moving to CASA (just keep calibrated target)
export target_dir=$(echo $target | tr '[:lower:]' '[:upper:]')
mkdir ../$target_dir
mv $target.$freq.cal ../$target_dir/.

print "DONE!"
