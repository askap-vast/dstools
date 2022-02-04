#!/bin/bash

source functions.sh

export proj_dir=$1
export pcode=C3369
export altpcode=$2
export process_dir=/import/ada1/jpri6587/phd/atca/$pcode
export data_dir=$(find $process_dir/raw_data/*/ | grep $proj_dir$)

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
	mkdir $proj_dir/
	cd $proj_dir

	export infiles=$(ls $data_dir/* | grep $pcode | tr '\n' ',')

	atlod in=$infiles out=$pcode.uv options=birdie,rfiflag,noauto,xycorr
	uvflag vis=$pcode.uv edge=40 flagval=flag
	uvsplit vis=$pcode.uv ;;

    [Nn]* )
	export reflag="n"
	print "Skipping data reload"
	cd $proj_dir ;;
esac

# Choose frequency
export freqs=$(ls | grep 00 | awk -F . '{print $2"."$3}' | sort | uniq)
print "Choose frequency / IF"
select f in $freqs; do
    if [[ "$REPLY" == stop ]]; then break; fi

    if [[ "$REPLY" == "" ]]; then
	print "'$REPLY' is not a valid choice"
	continue
    fi

    export freq=$f
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

# Restore pre-flagged versions
for vis in $pcal $scal; do
    cp -r ../$vis.$freq.flagged_backup $vis.$freq
done

# ------------------------------
# PRIMARY / BANDPASS CALIBRATION
# ------------------------------

prompt "Proceed with primary flagging?"
read reflag

case $reflag in

    [Yy]* )

	# Flag initial extreme values
	flag_extreme $pcal.$freq 200

	print "Running mfcal on $pcal.$freq"
	# Solve for initial bandpass on primary calibrator
	mfcal vis=$pcal.$freq interval=40.0 refant=1

	# Iterate between RFI flagging and bandpass solution until spectrum is clean
	while true; do

	    # Run blflag to clean up lower level RFI 
	    print "Running blflag on $pcal.$freq"
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav,nobase
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
	    blflag vis=$pcal.$freq device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav

	    # Solve for bandpass
	    print "Running mfcal on $pcal.$freq"
	    mfcal vis=$pcal.$freq interval=40.0 refant=1

	    # Automatic flagging
	    prompt "Run automatic flagging?"
	    read autoflag
	    case $autoflag in
		[Yy]* )
		    print "Running pgflag on $pcal.$freq"
		    pgflag vis=$pcal.$freq command="<b" device=/xs stokes=xx,yy,xy,yx
		    print "Running mfcal on $pcal.$freq"
		    mfcal vis=$pcal.$freq interval=40.0 refant=1
		    ;;
		[Nn]* )
		    : ;;
	    esac


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
	gpcal vis=$pcal.$freq interval=0.1 options=xyvary minants=3 nfbin=8 spec=2.1 refant=1;

	# Check primary calibration
	uvplt vis=$pcal.$freq stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs
	;;

    [Nn]* )
	print "Skipping primary gain calibration"
	;;
esac

# ----------------------------
# SECONDARY / GAIN CALIBRATION
# ----------------------------

# Copy calibration to secondary calibrator
gpcopy vis=$pcal.$freq out=$scal.$freq options=nocal;
flag_extreme $scal.$freq 200

prompt "Proceed with secondary flagging?"
read calibrate

case $calibrate in
	
	[Yy]* )

	    print "Flagging secondary calibrator"
	    # Iterate between gain solutions and RFI flagging
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
			;;
		    [Nn]* ) : ;;
		esac
		
		# Gain calibrations for secondary calibrator
		gpcal vis=$scal.$freq interval=0.1 options="xyvary,qusolve" minants=3 nfbin=4 refant=1;

		# Check secondary calibration
		# uvplt vis=$scal.$freq stokes=xx,yy axis=time,amp options=nofqav device=/xs
		# uvplt vis=$scal.$freq stokes=xx,yy axis=time,phase options=nofqav device=/xs
		uvplt vis=$scal.$freq stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs

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

print "Propagating flux scale to secondary calibrator"
# Transfer flux scale from primary to secondary
gpboot vis=$scal.$freq cal=$pcal.$freq;

print "Transferring calibration tables to science target"
# Transfer gain calibrations to target
gpcopy vis=$scal.$freq out=$target.$freq;
flag_extreme $target.$freq 20

print "Averaging 2 minute calibration samples"
# Average gain phase solutions over 2 minutes for better interpolation
gpaver vis=$target.$freq interval=2;

prompt "Proceed with target flagging?"
read calibrate
case $calibrate in
	
    [Yy]* )

	print "Flagging science target"
	while true; do

	    # Run blflag on science target prior to calibrating
	    print "Running blflag on $target.$freq"
	    blflag vis=$target.$freq device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav
	    blflag vis=$target.$freq device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav

	    # Run automatic flagging on secondary
	    prompt "Run automatic flagging?"
	    read autoflag
	    case $autoflag in
		[Yy]* )
		    print "Running pgflag on $target.$freq"
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
mkdir $proj_dir/miriad
mv $proj_dir/* $proj_dir/miriad/.
mv $proj_dir/miriad/$target.$freq.cal $proj_dir/.

print "DONE!"
