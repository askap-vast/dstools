#!/bin/bash

export proj_dir=$1
export pcode=C3369
export data_dir=/import/ada1/jpri6587/phd/atca/C3369/raw_data/OCTS2020/$proj_dir

# Load data
while true; do
    echo -e "\e[38;5;166mReload data? \e[38;5;196m(y/n)\e[39m"
    read reload
    case $reload in
	
	[Yy]* )
	    echo -e "\e[38;5;166mFlushing atlod files and reloading.\e[39m"
	    export reflag="y"
	    export reautoflag="y"
	    rm -rf $proj_dir/
	    mkdir $proj_dir/
	    cd $proj_dir
	    atlod in=`ls $data_dir/*.$pcode | tr '\n' ','` out=$pcode.uv options=birdie,rfiflag,noauto,xycorr
	    uvflag vis=$pcode.uv edge=40 flagval=flag
	    uvsplit vis=$pcode.uv
	    break;;
	
	[Nn]* )
	    export reflag="n"
	    echo -e "\e[38;5;166mSkipping data reload.\e[39m"
	    cd $proj_dir
	    break;;
    esac
done

# Choose from files which sources to use as calibrators and target
files=$(ls | grep 2100.1 | awk -F . '{print $1}')

declare -a orders
orders+=('pcal')
orders+=('scal')
orders+=('target')

declare -A sources
sources['pcal']="primary calibrator"
sources['scal']="secondary calibrator"
sources['target']="science target"

for src in "${orders[@]}"; do
    echo -e "\n\e[38;5;166mChoose ${sources[$src]}\e[39m"
    select file in $files
    do
	if [[ "$REPLY" == stop ]]; then break; fi

	if [[ "$REPLY" == "" ]]
	then
	    echo -e "\e[38;5;196m'$REPLY' \e[38;5;166mis not a valid choice.\e[39m"
	    continue
	fi

	export $src=$file

	break
    done
done

# --------
# FLAGGING
# --------

if [[ $reflag != 'y' ]]; then
    echo -e "\n\e[38;5;166mRe-do flagging? (y/n)\e[39m"
    read reflag;
else
    case $reflag in
	
	[Yy]* )
	    # Flag scan interrupted by wind stow
	    uvflag vis=$target.2100* select="time(21:38:30,22:30:00)" flagval=flag

	    # Flag whole 4-5 baseline
	    for vis in $pcal $scal $target; do
		uvflag vis=$vis.2100* "select=antenna(4)(5)" flagval=flag;
	    done
	    
	    # Initial automatic flagging
	    if [[ $reautoflag != 'y' ]]; then
		echo -e "\e[38;5;166mRe-run automatic flagging? (y/n)\n\e[39m"
		read reautoflag;
	    else 
		case $reautoflag in
		    [Yy]* )
			for vis in $pcal.2100* $scal.2100* $target.2100*;
			do echo -e "\n\e[38;5;166mRunning pgflag on $vis\e[39m"
			   pgflag vis=$vis command="<b" device=/xs stokes=xx,yy,xy,yx
			done
			break;;
		    [Nn]* )
			break;;
		esac
	    fi

	    # Iterate between bandpass solution and RFI flagging until spectrum is clean
	    while true; do

		# Solve for bandpass on primary calibrator
		for i in 1 2;
		do echo -e "\n\e[38;5;166mRunning mfcal on $i\e[39m"
		   mfcal vis=$pcal.2100.$i interval=40.0 refant=1
		done

		# Run blflag to clean up lower level RFI 
		for i in 1 2;
		do echo -e "\n\e[38;5;166mRunning blflag on $pcal.2100.$i\e[39m"
		   blflag vis=$pcal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
		   blflag vis=$pcal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav
		   blflag vis=$pcal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav
		done

		echo -e "\n\e[38;5;166mDo more interactive flagging? \e[38;5;196m(y/n)\e[39m"
		read flagmore
		case $flagmore in
		    [Yy]* ) true ;;
		    [Nn]* ) break ;;
		esac
		
	    done
	    
	    break;;
	
	[Nn]* )
	    echo -e "\e[38;5;166mSkipping flagging\e[39m"
	    break;;
    esac
fi

# -----------
# CALIBRATION
# -----------

echo -e "\n\e[38;5;166mProceed with calibration? \e[38;5;196m(y/n)\e[39m"
read calibrate

case $calibrate in
	
	[Yy]* )
	    echo -e "\n\e[38;5;166mPerforming gain calbration on primary\e[39m"
	    # Gain calibrations for primary
	    for i in 1 2;
	    do gpcal vis=$pcal.2100.$i interval=0.1 options=xyvary minants=3 nfbin=8 spec=2.1 refant=1;
	    done

	    # Check primary calibration
	    for i in 1 2; do
		uvplt vis=$pcal.2100.$i stokes=xx,yy axis=time,amp options=nofqav device=/xs
		uvplt vis=$pcal.2100.$i stokes=xx,yy axis=time,phase options=nofqav device=/xs
		uvplt vis=$pcal.2100.$i stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs
	    done
	    # Copy calibration to secondary calibrator
	    for i in 1 2;
	    do gpcopy vis=$pcal.2100.$i out=$scal.2100.$i;
	    done

	    echo -e "\n\e[38;5;166mFlagging secondary calibrator\e[39m"
	    # Iterate between gain solutions and RFI flagging
	    while true; do

		# Run blflag on secondary prior to calibrating
		for i in 1 2;
		do echo -e "\n\e[38;5;166mRunning blflag on $scal.2100.$i\e[39m"
		   blflag vis=$scal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
		   blflag vis=$scal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav
		   blflag vis=$scal.2100.$i device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav
		done
		
		# Gain calibrations for secondary calibrator
		for i in 1 2;
		do gpcal vis=$scal.2100.1  interval=0.1 options="xyvary,qusolve" minants=3 nfbin=4 refant=1;
		done

		# Check secondary calibration
		uvplt vis=$scal.2100.1 stokes=xx,yy axis=time,amp options=nofqav device=/xs
		uvplt vis=$scal.2100.1 stokes=xx,yy axis=time,phase options=nofqav device=/xs
		uvplt vis=$scal.2100.1 stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs

		echo -e "\n\e[38;5;166mDo more interactive flagging? \e[38;5;196m(y/n)\e[39m"
		read flagmore
		case $flagmore in
		    [Yy]* ) true ;;
		    [Nn]* ) break ;;
		esac
		
	    done

	    echo -e "\n\e[38;5;166mPropagating flux scale to secondary calibrator\e[39m"
	    # Transfer flux scale from primary to secondary
	    for i in 1 2;
	    do gpboot vis=$scal.2100.$i cal=$pcal.2100.$i;
	    done

	    echo -e "\n\e[38;5;166mTransferring calibration tables to science target\e[39m"
	    # Transfer gain calibrations to target
	    for i in 1 2;
	    do gpcopy vis=$scal.2100.$i out=$target.2100.$i;
	    done

	    echo -e "\n\e[38;5;166mAveraging 2 minute calibration samples\e[39m"
	    # Average gain phase solutions over 2 minutes for better interpolation
	    for i in 1 2;
	    do gpaver vis=$target.2100.$i interval=2;
	    done
	    
	    echo -e "\n\e[38;5;166mFlagging science target\e[39m"
	    while true; do
		
		# Run blflag on science target prior to calibrating
		for i in 1 2;
		do echo -e "\n\e[38;5;166mRunning blflag on $target.2100.$i\e[39m"
		   blflag vis=$target.2100.$i device=/xs stokes=xx,yy,xy,yx axis=chan,amp options=nofqav,nobase
		   blflag vis=$target.2100.$i device=/xs stokes=xx,yy,xy,yx axis=time,amp options=nofqav
		   blflag vis=$target.2100.$i device=/xs stokes=xx,yy,xy,yx axis=uvdistance,amp options=nofqav
		done

		# Check target calibration
		uvplt vis=$target.2100.1 stokes=xx,yy axis=time,amp options=nofqav device=/xs
		uvplt vis=$target.2100.1 stokes=xx,yy axis=time,phase options=nofqav device=/xs
		uvplt vis=$target.2100.1 stokes=xx,yy axis=real,imag options=nofqav,nobase,equal device=/xs

		echo -e "\n\e[38;5;166mDo more interactive flagging? \e[38;5;196m(y/n)\e[39m"
		read flagmore
		case $flagmore in
		    [Yy]* ) true ;;
		    [Nn]* ) break ;;
		esac
		
	    done
	    
	    echo -e "\n\e[38;5;166mApplying calibration to science target\e[39m"
	    # Apply calibrations to $target data
	    for i in 1 2;
	    do uvaver vis=$target.2100.$i out=$target.2100.cal.$i;
	       uvaver vis=$scal.2100.$i out=$scal.2100.cal.$i
	    done

	    break;;
	
	[Nn]* )
	    echo -e "\e[38;5;166mExiting prior to calibration\e[39m"
	    exit 0
	    break;;
    esac

