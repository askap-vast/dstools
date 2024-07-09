function flag_extreme {
    vis=$1
    amp=$2
    print "Flagging extreme values (amp > $amp) on $vis"
    uvflag vis=$vis select="amplitude($amp)" flagval=flag
}

function flag_channels {
    vis=$1
    band=$2
    print "Flagging known bad channels in band $band"
    case $band in
	[2100]* )

	    uvflag vis=$vis line=channel,60,435,1,1 flagval=flag
	    uvflag vis=$vis line=channel,5,1300,1,1 flagval=flag

	    ;;

	[5500]* )
	    ;;

	[9000]* )
	    ;;

	[17000]* )
	    ;;

    esac
}

function autoflag {

    vis=$1

    print "Running pgflag on $vis"
    pgflag vis=$vis command="<b" device=/xs stokes=xx,yy,xy,yx options=nodisp
    pgflag vis=$vis command="<b" device=/xs stokes=xx,yy,yx,xy options=nodisp
}

function cal_bandpass {

    vis=$1
    print "Running mfcal on $vis"

    mfcal vis=$vis interval=$mfinterval,$mfinterval,$bpinterval refant=$refant
}

function cal_flux {

    vis=$1

    gpcal vis=$vis interval=$gpinterval options=xyvary minants=3 nfbin=16 spec=$spec refant=$refant;

}

function print {
    msg=$1
    msgcolor="\e[38;5;166m"
    echo -e "$msgcolor$msg \e[39m"
}
function prompt {
    msg=$1
    promptcolor="\e[38;5;196m"
    print "$msg $promptcolor(y/n)"
}
