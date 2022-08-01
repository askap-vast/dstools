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
    esac
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
