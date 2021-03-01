function flag_extreme {
    vis=$1
    amp=$2
    print "Flagging extreme values (amp > $amp) on $vis"
    uvflag vis=$vis select="amplitude($amp)" flagval=flag
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
