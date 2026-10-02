{
        if (NR==1)            { printflag=0 }
	if ($0 ~ /NAMELIST-/) { printflag=1 }
	if (printflag)        { print $0    }
	if ($0 ~ /END-/)      { printflag=0 }
} 
