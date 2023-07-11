while read line; do ./testPhiANF $line; done <$1 >`basename $1`.out

