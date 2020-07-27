#!/bin/bash




# input parameter
if [ $# -ne 1 ]; then
	echo "please provide: <.non_homologsCQP>"
fi

#<.non_homologsCQP>
#0 1 0 5010000 0.31 C
#0 1 0 6010000 0.07 P
#0 1 0 9010000 0.64 C
#0 1 0 10010000 0.12 Q
#0 1 0 12010000 0.22 C



cat $1 | sed -n '/ C/p' | sed 's/ C//g' > $1.C
cat $1 | sed -n '/ [QP]/p' | sed 's/ [QP]//g' > $1.QP
echo "disected and written to $1.C and $1.QP"




