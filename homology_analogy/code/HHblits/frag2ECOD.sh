#!/bin/bash

# input parameter
if [ $# -ne 1 ]; then
	echo "please provide: <ONCE CALLED ONLY on OLT> "
	exit 1
fi

HHSEARCH="/ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/hhsearch"
DB="/ebio/abt1_share/toolkit_sync/databases/hh-suite/ECOD/ECOD_F70_20190225"




C=0
for VAR in 0 1 2 3 4 5 6 7 8 9; do
	rm a3m_$VAR/*frag2ecod* #remove old results
	for f in a3m_$VAR/frag_*.a3m; do
		echo "AT FRAGMENT $C"
		A3M=$f
		BASE=${f%.fasta}.frag2ecod

		echo $HHSEARCH -i $A3M -d $DB -norealign -o $BASE.hrr #changed to norealign on June3
		$HHSEARCH -i $A3M -d $DB -norealign -o $BASE.hrr

		tail -n +10 $BASE.hrr | sed -n '/>ECOD/p' > $BASE.ecodname
		tail -n +10 $BASE.hrr | head -n 250 | sed -n '/^[ 0-9]* ECOD/p' > $BASE.list

		C=$((C+1))
	done
done

#
#	echo a3m_$VAR/frag_*.fasta
#
#	BASE="test.frag2ecod"
#	
#	/ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/hhsearch -i test.a3m -d /ebio/abt1_share/toolkit_sync/databases/hh-suite/ECOD/ECOD_F70_20190225 -o $BASE.hrr
#	tail -n +10 $BASE.hrr | sed -n '/>ECOD/p' > $BASE.ecodname
#	tail -n +10 $BASE.hrr | head -n 250 | sed -n '/^[ 0-9]* ECOD/p' > $BASE.list
#
#	echo $VAR
#done
#
#
#for i in `seq 0 5`;
#do
#	echo "$i"
#done 


