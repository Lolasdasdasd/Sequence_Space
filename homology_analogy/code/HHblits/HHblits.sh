#!/bin/bash

# input parameter
if [ $# -ne 2 ]; then
	echo "please provide: <sample number> <start>"
	exit 1
fi

SAMPLE=$1
START=$2

P="/ebio/abt1/lweidmann/seq_space/final/homology/"

DB="/ebio/abt1_share/toolkit_sync/databases/hhblits/uniclust30"
export HHLIB="/ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/"

echo $HHLIB

DIR="$P/HHblits/a3m_${SAMPLE}"
mkdir $DIR

C=0
while read -r line; do
	let C++
	if (( $C <= 3)); then
	  continue
	fi
	#echo $C

	frag=$(echo $line | sed "s/\([A-Z]*\).*/\1/g" )
	pos=$(echo $line | sed "s/[A-Z]* \([0-9]*\)/\1/g" )

	if (( pos < START )); then continue; fi

	FASTA="$DIR/frag_${pos}.fasta"
	echo ">a" > $FASTA
	echo "$frag" >> $FASTA

	OUT="$DIR/frag_${pos}.a3m"
    echo "/ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/hhblits -i $FASTA -d $DB -oa3m $OUT -v 0 -cpu 4 -e 0.000001 -norealign"
    /ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/hhblits -i $FASTA -d $DB -oa3m $OUT -v 0 -cpu 4 -e 0.000001 -norealign

done < $P/sampling/sample$SAMPLE/fragments.csv

