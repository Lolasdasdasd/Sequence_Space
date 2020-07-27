#!/bin/bash

if [ $# -ne 1 ]; then
  echo "please provide: <homologs    / non_homologsCQP.C   >"
  echo "        format: 0 1 3000000 778010000 88.71 (B1 B2 P1 P2 Prob)"
  echo ""
  echo "  ---> outputs:  homologs.NW / non_homologsCQP.C.NW"
  exit 1
fi

#batch1 batch2 pos1 pos2 prob
#1 2 5010000 16020000 13.30
#1 4 5010000 5040000 40.47
#1 4 5010000 99040000 99.19
#1 4 9010000 54040000 19.79
#1 4 9010000 79040000 11.58
P="/ebio/abt1/lweidmann/seq_space/final/homology/"

#echo "# $1"
#echo "B1 B2 P1 P2 frag1 frag2 prob HD NWSn NWMn NWLn NWSg NWMg NWLg SH SWSn SWMn SWLn SWSg SWMg SWLg"
OUT=$1.NW
rm $OUT
touch $OUT

while read -r line; do
	B1=$(echo $line | sed 's/^\([0-9]*\).*/\1/')
	B2=$(echo $line | sed 's/^\([0-9]*\) \([0-9]*\).*/\2/')
	P1=$(echo $line | sed 's/^\([0-9]*\) \([0-9]*\) \([0-9]*\).*/\3/')
	P2=$(echo $line | sed 's/^\([0-9]*\) \([0-9]*\) \([0-9]*\) \([0-9]*\).*/\4/')
	PROB=$(echo $line | sed 's/^.* \([^ ]*\)$/\1/')


	FILE1="$P/sampling/sample${B1}/fragments.csv"
	FILE2="$P/sampling/sample${B2}/fragments.csv"

	FRAG1=$(cat $FILE1 | sed -n "/ ${P1}$/p" | sed "s/ .*//g")
	FRAG2=$(cat $FILE2 | sed -n "/ ${P2}$/p" | sed "s/ .*//g")

	S=$(/ebio/abt1/lweidmann/seq_space/final/homology/HHalign/align2seq $FRAG1 $FRAG2)
	#echo $B1 $B2 $P1 $P2 $FRAG1 $FRAG2 $PROB $S
	echo $B1 $B2 $P1 $P2 $FRAG1 $FRAG2 $PROB $S >> $OUT
done < $1
echo "NW of $1"
