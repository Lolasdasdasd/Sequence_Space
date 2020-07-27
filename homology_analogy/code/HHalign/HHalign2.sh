#!/bin/bash

# input parameter
if [ $# -ne 2 ]; then
	echo "please provide: <int [0-9]> <int [0-9]>"
	exit 1
fi

TMP1=$1
TMP2=$2

## fragment sampling for HHblits
## lf=100 / N samples =10 /  sample id = 0
#fragment position
#MYFVDNGNNYDASLNIALETYLVENRLVDEPILLFYINDPSIIVGRNQNTIEEVNQAYVEEKGIRVVRRMSGGGAVYHDRGNFSFCFIKDDDGSFRDFAS 0
#IVGYVGKKNAQDVILQGLEKLEYRGYDSAGIAIVEDGVIKSEKFKGRLAVLSDFLEANPIKGSLGIGHTRWATHGAPSDENSHPHLNADNTIAVVHNGII 4738320

#HHblits="/ebio/abt1_share/toolkit_dev/bioprogs/hhsuite/bin/hhalign"
DB="/ebio/abt1_share/toolkit_sync/databases/hhblits/uniclust30"
export HHLIB="/ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/"
MAX=100000
echo $HHLIB
HEADER=3
C=0
N=0
MAX=$3
P="/ebio/abt1/lweidmann/seq_space/final/homology/"


mkdir ali_${TMP1}
OUT="$P/HHalign/ali_${TMP1}/ali2_${TMP1}_${TMP2}.align"
echo "" > $OUT

while read -r line1; do
  let N++
  frag1=$(echo $line1 | sed "s/ .*$//" )
  pos1=$(echo $line1 | sed "s/^.* //" )
  if (( $N <= $HEADER)); then
    continue
  fi
  if (( $N > $MAX)); then
    break
  fi
  F1="$P/HHblits/a3m_${TMP1}/frag_${pos1}.a3m"
  if [ ! -f $F1 ]; then continue; fi
  M=0
  while read -r line2; do
	let M++
  	if (( $M <= $HEADER)); then
    	continue
  	fi
    if (( $M >= $N)); then
      break
    fi
    if (( $M > $MAX)); then
      break
    fi
	echo LINE1= $line1
	echo LINE2= $line2
	echo $N $M
    frag2=$(echo $line2 | sed "s/ .*$//" )
    pos2=$(echo $line2 | sed "s/^.* //" )
	F2="$P/HHblits/a3m_${TMP2}/frag_${pos2}.a3m"
  	if [ ! -f $F2 ]; then break; fi

	echo $frag1 $pos1 $frag2 $pos2
	let C++
    /ebio/abt1_share/small_projects/lweidmann/hh-suite/bin/hhalign -i $F1 -t $F2 -norealign >> $OUT
  done < $P/sampling/sample${TMP2}/fragments.csv
done < $P/sampling/sample${TMP1}/fragments.csv

echo "Total comparisons = $C"



