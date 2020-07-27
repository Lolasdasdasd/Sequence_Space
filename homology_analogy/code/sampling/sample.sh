#!/bin/bash

# ./frag_sampling
# please provide <path> <lf> <step> <offset>

STEP=1000000
POS_SAMP=100


for i in `seq 0 9`; do
	OFFSET=$(( i*$STEP/$POS_SAMP ))
	mkdir sample$i 
	echo "./frag_sampling ~/seq_space/final/hcM 100 $STEP $OFFSET > sample$i/fragments.csv 2> /dev/null"
	./frag_sampling ~/seq_space/final/hcM 100 $STEP $OFFSET > sample$i/fragments.csv 2> /dev/null
done 


