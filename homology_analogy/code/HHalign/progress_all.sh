#!/bin/bash

# input parameter

for V in 8; do
#for V in `seq 0 9`; do
#	bash progress_align.sh ali_${V}/ali_${V}_$(((V+1)%10)).align
	bash progress_align.sh ali_${V}/ali_${V}_$(((V+1)%10)).align
done

