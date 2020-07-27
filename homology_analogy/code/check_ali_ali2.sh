#!/bin/bash


for i in `seq 0 9`; do
	echo "*** Batch $i ***"
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*homologs.NW -tr | wc -l
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*C.NW -tr | wc -l
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*QP.NW -tr | wc -l
	echo "---"
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*homologs.shuf.NW -tr | wc -l
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*C.shuf.NW -tr | wc -l
	ls -l ~/seq_space/final/homology/HHalign/ali_$i/*QP.shuf.NW -tr | wc -l
done 


