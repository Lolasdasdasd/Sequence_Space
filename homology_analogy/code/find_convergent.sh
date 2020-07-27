#!/bin/bash

# input parameter
if [ $# -ne 3 ]; then
	echo "please provide: <batch1> <batch2> <BASE=ali/ali2>"
	exit 1
fi

BASE=$3

#./ecod_assign_non_homologs res/frag2ecod_results_all.txt.assign HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologs > HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP
#~/seq_space/final/homology/HHalign/disectCQP.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP
## ---> HHalign/ali_${1}/ali_${1}_${2}.align.non_homologsCQP.C without Class description C: "<Prob> C" 

## sequence bias
#echo "Sequence bias of homologs..."
#~/seq_space/final/homology/HHalign/HOMOLOGS2NW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.homologs
echo "Sequence bias of analogs..."
~/seq_space/final/homology/HHalign/HOMOLOGS2NW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.C
#echo "Sequence bias of unknowns..."
#~/seq_space/final/homology/HHalign/HOMOLOGS2NW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.QP


# D-model
echo "Composition bias of homologs..."
~/seq_space/final/homology/HHalign/DOMAIN_MODELNW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.homologs           > HHalign/ali_${1}/${BASE}_${1}_${2}.align.homologs.shuf.NW
echo "Composition bias of analogs..."
~/seq_space/final/homology/HHalign/DOMAIN_MODELNW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.C  > HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.C.shuf.NW
echo "Composition bias of unknown..."
~/seq_space/final/homology/HHalign/DOMAIN_MODELNW.sh HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.QP > HHalign/ali_${1}/${BASE}_${1}_${2}.align.non_homologsCQP.QP.shuf.NW
echo "...done"



