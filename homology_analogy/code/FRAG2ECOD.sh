#!/bin/bash


# input parameter
if [ $# -ne 2 ]; then
    echo "please provide: <Prob> <overlap>"
    exit 1
fi



for n in 0 1 2 3 4 5 6 7 8 9; do
	for f in ~/seq_space/final/homology/HHblits/a3m_$n/frag_*ecodname; do 
		echo ./frag2ecod ${f%.ecodname} $1 $2; 
	done > do_frag2ecod_${n}.sh
	RES="frag2ecod_results_${n}_${1}_${2}.txt"

	echo "bash do_frag2ecod_${n}.sh > $RES"
	bash do_frag2ecod_${n}.sh > $RES
	
	echo double domains?
	cat $RES | pcregrep -M '>>.*\n>>'
	
	echo "tested frags:"
	cat $RES | sed -n '/lweidmann/p' | wc -l
	
	
	cp $RES ${RES}.assign
	sed -n '/^[A-Z]/!p' ${RES}.assign | sed 's#^>>.* // ##' > ${RES}.assign2
	vi ${RES}.assign2  -c ':%s#\n\([1-9].*\)$#-\1#g' -c 'wq'
	echo "with domain:"
	sed -n '/-/p' ${RES}.assign2 | wc -l
	echo "with two domains"
	sed -n '/-[0-9\.]*-/p' ${RES}.assign2 | wc -l
	echo "with three domains ?!"
	sed -n '/-[0-9\.]*-[0-9\.]*-/p' ${RES}.assign2
	sed 's#\(^.*frag2ecod\)-\([^-]*\)-\([^-]*\)#\1-\2\n\1-\3#' ${RES}.assign2 > ${RES}.assign3
	sed -i 's#^.*a3m_\([0-9]*\)/frag_\([0-9]*\)[^-]*#\1 \2#' ${RES}.assign3
	sed 's#-# #' ${RES}.assign3 > ${RES}.assign
	
	rm ${RES}.assign2 ${RES}.assign3

done

