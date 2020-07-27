#!/bin/bash

C=0
for I in `seq 0 9`; do
	N=$(ls -l a3m_$I/*.a3m | wc -l)
	echo "$I) $N"
	C=$(( $C + $N ))
done
echo "-> total $C"
