#!/bin/bash
if [ $# -ne 0 ]; then
    echo "please provide: "
    exit 1
fi

rm homologs.txt
rm homologs.txt.tmp


#cat ali_*/*.align > all.align


for i in 0 1 2 3 4 5 6 7 8 ; do
	rm tmp${i}.align homologs.txt$i.tmp

	cp ali_$i/ali_${i}_$((i+1)).align tmp${i}.align
	bash progress_align.sh tmp${i}.align homologs.txt$i.tmp

	bash homologs2NW.sh homologs.txt$i.tmp > homologs$i.txt
done

i=9
rm tmp${i}.align homologs.txt$i.tmp
cp ali_$i/ali_${i}_0.align tmp${i}.align
bash progress_align.sh tmp${i}.align homologs.txt$i.tmp
bash homologs2NW.sh homologs.txt$i.tmp > homologs$i.txt

