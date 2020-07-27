#!/bin/bash

if [ $# -ne 1 ]; then
  echo "please provide: <IN.align>\n     (PROB-THRESHOLD is 90%)"
  exit 1
fi



IN=$1

F="align_process.tmp"
sed -n '/Probab=/p' $IN > $F


T=$(cat $F | wc -l)
echo "-----------------------------------------"
echo "Number of alignments: $T in $1"

H=$(sed -n '/Probab=[9][0-9]/p' $F | wc -l)
H2=$(sed -n '/Probab=100/p' $F | wc -l)
echo "Number of 90-100% homologous pairs: $H + $H2 (# of Prob = 100)"
python -c "print(100*($H.0 + $H2.0)/$T.0)"
python -c "print(100*($H.0 + $H2.0)/$T.0)" > $1.fraction_homologs

rm tmp123 tmp234
cat $IN | sed -n '/Command/,+7p' | sed -n '/^[CP]/p' > tmp123 
sed 's#^.*HHblits/a3m_\([0-9]*\)/frag_\([0-9]*\).*a3m_\([0-9]*\)/frag_\([0-9]*\).*#\1 \3 \2 \4#g' tmp123 > tmp234

sed 's#Probab=\([^ ]*\).*$#\1#g' tmp234 | sed -z 's#\([0-9]*\) \([0-9]*\)\n#\1 \2 #g' | sed -n '/ 9[0-9]\./p' > $1.homologs
sed 's#Probab=\([^ ]*\).*$#\1#g' tmp234 | sed -z 's#\([0-9]*\) \([0-9]*\)\n#\1 \2 #g' | sed -n '/ 100\./p' >> $1.homologs
#sed -n '/ [0-9]\{1,2\}\./p'

#                                                                                     ---v   removes 100.0% Prob     ---v removes 90-99.9%
sed 's#Probab=\([^ ]*\).*$#\1#g' tmp234 | sed -z 's#\([0-9]*\) \([0-9]*\)\n#\1 \2 #g' | sed -n '/ [0-9]\{1,2\}\./p' | sed -n '/9[0-9]\./!p' > $1.non_homologs
rm tmp123 tmp234

rm $F

