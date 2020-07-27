#!/bin/bash

FILE="Bacteria/species2id.table.unfiltered"




# ALL proteomes

NEWFILE=${FILE}.filter

cp $FILE $NEWFILE


# delete lines starting with:
sed -i '/^[a-z]'/d $NEWFILE
sed -i '/^Alpha'/d $NEWFILE
sed -i '/^Beta'/d $NEWFILE
sed -i '/^Gamma'/d $NEWFILE
sed -i '/^Gamma'/d $NEWFILE
sed -i '/^Candidatus'/d $NEWFILE
sed -i '/^Uncultured'/d $NEWFILE
sed -i '/[A-Z][a-z]* bacterium'/d $NEWFILE
sed -i '/-like'/d $NEWFILE

# delete single bacterium
#sed -i '/Crocosphaera watsonii'/d $NEWFILE
sed -i '/OM182'/d $NEWFILE
sed -i '/PVC'/d $NEWFILE
sed -i '/SAR'/d $NEWFILE

# search and replace
sed -i 's_\[__g' $NEWFILE
sed -i 's_\]__g' $NEWFILE
sed -i 's_"__g' $NEWFILE
sed -i "s_'__g" $NEWFILE
sed -i "s_(.*)__g" $NEWFILE
sed -i "s_\([A-Z][a-z]*\) [^\t]*_\1_g" $NEWFILE

sort $NEWFILE -o $NEWFILE


echo "number of fastas: " $(wc -l $FILE) " after filter: " $(wc -l $NEWFILE)
