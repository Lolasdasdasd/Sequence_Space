#!/bin/bash

# input parameter
if [ $# -ne 1 ]; then
	echo "cats thing of ali_X into ali_all <anything>"
	exit 1
fi

cat ali_[0-9]/*align.homologs.NW  > ali_all/all_homologs.NWf
cat ali_[0-9]/*homologs.shuf.NW  > ali_all/all_homologs.shuf.NWf

cat ali_[0-9]/*non_homologsCQP.C.NW  > ali_all/all_non_homologsCQP.C.NWf
cat ali_[0-9]/*non_homologsCQP.C.shuf.NW  > ali_all/all_non_homologsCQP.C.shuf.NWf

cat ali_[0-9]/*non_homologsCQP.QP.NW  > ali_all/all_non_homologsCQP.QP.NWf
cat ali_[0-9]/*non_homologsCQP.QP.shuf.NW  > ali_all/all_non_homologsCQP.QP.shuf.NWf
