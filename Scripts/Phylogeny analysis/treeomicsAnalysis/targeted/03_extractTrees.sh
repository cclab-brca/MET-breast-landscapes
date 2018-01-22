#!/bin/bash

DATADIR=~/Documents/Caldas/Data/dataProcessed/treeomics/output_targeted_selected/
RESDIR=~/Documents/Caldas/Data/dataResults/treeomics/targeted/

mkdir -p ${DATADIR}/extractedTrees/
mkdir -p ${RESDIR}

for CASE in 288 290 290exOvary 291 298 298main 308 315 323 328 330 DET52
do

    for TYPE in All Biopsies Fluids WESsamples
    do
	if [ -e ${DATADIR}/X${CASE}_${TYPE}/*.tex ]
	then

	    grep -E '\[|\]|% Acq'  ${DATADIR}/X${CASE}_${TYPE}/*.tex  | sed 's/\\node\[black,draw,text width=1\.09cm,inner sep=2pt,align=center\]//' | sed 's/\\edge node\[above, NavyBlue\]{\\footnotesize//' | sed 's/} node\[below, black!60\]{\\footnotesize//' | sed 's/\.\\small//' | sed 's/{[0-9]*\\%}//' | sed 's/ (.*)//' | sed 's/{//' | sed 's/};//' | sed 's/}//' | sed 's/\\Tree \[\.//' | sed 's/\\.*//' | sed 's/[ ]//g' | sed '/^[[:blank:]]*$/d'  > ${DATADIR}/extractedTrees/${CASE}_${TYPE}.txt
	    
	    (cd ${DATADIR}/X${CASE}_${TYPE}

	     ## run pdflatex
	     TEXFILE=$(ls *.tex)
	     pdflatex ${TEXFILE}

	     ## move treomics plot
	     FILE=$(ls X${CASE}*mlhtree_full.pdf)
	     cp ${FILE} ${RESDIR}/${CASE}_${TYPE}_treeomics_tree.pdf

	     ## move bayesian data table
	     FILE=$(ls mlh_mutnode_labels_*.txt)
	     cp ${FILE} ${RESDIR}/${CASE}_${TYPE}_treeomics_mlh_mutnode_labels.txt
	     
	     ## move bayesian data table plot
	     FILE=$(ls *_bayesian_data_table.png)
	     cp ${FILE} ${RESDIR}/${CASE}_${TYPE}_treeomics_bayesian_data_table.png	     
	     )
	fi
	
    done
    
done
