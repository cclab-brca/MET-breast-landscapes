
## 1. prepare treeomics input data
Rscript 01_prepData_WES.R

## 2. run treeomics to get posterior data table
python3.5 ${TREEOMICSDIR}/treeomics -u -r ${DATADIR}/input_WES_all/X${CASE}_altCount.txt -s ${DATADIR}/input_WES_all/X${CASE}_covCount.txt -o ${DATADIR}/output_WES_all/X${CASE}

## 3. run OncoNEM
Rscript 02_oncoNEMAnalysis.R



