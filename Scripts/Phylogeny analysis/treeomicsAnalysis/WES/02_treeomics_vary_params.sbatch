#!/bin/bash
#
#
#SBATCH --nodes 1 # one node
#SBATCH --ntasks 1 # one core
#SBATCH --mem 320000 # Memory request
#

source ~/.bash_profile

## WES All 
#/Users/ross01/.pyenv/versions/3.5.2/bin/python3.5 ${TREEOMICSDIR}/treeomics -u -r ${DATADIR}/input_WES_all/X${CASE}_altCount.txt -s ${DATADIR}/input_WES_all/X${CASE}_covCount.txt -o output_WES_all/X${CASE}

## WES selected
/Users/ross01/.pyenv/versions/3.5.2/bin/python3.5 ${TREEOMICSDIR}/treeomics -u -r ${DATADIR}/input_WES_selected/X${CASE}_altCount.txt -s ${DATADIR}/input_WES_selected/X${CASE}_covCount.txt -o output_WES_selected/X${CASE} -n ${NORMAL}


