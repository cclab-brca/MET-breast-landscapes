#!/bin/bash
#
#
#SBATCH --nodes 1 # one node
#SBATCH --ntasks 1 # one core
#SBATCH --mem 320000 # Memory request

source ~/.bash_profile

## targeted All 
#/Users/ross01/.pyenv/versions/3.5.2/bin/python3.5 ${TREEOMICSDIR}/treeomics -u -r ${DATADIR}/input_targeted_all/X${CASE}_altCount_${TYPE}.txt -s ${DATADIR}/input_targeted_all/X${CASE}_covCount_${TYPE}.txt -o output_targeted_all/X${CASE}_${TYPE} -n ${NORMAL}

## targeted selected
/Users/ross01/.pyenv/versions/3.5.2/bin/python3.5 ${TREEOMICSDIR}/treeomics -u -r ${DATADIR}/input_targeted_selected/X${CASE}_${TYPE}_altCount.txt -s ${DATADIR}/input_targeted_selected/X${CASE}_${TYPE}_covCount.txt -o output_targeted_selected/X${CASE}_${TYPE} -n ${NORMAL}




