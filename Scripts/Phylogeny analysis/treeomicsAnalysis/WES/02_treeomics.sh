#!/bin/bash

########################################################
## input files needed ##
# aligned BAM files

## programs needed
# treeomics + dependencies
TREEOMICSDIR=/Users/ross01/Programs/treeomics/src

## directories
DATADIR=/mnt/scratchb/fmlab/ross01/Caldas/treeomics/Data
RESULTSDIR=/mnt/scratchb/fmlab/ross01/Caldas/treeomics/

## export environment variables
export TREEOMICSDIR DATADIR

mkdir -p logfiles/out

while read CASE; do
    export CASE

	if [ "${CASE}" = "288" ]; then
		NORMAL="X024"
	elif [ "${CASE}" = "290" ]; then
		NORMAL="X004"
	elif [ "${CASE}" = "290exOvary" ]; then
		NORMAL="X004"
	elif [ "${CASE}" = "291" ]; then
		NORMAL="X005"
	elif [ "${CASE}" = "298" ]; then
		NORMAL="X023"
	elif [ "${CASE}" = "298main" ]; then
		NORMAL="X023"
	elif [ "${CASE}" = "298Brain" ]; then
		NORMAL="X023"
	elif [ "${CASE}" = "308" ]; then
		NORMAL="X024"
	elif [ "${CASE}" = "315" ]; then
		NORMAL="X022"
	elif [ "${CASE}" = "323" ]; then
		NORMAL="X007"
	elif [ "${CASE}" = "328" ]; then
		NORMAL="X010"
	elif [ "${CASE}" = "330" ]; then
		NORMAL="X007"
	elif [ "${CASE}" = "DET52" ]; then
		NORMAL="Xnormal"
	elif [ "${CASE}" = "DET52main" ]; then
		NORMAL="Xnormal"
	fi 

    export NORMAL

    sbatch -o logfiles/out/${CASE}.out \
	   -e logfiles/out/${CASE}.err \
	   --job-name=${CASE} \
	   01_treeomics_vary_params.sbatch
    sleep 1 # pause to be kind to the scheduler

done < allCases.txt
#
