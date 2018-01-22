#!/bin/bash

################################################################################
##### input files needed #####
## Output from 02_rerunASCAT.sh

##### programs and scripts needed #####
## prepareMEDICCinput.R
## MEDICC
## R package data.table

################################################################################
##### --- Adjust as needed --- #####
## Data directories
rawdatadir=~/Documents/Caldas/Data/dataRaw/
intdatadir=~/Documents/Caldas/Data/dataProcessed/
resultsdir=~/Documents/Caldas/Data/dataResults/medicc

## Program locations
mediccdir=~/Documents/Caldas/Programs/MEDICC/medicc
scriptdir=~/Documents/Caldas/Data/breastMetRepo/mediccAnalysis

################################################################################
## create directories
mkdir -p ${resultsdir}

prepareMedicc=true
runMedicc=true

## prepare medicc input
if [ ${prepareMedicc} = true ] ; then

    while read case; do

	Rscript ${scriptdir}/03_prepareMEDICCinput.R ${rawdatadir}/samples/medicc/samples_${case}.txt ${intdatadir}/ascat/ascatRefitted/${case} ${intdatadir}/medicc/mediccInput/${case}

    done <  ${rawdatadir}/samples/medicc/cases.txt
fi


## run medicc
source ${mediccdir}/../profile

if [ ${runMedicc} = true ] ; then

    while read case; do
	mkdir -p ${intdatadir}/medicc/mediccPhasing/${case}

	python ${mediccdir}/medicc_phase.py ${intdatadir}/medicc/mediccInput/${case}/desc.txt ${intdatadir}/medicc/mediccPhasing/${case} -v

	mkdir -p ${resultsdir}/${case}

	python ${mediccdir}/medicc.py  ${intdatadir}/medicc/mediccPhasing/${case}/phased_desc.txt ${resultsdir}/${case} -v -s

	rm /tmp/tmp*

    done < ${rawdatadir}/samples/medicc/cases.txt

fi

