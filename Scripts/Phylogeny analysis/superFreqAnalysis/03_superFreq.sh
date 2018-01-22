#!/bin/bash

########################################################
## input files needed ##
# metaDataFile
# preliminary vcf files
# captureRegsionFile
# reference
# 

## programs needed
# R
# R package superFreq and dependencies
# VEP
VEPCALL=~/Programs/ensembl-vep-release-87/vep.pl

## directories
RAWDATADIR=../dataRaw
ANNOTDIR=../dataAnnotation
INTDATADIR=../dataProcessed/superFreq
RESULTSDIR=../dataResults/superFreq

mkdir -p ${INTDATADIR}/referenceNormals/bam

while read NORMAL; do

    ln -s ${RAWDATADIR}/BAMs_WES/${NORMAL}.merged.realn.RG.bam ${INTDATADIR}/referenceNormals/bam/.
    ln -s ${RAWDATADIR}/BAMs_WES/${NORMAL}.merged.realn.RG.bam.bai ${INTDATADIR}/referenceNormals/bam/.

done < ${INTDATADIR}/allNormals.txt

while read CASE; do

    mkdir -p ${RESULTSDIR}/${CASE}/R
    mkdir -p ${RESULTSDIR}/${CASE}/plots

    ~/Programs/R-3.3.2/bin/Rscript superFreq.R ${RAWDATADIR}/samples/superfreq/metadata_${CASE}.txt ${INTDATADIR}/referenceNormals/ ${RESULTSDIR}/${CASE}/R ${RESULTSDIR}/${CASE}/plots ${ANNOTDIR}/captureRegions/nexterarapidcapture_exome_targetedregions.bed ${ANNOTDIR}/Reference/human_g1k_v37_decoy.fasta ${VEPCALL}

done < ${RAWDATADIR}/samples/superFreq/allCases.txt

