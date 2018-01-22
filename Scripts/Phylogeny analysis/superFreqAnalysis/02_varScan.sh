#!/bin/bash

########################################################
## input files needed ##
# aligned BAM files

## programs needed
# VarScan
VARSCANDIR=~/Programs/VarScan

## directories
RAWDATADIR=../dataRaw/
INTDATADIR=../dataProcessed/superFreq

mkdir -p ${INTDATADIR}/SNPs

while read SAMPLE; do

    java -jar ${VARSCANDIR}/VarScan.v2.4.3.jar mpileup2cns ${INTDATADIR}/pileups/${SAMPLE}.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.01 > ${INTDATADIR}/SNPs/${SAMPLE}.vcf

done < ${RAWDATADIR}/samples/superFreq/allSamples.txt

