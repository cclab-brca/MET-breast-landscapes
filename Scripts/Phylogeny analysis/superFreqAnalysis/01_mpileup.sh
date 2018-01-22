#!/bin/bash

########################################################
## input files needed ##
# aligned BAM files

## programs needed
# samtools
SAMTOOLSDIR=~/Programs/samtools-1.3.1/bin

## directories
RAWDATADIR=../dataRaw
ANNOTDIR=../dataAnnotation
INTDATADIR=../dataProcessed/superFreq

mkdir -p ${INTDATADIR}/pileups

while read SAMPLE; do

    BAM=$( ls ${RAWDATADIR}/BAMs_WES/${SAMPLE}*.bam)
    
    ${SAMTOOLSDIR}/samtools mpileup -d 10000 -q 1 -Q 15 -A -f ${ANNOTDIR}/Reference/human_g1k_v37_decoy.fasta ${BAM} > ${INTDATADIR}/pileups/${SAMPLE}.mpileup

done < ${RAWDATADIR}/samples/superFreq/allSamples.txt

#
