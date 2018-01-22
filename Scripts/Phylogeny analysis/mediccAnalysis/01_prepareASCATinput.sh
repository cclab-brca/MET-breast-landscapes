#!/bin/bash

################################################################################
##### input files needed #####
## METADATA
#OK ${rawdatadir}/samples/medicc/baseCases.txt ## file listing all cases (patient ids)
#OK ${rawdatadir}/samples/medicc/cases.txt ## file listing all cases including duplicate cases with subset of samples
#OK ${rawdatadir}/samples/medicc/samples_germline.txt ## file listing the germline sample for each case
#OK ${rawdatadir}/samples/ascat/samples_${case}.txt ## file listing all samples for each case including germline
## RAWDATA
#OK ${rawdatadir}/SNPs_normal/${n}.merged.realn.RG.VCF ## (germline vcf files)
#OK ${rawdatadir}/BAMs_WES/${sample}*.bam ## WES for all samples
#OK ${rawdatadir}/QDNAseq/ ## Directory containing QDNAseq files for each sample
## REFERENCE
#OK ${referencedir}/human_g1k_v37_decoy.fasta ## reference file

##### programs and scripts needed #####
# alleleCounter (part of ascatNGS)
# ./01_prepareASCATinput.R
# ./01_runASCAT_asmultipcf.R
# latest ASCAT R package version needs to be installed

################################################################################

makeLociFiles=true
runAlleleCounter=true
prepareSegData=true
runSeg=true


##### --- Adjust as needed --- #####
## Data directories
rawdatadir=~/Caldas/Data/dataRaw/
intdatadir=~/Caldas/Data/dataProcessed/ascat/
referencedir=~/Caldas/Data/dataAnnotation/Reference/
#outdatadir=~/Caldas/Data/finalAnalysis/results/
## Program locations
alleleCounter=~/Caldas/Programs/ascatNGS/opt/bin/alleleCounter
##### ------------------------ #####

## create directories
mkdir -p ${intdatadir}/loci
mkdir -p ${intdatadir}/alleleCount
mkdir -p ${intdatadir}/ascatInput
mkdir -p ${intdatadir}/ascatOutput

while read baseCase; do

    ## define normal sample
    n=$(grep ${baseCase} ${rawdatadir}/samples/medicc/samples_germline.txt)
    echo ${baseCase}
    echo ${n}

    if [ ${makeLociFiles} = true ] ; then
	SNPFILE=$(ls ${rawdatadir}/SNPs_normal/ | grep ${n})
	## create loci files for allele counter and ASCAT
	cut -f 1,2 ${rawdatadir}/SNPs_normal/${SNPFILE} | grep -vP '#' | grep -E '^([0-9]|X|Y)'  > ${intdatadir}/loci/${baseCase}_loci_file.txt
	awk -v OFS='\t' '{print "SNP_" NR, $0}' ${intdatadir}/loci/${baseCase}_loci_file.txt > ${intdatadir}/loci/${baseCase}_SnpPositions.tsv
	echo -e 'Probe\tChr\tPosition' | cat - ${intdatadir}/loci/${baseCase}_SnpPositions.tsv > temp && mv temp ${intdatadir}/loci/${baseCase}_SnpPositions.tsv
    fi
	
    ################################################################################
    ## run ascatNGS allele count step
    
    if [ ${runAlleleCounter} = true ] ; then
	
	## run allele counter
	while read sample; do
	    
       	    echo "Running alleleCounter for sample ${sample}"
	    
	    bam=$(ls ${rawdatadir}/BAMs_WES/${sample}*.bam)
	    
	    if [ ! -e  ${intdatadir}/alleleCount/${sample}.count ]; then
		${alleleCounter} -b ${bam} -o ${intdatadir}/alleleCount/${sample}.count -l ${intdatadir}/loci/${baseCase}_loci_file.txt -r ${referencedir}/human_g1k_v37_decoy.fasta -m 20
	    fi
	    
	done < ${rawdatadir}/samples/medicc/samples_${baseCase}.txt

    fi

    ################################################################################
    ## prepare data for segmentation

    cases=$(grep ${baseCase} ${rawdatadir}/samples/medicc/cases.txt)

    for case in ${cases}; do
    
    if [ ${prepareSegData} = true ] ; then
	Rscript ./01_prepareASCATinput.R ${baseCase} ${case} ${rawdatadir}/samples/ascat/samples_${case}.txt ${n} ${intdatadir} ${rawdatadir}/QDNAseq ${intdatadir}/loci/${baseCase}_SnpPositions.tsv "c('${intdatadir}/alleleCount/',sampleID,'.count')" ${intdatadir}/alleleCount/${n}.count 24
    fi

    ################################################################################
    ## run multi sample segmentation and ascat
    
    if [ ${runSeg} = true ] ;  then	
	Rscript ./01_runASCAT_asmultipcf.R ${intdatadir}/ascatInput/${case} ${intdatadir}/ascatOutput/${case} ${case} ${rawdatadir}/samples/ascat/samples_${case}.txt ${n} "NULL" XX 24
    fi

    done

done < ${rawdatadir}/samples/medicc/baseCases.txt

################################################################################
## manually check ascat output and rerun with manually selected ploidy and purity estimates if necessary


