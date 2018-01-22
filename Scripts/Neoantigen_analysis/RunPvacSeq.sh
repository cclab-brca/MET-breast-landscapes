#!/bin/bash

#SBATCH -J pvac			      # job name
#SBATCH -a 1-9%3			  # job array
#SBATCH -N 1                      # number of nodes
#SBATCH -n 1                     # number of min/max cores
#SBATCH --mem 10000                # memory pool for all cores
#SBATCH -o pvac.%a.out        # STDOUT
#SBATCH -e pvac.%a.err        # STDERR

# Runs pVacSeq pipeline

pvacseqDir=/scratchb/cclab/autopsies/pvacseq
polysolverDir=/scratchb/cclab/autopsies/pvacseq/hla

splitVepVCFDir=$pvacseqDir/vcf-vep-split
logDir=$pvacseqDir/logs


INDEX=$SLURM_ARRAY_TASK_ID
VCF_PATH=$(ls $splitVepVCFDir/*vcf | sed -n "$INDEX"p)
VCF_FILENAME=${VCF_PATH##*/}
sampleName=${VCF_FILENAME%.vcf*}
logFile=$logDir/$sampleName.log.txt

outDir=$pvacseqDir/output/$sampleName
mkdir -p $outDir

echo "##### RUNNING ON FILE: "$sampleName > $logFile; echo "" >> $logFile

tumourPolySolverPath="$polysolverDir/$sampleName".winners.hla.txt
patient=($(echo $sampleName | sed -e 's/_/\n/g'))
normalPolySolverPath=$(ls $polysolverDir/${patient[0]}*BC*txt)
    
POLYSOLVER_ALLELES=$(sed 's/\t/\n/g'  $normalPolySolverPath | grep hla_ | sort | uniq)

#extract and reformat polysolver results from sample name file
POLYSOLVER_STRING=""
    
for POLYSOLVER_ALLELE in $POLYSOLVER_ALLELES; do
 FORMAT_ALLELE1=$(echo $POLYSOLVER_ALLELE | cut -d'_' -f1-4 | tr [a-z] [A-Z] | sed 's/_/-/' | sed 's/_/*/' | sed 's/_/:/')
 v3=$(grep -Fc $FORMAT_ALLELE1 <<< "$POLYSOLVER_STRING")
 if  [ "$v3" -eq "0" ]; then
 POLYSOLVER_STRING=$POLYSOLVER_STRING","$FORMAT_ALLELE1
 fi
done

POLYSOLVER_PVACSEQ=$(echo $POLYSOLVER_STRING | sed 's/,//')

echo "##### POLYSOLVER ALLELES TO RUN ON: "$POLYSOLVER_PVACSEQ >> $logFile; echo "" >> $logFile
source activate py35

pvacseq run \
-e=8,9,10,11 \
--netmhc-stab \
--iedb-install-directory /tmp \
$VCF_PATH $sampleName $POLYSOLVER_PVACSEQ \
NetMHC SMM NetMHCpan \
$outDir >> $logFile

source deactivate py35

