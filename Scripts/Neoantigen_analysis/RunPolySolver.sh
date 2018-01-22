#!/bin/bash

#SBATCH -J polysolver			      # job name
#SBATCH -a 1-10%10				  # job array
#SBATCH -N 1                      # number of nodes
#SBATCH -n 8                     # number of min/max cores
#SBATCH --mem 40000                # memory pool for all cores
#SBATCH -o polysolver.%j.%a.out        # STDOUT
#SBATCH -e polysolver.%j.%a.err        # STDERR
 

# Runs PolySolver for HLA detection

INDEX=$SLURM_ARRAY_TASK_ID

# load preconfigured variables
source /scratcha/cclab/sammut01/software/polysolver/scripts/config.bash

wd=/scratchb/cclab/autopsies
bamDir=$wd/bam
outputDir=$wd/polysolver/results
sortTempDir=$wd/polysolver/sortTemp

mkdir -p $outputDir

bamPath=$(ls $bamDir/*m | sed -n "$INDEX"p)
fileName=${bamPath##*/}
sampleName=${fileName%.bam*}

tempDir=$wd/polysolver/$sampleName
TMP_DIR=$sortTempDir/$sampleName

/software/polysolver/scripts/shell_call_hla_type $bamDir/$fileName Unknown 1 hg19 STDFQ 0 $tempDir

mv $tempDir/winners.hla.txt $outputDir/$sampleName.winners.hla.txt
