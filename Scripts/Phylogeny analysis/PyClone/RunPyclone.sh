#!/bin/bash

# Title:	Run Pyclone 0.13

#SBATCH -J Py			      # job name
#SBATCH -a 1-10%10				  # job array
#SBATCH -N 1                      # number of nodes
#SBATCH -n 1                      # number of min/max cores
#SBATCH --mem 10000                # memory pool for all cores
#SBATCH -o pyclone.%a.out        # STDOUT
#SBATCH -e pyclone.%a.err        # STDERR

set -e


pycloneDir="/projects/pyclone"

cellularityFile=$pycloneDir/data/ascat-final-cellularity-mod.txt

patientID=$(ls *tsv | cut -d'-' -f1| uniq| sed -n "$SLURM_ARRAY_TASK_ID"p )
burnin=20000
cd $pycloneDir
patientID=$(ls *tsv | cut -d'_' -f1| uniq| sed -n "$SLURM_ARRAY_TASK_ID"p )
tsvFiles=$(ls $patientID*)
mkdir -p $pycloneDir/workspace/$patientID

#tumPurity=$(grep $patientID $cellularityFile | cut -f2 | xargs)
tumPurity=$(grep "$patientID" $cellularityFile | cut -f2 | xargs)

PyClone run_analysis_pipeline --in_files $tsvFiles --working_dir $pycloneDir/workspace/$patientID \
--tumour_contents $tumPurity \
--density pyclone_beta_binomial \
--num_iters 40000 \
--prior parental_copy_number \
--burnin $burnin \
--min_cluster_size 3
