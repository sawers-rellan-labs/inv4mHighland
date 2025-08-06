#!/bin/tcsh

# LSF Batch Submission Script for GO Annotation Building
# Usage: bsub < queue_build_go_annotation.sh
# 
# This script can be submitted directly to LSF for GO annotation processing

#BSUB -J build_go_annotation
#BSUB -q sara
#BSUB -n 1
#BSUB -W 2:00
#BSUB -R "rusage[mem=4GB]"
#BSUB -o %J.stdout
#BSUB -e %J.stderr

# Activate R environment and run script
conda activate /usr/local/usrapps/maize/sorghum/conda/envs/r_env
Rscript get_full_GOannotation.R