#!/bin/tcsh

# LSF Batch Submission Script for Crow 2020 RNA-seq Reanalysis
# Usage: bsub < q_crow_reanalysis.sh
# OR: ./q_crow_reanalysis.sh [SAMPLE_LIST_FILE] [OPTIONS]
# 
# This script submits the Crow 2020 reanalysis pipeline to LSF

#BSUB -J crow_reanalysis
#BSUB -q sara
#BSUB -n 8
#BSUB -W 12:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o crow_reanalysis_%J.stdout
#BSUB -e crow_reanalysis_%J.stderr

# Change to the crow reanalysis directory
cd ../scripts/crow_reanalysis

# Run the comprehensive RNA-seq processing pipeline
# Modify arguments as needed - currently processes all samples
bash process_crow2020_rnaseq.sh --threads 8 --bootstrap 100

echo "Crow 2020 reanalysis pipeline completed"