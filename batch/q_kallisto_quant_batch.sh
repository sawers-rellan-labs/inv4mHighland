#!/bin/tcsh

# LSF Batch Submission Script for Kallisto RNA-seq Quantification
# Usage: ./q_kallisto_quant_batch.sh SAMPLE_LIST_FILE
# 
# This script submits Kallisto quantification jobs for multiple samples using LSF bsub
# Each sample gets submitted as a separate job with appropriate resource allocation

set sampleList=$1

# Validate input file
if ( ! -f "$sampleList" ) then
    echo "ERROR: Sample list file not found: $sampleList"
    echo "Usage: $0 SAMPLE_LIST_FILE"
    exit 1
endif

echo "Submitting Kallisto quantification jobs for samples in: $sampleList"

# Loop over the sample names in the sample list file
foreach sample (`cat ${sampleList}`)
    echo "Submitting job for sample: $sample"
    
    # LSF parameters for Kallisto quantification (moderate resource requirements)
    set par="-q sara -n 8 -W 2:00 -o %J.stdout -e %J.stderr"
    
    # Submit the bsub command to the job scheduler
    bsub $par -R 'span[hosts=1]' -R 'rusage[mem=8GB]' ../scripts/genomics_pipelines/rnaseq_processing/kallisto_quantify_sample.sh ${sample}
end

echo "All Kallisto quantification jobs submitted"