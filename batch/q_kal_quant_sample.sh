#!/bin/tcsh

# LSF Batch Submission Script for Kallisto RNA-seq Quantification (Legacy)
# Usage: ./q_kal_quant_sample.sh SAMPLE_LIST_FILE
# 
# Note: This script calls ./quant_sample.sh which should be in the same directory
# For new genomics pipeline, use q_kallisto_quant_batch.sh instead

set sampleList=$1
set index=Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx
set transcripts=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa

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
    
    # LSF parameters for Kallisto quantification - updated for consistency
    set par="-q sara -n 8 -W 2:00 -o %J.stdout -e %J.stderr"
    
    # Submit the bsub command to the job scheduler
    bsub $par -R 'rusage[mem=4GB]' ./quant_sample.sh ${sample}
end

echo "All Kallisto quantification jobs submitted"
