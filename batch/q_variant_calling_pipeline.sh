#!/bin/tcsh

# LSF Batch Submission Script for GATK Variant Calling Pipeline
# Usage: ./q_variant_calling_pipeline.sh SAMPLE_LIST_FILE
# 
# This script submits the complete variant calling pipeline jobs for multiple samples
# Pipeline includes: alignment -> add read groups -> mark duplicates -> BQSR -> call variants

set sampleList=$1

# Validate input file
if ( ! -f "$sampleList" ) then
    echo "ERROR: Sample list file not found: $sampleList"
    echo "Usage: $0 SAMPLE_LIST_FILE"
    exit 1
endif

echo "Submitting variant calling pipeline jobs for samples in: $sampleList"

# Loop over the sample names in the sample list file
foreach sample (`cat ${sampleList}`)
    echo "Processing sample: $sample"
    
    # 1. Submit alignment job
    echo "  Submitting alignment job for: $sample"
    set align_par="-q sara -n 16 -W 6:00 -o align_${sample}_%J.stdout -e align_${sample}_%J.stderr"
    set align_job=`bsub $align_par -R 'span[hosts=1]' -R 'rusage[mem=32GB]' ../scripts/genomics_pipelines/variant_calling/align_sample.sh ${sample} | grep -o '<[0-9]*>'`
    
    # 2. Submit add read group job (depends on alignment)
    echo "  Submitting add read group job for: $sample"
    set rg_par="-q sara -n 1 -W 1:00 -o rg_${sample}_%J.stdout -e rg_${sample}_%J.stderr"
    set rg_job=`bsub $rg_par -R 'rusage[mem=4GB]' -w "done($align_job)" ../scripts/genomics_pipelines/variant_calling/add_read_group.sh ${sample} | grep -o '<[0-9]*>'`
    
    # 3. Submit mark duplicates job (depends on read group)
    echo "  Submitting mark duplicates job for: $sample"
    set md_par="-q sara -n 1 -W 2:00 -o md_${sample}_%J.stdout -e md_${sample}_%J.stderr"
    set md_job=`bsub $md_par -R 'rusage[mem=8GB]' -w "done($rg_job)" ../scripts/genomics_pipelines/variant_calling/mark_duplicates.sh ${sample} | grep -o '<[0-9]*>'`
    
    # 4. Submit variant calling job (depends on mark duplicates)
    echo "  Submitting variant calling job for: $sample"
    set vc_par="-q sara -n 1 -W 4:00 -o vc_${sample}_%J.stdout -e vc_${sample}_%J.stderr"
    bsub $vc_par -R 'rusage[mem=6GB]' -w "done($md_job)" ../scripts/genomics_pipelines/variant_calling/call_variants.sh ${sample}
    
    echo "  Pipeline jobs submitted for sample: $sample"
end

echo "All variant calling pipeline jobs submitted with proper dependencies"