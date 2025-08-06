#!/bin/tcsh

# LSF Batch Submission Script for Structural Variants Detection
# Usage: ./q_structural_variants.sh SAMPLE_LIST_FILE
# 
# This script submits AnchorWave structural variant detection jobs for multiple samples

set sampleList=$1

# Validate input file
if ( ! -f "$sampleList" ) then
    echo "ERROR: Sample list file not found: $sampleList"
    echo "Usage: $0 SAMPLE_LIST_FILE"
    exit 1
endif

echo "Submitting structural variant detection jobs for samples in: $sampleList"

# Loop over the sample names in the sample list file
foreach sample (`cat ${sampleList}`)
    echo "Submitting structural variant jobs for sample: $sample"
    
    # Submit AnchorWave part 1 job (alignment)
    echo "  Submitting AnchorWave part 1 for: $sample"
    set aw1_par="-q sara -n 8 -W 6:00 -o aw1_${sample}_%J.stdout -e aw1_${sample}_%J.stderr"
    set aw1_job=`bsub $aw1_par -R 'span[hosts=1]' -R 'rusage[mem=16GB]' ../scripts/genomics_pipelines/structural_variants/anchorwave_part1.sh ${sample} | grep -o '<[0-9]*>'`
    
    # Submit AnchorWave part 2 job (structural variant detection)
    echo "  Submitting AnchorWave part 2 for: $sample"
    set aw2_par="-q sara -n 4 -W 4:00 -o aw2_${sample}_%J.stdout -e aw2_${sample}_%J.stderr"
    set aw2_job=`bsub $aw2_par -R 'rusage[mem=8GB]' -w "done($aw1_job)" ../scripts/genomics_pipelines/structural_variants/anchorwave_part2.sh ${sample} | grep -o '<[0-9]*>'`
    
    # Submit breakpoint plotting job (depends on AnchorWave part 2)
    echo "  Submitting breakpoint plotting for: $sample"
    set plot_par="-q sara -n 1 -W 1:00 -o plot_${sample}_%J.stdout -e plot_${sample}_%J.stderr"
    bsub $plot_par -R 'rusage[mem=4GB]' -w "done($aw2_job)" ../scripts/genomics_pipelines/structural_variants/make_breakpoint_dotplots.sh ${sample}
    
    echo "  Structural variant pipeline jobs submitted for sample: $sample"
end

echo "All structural variant detection jobs submitted with proper dependencies"