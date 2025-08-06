#!/bin/tcsh

set $1

# Specify your input file (replace 'your_file.txt' with your actual file name)
set picard="/usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar"

set inputDir="split"
set outputDir="MQ"

set input_suffix="_split.bam"
set output_suffix="_qual.pdf"
set txt_suffix="_qual.txt"

set bamin=${inputDir}/${field2}${input_suffix}
set pdfout=${outputDir}/${field2}_${inputDir}${output_suffix}
set txtout=${outputDir}/${field2}_${inputDir}${txt_suffix}


# Create the output directory
mkdir -p ${outputDir}

# Run  picard add readgroup

java -Xms4g -jar $picard QualityScoreDistribution \
      -I $bamin\
      -O $txtout \
      -CHART $pdfout


java -Xms4g -jar $picard QualityScoreDistribution \
      -I recal/4131_recal.bam\
      -O MQ/4131_recal.txt \
      -CHART MQ/4131_recal.pdf
