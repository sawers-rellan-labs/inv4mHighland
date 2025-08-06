#!/usr/bin/tcsh

# Get a command line arguments
set field2=$1 
set field3=$2

# Load required modules


set PICARD="/usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar"

set inputDir="with_rg"
set outputDir="dedup"

set rg_suffix="_rg_added_sorted.bam"
set dedup_suffix="_rg_added_dedup.bam"
set metrics_suffix=".metrics.txt"

set bamin=${inputDir}/${field2}${rg_suffix}
set bamout=${outputDir}/${field2}${dedup_suffix}
set metricsout=${outputDir}/${field2}${metrics_suffix}


# Create the output directory
mkdir -p ${outputDir}

# Run  picard add readgroup

java -jar ${PICARD} MarkDuplicates \
   -I ${bamin} \
   -O ${bamout} \
   -M ${metricsout} \
   -CREATE_INDEX true \
   -VALIDATION_STRINGENCY SILENT
