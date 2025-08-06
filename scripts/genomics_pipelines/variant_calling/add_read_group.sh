#!/usr/bin/tcsh

# Get a command line arguments
set field2=$1 
set field3=$2

# Load required modules


set PICARD="/usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar"

set inputDir="aln_out"
set outputDir="with_rg"

set bam_suffix="Aligned.sortedByCoord.out.bam"
set rg_suffix="_rg_added_sorted.bam"

set bamin=${inputDir}/${field2}${bam_suffix}
set bamout=${outputDir}/${field2}${rg_suffix}


# Create the output directory
mkdir -p ${outputDir}

# Run  picard add readgroup

java -jar ${PICARD} AddOrReplaceReadGroups \
   -I  ${bamin} \
   -O  ${bamout}\
   -SO coordinate \
   -RGID 1 \
   -RGLB ${field2} \
   -RGPL illumina \
   -RGPU L002 \
   -RGSM ${field3}
