#!/usr/bin/tcsh

# Get a command line arguments
set field2=$1
set field3=$2

# Load required modules
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"


# Set the path to the genome reference
set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

set inputDir="dedup"
set outputDir="split"

set split_suffix="_split.bam"
set dedup_suffix="_rg_added_dedup.bam"

set bamin=${inputDir}/${field2}${dedup_suffix}
set bamout=${outputDir}/${field2}${split_suffix}


# Create the output directory
mkdir -p ${outputDir}


java -jar $GATK SplitNCigarReads \
   --use-original-qualities \
   -R $genome \
   -I $bamin \
   -O $bamout


# In GATK4 the following command doesn't work
# It seems that it automatically rescales MAP quality of single reads
# https://github.com/bcbio/bcbio-nextgen/issues/2163

# java -jar $GATK -T SplitNCigarReads \
#    -R $genome \
#    -I $bamin -o $bamout \
#    -rf ReassignOneMappingQuality \
#    -RMQF 255 \
#    -RMQT 60 \
#    -U ALLOW_N_CIGAR_READS \
#    -fixMisencodedQuals


# see 
# https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl