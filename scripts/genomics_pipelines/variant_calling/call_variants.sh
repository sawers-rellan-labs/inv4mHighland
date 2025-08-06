#!/usr/bin/tcsh

# Get a command line arguments
set sampleid=$1

# Load required modules
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"


# Set the path to the genome reference
set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

set inputDir="recal"
set outputDir="hc"

set bam_suffix="_recal.bam"
set vcf_suffix=".vcf.gz"

set bamin=${inputDir}/${sampleid}${bam_suffix}
set vcfout=${outputDir}/${sampleid}${gvcf_suffix}


# Create the output directory
mkdir -p ${outputDir}


# GATK4: call with haplotypecaller

java -Xmx6g -jar $GATK HaplotypeCaller \
    -R ${genome} \
    -I ${bamin} \
    -O ${vcfout} \
    -dont-use-soft-clipped-bases \
    --stand-call-conf 20.0 
