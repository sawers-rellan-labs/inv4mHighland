#!/usr/bin/tcsh

# Get a command line arguments
set sampleid=$1

# Load required modules
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"


# Set the path to the genome reference
set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

set recalDir="recal"

set arg_suffix=".list"
set bam_suffix="_recal.bam"

set argfile=${recalDir}/${sampleid}${arg_suffix}
set bamout=${recalDir}/${sampleid}${bam_suffix}
echo $bamout

# index vcf beforehand
# java -jar $GATK IndexFeatureFile \
#    -I /rsstu/users/r/rrellan/sara/ref/zea_mays.vcf



# GATK4: ApplyBQSR
java -jar $GATK  ApplyBQSR \
   -R $genome  \
   --arguments_file $argfile \
   --bqsr-recal-file $recalDir/recal.table \
   -O $bamout

