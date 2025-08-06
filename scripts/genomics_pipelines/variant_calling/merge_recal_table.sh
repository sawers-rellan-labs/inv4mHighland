#!/usr/bin/tcsh
#BSUB -q sara
#BSUB -J recalibrate
#BSUB -W 12:00
#BSUB -n 2
#BSUB -R "rusage[mem=2GB]"
#BSUB -o %J.stdout
#BSUB -e %J.stderr

#activate environment

conda activate /usr/local/usrapps/maize/sorghum/conda/envs/r_env


# Set GATK path
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"


java  -jar $GATK  GatherBQSRReports \
   -I recal/chr1_recal.table \
   -I recal/chr2_recal.table \
   -I recal/chr3_recal.table \
   -I recal/chr4_recal.table \
   -I recal/chr5_recal.table \
   -I recal/chr6_recal.table \
   -I recal/chr7_recal.table \
   -I recal/chr8_recal.table \
   -I recal/chr9_recal.table \
   -I recal/chr10_recal.table \
   -O recal/recal.table 