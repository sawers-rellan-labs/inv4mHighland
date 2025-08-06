#!/usr/bin/tcsh

chr=$1

#activate environment

conda activate /usr/local/usrapps/maize/sorghum/conda/envs/r_env


# Set GATK path
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"


# Set the path to the genome reference
set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
set knownSites="/rsstu/users/r/rrellan/sara/ref/zea_mays.vcf"

set inputDir="split"
set outputDir="recal"

set input_suffix="_rg_added_dedup.bam"
set recal_suffix="_recal.table"



# Create the output directory
mkdir -p ${outputDir}


java -jar $GATK  BaseRecalibrator \
   -R $genome \
   --known-sites  $knownSites \
   -I $inputDir/L02_S2$input_suffix \
   -I $inputDir/L03_S3$input_suffix \
   -I $inputDir/L04_S4$input_suffix \
   -I $inputDir/L05_S5$input_suffix \
   -I $inputDir/L06_S6$input_suffix \
   -I $inputDir/L07_S7$input_suffix \
   -I $inputDir/L08_S8$input_suffix \
   -I $inputDir/L09_S9$input_suffix \
   -I $inputDir/L10_S10$input_suffix \
   -I $inputDir/L11_S11$input_suffix \
   -I $inputDir/L12_S12$input_suffix \
   -I $inputDir/L13_S13$input_suffix \
   -I $inputDir/L14_S14$input_suffix \
   -I $inputDir/L15_S15$input_suffix \
   -I $inputDir/L16_S16$input_suffix \
   -I $inputDir/L17_S17$input_suffix \
   -I $inputDir/L18_S18$input_suffix \
   -I $inputDir/L19_S19$input_suffix \
   -I $inputDir/L20_S20$input_suffix \
   -I $inputDir/L21_S21$input_suffix \
   -I $inputDir/L22_S22$input_suffix \
   -I $inputDir/L23_S23$input_suffix \
   -I $inputDir/L24_S24$input_suffix \
   -I $inputDir/L25_S25$input_suffix \
   -I $inputDir/L26_S26$input_suffix \
   -I $inputDir/L27_S27$input_suffix \
   -I $inputDir/L28_S28$input_suffix \
   -I $inputDir/L29_S29$input_suffix \
   -I $inputDir/L30_S30$input_suffix \
   -I $inputDir/L31_S31$input_suffix \
   -I $inputDir/L32_S32$input_suffix \
   -I $inputDir/R33_S33$input_suffix \
   -I $inputDir/R34_S34$input_suffix \
   -I $inputDir/R35_S35$input_suffix \
   -I $inputDir/R36_S36$input_suffix \
   -I $inputDir/R37_S37$input_suffix \
   -I $inputDir/R38_S38$input_suffix \
   -I $inputDir/R39_S39$input_suffix \
   -I $inputDir/R40_S40$input_suffix \
   -I $inputDir/R45_S41$input_suffix \
   -I $inputDir/R46_S42$input_suffix \
   -I $inputDir/R47_S43$input_suffix \
   -I $inputDir/R48_S44$input_suffix \
   -I $inputDir/L49_S45$input_suffix \
   -I $inputDir/L50_S46$input_suffix \
   -I $inputDir/L51_S47$input_suffix \
   -I $inputDir/L52_S48$input_suffix \
   -I $inputDir/L53_S49$input_suffix \
   -I $inputDir/L54_S50$input_suffix \
   -I $inputDir/L55_S51$input_suffix \
   -I $inputDir/L56_S52$input_suffix \
   -I $inputDir/L57_S53$input_suffix \
   -I $inputDir/L58_S54$input_suffix \
   -I $inputDir/L59_S55$input_suffix \
   -I $inputDir/L60_S56$input_suffix \
   -I $inputDir/L61_S57$input_suffix \
   -I $inputDir/L62_S58$input_suffix \
   -I $inputDir/L63_S59$input_suffix \
   -I $inputDir/L64_S60$input_suffix \
   -O $outputDir/${chr}_recal.table

conda deactivate