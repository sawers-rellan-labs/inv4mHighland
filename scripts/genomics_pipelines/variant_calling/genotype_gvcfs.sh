
#!/usr/bin/tcsh
conda activate /usr/local/usrapps/maize/sorghum/conda/envs/r_env


set outputDir="vcf"
set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk"

# Create the output directory
mkdir -p ${outputDir}



$GATK --java-options "-Xmx6g" GenotypeGVCFs \
   -R $genome \
   -V gvcf/inv4m_PSU.g.vcf.gz \
   --force-output-intervals v5NAM.genes_interval_list \
   -O $outputDir/inv4m_PSU.vcf.gz