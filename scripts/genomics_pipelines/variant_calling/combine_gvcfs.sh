
#!/usr/bin/tcsh
conda activate /usr/local/usrapps/maize/sorghum/conda/envs/r_env



set outputDir="gvcf"

set genome="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
set GATK="/usr/local/usrapps/maize/gatk-4.5.0.0/gatk"


# Create the output directory
mkdir -p ${outputDir}

$GATK --java-options "-Xms6g" CombineGVCFs \
   -R $genome \
   -V    gvcf/3056.g.vcf.gz \
   -V    gvcf/3059.g.vcf.gz \
   -V    gvcf/3080.g.vcf.gz \
   -V    gvcf/3083.g.vcf.gz \
   -V    gvcf/3092.g.vcf.gz \
   -V    gvcf/3095.g.vcf.gz \
   -V    gvcf/3113.g.vcf.gz \
   -V    gvcf/3131.g.vcf.gz \
   -V    gvcf/4080.g.vcf.gz \
   -V    gvcf/4083.g.vcf.gz \
   -V    gvcf/4107.g.vcf.gz \
   -V    gvcf/4113.g.vcf.gz \
   -V    gvcf/4122.g.vcf.gz \
   -V    gvcf/4125.g.vcf.gz \
   -V    gvcf/4131.g.vcf.gz \
   -O    gvcf/inv4m_PSU.g.vcf.gz


