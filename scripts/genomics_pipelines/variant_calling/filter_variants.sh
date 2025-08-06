

$GATK --java-options "-Xmx6g" SelectVariants \
   -V $outputDir/inv4m_PSU.vcf.gz \
   -select-type SNP \
   --restrict-alleles-to BIALLELIC \
   -O $outputDir/inv4m_PSU_biallelic_snps.vcf.gz





$GATK --java-options "-Xmx6g" VariantFiltration \
    -V $outputDir/inv4m_PSU_biallelic_snps.vcf.gz \
    --window 35 \
    --cluster 3 \
    -filter "QD < 2.0" -filter-name "QD2" \
    -filter "FS > 30.0" -filter-name "FS30" \
    -filter "SOR > 3.0" -filter-name "SOR3" \
    -filter "MQ < 40.0" -filter-name "MQ40" \
    -O $outputDir/inv4m_PSU_filtered_snps.vcf.gz

bcftools view --apply-filters .,PASS vcf/inv4m_PSU_filtered_snps.vcf.gz | bgzip -c > vcf/inv4m_PSU_snps_masked.vcf.gz


# TASSEL filter
# Filter to just chromosomes 1-10
# Filter minimum allele count to 12
# Filter taxa with more than 50% missing data
# check matrix missing %
# Impute the markers by KNN
# check matrix missing %


