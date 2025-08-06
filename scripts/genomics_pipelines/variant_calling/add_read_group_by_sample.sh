#!/bin/tcsh

# Specify your input file (replace 'your_file.txt' with your actual file name)
set picard="/usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar"
set input_file=$1
set bam_suffix="Aligned.sortedByCoord.out.bam"
set rg_suffix="_rg_added_sorted.bam"

set aln_dir="aln_out"

# Loop through each line in the file
# (I feel dirty)
foreach line ("`cat $input_file`")
    echo $line
    # Extract the fields 2 and 3
    set field2 = `echo $line | awk '{print $2}'`
    echo "field2" $field2
    set field3 = `echo $line | awk '{print $3}'`
    echo "field3" $field3
    set bamin={$aln_dir}/{$field2}{$bam_suffix}
    set bamout={$aln_dir}/{$field2}{$rg_suffix}

    # Execute the desired command with the extracted fields
   java -jar picard AddOrReplaceReadGroups I=$bamin O=$bamout SO=coordinate RGID=1 RGLB=$field2 RGPL=illumina RGPU=L002 RGSM=$field3
end

AddOrReplaceReadGroups -I aln_out/L01_S1_L002Aligned.sortedByCoord.out.bam -O aln_out/L01_S1_L002_rg_added_sorted.bam -SO coordinate -RGID 1 -RGLB L01_S1_L002 -RGPL illumina -RGPU L002 -RGSM 3056


java -jar /usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar AddOrReplaceReadGroups I=aln_out/L01_S1_L002Aligned.sortedByCoord.out.bam O=aln_out/L01_S1_L002_rg_added_sorted.bam SO=coordinate RGID=1 RGLB=L01_S1_L002 RGPL=illumina RGPU=L002 RGSM=3056


java -jar picard AddOrReplaceReadGroups I=output.sam O=rg_added_sorted.bam SO=coordinate RGID=ID_NAME RGLB=library RGPL=illumina RGPU=identifier RGSM=sample_name

Mark duplicate reads.

java -jar /usr/local/share/picard-tools/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

Identify and split Cigar N Reads and reassign quality scores.

java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T SplitNCigarReads -R /path/to/genome/fasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -fixMisencodedQuals

Perform BaseRecalibration.
Calibration files can be found here for hg38 ​Recalibration Files
NOTE: Calibration Files are only available for a few genomes (Human, Mouse, etc).
java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T BaseRecalibrator -R /path/to/genome/fasta -I dedupped.bam -knownSites  /path/to/calibration/files -o recalibration.table

java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T PrintReads  -R /path/to/genome/fasta -I dedupped.bam -BQSR recalibration.table -o recalibrated.bam

3 - Call and Filter Variants.

java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T HaplotypeCaller  -R /path/to/genome/fasta -I  recalibrated.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o Variants_called.vcf

java -jar /usr/local/gatk3/GenomeAnalysisTK.jar -T VariantFiltration -R /path/to/genome/fasta -V  Variants_called.vcf -window 35 -cluster 3 -filterName Filter -filter "QD < 2.0" -filterName Filter -filter "FS > 30.0" -o Filtered_variants_called.vcf

