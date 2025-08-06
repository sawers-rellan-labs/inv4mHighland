#!/usr/bin/csh

# Load required modules
# module load STAR/2.7.7a-foss-2020b

set STAR=/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR

# Set the path to the genome reference directory
set genomeDir=/rsstu/users/r/rrellan/sara/ref/NAM5_CHR

# Set the path to the annotation file
set annotation=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf

# Set the path to the fastq files directory
set fastqDir=/rsstu/users/r/rrellan/sara/RNAseq/RellanAlvarez


# Set the path to the output directory
# this only writes Log.out 
set outputdir=counts_out
set bamdir=aln_out
# Get the sample name as a command line argument
set sample=$1

set bamfile=$bamdir/${sample}Aligned.sortedByCoord.out.bam

# Set the paths to the fastq files for the current sample
set read1=${fastqDir}/${sample}${read1Suffix}
set read2=${fastqDir}/${sample}${read2Suffix}

# Set the output file name for the current sample
set outputFile=${outputDir}/${sample}

# Create the output directory
mkdir -p ${outputDir}

# Run STARcount on the current sample's bam files
${STAR} \
     --quantMode GeneCounts \
     --runMode inputAlignmentsFromBAM \
     --inputBAMfile ${bamfile}\
     --runThreadN 8 \
     # --genomeDir ${genomeDir} \
     # --sjdbGTFfile ${annotation} \
     # --readFilesCommand zcat \
     # --readFilesIn ${read1} ${read2} \
     # --outSAMtype BAM SortedByCoordinate \
     # --outFileNamePrefix ${outputFile} \
     # --outReadsUnmapped Fastx