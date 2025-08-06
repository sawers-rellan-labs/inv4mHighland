#!/bin/csh
#BSUB -J star_index
#BSUB -q sara
#BSUB -n 24
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=128GB]"
#BSUB -W 24:00
#BSUB -o %J.out
#BSUB -e %J.err

# Load required modules
#module load STAR/2.7.7a-foss-2020b

set STAR=/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR

# Set the path to the genome reference directory
# in this directory the index SAindex will be written
set genomeDir=/rsstu/users/r/rrellan/sara/ref/NAM5_CHR

# Set the path to the annotation file
set annotation=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf

# Set the path to the output directory
# this only writes Log.out 
set outputDir=aln_out

# Set the number of threads
set threads=24



# Create the output directory
mkdir -p ${outputDir}

# Create the STAR index
${STAR} \
  --runThreadN ${threads} \
  --runMode genomeGenerate \
  --genomeDir ${genomeDir} \
  --genomeFastaFiles ${genomeDir}/*.fa \
  --sjdbGTFfile ${annotation} \
  --sjdbOverhang 100 \
  --outFileNamePrefix ${outputDir}/


# This only needs 30 GB (28GB used max) ~ 1 hour indexing
# In outputDir it writes only a log
#
# -bash-4.2$ ls -aFlhtr aln_out/
# total 897K
# drwxr-sr-x. 3 frodrig4 maize 4.0K Apr  5 15:55 ../
# drwxr-sr-x. 2 frodrig4 maize 4.0K Apr  5 16:08 ./
# -rw-r--r--. 1 frodrig4 maize 864K Apr  5 16:08 Log.out
#
# In genomeDir it writes the index SAindex
#
# Writing 1565873491 bytes into /rsstu/users/r/rrellan/sara/ref/NAM5_CHR//SAindex ; empty space on disk = 7344799875072 bytes ... done
# Apr 05 16:08:01 ..... finished successfully
# DONE: Genome generation, EXITING