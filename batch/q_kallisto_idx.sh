#!/usr/bin/csh
#BSUB -J star_index
#BSUB -n 1
#BSUB -R "rusage[mem=2GB]"
#BSUB -W 00:20:00
#BSUB -o %J.stdout
#BSUB -e %J.stderr


conda activate /usr/local/usrapps/maize/sorghum/conda/envs/rnaseq

set index=Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx
set transcripts=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa


kallisto index -i ${index} ${transcripts}