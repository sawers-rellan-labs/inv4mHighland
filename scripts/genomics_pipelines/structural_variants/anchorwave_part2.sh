#!/bin/bash
source ~/.bashrc
conda activate /usr/local/usrapps/maize/frodrig/conda/envs/mcscan


anchorwave genoAli -t 2 \
  -i Zm-B73-REFERENCE-NAM-5.0_chr4.gff3 \
  -as B73_cds.fa \
  -r Zm-B73-REFERENCE-NAM-5.0_chr4.fa \
  -a  PT_B73_cds.sam \
  -ar ref.sam \
  -s Zm-PT-REFERENCE-HiLo-1.0_chr4.fa \
  -n anchors \
  -o PT_anchorwave.maf \
  -f PT_anchorwave.f.maf \
  -IV

conda deactivate


# #!/bin/bash
# source ~/.bashrc
# conda activate /usr/local/usrapps/maize/frodrig/conda/envs/mcscan
# 
# minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.fa B73_cds.fa > TIL18_B73_cds.sam
# 
# anchorwave genoAli -t 2 \
#   -i Zm-B73-REFERENCE-NAM-5.0_chr4.gff3 \
#   -as B73_cds.fa \
#   -r Zm-B73-REFERENCE-NAM-5.0_chr4.fa \
#   -a TIL18_B73_cds.sam \
#   -ar ref.sam \
#   -s Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.fa \
#   -n TIL18.anchors \
#   -o TIL18_anchorwave.maf \
#   -f TIL18_anchorwave.f.maf \
#   -IV
# 
# conda deactivate
