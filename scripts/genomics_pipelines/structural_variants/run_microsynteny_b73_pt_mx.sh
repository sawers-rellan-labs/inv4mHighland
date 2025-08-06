# format gff
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.gff3 -o PT.bed
python -m jcvi.formats.gff bed --type=transcript --key=transcript_id --Parent=gene_id --primary_only Zm-B73-REFERENCE-NAM-5.0.59_plus_jmj.gff3  -o B73.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3 -o mexicana.bed


# format fasta
python -m jcvi.formats.fasta format Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.cds.fa PT.cds
python -m jcvi.formats.fasta format Zm-B73-REFERENCE-NAM-5.0.59.cds_plus_jmj.fasta B73.cds
python -m jcvi.formats.fasta format Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.cds.fa mexicana.cds



# Pairwise synteny search
ls *.???
python -m jcvi.compara.catalog ortholog B73 mexicana --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog B73 PT --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog PT mexicana --cscore=.99 --no_strip_names

ls B73.mexicana.*

# Pairwise synteny visualization
python -m jcvi.graphics.dotplot B73.mexicana.anchors
python -m jcvi.graphics.dotplot B73.PT.anchors
python -m jcvi.graphics.dotplot PT.mexicana.anchors


python -m jcvi.compara.synteny screen --minspan=30 --simple B73.mexicana.anchors B73.mexicana.anchors.new 
python -m jcvi.compara.synteny screen --minspan=30 --simple B73.PT.anchors B73.PT.anchors.new 



python -m jcvi.formats.bed merge B73.bed PT.bed mexicana.bed -o B73.PT.mexicana.bed


# 
# python -m jcvi.graphics.synteny blocks B73.mexicana.bed blocks.layout
# 
# 
# #synteny 1:1
# python -m jcvi.compara.synteny depth --histogram B73.mexicana.anchors
# 
# 
# # Macrosynteny visualization
# #seqids file
# Chr4
# chr4
# 
# #layout file
# 
# # y, xstart, xend, rotation, color, label, va,  bed
#  .6,     .1,    .8,       0,      , B73, top, B73.bed
#  .4,     .1,    .8,       0,      , mexicana, top, mexicana.bed
# # edges
# e, 0, 1, B73.mexicana.anchors.simple
# 
# 
# 
# 
# python -m jcvi.graphics.karyotype seqids layout
# 

# Microsynteny visualization

python -m jcvi.compara.synteny mcscan B73.bed B73.mexicana.lifted.anchors --iter=1 -o B73.mexicana.i1.blocks


grep -A 25 -B 25 Zm00001d051961 B73.mexicana.i1.blocks > blocks

# blocks.layout file

# x,   y, rotation,   ha,     va,   color, ratio,            label
0.5, 0.6,        0, left, center,       m,     1,       B73 chr4
0.5, 0.4,        0, left, center, #fc8d62,     1, mexicana chr4
# edges
e, 0, 1




python -m jcvi.graphics.synteny blocks B73_mexicana.bed blocks.layout --glyphcolor=orthogroup



python -m jcvi.graphics.synteny blocks  B73_mexicana.bed blocks.layout --genelabelsize=5 --genelabels Zm00001eb191790_T013,Zm00001eb191820_T001,Zx00002aa015800_T001


# add PT

python -m jcvi.compara.synteny mcscan B73.bed B73.PT.lifted.anchors --iter=1 -o B73.PT.i1.blocks


# merge blocks
python -m jcvi.formats.base join B73.PT.i1.blocks B73.mexicana.i1.blocks  --noheader | cut -f1,2,4,6 > B73.blocks

grep -A 25 -B 25 Zm00001d051961 B73.blocks > blocks2

# blocks2.layout

# x,   y, rotation, ha, va, color, ratio, label
0.5, 0.3, 0, left, center, gold, 1, B73
0.5, 0.4, 0, left, center, indigo,  1, PT
0.5, 0.5, 0, left, center, indigo,    1, TIL18
# edges
e, 0, 1
e, 1, 2

python -m jcvi.formats.bed merge B73.bed PT.bed mexicana.bed  -o B73_PT_mexicana.bed
python -m jcvi.graphics.synteny blocks2 B73_PT_mexicana.bed blocks2.layout --glyphstyle=arrow  --genelabelsize=5
python -m jcvi.graphics.synteny blocks2 B73_PT_mexicana.bed blocks2.layout --glyphstyle=arrow  --genelabelsize=5 --genelabels Zm00001d051995_T001,Zx00002aa015800_T001,Zm00109aa017906_T001


