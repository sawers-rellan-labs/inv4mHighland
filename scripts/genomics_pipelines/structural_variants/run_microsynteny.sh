# format gff
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3 -o B73_v4.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o B73_v5.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3 -o mexicana.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.gff3 -o parviglumis.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only  Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.gff3 -o PT.bed

# format fasta
python -m jcvi.formats.fasta format Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.cds.fa B73.cds
python -m jcvi.formats.fasta format Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical.cds.fa B73_v5.cds
python -m jcvi.formats.fasta format Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.cds.fa mexicana.cds
python -m jcvi.formats.fasta format Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.cds.fa parviglumis.cds
python -m jcvi.formats.fasta format Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.cds.fa PT.cds


# Pairwise synteny search
ls *.???
python -m jcvi.compara.catalog ortholog B73_v5 mexicana --cscore=.99 --no_strip_names

ls B73.mexicana.*

# Pairwise synteny visualization
python -m jcvi.graphics.dotplot B73.mexicana.anchors



#synteny 1:1
python -m jcvi.compara.synteny depth --histogram B73.mexicana.anchors


# Macrosynteny visualization
#seqids file
Chr4
chr4

#layout file

# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , B73, top, B73.bed
 .4,     .1,    .8,       0,      , mexicana, top, mexicana.bed
# edges
e, 0, 1, B73.mexicana.anchors.simple


python -m jcvi.compara.synteny screen --minspan=30 --simple B73.mexicana.anchors B73.mexicana.anchors.new 

python -m jcvi.graphics.karyotype seqids layout


# Microsynteny visualization

python -m jcvi.compara.synteny mcscan B73.bed B73.mexicana.lifted.anchors --iter=1 -o B73.mexicana.i1.blocks


grep -A 25 -B 25 Zm00001eb191790 B73.mexicana.i1.blocks > blocks

# blocks.layout file

# x,   y, rotation,   ha,     va,   color, ratio,            label
0.5, 0.6,        0, left, center,       m,     1,       B73 chr4
0.5, 0.4,        0, left, center, #fc8d62,     1, mexicana chr4
# edges
e, 0, 1


python -m jcvi.formats.bed merge B73.bed mexicana.bed -o B73_mexicana.bed
python -m jcvi.graphics.synteny blocks B73_mexicana.bed blocks.layout

python -m jcvi.graphics.synteny blocks B73_mexicana.bed blocks.layout --glyphcolor=orthogroup



python -m jcvi.graphics.synteny blocks  B73_mexicana.bed blocks.layout --genelabelsize=5 --genelabels Zm00001eb191790_T013,Zm00001eb191820_T001,Zx00002aa015800_T001


# add parviglumis

python -m jcvi.compara.catalog ortholog B73 parviglumis --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny mcscan B73.bed B73.parviglumis.lifted.anchors --iter=1 -o B73.parviglumis.i1.blocks


# merge blocks
python -m jcvi.formats.base join B73.mexicana.i1.blocks B73.parviglumis.i1.blocks --noheader | cut -f1,2,4,6 > B73.blocks

grep -A 25 -B 25 Zm00001eb191790 B73.blocks > blocks2

# blocks2.layout

# x,   y, rotation,     ha,     va, color, ratio,            label
0.5, 0.2,        0, left, center,      ,    1, B73 Chr4
0.5, 0.6,        0, left, center,      ,     1, parviglumis Chr4
0.5, 0.4,        0, left, center,      ,    1, mexicana Chr4
# edges
e, 0, 1
e, 1, 2

python -m jcvi.formats.bed merge B73.bed mexicana.bed parviglumis.bed -o B73_mexicana_parviglumis.bed
python -m jcvi.graphics.synteny blocks2 B73_mexicana_parviglumis.bed blocks2.layout --glyphstyle=arrow  --genelabelsize=5 --genelabels Zm00001eb191790_T013,Zm00001eb191820_T001,Zx00002aa015800_T001,Zv00001aa016181_T001



# Macrosynteny visualization

python -m jcvi.compara.synteny screen --minspan=30 --simple B73.parviglumis.anchors B73.parviglumis.new
ls B73.mexicana.*



#synteny 1:1
python -m jcvi.compara.synteny depth --histogram B73.mexicana.anchors


#seqids file
Chr4
chr4
chr4
PT04

#layout2 file

# y, xstart, xend, rotation, color, label, va, bed
 .3,     .2,    .9,       0, gold, B73, center, B73.bed
 .5,     .2,    .9,       0, indigo, mex, center, mexicana.bed
 .2,     .2,    .9,       0, gold, parv, center, parviglumis.bed
 .4,     .2,    .9,       0, indigo, PT, center, PT.bed
# edges
e, 0, 2, B73.parviglumis.anchors.simple
e, 0, 3, B73.PT.anchors.simple
e, 3, 1, PT.mexicana.anchors.simple



python -m jcvi.compara.synteny screen --minspan=30 --simple B73.parviglumis.anchors B73.parviglumis.anchors.new 


# add PT
python -m jcvi.compara.catalog ortholog B73 PT --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple B73.PT.anchors B73.PT.anchors.new 

python -m jcvi.compara.catalog ortholog PT mexicana --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple PT.mexicana.anchors PT.mexicana.anchors.new 

# search inversion
grep -A 150 -B 150 Zm00001eb191820  B73.PT.anchors| cut -f1 |xargs -I {} grep {} B73.PT.anchors.simple

# this is the inversion add g* to make it green
# Zm00001eb190470_T001	Zm00001eb194780_T001	Zm00109aa017631_T001	Zm00109aa018009_T001	379	-

perl -i -pe 's/^Zm00001eb190470_T001/indigo\*Zm00001eb190470_T001/' B73.PT.anchors.simple


python -m jcvi.graphics.karyotype   --basepair  --keep-chrlabels seqids2 layout2

# merge blocks
python -m jcvi.formats.base join B73.blocks B73.PT.i1.blocks --noheader | cut -f1,2,3,5,6 > with_PT.blocks

# search jmj
grep -A 25 -B 25 Zm00001eb191790 with_PT.blocks > blocks2

# blocks2.layout

# x,   y, rotation, ha, va, color, ratio, label
0.5, 0.3, 0, left, center, gold, 1, B73
0.5, 0.5, 0, left, center, indigo, 1, TIL18
0.5, 0.2, 0, left, center, gold, 1, TIL01
0.5, 0.4, 0, left, center, indigo, 1, PT 
# edges
e, 1, 3
e, 3, 0
e, 0, 2

python -m jcvi.formats.bed merge B73.bed mexicana.bed parviglumis.bed PT.bed -o B73_mexicana_parviglumis_PT.bed
python -m jcvi.graphics.synteny blocks2 B73_mexicana_parviglumis_PT.bed blocks2.layout --glyphstyle=arrow  --genelabelsize=5 --genelabels Zm00001eb191790_T013,Zm00001eb191820_T001,Zx00002aa015800_T001,Zv00001aa016181_T001,Zm00109aa017906_T001


# blocks3.layout
# x,   y, rotation, ha, va, color, ratio, label
0.5, 0.3, 0, left, center, gold, 1, B73
0.5, 0.4, 0, left, center, indigo,  1, PT
0.5, 0.5, 0, left, center, indigo,    1, TIL18
# edges
e, 0, 1
e, 1, 2


python -m jcvi.formats.base join B73.PT.i1.blocks B73.mexicana.i1.blocks --noheader | cut -f1,2,4,6 > B73.PT.mexicana.blocks


grep -A 25 -B 25 Zm00001d051961 B73.PT.mexicana.blocks > blocks3

#grep -A 25 -B 25 Zm00001eb191790 B73.PT.mexicana.blocks > blocks3

python -m jcvi.formats.bed merge B73.bed PT.bed mexicana.bed -o B73_PT_mexicana.bed

python -m jcvi.graphics.synteny blocks3 B73_PT_mexicana.bed blocks3.layout --glyphstyle=arrow  
python -m jcvi.graphics.synteny blocks3 B73_PT_mexicana.bed blocks3.layout --glyphstyle=arrow  --genelabelsize=5 --genelabels Zm00001eb191790_T013,Zm00001eb191820_T001,Zx00002aa015800_T001,Zm00109aa017906_T001



#ra1
grep Zm00001eb312340 with_PT.blocks
grep -A 25 -B 25 Zm00001eb312340 with_PT.blocks > blocks.ra1
grep Zm00001eb312340 blocks.ra1



python -m jcvi.graphics.synteny blocks.ra1 B73_mexicana_parviglumis_PT.bed blocks.ra1.layout --glyphstyle=arrow  --genelabelsize=5 --genelabels Zm00001eb312340_T001,Zx00002aa025486_T001,Zv00001aa026720_T001,Zm00109aa029061_T001,Zm00109aa029065_T001

#seqids file
chr4
chr4
chr4

#layout2 file

# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , B73 Chr4, top, B73.bed
 .4,     .1,    .8,       0,      , parviglumis Chr4, top, parviglumis.bed
 .2,     .1,    .8,       0,      , mexicana Chr4, top, mexicana.bed
# edges
e, 0, 1, B73.parviglumis.anchors.simple
e, 1, 2, B73.mexicana.anchors.simple


python -m jcvi.compara.synteny screen --minspan=30 --simple B73.PT.anchors B73.PT.anchors.new 

python -m jcvi.graphics.karyotype seqids layout2
