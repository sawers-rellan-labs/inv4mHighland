blastdbcmd -db ../synteny/ref/B73/Zm-B73-REFERENCE-NAM-5.0.fa -entry chr4 > Zm-B73-REFERENCE-NAM-5.0_chr4.fa
blastdbcmd -db ../synteny/ref/mexicana/Zx-TIL18-REFERENCE-PanAnd-1.0.fa -entry chr4 > Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.fa
blastdbcmd -db ../synteny/ref/PT/Zm-PT-REFERENCE-HiLo-1.0.fa -entry PT04 > Zm-PT-REFERENCE-HiLo-1.0_chr4.fa
perl -i -pe 's/^>PT04/>chr4/'  Zm-PT-REFERENCE-HiLo-1.0_chr4.fa

grep chr4 ../microsynteny/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 >Zm-B73-REFERENCE-NAM-5.0_chr4.gff3
grep chr4 ../microsynteny/Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3 >Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.gff3
grep PT04 ../microsynteny/Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.gff3 > Zm-PT-REFERENCE-HiLo-1.0_chr4.gff3
perl -i -pe 's/^PT04/chr4/' Zm-PT-REFERENCE-HiLo-1.0_chr4.gff3

anchorwave gff2seq -i Zm-B73-REFERENCE-NAM-5.0_chr4.gff3 -r Zm-B73-REFERENCE-NAM-5.0_chr4.fa -o B73_cds.fa
anchorwave gff2seq -i Zm-PT-REFERENCE-HiLo-1.0_chr4.gff3 -r Zm-PT-REFERENCE-HiLo-1.0_chr4.fa -o PT_cds.fa
anchorwave gff2seq -i Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.gff3 -r Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.fa -o TIL18_cds.fa


minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-B73-REFERENCE-NAM-5.0_chr4.fa B73_cds.fa > ref.sam
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-PT-REFERENCE-HiLo-1.0_chr4.fa B73_cds.fa > PT_B73_cds.sam
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zx-TIL18-REFERENCE-PanAnd-1.0_chr4.fa B73_cds.fa > TIL18_B73_cds.sam

 