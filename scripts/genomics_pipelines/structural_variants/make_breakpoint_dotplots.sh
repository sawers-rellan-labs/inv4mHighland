# make breakpoint self-dotplots

# B73 

# upstream breakpoint
blastdbcmd -db ../synteny/ref/B73/Zm-B73-REFERENCE-NAM-5.0.fa -entry chr4  -range 172561959-172882309 > B73_brkpt_up.fasta
lastdb B73_brkpt_up B73_brkpt_up.fasta
lastal -f TAB  B73_brkpt_up B73_brkpt_up.fasta > B73_brkpt_up.last
last-dotplot -x 200 -y 200 B73_brkpt_up.last B73_brkpt_up.png

# downstream breakpoint
blastdbcmd -db ../synteny/ref/B73/Zm-B73-REFERENCE-NAM-5.0.fa -entry chr4  -range 188131461-188220418 > B73_brkpt_down.fasta
lastdb B73_brkpt_down B73_brkpt_down.fasta
lastal -f TAB  B73_brkpt_down B73_brkpt_down.fasta > B73_brkpt_down.last
last-dotplot -x 200 -y 200 B73_brkpt_down.last B73_brkpt_down.png

# PT

# upstream breakpoint
blastdbcmd -db ../synteny/ref/PT/Zm-PT-REFERENCE-HiLo-1.0.fa -entry PT04  -range 173369064-173486186 > PT_brkpt_up.fasta
lastdb PT_brkpt_up PT_brkpt_up.fasta
lastal -f TAB  PT_brkpt_up PT_brkpt_up.fasta > PT_brkpt_up.last
last-dotplot -x 200 -y 200 PT_brkpt_up.last PT_brkpt_up.png

# downstream breakpoint
blastdbcmd -db ../synteny/ref/PT/Zm-PT-REFERENCE-HiLo-1.0.fa -entry PT04  -range 186925483-187092654 > PT_brkpt_down.fasta
lastdb PT_brkpt_down PT_brkpt_down.fasta
lastal -f TAB  PT_brkpt_down PT_brkpt_down.fasta > PT_brkpt_down.last
last-dotplot -x 200 -y 200 PT_brkpt_down.last PT_brkpt_down.png


# TIL18

# upstream breakpoint
blastdbcmd -db ../synteny/ref/mexicana/Zx-TIL18-REFERENCE-PanAnd-1.0.fa -entry chr4  -range 180269950-180365316 > TIL18_brkpt_up.fasta
lastdb TIL18_brkpt_up TIL18_brkpt_up.fasta
lastal -f TAB  TIL18_brkpt_up TIL18_brkpt_up.fasta > TIL18_brkpt_up.last
last-dotplot -x 200 -y 200 TIL18_brkpt_up.last TIL18_brkpt_up.png

# downstream breakpoint
blastdbcmd -db ../synteny/ref/mexicana/Zx-TIL18-REFERENCE-PanAnd-1.0.fa -entry chr4  -range 193570651-193734606 > TIL18_brkpt_down.fasta
lastdb TIL18_brkpt_down TIL18_brkpt_down.fasta
lastal -f TAB  TIL18_brkpt_down TIL18_brkpt_down.fasta > TIL18_brkpt_down.last
last-dotplot -x 200 -y 200 TIL18_brkpt_down.last TIL18_brkpt_down.png

