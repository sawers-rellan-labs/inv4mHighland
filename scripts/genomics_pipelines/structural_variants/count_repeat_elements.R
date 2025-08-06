# extract tandem array connected components
library(dplyr)
library(rtracklayer)
reps <- rtracklayer::import("~/Desktop/Zm-B73-REFERENCE-NAM-5.0.TE.gff3")
class(reps)

b73_up_reps <- reps[seqnames(reps) == "chr4" & start(reps) < 172882309 & end(reps) > 172561959 ]
table(b73_up_reps$Name) %>% sort()
table(b73_up_reps$Name) %>% prop.table() %>% sort()
length(b73_up_reps$Name)
# upstream  636, knob180,  438 (69%), TR-1 29 (5%), 

b73_down_reps <- reps[seqnames(reps) == "chr4" & start(reps) < 188220418  & end(reps) > 188131461 ]
table(b73_down_reps$Name) %>% sort()
table(b73_down_reps$Name) %>% sort() %>% prop.table()
length(b73_down_reps$Name)

# downstream  81, knob180, 37 (46%)

# upstream  636, knob180,  438 (69%), TR-1 29 (5%)
# downstream  81, knob180, 37 (46%)



b73_down_reps <- reps[seqnames(reps) == "chr4" & start(reps) < 172882309 & end(reps) > 172561959 ]
table(b73_up_reps$Name) %>% sort()


lastcols <- c("score","name1","start1","alnSize1","strand1","seqSize1","name2","start2","alnSize2","strand2","seqSize2","blocks")

last  <- read.table("~/Desktop/B73_brkpt_up.last", sep= "\t", header = FALSE) %>% arrange(V7)

colnames(knob_mex) <- blastcols
knob_mex$bitscaled <- scales::rescale(knob_mex$bitscore, to=c(0,1))
# TR-1_b73.blast      TR-1_hsu_2002.fasta TR-1_pt.blast       TR-1_til18.blast
tr1_mex <- read.table("~/Desktop/TR-1_til18.blast", sep= "\t", header = FALSE)
colnames(tr1_mex) <- blastcols 
tr1_mex$bitscaled <- scales::rescale(tr1_mex$bitscore, to=c(0,1))


to_rle <- inv4m_knobs %>%
  filter(qstart > 172000000, qend <200000000) %>%
  arrange(qstart) %>%
  mutate(qmid = round((qstart+qend)/2)) %>% 
  mutate( d = qmid- lag(qmid,default = 0)) %>%
  mutate(is_linked= if_else(d<100000,1,0))