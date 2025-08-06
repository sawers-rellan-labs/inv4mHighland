library(dplyr)

setwd("/rsstu/users/r/rrellan/DOE_CAREER/RNAseq/simreadsPT")
last <- read.table(file='/rsstu/users/r/rrellan/DOE_CAREER/inv4m_microsynteny/B73_v5.PT.last.filtered',sep='\t', header=FALSE)

last$B73_gene <- gsub("_T.*","",last$V1,perl=TRUE)


expressed <-read.csv("predictor_effects.csv") %>%
  filter(predictor=="leaf_tissue") %>% select(B73_gene="gene")

colnames(expressed)[1] <- "B73_gene"
b73_pt <- inner_join(expressed,last,relationship = "many-to-many") 

pt_genes <-last$V2 %>% sort %>%unique()
length(pt_genes)
cat(pt_genes,file="PT_gene_list", sep="\n")
