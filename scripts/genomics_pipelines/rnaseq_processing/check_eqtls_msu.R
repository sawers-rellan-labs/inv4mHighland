trans <- read.table(file="~/Desktop/trans_L3Tip_eQTL_results.txt", sep="\t", header= TRUE)

quartz()
hist(trans$snp_distance_to_TSS %>% log10())

trans%>%
  filter(gene_id %in% targets$v4 )  %>% pull(gene_id) %>% table()
  

targets%>%
  filter(v4=="Zm00001d038675")
