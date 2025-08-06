# DEG characterization.
library(dplyr)
library(ggplot)
library(scales)
source("/Users/fvrodriguez/Projects/NCSU/06_GEA/association/sorghum/fastman/fastman.R")

myGFF <- "/System/Volumes/Data/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gff3"
newGFF <- rtracklayer::import(myGFF)

genes = as.data.frame(newGFF)  %>%
        dplyr::filter(type=="gene") 
genes$start
colnames(genes)



#### P effect ####

# DEG for GO terms
# by_leaf with local false sign rate
# lfsr <- as.data.frame(get_lfsr(results$`2`$m.Vem_ed))

# trt_deg <- names(get_significant_results(results$`2`$m.Vem_ed))

# de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)
# write.table(trt_deg, file= "P_trt_DEG.txt",
#             row.names = FALSE, col.names = FALSE,
#             quote= FALSE)


# FDR from limma
p <- r$`2`$adj.P.Val
names(p) <- row.names(r$`2`)

maxnegLogP <- -log10(min(p[p>0]))

names(p)[p<0.05]

sum(p<0.05)


to_plot <- data.frame( gene_id = names(p), P = p) %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)

to_plot %>%
  filter(CHR==4 & BP > 155e6 &  BP < 195e6 & P<0.05) %>%
  arrange(P) %>% head()


quartz(height = 4, width = 12)
par(mar=c(5,6,4,1)+.1)
fastman (to_plot, maxP = maxnegLogP,
         ylab=expression(-log[10](italic(FDR))),
         cex.axis = 2, cex.lab=2, main = "DEG Phosphorus Treatment",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))


## inv4m effect ####
p <- r$`3`$adj.P.Val
names(p) <- row.names(r$`3`)

maxnegLogP <- -log10(min(p[p>0]))

names(p)[p<0.05]

sum(p<0.05)


to_plot <- data.frame( gene_id = names(p), P = p) %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)

to_plot %>%
  filter(CHR==4 & BP > 173905023 &  BP < 187937996 & P<0.05) %>%
  arrange(P) %>% nrow()

to_plot %>%
  filter(CHR==4 & BP > 155e6 &  BP < 195e6 & P<0.05) %>%
  arrange(P) %>% nrow()

to_plot %>% arrange(P) %>% head()


quartz(height = 4, width = 12)
par(mar=c(5,6,4,1)+.1)
fastman (to_plot, maxP = maxnegLogP,
         ylab=expression(-log[10](italic(FDR))),
         cex.axis = 2, cex.lab=2, main = "Inv4m effect",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))


quartz(height = 4, width = 10)
to_plot %>%
  filter(CHR==4, BP > 150e6 & BP < 200e6) %>%
  ggplot(aes( x= BP, y= -log10(P)))+
  geom_vline(xintercept = 173905023) +
  geom_vline(xintercept = 187937996) +
  geom_point(col= "#60a500") +
  xlab("Chr4")+
  scale_x_continuous(labels =  unit_format(unit = "",scale = 1e-6))+
  ggpubr::theme_classic2(base_size = 30)

# genotype x treatment effect

# trt_deg <- names(get_significant_results(results$`4`$m.Vem_ed))
# de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)

#:) this look good
# lfsr <- as.data.frame(get_lfsr(results$`4`$m.Vem_ed))
# p <- pmin( lfsr$V1,lfsr$V2,lfsr$V3,lfsr$V4)

p <- r$`4`$adj.P.Val
names(p) <- row.names(r$`4`)
maxnegLogP <- -log10(min(p[p>0]))



to_plot <- data.frame( gene_id = names(p), P = p) %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)


# to_plot <- as.data.frame(get_lfsr(results$`4`$m.Vem_ed)) %>%
#   mutate(P=pmin(V1,V2,V3,V4)) %>%
#   tibble::rownames_to_column("gene_id") %>%
#   inner_join(genes %>%
#                dplyr::select(gene_id, CHR = "seqnames",BP = "start")
#   ) %>%
#   dplyr::select(gene_id, CHR, BP, P)


to_plot %>%
  filter(CHR==4 & BP > 155e6 &  BP < 195e6 & P<0.05) %>%
  arrange(P)

quartz(height = 4, width = 12)
par(mar=c(5,6,4,1)+.1)
fastman (to_plot, maxP = maxnegLogP,
         ylab=expression(-log[10](italic(FDR))),
         cex.axis = 1.5, cex.lab=1.5, main = "P x Inv4m effect",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))

quartz(height = 4, width = 10)
to_plot %>%
  filter(CHR==4, BP > 150e6 & BP < 200e6) %>%
  ggplot(aes( x= BP, y= -log10(P)))+
  geom_vline(xintercept = 173905023) +
  geom_vline(xintercept = 187937996) +
  geom_point(col= "#60a500") +
  xlab("Chr4")+
  scale_x_continuous(labels =  unit_format(unit = "",scale = 1e-6))+
  ggpubr::theme_classic2(base_size = 30)

