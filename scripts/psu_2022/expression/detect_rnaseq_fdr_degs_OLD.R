# Get list of differentially expressed genes for
# GWAS plotting
# Gene Ontology analysis

# Load required libraries
library(edgeR)      # Bioconductor: differential expression analysis
library(limma)      # Bioconductor: linear models for microarray/RNA-seq data
library(rtracklayer) # Bioconductor: genomic annotation handling
library(GenomicRanges) # Bioconductor: genomic ranges operations
library(dplyr)      # CRAN: data manipulation
library(ggplot2)    # CRAN: plotting
library(ggpubr)     # CRAN: publication ready plots
library(ggtext)     # CRAN: formatted text in plots


# Get expression data
counts <- read.csv('data/inv4mRNAseq_gene_sample_exp.csv')
sampleInfo <- read.csv('data/PSU-PHO22_Metadata.csv') 

tag <- sampleInfo$side_tag
names(tag) <-sampleInfo$library

genes <- data.frame(gene = counts[,2])
counts <- as.matrix(counts[,-c(1:2)])

rownames(counts) <- genes$gene

sampleNames <- tag[colnames(counts)]
colnames(counts) <- sampleNames 
sampleNames %in% sampleInfo$side_tag
sampleInfo <- sampleInfo[match(sampleNames,sampleInfo$side_tag),]

y <- DGEList(counts = counts,samples = sampleInfo)
y$group <- interaction(y$samples$Treatment,y$samples$Genotype)
keep <- filterByExpr(y,group = y$group)
y_filtered <- y[keep,]


#  There are 17 low quality libraries

# mds <- plotMDS(y_filtered,pch=21,plot = FALSE)
# plot(mds$x,mds$y,pch=21,bg = (as.factor(y_filtered$samples$lowCount)==T)+1, 
#      main = "Expression MDS")

# Filter out low quality libraries by library read count 

y_filtered$samples$lowCount <- y_filtered$samples$lib.size < 2e7

y_filtered_bySample <- y_filtered[,!y_filtered$samples$lowCount]

y_filtered_bySample$samples
with(y_filtered_bySample$samples,
     table(Genotype))

table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$leaf_tissue)
table(y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

# plot MDS

mds <- plotMDS(y_filtered_bySample,pch=21,
                label = y_filtered_bySample$samples$side_tag,
                plot=FALSE)


mds2 <- plotMDS(y_filtered_bySample,pch=21,
               label = y_filtered_bySample$samples$side_tag,
               dim.plot = c(3,4),
               plot=FALSE)

# The third dimension in the MDS separates by 
# genotype!!!

d<-y_filtered_bySample$samples
d$dim1 <- mds$x
d$dim2 <- mds$y
d$dim3 <- mds2$x
d$dim4 <- mds2$y
d$Treatment<- factor(to_plot$Treatment)
levels(d$Treatment) <- c("+P","-P")
d$Genotype <- factor(to_plot$Genotype)
d$RNA_Plant <- factor(d$RNA_Plant)

cor(d$dim3,d$Genotype %>% as.numeric)

labels<-c("CTRL", "*Inv4m*")
names(labels) <- c("CTRL","INV4")

quartz()  
d %>%
  ggplot(aes(x=dim3,y=dim4, fill=Genotype, shape = Treatment)) + 
  xlab(paste0("dim3 ","(",round(100*mds2$var.explained[3],),"%)")) +
  ylab(paste0("dim4 ","(",round(100*mds2$var.explained[4],),"%)"))+
  geom_point(size=4) +
  scale_fill_viridis_d(direction = -1, labels = labels)+
  scale_shape_manual(values=c(24,21))+
  guides(   shape = "none",
            fill = guide_legend( title = "Genotype", order = 2,
                                 #labels = c("CTRL", expression(italic(Inv4m))),
                                 override.aes = list(geom = "point", shape = 22,size=7,reverse =TRUE))) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position = c(0.89,0.9),
        legend.text = ggtext::element_markdown(),
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "line")) 




pdf(file = "~/Desktop/inv4m_expression_MDS_by_var.pdf")
ggplot(d,aes(x=dim1,y=dim2)) + 
  geom_point(aes(color = Treatment)) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
      legend.spacing =  unit(0,"line"),
      legend.box.spacing =unit(0, "in"),
      legend.position = c(0.75,0.2))

ggplot(d,aes(x=dim1,y=dim2)) + 
  geom_point(aes(color = row,shape = Treatment)) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.2))

ggplot(d,aes(x=dim1,y=dim2)) + 
  geom_point(aes(color = decimal_time )) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.2))
ggplot(d,aes(x=dim1,y=dim2)) + 
  geom_point(aes(color = COLLECTOR)) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.2))

ggplot(d,aes(x=dim1,y=dim2)) + 
  geom_point(aes(color = Genotype)) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.2))
dev.off()

png(file = "~/Desktop/inv4m_expression_MDS.png",
    width = 7,
    height = 7,
    units= "in",
    res = 300)

quartz()
d %>%
  mutate(leaf =factor(leaf_tissue)) %>%
  ggplot(aes(x=dim1,y=dim2)) + 
  xlab(paste0("dim1 ","(",round(100*mds$var.explained[1],),"%)")) +
  ylab(paste0("dim2 ","(",round(100*mds$var.explained[2],),"%)"))+
  geom_point(aes(fill = leaf,shape = Treatment),size=4) +
  scale_fill_viridis_d()+
  scale_shape_manual(values=c(24,21))+
  guides(   shape = guide_legend(title = "Treatment",order = 1,
                                 override.aes = list(size=7)),
            fill = guide_legend( title = "Leaf", order = 2, 
                                 override.aes = list(geom = "point", shape = 22,size=7))) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.17))
dev.off()

cor(d$dim1,d$Treatment %>% as.numeric)
cor(d$dim2,d$leaf_tissue)
cor(d$dim3,d$Genotype %>% as.numeric)

cor.pvalues <- c(
  cor.test(d$dim1,d$Treatment %>% as.numeric)$p.value,
  cor.test(d$dim2,d$leaf_tissue)$p.value,
  cor.test(d$dim3,d$Genotype %>% as.numeric)$p.value
)

p.adjust(cor.pvalues)

# Normalize
y_filtered_bySample <- calcNormFactors(y_filtered_bySample)


# specify  reduced model

design <- model.matrix(~ Plot_Column + Plot_Row + leaf_tissue + Treatment*Genotype,d)

voomR <- voom(y_filtered_bySample,design=design,plot= FALSE)

saveRDS(voomR$E,file="~/Desktop/normalized_expression_logCPM.rda")

saveRDS(voomR,file="~/Desktop/normalized_expression_voom_object.rda")


fit = lmFit(voomR)

ebfit <- eBayes(fit, robust=TRUE)

colSums(abs(decideTests(ebfit)))


# Variance partition analysis
# output is distribution of variance partition over all the genes
#
# library(BiocParallel)
# library(variancePartition)
# 
# # Define formula
# 
# form <- ~ Plot_Row  + Plot_Column + leaf_tissue +(1 | Treatment)+  (1 | Genotype) +  (1 | Treatment:Genotype)
# 
# # variancePartition seamlessly deals with the result of voom()
# # by default, it seamlessly models the precision weights
# # This can be turned off with useWeights=FALSE
# 
# # I must send this over to the server
# # lmer option verbose=TRUE  passed through ... is not working, 
# # it shows no output
# 
# varPart <- fitExtractVarPartModel(voomR, form, d, REML=TRUE,verbose = 2L)
# 

# Fit bayesian mixed effect model 

coef_cols <- ebfit$coefficients %>% colnames()


# Coefficients of interest in the model
of_interest <-  c("leaf_tissue","Treatment-P","GenotypeINV4","Treatment-P:GenotypeINV4")

# of_interest <-  c("leaf_tissue","Treatment-P","GenotypeINV4")

is_interesting <- which(coef_cols %in% of_interest )

# Collect results

results <-list()

for(x in of_interest){
  r <- cbind(topTable(ebfit,coef = x,sort.by = 'none',n=Inf),
             data.frame(predictor=x)) %>%
    tibble::rownames_to_column("Response")
  cr <- qt(0.975, ebfit$df.residual + ebfit$df.prior) * ebfit$stdev.unscaled[,x] * sqrt(ebfit$s2.post)
  # Calculating  confidence interval
  r$upper  <- r$logFC +cr
  r$lower  <- r$logFC -cr
  results[[x]]<-r
}


# Add gene symbol (locus name)

gene_symbol <- read.table("/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My\ Drive/repos/inv4mRNA/data/gene_symbol.tab",quote="",header=TRUE, sep ="\t", na.strings = "")
head(gene_symbol )

pannzer <-read.table("/Users/fvrodriguez/Desktop/PANNZER/PANNZER_DESC.tab",quote="",header=TRUE, sep ="\t", na.strings = "") %>%
        group_by(gene_model) %>%
        dplyr::slice(1) %>%
        dplyr::select(gene_model, desc)
     
gene_pannzer <- gene_symbol %>%
  left_join(pannzer)
  
effect_order <-  c("leaf_tissue","Treatment-P","GenotypeINV4","Treatment-P:GenotypeINV4")
# effect_order <-  c("leaf_tissue","Treatment-P","Treatment-P:GenotypeINV4")

# get gene coordinates ----
# get inversion coordinates

# I think I should leave this info in a separate source file
# so I can call the inversion coordinates as a function


v5_GFF <- "~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
v5<- rtracklayer::import(v5_GFF )  %>%
  subset(type=="gene"  & seqnames %in% 1:10) 
B73 <- v5
genes <- as.data.frame(B73) 
nrow(genes)

genes$ID <- gsub("gene:","",genes$ID)

# Extract these coordinates from the MCScan data

inv4m_start <- genes[genes$ID=="Zm00001eb190470","start"]
inv4m_start
# 172883675

inv4m_end <- genes[genes$ID=="Zm00001eb194800","end"]
#TIL18 Zx00002aa015905
inv4m_end
# 188132113


inv4m_gene_ids <- genes %>%
  filter(
    seqnames==4,
    start >= inv4m_start,
    end <= inv4m_end
  ) %>% pull(ID)

length(inv4m_gene_ids)

# checked manually in the RNAseq genotype table v5
introgression_start <- 157012149
introgression_end <- 195900523

shared_introgression_gene_ids<-  genes %>%
  filter(
    seqnames==4,
    start >= introgression_start,
    end <= introgression_end
  ) %>% pull(ID)

length(shared_introgression_gene_ids)

flanking_introgression_gene_ids <- shared_introgression_gene_ids[!(shared_introgression_gene_ids %in% inv4m_gene_ids)]

length(flanking_introgression_gene_ids)


# Make  a single dataframe  for size effects

# the conditions here are for a very specific type of analysis
# but I will leave them as default for now

effects <- results %>% dplyr::bind_rows() %>%
  mutate(predictor = factor(predictor, levels=effect_order))%>%
  mutate(P = adj.P.Val) %>%
  mutate(neglogP = -log10(adj.P.Val)) %>%
  mutate(is_significant = adj.P.Val<0.05) %>%
  mutate(upregulated = case_when(
    (predictor !="leaf_tissue" ) &   (logFC > 2)    & is_significant ~ TRUE,
    (predictor =="leaf_tissue" ) &   (logFC > 0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(downregulated = case_when(
     (predictor !="leaf_tissue" ) & (logFC < -2)   & is_significant ~ TRUE,
     (predictor =="leaf_tissue" ) & (logFC < -0.7) & is_significant ~ TRUE,
     .default = FALSE))  %>%
  mutate(regulation = case_when(
      is_significant & upregulated ~ "Upregulated",
      is_significant & downregulated ~ "Downregulated",
      .default = "Unregulated"
    )) %>%
  mutate(is_DEG= is_significant & regulation != "Unregulated") %>%
  left_join(gene_pannzer, by =c(Response="gene_model"), relationship = "many-to-many") %>%
  mutate(desc_merged = coalesce(locus_name, desc)) %>%
  dplyr::select(predictor,regulation,Response,locus_symbol, desc_merged,  everything())  %>%
  dplyr::rename(gene="Response") %>%
  inner_join(
    genes %>% as.data.frame() %>%
      dplyr::select(gene=gene_id, 
                    CHR=seqnames,
                    BP= start) %>%
      mutate(CHR=CHR %>% as.character() %>% as.integer())
  ) 
colnames(effects)


# quartz()
# boxplot(data=effects,logFC~regulation)

# leave this function to another file so I source it

add_malahanobis_outliers <- function(
    data = NULL,
    distance_quantile = 0.05,
    FDR = 0.05
){
  data <- lapply(
    split(data, factor(data$predictor)),
    function(per_predictor){
      bivariate <- per_predictor %>% dplyr::select(logFC, neglogP)
      per_predictor$mahalanobis <- mahalanobis(
        x =bivariate,
        #center = c(0.0),
        center = colMeans( bivariate) ,
        cov = cov(bivariate))
      per_predictor
    }
  ) %>% dplyr::bind_rows()
  
  cutoff <- qchisq(p = 1-distance_quantile, df = 2)
  
  data$is_mh_outlier <- (data$adj.P.Val < FDR) & (data$mahalanobis > cutoff)
  
  data %>%
    ungroup%>%
    filter() %>%
    group_by(predictor,regulation) %>%
    arrange(-mahalanobis, .by_group = TRUE) %>%
    ungroup()
}

# effects$odistance <- sqrt(effects$logFC^2+ effects$neglogP^2)

# 5%
# 5% Malahanobis outliers
# add column  for malahanobis outliers



# 38.9




#  P-value from limma -----


# effects$neglogP <- -log10(effects$P) 
#   
# bonferroni_threshold <- 0.05/nrow(effects %>% filter(predictor=="GenotypeINV4"))
# 
# effects$is_significant <- effects$P < bonferroni_threshold
# 
# table(effects$predictor, effects$is_significant)
# 
# # Difference in effect size with QQplot and Kolmogorov-Smirnov test ######
# # no significant differece!
# # effects$is_regulated <-  effects$regulation != "Unregulated
# 
# effects$in_shared <- effects$CHR == 4  & (effects$BP >= 157012149) & (effects$BP <= 195900523)
# 
# effects$in_inv4m <-   effects$CHR == 4 & (effects$BP >= inv4m_start) & (effects$BP <= inv4m_end)
# 
# 
# effects$in_flanking <-  effects$in_shared & !effects$in_inv4m
# 
# sum(effects$in_flanking)
# levels(factor(effects$predictor))
# 
# effects <- effects %>%
#   mutate(region = case_when(
#     in_inv4m ~ "inv4m",
#     in_flanking ~ "flanking",
#     .default= "outside")
#   )
# 
# effect_inv4m <- effects %>%
#   dplyr::filter(predictor=="GenotypeINV4") %>%
#   arrange(CHR,BP) %>%
#   droplevels()
# 
# contingency <- with(
#   effect_inv4m %>% 
#     droplevels(),
#   table(regulation, region)
# )
# 
# contingency 
# 
# chisq.test(contingency)
# 
# contingency <- with(
#   effect_inv4m %>% 
#     filter(regulation !="Unregulated") %>% 
#     droplevels(),
#   table(regulation, region)
# )
# chisq.test(contingency)
# 
# contingency
# 
# contingency <- with(
#   effect_inv4m %>% 
#     filter(regulation !="Unregulated" &  region !="outside") %>% 
#     droplevels(),
#   table(regulation, region)
# )
# chisq.test(contingency)
# 
# contingency
# 
# effects %>%
#   filter(is_significant & regulation != "Unregulated", predictor=="Treatment-P") %>%
#   group_by(predictor,regulation) %>% dplyr::slice(1:10) %>% print(n=200)


to_table <- add_malahanobis_outliers(effects) %>%
  filter(is_significant & regulation != "Unregulated") %>%
  group_by(predictor,regulation) %>% dplyr::slice(1:10) %>%
#  filter(!is.na(locus_name)) %>%
 dplyr::select(predictor, gene, locus_symbol,desc_merged,logFC,neglogP,mahalanobis) %>% 
  arrange(regulation, -mahalanobis,.by_group = TRUE) %>%
  dplyr::slice(1:10) %>%
  arrange(regulation, -neglogP,.by_group = TRUE) %>%
  dplyr::select(-mahalanobis) %>% print(n=200)

to_table$predictor <- factor(to_table$predictor)
levels(to_table$predictor) <- c("Leaf stage","-P","invfour")


DEGs <- effects %>%  
  group_by(predictor,regulation) %>% 
  filter(is_significant) %>%
  mutate(in_cis = gene %in% shared_introgression_gene_ids,
         in_trans = ! in_cis,
         in_Inv4m = gene %in% inv4m_gene_ids) %>%
  dplyr::select(predictor,regulation,gene, locus_symbol,description= desc_merged,logFC,neglogP, in_cis:in_Inv4m) %>% 
  arrange(regulation, -neglogP,.by_group = TRUE) %>%
  droplevels()

# Inv4m upregualtion or downregulation DEGs

ud <- with(
  DEGs %>% 
    filter(predictor== "GenotypeINV4" & in_cis) %>%
    mutate(FDR_regulation = case_when(
      logFC > 0 ~ "Up",
      logFC < 0 ~"Down",
      .default = NA)
      ) %>% 
    group_by(gene) %>%
    dplyr::slice(1),
      table(in_Inv4m,FDR_regulation)
    )

ud 
fisher.test(ud)


top_DEGs <- effects %>%  
  group_by(predictor,regulation) %>% 
  filter(is_significant, abs(logFC) >2) %>%
  mutate(in_cis = gene %in% shared_introgression_gene_ids,
         in_trans = ! in_cis,
         in_Inv4m = gene %in% inv4m_gene_ids)%>%
  dplyr::select(predictor,regulation,gene, locus_symbol,description= desc_merged,logFC,neglogP, in_cis:in_Inv4m) %>% 
  arrange(regulation, -neglogP,.by_group = TRUE) %>%
  droplevels()

top_DEGs %>% filter(predictor=="GenotypeINV4" , in_trans==TRUE) %>% pull(gene) %>%cat()

nrow(top_DEGs)
with(top_DEGs %>% filter(predictor!="leaf_tissue"),
     table(predictor,regulation)
)

# Inv4m upregualtion or downregulation topDEGs

ud <- with(
  top_DEGs %>% 
    filter(predictor== "GenotypeINV4" & in_cis) %>%
    mutate(FDR_regulation = case_when(
      logFC > 0 ~ "Up",
      logFC < 0 ~"Down",
      .default = NA)
    ) %>% 
    group_by(gene) %>%
    dplyr::slice(1),
  table(in_Inv4m,FDR_regulation)
)

ud

fisher.test(ud)

# There is  no over representation of upregulated or downregulated genes
# in the inversion relative to the flanking regions


write.csv( DEGs %>% filter(predictor=="Treatment-P") %>%
             mutate(),
           file= "~/Desktop/PSU_phosphorus_DEGs.csv",row.names=FALSE)

DEGs %>% filter(predictor=="GenotypeINV4")  print(n=200)


write.csv( DEGs %>% filter(predictor=="GenotypeINV4"),
           file= "~/Desktop/PSU_Inv4m_DEGs.csv",row.names=FALSE)


write.csv( add_malahanobis_outliers(effects),
           file= "~/Desktop/predictor_effects.csv",row.names=FALSE)


