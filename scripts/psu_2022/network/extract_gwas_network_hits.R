library(dplyr)
library(rtracklayer)
library(GenomicRanges)

v4_GFF <- "/Users/fvrodriguez/ref/zea/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3"
v4<- rtracklayer::import(v4_GFF)  %>%
  subset(type=="gene"  & seqnames %in% 1:10) + 5000

v5_GFF <- "~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
v5<- rtracklayer::import(v5_GFF )  %>%
  subset(type=="gene"  & seqnames %in% 1:10) + 5000

gene_symbol <- read.table("/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My\ Drive/repos/inv4mRNA/data/gene_symbol.tab",quote="",header=TRUE, sep ="\t", na.strings = "")
head(gene_symbol )

pannzer <-read.table("/Users/fvrodriguez/Desktop/PANNZER/PANNZER_DESC.tab",quote="",header=TRUE, sep ="\t", na.strings = "") %>%
  group_by(gene_model) %>%
  dplyr::slice(1) %>%
  dplyr::select(gene_model, desc)

gene_pannzer <- gene_symbol %>%
  left_join(pannzer)


# Get introgression coordinates in v5 ----

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
inv4m_end
# 188132113


inv4m_gene_ids <- genes %>%
  filter(
    seqnames==4,
    start >= inv4m_start,
    end <= inv4m_end
  ) %>% pull(ID)

length(inv4m_gene_ids)

# checked manually in the RNAseq geontype table
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

# sum(!rownames(voomR$E) %in% no_chr & !rownames(voomR$E) %in% shared_introgression_gene_ids)


# gene_pannzer 

v3_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v3_to_B73v5.tsv", sep = "\t") %>%
  dplyr::rename(v3="V1",v5="V2") %>%
  tidyr::separate_longer_delim(cols = v5,delim=",")

v4_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v4_to_B73v5.tsv", sep = "\t") %>%
  dplyr::rename(v4="V1",v5="V2") %>%
  tidyr::separate_longer_delim(cols = v5,delim=",")
grassius <- read.csv("~/Desktop/grassius_annotation_v5.csv")

# load effects


effects <- read.csv("~/Desktop/predictor_effects.csv") %>%
  dplyr::filter(predictor=="GenotypeINV4") %>%
  mutate(is_significant = adj.P.Val <0.05,
         is_topDEG = is_significant & abs(logFC) > 2,
         in_shared = gene %in% shared_introgression_gene_ids,
         in_Inv4m = gene %in% inv4m_gene_ids,
         in_flanking = gene %in% flanking_introgression_gene_ids,
         is_outside = !in_shared,
         is_TF = gene %in% grassius$gene.ID
  ) %>%
  group_by(gene) %>%
  arrange(adj.P.Val) %>%
  dplyr::slice(1) %>%  # This slice is necessary because I have a many to many relationship
  ungroup()

sum((effects$gene %>% sort %>% unique()) %in% flanking_introgression_gene_ids)

effects$in_inv4m <-NULL 



#load GWAS atlas
# wget ftp://download.big.ac.cn/gwas_atlas/gwas_association_result_for_maize.txt.gz

atlas <- read.table("~/Desktop/gwas_association_result_for_maize.txt", na.strings = c("-",""," "), head = FALSE, sep="\t", quote="")
atlas_cols <-
  c(
    "idx","chr","pos","assembly", "id1","assay_type","location","year","condition","germplasm","n","tissue","trait","TOID",
    "trait_type","genotype_tech","stat_model","native_pos","alt","num1","p.value","effect.size","native_gene", "gene_name",
    "note3","PUBMID","journal","title","first_author","note4","date"
  )

length(atlas_cols)

colnames(atlas) <- atlas_cols
atlas$p.value <- as.numeric(atlas$p.value)


# ouput columns:
# predictor regulation gene locus_symbol description logFC neglogP in_cis in_trans in_Inv4m

# Inv4m DEGs -----

# effects <- read.csv("~/Desktop/predictor_effects.csv")

effects %>%
  filter(gene=="Zm00001eb070790")

# Zinc Finger in inv4m from Vlad's TWAS
effects %>%
  filter(gene=="Zm00001eb193240") %>% tibble()
# p = 2.21e-5 but logFC = -0.833

Inv4m_DEGs <- effects%>%
  mutate(P = adj.P.Val) %>%
  mutate(neglogP=-log10(P)) %>%
  filter(predictor=="GenotypeINV4",
         P < 0.05) %>% 
  mutate(label_long = locus_symbol) %>%
  mutate(label_long = case_when(
    is.na(locus_name) ~ NA,
    .default = label_long)
  ) %>%
  mutate(label_short = coalesce(locus_symbol, gene)) %>%
  mutate(label_long = coalesce(label_long, desc)) %>%
  mutate(label_long = coalesce(label_long, gene)) %>%
  arrange(P) %>%
  mutate(is_TF = gene %in%  grassius$gene.ID ) %>%
  filter(!is.na(logFC)) %>% arrange(-neglogP) %>%
  dplyr::select(predictor,regulation,gene,symbol=locus_symbol,desc, label_long,label_short,logFC,neglogP,in_Inv4m,in_flanking,is_outside,is_TF) %>%
  tibble::tibble()

nrow(Inv4m_DEGs)


write.csv(Inv4m_DEGs, 
          file="~/Desktop/pheno_candidates/Inv4m_DEGs_ALL.csv")

write.csv(Inv4m_DEGs  %>%
            dplyr::filter(abs(logFC) >2 ), 
          file="~/Desktop/pheno_candidates/Inv4m_DEGs_TOP.csv")

Inv4m_DEGs  %>%
  dplyr::filter(abs(logFC) >2 ) %>% arrange(-abs(logFC)) %>% print(n=200)

# Flowering time candidates.-----

# Supplementary table 6 from
# Genomic insights into historical improvement of heterotic groups during modern hybrid maize breeding
# https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.16260?

li2022 <- read.table("~/Desktop/li2022_FT.list", header=FALSE) %>%
  dplyr::rename(v4= "V1") %>%
  left_join(v4_v5 %>% dplyr::select(v4,gene=v5)) %>%
  filter(!is.na(gene) & gene!="") %>%
  arrange(gene) %>% distinct()

li2022$gene %>% sort() %>% unique() %>% length()

# Molecular Parallelism Underlies Convergent Highland  Adaptation of Maize Landraces
liwang <- read.csv("/Users/fvrodriguez/Library/CloudStorage/GoogleDrive-frodrig4@ncsu.edu/My\ Drive/floweringTime_genes_liwang.csv",
                   na.strings = "") %>% 
  dplyr::select(v3=gene, MH:sum) %>%
  left_join(v3_v5 %>% dplyr::select(v3,gene=v5)) %>%
  filter(gene !="")

liwang$gene[liwang$gene==""] <-NA
liwang[liwang$gene=="Zm00001eb060540",]
liwang$gene %>% sort %>% unique %>% length

# These GWAS SNPs are in NAM v5 coordinates
# tibbs

plasticity<-read.csv("~/Desktop/Supplemental_Table_S3.csv", skip=1)
colnames(plasticity)
plasticity_long <- plasticity %>%
  dplyr::select(seqnames=Chr,ID=Marker, start=bp,end=bp,DTS_slope,DTS_intcp,DTA_slope,DTA_intcp) %>%
  tidyr::pivot_longer(cols=c("DTS_slope","DTS_intcp","DTA_slope","DTA_intcp"), values_to = "P", names_to = "trait") %>%
  filter(!is.na(P)) %>%
  mutate(neglogP = -log10(P)) %>%
  arrange(P)



FT <- with(plasticity_long,
           GRanges(seqnames=seqnames,
                   ID=ID,
                   ranges=IRanges(start=start, end=end),
                   P=P,
                   neglogP=neglogP,
                   strand = "+")
)

FT_df <- as.data.frame(FT)

genes <- as.data.frame(v5) 

genes$ID <- gsub("gene:","",genes$ID)

olap <-findOverlaps(FT,v5)

tibbs <- data.frame(gene=genes$ID[subjectHits(olap)], 
                    ID= FT_df$ID[queryHits(olap)],
                    P= FT_df$P[queryHits(olap)]) %>%
   group_by(gene) %>%
  arrange(P) %>%
  dplyr::slice(1) %>%
  arrange(P)
tibbs$gene[tibbs$gene==""] <- NA

tibbs_candidates <- data.frame(
  gene = genes$ID[subjectHits(olap)], 
  p.value = FT$P[queryHits(olap)]) %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::distinct() %>%
  group_by(gene) %>%
  arrange(p.value, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(p.value)

# GWAS ATLAS with flowering time PPTO:0000155-----

# the datas atlas is in v4 coordinates

atlas %>% filter(TOID == "PPTO:0000155") %>% pull(trait)
with(atlas %>% filter(TOID == "PPTO:0000155"),
     table(assembly))


FT <- with(atlas %>% filter(TOID == "PPTO:0000155"),
           GRanges(seqnames=chr,
                   ID=idx,
                   ranges=IRanges(start=pos, end=pos),
                   strand = "+")
)


genes <- as.data.frame(v4) 

genes$ID <- gsub("gene:","",genes$ID)

olap <-findOverlaps(FT,v4)

# Genic SNPs
olap@from %>% sort() %>% unique() %>% length()

# number of candidate genes
olap@to %>% sort() %>% unique() %>% length()
atlas[atlas$idx=="28622",]

atlas_candidates <- cbind(
  data.frame(gene=genes$ID[subjectHits(olap)]),
  
  atlas[FT$ID[queryHits(olap)], c("idx","chr","pos","TOID","p.value","PUBMID","first_author")]
) %>%
  group_by(gene,TOID) %>%
  arrange(gene,TOID,p.value) %>%
  dplyr::slice(1) %>%
  filter(TOID=="PPTO:0000155") %>%
  ungroup() %>%
  arrange(p.value)  %>% 
  left_join(v4_v5 %>% dplyr::select(v4,gene=v5), relationship = "many-to-many") %>%
  left_join(gene_pannzer, by=c(gene="gene_model"))
atlas_candidates$gene[atlas_candidates$gene==""] <- NA



FT_gwas <- rbind(
  atlas_candidates %>% dplyr::select(gene,p.value,),
  tibbs_candidates
) %>%
  filter(gene!="") %>%
  filter(!is.na(gene)) %>%
  filter(!is.na(p.value)) %>%
  group_by(gene) %>%
  dplyr::arrange(p.value) %>%
  dplyr::slice(1) %>%
  mutate(pheno="Flowering Time") %>%
  dplyr::arrange(p.value) 
FT_gwas

atlas_candidates$gene %>% sort() %>% unique() %>% length()

FT_candidates <- c(liwang$gene[!is.na(liwang$gene)],
                   li2022$gene[!is.na(li2022$gene)],
                   tibbs_candidates$gene[!is.na(tibbs_candidates$gene)],
                   atlas_candidates$gene[!is.na(atlas_candidates$gene)])%>%
                 sort() %>% unique()

length(FT_candidates)
"Zm00001eb060540" %in% FT_candidates 
"Zm00001eb060540" %in% liwang$gene
"Zm00001eb060540" %in% li2022$gene


FT_table  <- data.frame(gene=FT_candidates) %>%
  left_join(gene_pannzer, by=c(gene="gene_model")) %>%
  mutate(label_long = locus_symbol) %>%
  mutate(label_long = case_when(
    is.na(locus_name) ~ NA,
    .default = label_long)
  ) %>%
  mutate(label_short = coalesce(locus_symbol, gene)) %>%
  mutate(label_long = coalesce(label_long, desc)) %>%
  mutate(label_long = coalesce(label_long, gene)) %>%
  left_join( effects %>% filter(predictor=="GenotypeINV4") %>%
               dplyr::select( predictor,regulation,gene,logFC,P=adj.P.Val,in_Inv4m,in_flanking,is_outside, is_TF)
  ) %>% arrange(P) %>%
  distinct() %>%
  mutate(neglogP=-log10(P)) %>%
  dplyr::filter(!is.na(logFC)) %>% arrange(-neglogP) %>%
  dplyr::filter(10^-neglogP <0.05) %>% 
  dplyr::select(predictor,regulation,gene,symbol=locus_symbol,desc, label_long,label_short,logFC,neglogP,in_Inv4m,in_flanking,is_outside, is_TF) %>%
  tibble::tibble() %>%
  ungroup() %>%
  arrange(-abs(logFC)) %>% print (n=100)  
  
length(FT_candidates)

quartz()
FT_table %>%
  filter(!in_flanking) %>%arrange(-neglogP) %>%  print (n=100)  %>% pull(neglogP) %>% plot()
  
  
write.csv(FT_table, file="~/Desktop/FT_candidates_inv4m_perturbation_DEGs.csv")

write.csv(FT_table  %>%
            filter(abs(logFC)>2) %>% arrange(10^-neglogP), 
          file="~/Desktop/FT_candidates_inv4m_perturbation_TopDEGs.csv")

##  FT Candididate enrichment ----

FT_table  %>%
  filter(abs(logFC)>2) %>% arrange(10^-neglogP)

effects$is_FT_candidate <- effects$gene %in% FT_candidates

FT_topDEG <- with(effects %>% group_by(gene) %>% arrange(-neglogP) %>% dplyr::slice(1),
          table(effects$is_topDEG, is_FT_candidate)
)

FT_topDEG
fisher.test(FT_topDEG)

FT_DEG <- with(effects %>% group_by(gene) %>% arrange(-neglogP) %>% dplyr::slice(1),
     table(is_significant, is_FT_candidate)
)

FT_DEG

fisher.test(FT_DEG)

# A hypergeometic test is a one tailed  fisher exact test
fisher.test(FT_DEG, alternative = "greater")

# qnorm(p=1-0.001706)


# Plant Height candidates.-----


plasticity<-read.csv("~/Desktop/Supplemental_Table_S3.csv", skip=1)

plasticity_long <- plasticity %>%
  dplyr::select(seqnames=Chr,ID=Marker, start=bp,end=bp,PH_slope,PH_intcp) %>%
  tidyr::pivot_longer(cols=c("PH_slope","PH_intcp"), values_to = "P", names_to = "trait") %>%
  filter(!is.na(P)) %>%
  mutate(neglogP = -log10(P)) 


PH <- with(plasticity_long,
GRanges(seqnames=seqnames,
        ID=ID,
        ranges=IRanges(start=start, end=end),
        p.value=P,
        strand = "+")
)

 

olap <-findOverlaps(PH,v5)
tibbs_candidates <- data.frame(
  gene = genes$ID[subjectHits(olap)], 
  p.value = PH$p.value[queryHits(olap)])



tibbs <- tibbs_candidates$gene[!is.na(tibbs_candidates$gene)] %>% sort() %>% unique()
length(tibbs)


# head(atlas)
# tail(atlas)
# 
# table(atlas$alt) 
# 
# quartz()
# hist(-log10(atlas$p.value))

# with(atlas %>% filter(TOID == "PPTO:0000126"),
#      table(assembly))

PH <- with(atlas %>% filter(TOID  == "PPTO:0000126"),
           GRanges(seqnames=chr,
                   ID=idx,
                   ranges=IRanges(start=pos, end=pos),
                   strand = "+")
)



genes <- as.data.frame(v4) 

genes$ID <- gsub("gene:","",genes$ID)

olap <-findOverlaps(PH,v4)

atlas_candidates <- cbind(
data.frame(gene=genes$ID[subjectHits(olap)]),
atlas[PH$ID[queryHits(olap)], c("idx","chr","pos","TOID","p.value")]
) %>%
  group_by(gene,TOID) %>%
  arrange(gene,TOID,p.value) %>%
  dplyr::slice(1) %>%
  filter(TOID=="PPTO:0000126") %>%
  ungroup() %>%
  arrange(p.value)  %>% 
  rename(v4="gene") %>%
  left_join(v4_v5 %>% dplyr::select(v4,gene=v5), relationship = "many-to-many") %>%
  left_join(gene_pannzer, by=c(gene="gene_model"))


PH_gwas <- rbind(
  atlas_candidates %>% dplyr::select(gene,p.value),
  tibbs_candidates
) %>%
  filter(gene!="") %>%
  group_by(gene) %>%
  dplyr::arrange(p.value) %>%
  dplyr::slice(1) %>%
  mutate(pheno="Plant Height")

atlas_candidates$gene %>% sort() %>% unique() %>% length()



tail(atlas_candidates) 

li2022 <- data.frame(v4 = c("Zm00001d010396","Zm00001d011069","Zm00001d029657","Zm00001d029656","Zm00001d029658","Zm00001d029659","Zm00001d029758","Zm00001d029884","Zm00001d032747","Zm00001d032776","Zm00001d032788","Zm00001d032790","Zm00001d032791","Zm00001d033204","Zm00001d033206","Zm00001d033205","Zm00001d033223","Zm00001d033224","Zm00001d033225","Zm00001d033226","Zm00001d033228","Zm00001d033225","Zm00001d033227")) %>%
   inner_join(v4_v5) %>% pull(v5)

length(li2022)

PH_candidates <- c(tibbs,
                li2022,
                atlas_candidates$gene) %>% sort() %>% unique()

length(PH_candidates)

PH_genes <-data.frame( gene= PH_candidates) 

PH_genes

intersection <- effects %>% dplyr::select(predictor,regulation,gene,symbol=locus_symbol,logFC,P,neglogP, in_Inv4m,in_flanking,is_outside,is_TF) %>%
  filter(predictor=="GenotypeINV4") %>%
  inner_join(PH_genes ) %>%
    left_join(gene_pannzer%>% rename(gene="gene_model") , relationship = "many-to-many") %>%
  arrange(P) 



intersection$label_long <- intersection$locus_symbol
intersection$label_long[is.na(intersection$locus_name)] <- NA



PH_candidates <- intersection %>%
  mutate(label_short = coalesce(locus_symbol, gene)) %>%
  mutate(label_long = coalesce(label_long, desc)) %>%
  mutate(label_long = coalesce(label_long, gene)) %>%
  dplyr::filter(!is.na(logFC)) %>% arrange(-neglogP) %>%
  dplyr::select(
    predictor,regulation,gene,symbol=locus_symbol,desc, 
    label_long,label_short,
    logFC,neglogP,in_Inv4m,in_flanking,is_outside, is_TF) %>%
  ungroup() %>%
  arrange(-abs(logFC)) %>%
  tibble() %>% print(n=200)

nrow(PH_candidates)

effects$is_PH_candidate <- effects$gene %in% PH_genes$gene
  
## Plant Height candidate enrichment -----

PH_topDEG <- with(effects %>% group_by(gene) %>% arrange(-neglogP) %>% dplyr::slice(1),
                  table(effects$is_topDEG, is_PH_candidate)
)

PH_topDEG
fisher.test(PH_topDEG)


PH_DEG <- with(effects %>% group_by(gene) %>% arrange(-neglogP) %>% dplyr::slice(1),
               table(is_significant, is_PH_candidate)
)

PH_DEG
fisher.test(PH_DEG)


write.csv(PH_candidates,
          file="~/Desktop/PH_candidates_inv4m_perturbation_ALL.csv")

write.csv(PH_candidates  %>%
            dplyr::filter(10^-neglogP <0.05), 
          file="~/Desktop/PH_candidates_inv4m_perturbation_DEGs.csv")

gwas <-rbind(FT_gwas,PH_gwas)
  tail(FT_gwas)

out <- rbind(
FT_table %>%
  dplyr::filter(10^-neglogP <0.05 & abs(logFC)>2) %>% 
  mutate(pheno = "Flowering Time"),
PH_candidates %>%
  dplyr::filter(10^-neglogP <0.05 & abs(logFC)>2) %>% 
  mutate(pheno = "Plant Height")
) %>% 
  ungroup() %>%
  dplyr::select(pheno,  everything(), ) %>%
  dplyr::select(-predictor, -regulation, -starts_with("label"), -is_TF) %>% 
  left_join(gwas) %>%
  arrange(pheno,-neglogP) %>%
  print(n=200)

tibbs[tibbs %in% out$gene]


atlas_candidates %>%
  dplyr::filter(gene %in% out$gene) 

effects[effects$gene=="Zm00001eb011450",]

FT_table[FT_table$gene=="Zm00001eb011450",] 

write.csv(out, 
          file="~/Desktop/inv4m_TopDEGs_by_phenotype.csv")

