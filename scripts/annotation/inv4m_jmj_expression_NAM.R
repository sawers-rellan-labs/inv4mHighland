library(dplyr)
library(ggplot)
library(ggcorrplot)

setwd("~/Desktop/inv4m")


# Make phenotype table #########################################################
# Summer - winter data Hung et al 2012

# https://datacommons.cyverse.org/browse/iplant/home/shared/panzea/phenotypes/Hung_etal_2012_PNAS_data-120523.zip
gddtf <- read.csv("Dataset_S2.Photoperiod_paper_data_supplement_line_BLUPs.csv") %>%
  dplyr::filter(pop==27) %>%
  dplyr::select( Taxa = geno_code, gdd_anth_long, gdd_anth_short,gdd_anth_photo_resp
                 , gdd_silk_long, gdd_silk_short, gdd_silk_photo_resp)

pfeiffer <- read.csv("~/Desktop/Peiffer2014Genetics_blupPhenos20150325.csv")

# https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Kremling_Nature3RNASeq282_March2018
xpdir <- "TASSEL_fmt_expression_w_covariates"
tissue_file <- list.files(xpdir,pattern = "TASSEL.*Dups.txt") # This dataset is huge add to .gitignorre
tissue <-  gsub(".*and_Altnames_","",tissue_file)
tissue <-  gsub("_.*","",tissue)
tissue  <- tissue[1] 
names(tissue_file) <- tissue[1]

# slow step you don't wanna repeat it
# RNA <-read.table(file=file.path(xpdir,tissue_file["LMAD"]), skip =2, header = TRUE)
RNA$Taxa <- gsub("282set_","",RNA$Taxa)

pheno <-gddtf %>%
  inner_join(
    pfeiffer %>% filter(Panel=="ASSO") %>%
    dplyr::select(-(Panel:Z_Num), Taxa="Family_Inbred_Name")
    ) 
             
PHENO_RNA <- pheno %>%
  dplyr::inner_join(RNA) %>%
  dplyr::select(all_of(colnames(pheno)),everything())

PHENO_RNA[1:10,1:20]
# colnames(PHENO_RNA) <- gsub("gdd_","",colnames(PHENO_RNA))


TWAS <- read.csv("~/Desktop/inv4m/VLAD_TWAS_hits.csv", na.strings = "") %>%
  dplyr::mutate(label=coalesce(symbol,gene))

NK_candidates <- TWAS$gene[TWAS$NE_significant]
names(NK_candidates) <- TWAS$label[TWAS$NE_significant]
NK_candidates  <- NK_candidates[!is.na(NK_candidates)]
jmj <- c(
  jmj4 = "Zm00001eb191820",
  jmj6 = "Zm00001eb191790"
  # jmj21 = "Zm00001eb190750"
)

TWAS_candidates <- read.csv("TWAS_candidates_AGPv3.csv") %>%
  filter(v5 %in% c(jmj,NK_candidates))




of_interest <- TWAS_candidates$v3[TWAS_candidates$v3 %in% colnames(PHENO_RNA )] 

length(of_interest)

names(of_interest) <- TWAS_candidates$symbol[TWAS_candidates$v3 %in% colnames(PHENO_RNA )] 

PHENO_RNA[,c("DTS","DTA")]


xp_ft <- PHENO_RNA[,c("DTS","DTA",of_interest)] %>% data.matrix()
 
ncol(xp_ft)
colnames(xp_ft) <- c("DTS","DTA",names(of_interest ))





quartz()
hist(xp_ft[,-1])

round(cor(xp_ft, use = "pairwise.complete.obs"),2)

corr <- round(cor(xp_ft, use = "pairwise.complete.obs"),2)

order <- order(-corr[,"DTS"])

by_cor_pheno <- corr[order,order %>% rev()]


quartz()
ggcorrplot(by_cor_pheno )


TWAS_cor_plot <- ggcorrplot(by_cor_pheno) +
  theme(axis.text = element_text(face="italic"))





PHENO_RNA <- gddtf %>%
  dplyr::inner_join(RNA) %>%
  dplyr::select(Taxa, gdd_anth_long, gdd_anth_short, gdd_anth_photo_resp, gdd_silk_long, gdd_silk_short, gdd_silk_photo_resp,everything())

