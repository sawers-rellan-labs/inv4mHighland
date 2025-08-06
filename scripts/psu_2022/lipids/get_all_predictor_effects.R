library(dplyr)
library(janitor)
library(ggplot2)
library(edgeR)
library(limma)
library(ashr)
library(mashr)


pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")

lipid_meta <- read.csv("data/lipid_metadata.csv")

# filter internal standards

internal_standards<- c(
  "CUDA",
  "Cholesterol",
  "CHOLESTEROL_D7_H20",
  "C17_CERAMIDE_H2O",
  "DG_12_0",
  "DG_18_1",
  "DG_18_1_2_0_0_0",
  "LPC_17_0",
  "LPE_17_1",
  "PC_25_0",
  "PE_34_0",
  "PG_34_0",
  "PG_17_0_17_0",
  "PG_17_0",
  "SM_35_1",
  "Sphingosine_17_1",
  "Sphingosine_D_17_1",
  "SPHINGOSINE_D_17_1",
  "TG_17_0",
  "TG_17_0_17_1_17_0_D5"
)


lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) 

to_exclude <-internal_standards[internal_standards %in% colnames(raw)]

raw  <- raw %>%  dplyr::select(-all_of(to_exclude)) %>%
  filter(!grepl("Methanol",sampleID)) %>%
  filter(!grepl("Pool",sampleID))
nrow(raw)
ncol(raw)


lipid_meta <- read.csv("data/lipid_metadata.csv") 

ms_order<- read.csv("data/PSU-PHO22 _ms_order.csv")
ms_order$tube <- gsub("L|R","S",ms_order$top_tag, perl =TRUE)


psu <- read.csv( "data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv") %>%
  rename(Genotype = "Who.What", row = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

metadata <- read.csv('data/PSU-PHO22_Metadata.csv')

metadata$tube <- gsub("L|R","S",metadata$tube, perl =TRUE)
colnames(metadata)


sampleInfo <- psu  %>% 
  dplyr::select(row,Rep,Plot_Row, Plot_Column,DTA,DTS) %>%
  dplyr::inner_join(metadata) %>% filter(!is.na(DTA)) %>%
  dplyr::right_join(ms_order) %>% distinct()
nrow(sampleInfo)
colnames(sampleInfo)


sampleInfo$leaf_node <- 2*sampleInfo$leaf_tissue -1 
sampleInfo$leaf_group <-NA
sampleInfo$leaf_group[sampleInfo$leaf_tissue <3] <- "apical"
sampleInfo$leaf_group[sampleInfo$leaf_tissue >2] <- "basal"





sampleInfo$rowid  <-  sampleInfo$row

raw$tube <- gsub(".*_|.raw","",raw$sampleID, perl =TRUE)
raw$tube <- gsub("L|R","S",raw$tube, perl =TRUE)
colnames(raw)


pheno <- sampleInfo %>%
  dplyr::select(tube,row,Rep,Plot:Plot_Row,Treatment:leaf_tissue,leaf_node,leaf_group) %>%
  inner_join(
    raw %>%
      dplyr::select(tube,DGDG_34_0:TG_58_5)) %>% 
  mutate(block=as.factor(Rep)) %>%
  droplevels() %>%
  tidyr::pivot_wider(names_from = Rep, values_from = Rep, names_prefix = "block_",
                     values_fn = function(x) 1, values_fill = 0)  %>%
  dplyr::select(tube:Plot_Row,block,block_6:block_12,everything()) %>%
  mutate(Treatment = factor(Treatment,levels=c("Low_P","High_P"))) %>% unique()

nrow(pheno)
colnames(pheno)




# MDS with raw data -----

colnames(pheno)

m <- pheno %>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() 

m[is.na(m)] <- 0

dim(m)
diag(m)
rownames(m) <- pheno$tube
# 
# n_factor <- mean(m)/rowSums(m)
# dim(m)
# length(n_factor)
# 
# norm_counts <- round(diag(n_factor) %*% m,0) %>% t()
# colnames(norm_counts) <- pheno$tube
# 
# norm <-cbind(
#   data.frame(tube=rownames(m)),
#   t(norm_counts)) # maximum_signal/total_lipid_sum persample
# 
# 
# counts <-  norm %>%
#   dplyr::select(DGDG_34_0:TG_58_5) %>%
#   as.matrix() %>% t()

counts <- t(m)
colnames(counts) <- pheno$tube


quartz()
hist(counts)


genes = data.frame(gene= rownames(counts))

sampleNames = pheno$tube
sampleNames %in% sampleInfo$tube
sampleInfo = sampleInfo[match(sampleNames,sampleInfo$tube),]

y = DGEList(counts = counts,samples = sampleInfo)
y$group = interaction(y$samples$Treatment,y$samples$Genotype)
d <- y$samples

# Modelling ---

design = model.matrix(~ injection_order + Plot_Column + Plot_Row + leaf_node+Treatment*Genotype,d)
voomR = voom(y,design=design,plot=T)

fit = lmFit(voomR)
ebfit = eBayes(fit)



quartz()
hist(ebfit$p.value)

ebfit$p.value

quartz()
hist(ebfit$p.value[,"injection_order"])

quartz()
hist(ebfit$p.value[,"leaf_node"])



quartz()
hist(ebfit$p.value[,"TreatmentLow_P"])

quartz()
hist(ebfit$p.value[,"GenotypeINV4"])

head(ebfit$coefficients)


ebfit$coefficients

nrow(ebfit$coefficients)
coef_cols <- ebfit$coefficients %>% colnames()
of_interest <-  c("leaf_node","TreatmentLow_P",
                  "GenotypeINV4",
                  "TreatmentLow_P:GenotypeINV4")
is_interesting <- which(coef_cols %in% of_interest )

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

effects <- results %>% dplyr::bind_rows() %>%
  mutate(predictor=factor(predictor, levels=of_interest)) %>%
  inner_join(lipid_meta, by=c(Response="col_name"))  %>%
  arrange(head_group,predictor,Response) 



# plot effects ----
leaf_nodes <- c("node1","node3","node5","node7")
names(leaf_nodes) <- c("leaf1","leaf2","leaf3","leaf4")

to_plot <- effects %>%
  filter(adj.P.Val< 0.05)   %>%
  mutate(sign=sign(logFC) %>% as.factor()) %>%
  mutate(head_group = factor(head_group, levels=c("phospholipid","glycolipid","neutral")) )%>%
  mutate(class = forcats::fct_reorder(class,-abs(logFC)))  %>%
  arrange(head_group,class,-abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  %>%
  mutate(class= forcats::fct_reorder(class,.r)) 
levels(to_plot$sign)

pd = position_dodge(0.4)

quartz()
to_plot  %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              color = head_group,
              fill  = head_group)) +
  ggtitle("Differential lipid analysis joint leaf analysis \n multiple hypothesis corrected FDR < 0.05 ") +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 2)+
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides( color = guide_legend(reverse = TRUE, title =""),
          shape = guide_legend(title =""),
          fill="none") +
  facet_grid(class~predictor, scales="free", space="free_y")  +
  scale_color_manual(values = c("red","blue","orange"))+
  scale_fill_manual(values = c("red","blue","orange"))+
  scale_shape_manual(values=c(17,25),
                     labels= c( "leaf 3","leaf 4")
  ) +
  scale_y_discrete(limits=rev) +
  theme_light(base_size =12) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        strip.text.y= element_text(color = "white",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


to_plot <- effects %>%
  #mutate(predictor=factor(predictor, levels= c("leaf","-P","Inv4m","-P:Inv4m"))) %>%
  filter(P.Value< 0.05) %>%  
  mutate(head_group = factor(head_group, levels=c("phospholipid","glycolipid","neutral")) )%>%
  mutate(class = forcats::fct_reorder(class,-abs(logFC)))  %>%
  arrange(head_group,class,-abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  %>%
  mutate(class= forcats::fct_reorder(class,.r)) 



pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              col = head_group)) +
  ggtitle("p-value < 0.05 (no multiple hypothesis adjustment)") +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free", ncol =3)  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("red","blue","orange"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))

