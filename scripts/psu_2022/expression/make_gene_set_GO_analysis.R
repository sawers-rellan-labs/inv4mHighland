library(dplyr)
library(ggplot2)


effect_order <-  c("Leaf","-P","inv4m")


effects <- read.csv("~/Desktop/DEGs_FDR_for_table.csv") %>%
  mutate(neglogP = -log10(adj.P.Val)) %>%
  # mutate(is_significant = adj.P.Val<0.05) %>%
  # mutate(upregulated = case_when(
  #   (logFC > 0)    & effects$is_mh_outlier ~ TRUE,
  #   .default = FALSE)) %>%
  # mutate(downregulated = case_when(
  #   (logFC < 0)    & effects$is_mh_outlier ~ TRUE,
  #   .default = FALSE))  %>%
##########
  # mutate(is_significant = adj.P.Val<0.05) %>%
  # mutate(upregulated = case_when(
  #   (logFC > 0)    & is_significant ~ TRUE,
  #   .default = FALSE)) %>%
  # mutate(downregulated = case_when(
  #   (logFC < 0)    & is_significant ~ TRUE,
  #   .default = FALSE))  %>%

#######
  mutate(is_significant = adj.P.Val<0.05) %>%
  mutate(upregulated = case_when(
    (predictor !="Leaf" ) &   (logFC >= 2)    & is_significant ~ TRUE,
    (predictor =="Leaf" ) &   (logFC >= 0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(dowbregulated = case_when(
    (predictor !="Leaf" ) &   (logFC < 2)    & is_significant ~ TRUE,
    (predictor =="Leaf" ) &   (logFC < 0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(regulation = case_when(
    is_significant & upregulated ~ "Upregulated",
    is_significant & downregulated ~ "Downregulated",
    .default = "Unregulated"
  )) %>%
  mutate(is_DEG= is_significant & regulation != "Unregulated") %>%
  arrange(predictor,adj.P.Val)

effects %>%
  filter(regulation =="Upregulated" & predictor=="-P")

with(effects,
     table(predictor,regulation)
)


# colnames(effects)
# quartz()
# effects %>%
#   group_by(predictor) %>%
#   ggplot( aes(x=logFC, group= predictor))  + 
#   geom_histogram() +
#   facet_wrap(.~predictor, scales = "free_y")

nrow(effects)
# Get differential gene expression lists
# as: (up,down,) x (leaf,-P,inv4m,inv4m:-P) 


# All significant.

# 
# get_outlier_label <- function(data=data, outliers=NULL){
#   in_outliers <- data$predictor == outliers$predictor[1] & data$gene %in% outliers$gene
#   labels <-  which(in_outliers) 
#   names(labels) <- data$locus_symbol[in_outliers]
#   labels
# }
# 
# 
# outlier_list <- lapply(levels(effects$predictor), function(x){
#   get_mahalanobis_outliers(data=effects, condition=x)
# }
# ) 

# effects$label <-""
# 
# effects$label[outlier_labels] <- names(outlier_labels)
# 
# effects$label[effects$label==""] <- NA
# 
# outlier_list[[1]]$all %>% select(-mahalanobis)



library(rtracklayer)

B73v5 <- rtracklayer::import("~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3") %>%
 subset(type=="gene")

# B73v5 <- read.table("/System/Volumes/Data/ref/zea/Zm00001eb.1.fulldata.txt", header= TRUE,sep= "\t", quote="", fill = TRUE)
#nrow(B73v5)
colnames(B73v5)
# lm(data=effect_by_region %>% filter(regulation=="up"), logFC~region) %>% summary()

v4_v5 <-read.csv("data/Results-Zea_mays_Tools_IDMapper_.csv") %>% dplyr::select(-Releases)
v3_v4 <-read.csv("data/v3_v4_2.csv")
v4_v5_2 <-read.csv("data/v4_v5_2.csv")  %>% dplyr::select(-Releases)




levels(factor(effects$predictor))

effects %>%
  filter(predictor =="inv4m"& is_significant) %>% 
  pull(gene) %>%
  cat(., file = "~/Desktop/inv4_bg.gene_list")


# Send this over to the server
# library(AnnotationHub)
# hub <- AnnotationHub()
# org.Zm.eg.db <- hub[["AH114308"]]
# gomap <-buildGOmap(goa)


library(clusterProfiler )

source("~/Desktop/GOMAP_maize_B73_NAM5/go_functions.R", chdir=TRUE)

# interaction(levels(effects$predictor),c(levels(effects$regulation),"all")) %>% levels()
# 
# go_results <- list()
# 
# 
# 
# for(condition in levels(effects$predictor)){
#   r <- effects %>% 
#     filter(predictor==condition) 
#   
#   universe <- r$gene
#   
#   go_list_name <- paste(condition,"all",sep=".")
#   
#   regulated <- r %>% 
#     filter(is_DEG ) %>%
#     pull(gene)
#   print(go_list_name)
#   print(length(regulated))
#   
#   go_results[[go_list_name]] <- ego_analysis(regulated, universe)
#   go_list_name <- NULL
#   regulated <- NULL
#   
#   for(reg  in c("Upregulated","Downregulated")){
#     go_list_name <- paste(condition,reg,sep=".")
#     regulated <- r %>% filter( regulation==reg ) %>%
#       pull(gene)
#     print(go_list_name)
#     print(length(regulated))
#     
#     if(length(regulated)>1){
#       go_results[[go_list_name]] <- ego_analysis(regulated, universe)
#     }
#   }
# }


DEGs <- list(
  Leaf.Up = effects  %>% filter(predictor=="Leaf" & regulation=="Upregulated") %>% pull(gene),
  Leaf.Down = effects  %>% filter(predictor=="Leaf" & regulation=="Downregulated") %>% pull(gene),
  `-P.Up` = effects  %>% filter(predictor=="-P" & regulation=="Upregulated") %>% pull(gene),
  `-P.Down` = effects  %>% filter(predictor=="-P" & regulation=="Downregulated") %>% pull(gene),
  `inv4m.Up` = effects  %>% filter(predictor=="inv4m" & regulation=="Upregulated") %>% pull(gene),
  `inv4m.Down` = effects  %>% filter(predictor=="inv4m" & regulation=="Downregulated") %>% pull(gene)
)


DEGs <- list(
  Leaf.Up = effects  %>% filter(predictor=="Leaf" & is_significant & logFC > 0) %>% pull(gene),
  Leaf.Down = effects  %>% filter(predictor=="Leaf" & is_significant & logFC < 0) %>% pull(gene),
  `-P.Up` = effects  %>% filter(predictor=="-P" & is_significant & logFC > 0) %>% pull(gene),
  `-P.Down` = effects  %>% filter(predictor=="-P" & is_significant & logFC < 0) %>% pull(gene),
  `inv4m.Up` = effects  %>% filter(predictor=="inv4m" & is_significant & logFC > 0) %>% pull(gene),
  `inv4m.Down` = effects  %>% filter(predictor=="inv4m" & is_significant & logFC < 0) %>% pull(gene)
)


library(clusterProfiler)

# Objects needed for GO analysis
TERM2NAME <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2NAME.rds")
# TERM2GENE <-  readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024.rds") 
# TERM2GENE <-  readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024_full.rds") # Fattel_2024 annotation with all ancestor names
# TERM2GENE <-  readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE.rds")

colnames(TERM2GENE)
# "GO"   "gene"
# so this is. the panzer annotation

# TERM2GENE <- read.table("/Users/fvrodriguez/Desktop/B73_GO.out", sep ="\t", quote="", header=TRUE) %>%
#   mutate(GO= paste0("GO:",sprintf("%07d", goid)), gene= gsub("_.*","",qpid, perl=TRUE)) %>%
#   dplyr:::select(GO,gene) %>%
#   arrange(gene,GO) %>%
#   distinct()
# saveRDS(TERM2GENE,"~/Desktop/TERM2GENE_PANZER_2022.rds")
# TERM2GENE$gene %>% sort() %>% unique() %>% length()

#PANZER ANNPOTATIONS
TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/TERM2GENE_PANZER_2022.rds")
TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/TERM2GENE_PANZER_2022_full.rds")


#FASSO ANNOTATIONS
id_map<- read.csv("/Users/fvrodriguez/Desktop/maize_FASSO_annotations/Protein_Transcript_Gene.csv",header=FALSE)
colnames(id_map)[c(2,3)] <- c("qpid", "gene")
head(id_map)
uniprot_goa <- read.table("/Users/fvrodriguez/Desktop/maize_FASSO_annotations/maize_GO.out", header=TRUE, fill = TRUE,sep="\t", row.names = NULL)
uniprot_goa$GO<- sprintf("GO:%07d",uniprot_goa$goid)

TERM2GENE <- id_map %>% inner_join(uniprot_goa, relationship = "many-to-many") %>% select(GO,gene)

#saveRDS(TERM2GENE,"~/Desktop/TERM2GENE_FASSO_PANZER_2022.rds")
#gomap <-buildGOmap(TERM2GENE)

#saveRDS(gomap, file="~/Desktop/TERM2GENE_FASSO_PANZER_2022_full.rds")

TERM2GENE <-gomap

library(clusterProfiler)

nrow(TERM2GENE)
head(TERM2GENE)
#TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024.rds")

# Gather all Terms for the GO_search function
TERM2NAME_ALL <- do.call("rbind", TERM2NAME)

with(effects,
     table(predictor)
)

with(effects,
     table(predictor,is_DEG)
)

TERM2GENE_MF <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$MF$go_id)
TERM2GENE_BP <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$BP$go_id)
TERM2GENE_CC <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$CC$go_id)

universe <- effects  %>% filter(predictor=="Leaf") %>% pull(gene)

length(universe)

ego_BP <- function(x){ 
  n_genes <- length(x)
  print(paste("calculating enrichment for:",n_genes,"genes"))
  enricher(
    gene=x,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe=universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_BP,
    TERM2NAME = TERM2NAME$BP
  )}


ego_MF <- function(x){ 
  n_genes <- length(x)
  print(paste("calculating enrichment for:",n_genes,"genes"))
  enricher(
    gene=x,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe=universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_MF,
    TERM2NAME = TERM2NAME$MF
  )}

ego_CC <- function(x){ 
  n_genes <- length(x)
  print(paste("calculating enrichment for:",n_genes,"genes"))
  enricher(
    gene=x,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe=universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_CC,
    TERM2NAME = TERM2NAME$CC
  )}

compGO_CC <- compareCluster(geneCluster   = DEGs,
                            fun           = "ego_CC")
quartz()
dotplot(compGO_CC)



compGO_MF <- compareCluster(geneCluster   = DEGs,
                         fun           = "ego_MF")

quartz()
dotplot(compGO_MF)

compGO_BP <- compareCluster(geneCluster   = DEGs,
                         fun           = "ego_BP")
#write.csv(compGO_BP,file="~/Desktop/inv4mP_FDR_ALL_goa_PANZER.csv")


quartz()
dotplot(compGO_BP)

compGO_BP@compareClusterResult$Cluster %>% sort %>% unique()

compGO_BP@compareClusterResult %>% select(-geneID) %>% tidyr::tibble() %>%
  arrange(Cluster,pvalue) %>% print(n=250)

# curated_goa <-read.csv(file="~/Desktop/inv4mP_FDR_ALL_goa_curated.csv")
curated_goa <-NULL
curated_goa <-read.csv(file="~/Desktop/inv4mP_FDR_ALL_goa_PANZER_curated.csv")

slim <- curated_goa %>%
  dplyr::filter(to_slim2==TRUE) %>%
  #dplyr::select(ontology,ID,Description) %>%
  dplyr::select(ID,Description)  %>%
  dplyr::arrange(ID) %>%
  dplyr::distinct() 

slim %>%dplyr::filter(ID=="GO:0019684")  

slim_comp <- compGO_BP
slim_comp@compareClusterResult <- slim_comp@compareClusterResult %>%
  filter(ID %in% slim$ID)


slim_comp@compareClusterResult %>%
  dplyr::filter(ID=="GO:0019684")  


cluster_order <- c("Leaf.Down","Leaf.Up","-P.Down","-P.Up","inv4m.Down","inv4m.Up")

# slect the minimum (p_value per term)
sorted_by_p <- slim_comp@compareClusterResult  %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  dplyr::slice(1:8) %>%
  mutate(label.order = 1:n()) %>%
  mutate(Description= forcats::fct_reorder(Description,label.order))

top_terms <- levels(sorted_by_p$Description)

# to_bold <- c("GO:0015706","GO:0016036","GO:0010073","GO:0006032","GO:0009813","GO:0006479","GO:0019374")

to_plot <-compGO_BP@compareClusterResult  %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  filter(Description %in% top_terms) %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  mutate(Description = factor(Description, levels = top_terms)) %>% 
  mutate(label_order = 1:n()) %>%
  droplevels() %>%
  mutate(Description= forcats::fct_reorder(Description,label_order))

to_plot$Cluster %>% levels()
colnames(effects)
colnames(TERM2GENE)

gene_symbol <- effects %>% 
  mutate(Cluster = case_when(
    regulation == "Upregulated" ~ paste(predictor,"Up", sep="."),
    regulation == "Downregulated" ~ paste(predictor,"Down",sep=".")
  )) %>%
  dplyr::select(Cluster,gene, SYMBOL="locus_symbol",
                  locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
  filter(!is.na(SYMBOL) & !is.na(locus_name) & !is.na(Cluster)) %>%
  mutate(neglogPG=-log10(P))  %>%
  inner_join(TERM2GENE, relationship =  "many-to-many") %>%
  rename(geneID="gene", ID="GO")

gene_symbol %>%
  filter(ID=="GO:0015706")

TERM2GENE %>%
  filter(gene=="Zm00001eb126380")

  

term2symbol <- to_plot %>%   
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  inner_join(gene_symbol) %>%
  dplyr::select(Cluster,geneID,SYMBOL,ID, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG) %>%
  distinct() %>%
  group_by(Cluster,ID) %>%
  arrange(Cluster,ID,-mahalanobis) %>%
  slice(1)

term2symbol %>% print(n=150)



# term2symbol %>%
#   filter(SYMBOL=="spx4")
# term2symbol %>%
#   filter(ID=="GO:0019374")
# 
# term2symbol %>%
#   filter(ID=="GO:0019374") %>% pull(geneID)
# 
# galacto_genes <- term2symbol %>%
#     filter(ID=="GO:0019374") %>% pull(geneID)



# labels <- c(
#   "Leaf.Up" = expression(""%up% Leaf), 
#   "Leaf.Down" = expression("" %down% Leaf),
#   "-P.Up" = expression("" %up% -P),
#   "inv4m.Up" = expression("" %up% italic(Inv4m))
# )
# 

labels <- c(
  "Leaf.Up" = "Leaf Up", 
  "Leaf.Down" = "Leaf Down",
  "-P.Up" = "-P\nUp",
  # "inv4m.Up" = expression(atop("Up",italic(Inv4m),))
  "inv4m.Up" = "Inv4m\nUp"
)

labels <- c(
  "Leaf.Up" = "Leaf Up", 
  "Leaf.Down" = "Leaf Down",
  "-P.Up" = "-P Up",
  "-P.Down" = "-P Down",
  # "inv4m.Up" = expression(atop("Up",italic(Inv4m),))
  "inv4m.Up" = "Inv4m Up"
)


bg <- to_plot %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, fill= neglogP)) +
  ggtitle("Biological Process Annotation of DEGs ") +
  scale_fill_continuous(low="dodgerblue", high="tomato", name = "-log10(FDR)",
                        guide=guide_colorbar())+ 
                        # limits=c(1.3, 8.3), 
                        # breaks=2*(1:4),
                        # labels=2*(1:4)) +
  scale_y_discrete(labels = scales::wrap_format(50),
                   limits = rev(levels(to_plot$Description))) +
  scale_x_discrete(labels=labels) +
  geom_point(shape=21)  +
  scale_size(range = c(2, 8))+
  ylab(NULL) +  xlab(NULL) +
  DOSE::theme_dose() +
  theme(axis.text.x = element_text(vjust=0.5))
        # xis.text.y =ggtext::element_markdown())


term2symbol$neglogP <-  -log10(term2symbol$p.adjust)
  
quartz(height=6, width = 12)
bg + # data point size
geom_text(
    data=term2symbol,  position=position_nudge(0.1), hjust=0,vjust=0.5,
    fontface="italic",
    mapping =aes(x=Cluster,y =Description, label =SYMBOL),
    inherit.aes = FALSE
    )

quartz()
dotplot(slim_comp, color = "p.adjust", 
        showCategory =9,
        size="count", label_format=50) +
  theme(axis.text.x = element_text(angle=90,hjust = 0.5))



# lapply(go_results, function(list_ego_results){ sum(x@result$p.adjust < 0.05)}
#        )


ego_df <- lapply(names(go_results),FUN=function(set){
  ego_list <- go_results[[set]]
  
  lapply(names(ego_list), FUN=function(ontology){
    df_out <- as.data.frame(ego_list[[ontology]])
    if(nrow(df_out)>0){
      prefix <- paste0(ontology,".")
      colnames(df_out) <- gsub(prefix ,"",colnames(df_out))
      df_out$ontology <- ontology
      df_out$set <- set
      df_out %>%
        dplyr::filter(p.adjust <0.05) %>%
        dplyr::select(set,ontology,ID,Description,p.adjust,FoldEnrich,Count)
    }
  }) %>% dplyr::bind_rows() 
  
}) %>% dplyr::bind_rows()


ego_df$ontology <-factor(ego_df$ontology)
levels(ego_df$ontology) <- c("BP","CC","MF")
ego_df$ontology <-factor(ego_df$ontology)
ego_df <- ego_df %>% tidyr::separate_wider_delim(set, ".", names = c("predictor","reg_set"), cols_remove = FALSE)

predictor_order <-c("Leaf","-P","inv4m")
ego_df$predictor <- factor(ego_df$predictor,levels= predictor_order)

reg_set_order<- c("Upregulated","Downregulated","all")
ego_df$reg_set <- factor(ego_df$reg_set,levels= reg_set_order)


with(ego_df,
     table(set,ontology)
)


write.csv(ego_df, file="~/Desktop/inv4mP_goa.csv")

curation <- read.csv(file="~/Desktop/inv4mP_goa_curated.csv") %>%
            dplyr::select(ID,to_slim1,to_slim2) %>%
            dplyr::distinct()
nrow(curation)            

for_curation <- ego_df %>% dplyr::select(-(predictor:set)) %>%
  group_by(ID) %>%
  arrange(ID,p.adjust) %>%
  slice(1) %>%
  dplyr::distinct() %>%
  left_join(curation)

nrow(for_curation)

for_curation$to_slim1[is.na(for_curation$to_slim1)] <- FALSE
for_curation$to_slim2[is.na(for_curation$to_slim2)] <- FALSE
write.csv(for_curation, file="~/Desktop/inv4mP_goa_forcuration.csv")


enrichment <- go_results$`-P.Upregulated`$ego_BP

simple_idx <- which(enrichment@result$ID %in% slim$ID)
names(simple_idx) <- enrichment@result$ID[simple_idx]
simple_enrichment <- subset_enrichResult(enrichment,simple_idx)
enrichment@result$Description[simple_idx]

names(simple_idx) <- enrichment@result$Description[simple_idx]
enrichment@result$Description[simple_idx]

write.csv(enrichment@result %>%select(-geneID), file="~/Desktop/-P.Upregulated_goa.csv")



length(simple_idx)


quartz()
dotplot(simple_enrichment ,
        #orderBy ="p.adjust",
        showCategory = names(simple_idx)[1:10],
        color= "p.adjust",
        size="FoldEnrich",
        decreasing =FALSE,
        x='p.adjust',
        title="-P upregulated: Biological Process") +
  scale_x_log10()

title="-P upregulated: Biological Process"

quartz()
simple_enrichment@result %>%
  dplyr::mutate(Description=factor(Description)) %>%
  dplyr::mutate(Description=forcats::fct_reorder(Description,-p.adjust)) %>%
  slice(1:10) %>%
  ggplot( aes(x=p.adjust, y=Description, size=Count, color=FoldEnrich)) +
  geom_point() +
  scale_x_log10() +
  scale_y_discrete(labels = scales::wrap_format(30))+
  scale_color_continuous(low="blue", high="red", name = "Fold\nEnrinchment",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) + ggtitle(title) + theme_dose() 



enrichment <- go_results$`Leaf.Downregulated`$ego_BP

simple_idx <- which(enrichment@result$ID %in% slim$ID)
names(simple_idx) <- enrichment@result$ID[simple_idx]
simple_enrichment <- subset_enrichResult(enrichment,simple_idx)
enrichment@result$Description[simple_idx]

names(simple_idx) <- enrichment@result$Description[simple_idx]
enrichment@result$Description[simple_idx]

length(simple_idx)

quartz()
dotplot(simple_enrichment ,
        #orderBy ="p.adjust",
        showCategory = names(simple_idx)[1:10],
        size="FoldEnrich",
        decreasing =FALSE,
        x='p.adjust',
        title="Leaf Downregulated: Biological Process") +
  scale_x_log10()



title="Leaf Downregulated: Biological Process"



quartz()
simple_enrichment@result %>%
  dplyr::mutate(Description=factor(Description)) %>%
  dplyr::mutate(Description=forcats::fct_reorder(Description,-p.adjust)) %>%
  slice(1:10) %>%
  ggplot( aes(x=p.adjust, y=Description, size=Count, color=FoldEnrich)) +
  geom_point() +
  scale_x_continuous(
    trans  = scales::compose_trans("log10", "reverse")
  ) +
  scale_y_discrete(labels = scales::wrap_format(30))+
  scale_color_continuous(low="dodgerblue", high="tomato", name = "Fold\nEnrinchment",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) + ggtitle(title) + theme_dose() 

quartz()
dotplot(simple_enrichment , x="FoldEnrich",showCategory=10, title="-P upregulated: Biological Process") 

quartz()
dotplot(go_results$`Leaf.Downregulated`$ego_BP, x="FoldEnrich",showCategory=10, title="Leaf Downregulated: Biological Process")

quartz()
dotplot(go_results$`Leaf.Upregulated`$ego_BP, x="FoldEnrich",showCategory=10, title="Leaf Downregulated: Biological Process")


ego_df %>%
  group_by(set,ontology,ID) %>%
  arrange(p.adjust) %>%
  filter(set=="inv4m.Upregulated")


library(AnnotationHub)
library(GOSemSim)
library(igraph)

hub <- AnnotationHub()
class(hub)

org.Zm.eg.db <- hub[["AH114308"]]

ZmGO <- godata(org.Zm.eg.db, ont="BP")


go_ids <- unique(ego_df$ID[ego_df$ontology=="BP"] %>% sort())

mat <- GOSemSim::termSim(go_ids,go_ids,ZmGO,method ="Rel")
mat[is.na(mat)] <- 0

quartz()
df = simplifyGO(mat)

go_set <- ego_df[ego_df$ontology=="BP",]
go_set$ID
go_set %>%
  dplyr::inner_join(df %>% dplyr::rename(ID="id")) %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,cluster) %>%
  arrange(predictor,reg_set,ontology,cluster,-FoldEnrich) %>%
  slice(1) %>%   ungroup() %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  dplyr::select(-predictor,-reg_set,-Count)  %>%
  print(n=200)


# semantic distance



semsim <-  mgoSim(go_ids, go_ids, semData=ZmGO, combine=NULL,measure = "Rel")
hist(semsim)


library(igraph)

adjacency <- semsim 
adjacency[adjacency >= 0.5 ] <-1
adjacency[adjacency< 0.5]  <-0
adjacency[is.na(adjacency)]  <- 0

G <- graph_from_adjacency_matrix(adjacency,mode = "undirected")

cc<- components(G)

go_set$semcluster <- cc$membership[go_set$ID]

# library(rrvgo)

colnames(go_set)




go_set <- go_set %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology) %>%
  arrange(predictor,reg_set,ontology, p.adjust) 

runs <- rle(unname(go_set$semcluster))

go_set$run <- rep(1:length(runs$lengths), runs$lengths)


leaf_down <- go_set %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,run) %>%
  arrange(predictor,reg_set,ontology,run, -FoldEnrich) %>%
  slice(1) %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  ungroup() %>%
  dplyr::select(-semcluster,-predictor,-reg_set,-Count) %>% 
  dplyr::filter(set=="Leaf.Downregulated") %>% print(n=200)


P_up <- go_set %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,run) %>%
  arrange(predictor,reg_set,ontology,run, -FoldEnrich) %>%
  slice(1) %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  ungroup() %>%
  dplyr::select(-semcluster,-predictor,-reg_set,-Count)  %>%
  dplyr::filter(set=="-P.Upregulated") %>% print(n=200)


go_set %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,run) %>%
  arrange(predictor,reg_set,ontology,run,-FoldEnrich) %>%
  slice(1) %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  ungroup() %>%
  dplyr::select(-semcluster,-predictor,-reg_set,-Count) %>% print(n=200) %>%
  dplyr::filter(set=="inv4m.Upregulated")


go_set %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,run) %>%
  arrange(predictor,reg_set,ontology,run,-FoldEnrich) %>%
  slice(1) %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  ungroup() %>%
  dplyr::select(-semcluster,-predictor,-reg_set,-Count)  %>%
  dplyr::filter(set=="inv4m.Upregulated") %>% print(n=200)






ego_df %>%
  dplyr::filter(set=="-P.Upregulated") %>%
  ungroup() %>%
  group_by(predictor,reg_set,ontology,run) %>%
  arrange(predictor,reg_set,ontology,run,-p.adjust) %>%
  slice(1) %>%
  arrange(predictor,reg_set,ontology,p.adjust) %>%
  ungroup() %>%
  dplyr::select(-semcluster,-predictor,-reg_set,-Count) %>% print(n=100)






r_go_set %>% 
  dplyr::filter(predictor=="-P") %>%
  dplyr::select(-predictor,-reg_set)  %>%
  dplyr::group_by(set,run,semcluster) %>%
  #  dplyr::arrange(set,run,semcluster,p.adjust) %>%
  dplyr::arrange(run,semcluster,-FoldEnrich) %>%
  dplyr::slice(1) %>%
  dplyr::arrange(set,p.adjust) %>% print(n=100)




r_go_set <- ego_df %>% dplyr::filter(set=="-P.all")

simMatrix <- calculateSimMatrix(r_go_set$ID,
                                org.Zm.eg.db,
                                semdata = ZmGO,
                                ont="BP",
                                method="Wang")

#scores <- setNames(-log10(r_go_set$p.adjust), r_go_set$ID)
scores <- setNames(r_go_set$FoldEnrich, r_go_set$ID)
reduced_terms <- reduceSimMatrix(simMatrix,
                                 scores = scores,
                                 threshold= 0.5,
                                 orgdb=org.Zm.eg.db)

runs <- rle(unname(r_go_set$semcluster))
r_go_set$run <- rep(1:length(runs$lengths), runs$lengths)

r_go_set %>% 
  dplyr::filter(set=="-P.all") %>%
  dplyr::filter(ID %in% reduced_terms$go) %>%
  dplyr::select(-predictor,-reg_set)  %>%
  dplyr::group_by(set,semcluster) %>%
  #  dplyr::arrange(set,run,semcluster,p.adjust) %>%
  dplyr::arrange(set,run,semcluster,-FoldEnrich) %>%
  dplyr::slice(1) %>%
  dplyr::arrange(set,p.adjust)



r_go_set <- ego_df %>% dplyr::filter(set=="Leaf.down")
r_go_set %>% dplyr::filter(set=="Leaf.down")

simMatrix <- calculateSimMatrix(r_go_set$ID,
                                org.Zm.eg.db,
                                semdata = ZmGO,
                                ont="BP",
                                method="Wang")

# scores <- setNames(-log10(r_go_set$p.adjust), r_go_set$ID)
scores <- setNames(r_go_set$FoldEnrich, r_go_set$ID)
reduced_terms <- reduceSimMatrix(simMatrix,
                                 scores = scores,
                                 threshold= 0.5,
                                 orgdb=org.Zm.eg.db)

runs <- rle(unname(r_go_set$semcluster))
r_go_set$run <- rep(1:length(runs$lengths), runs$lengths)


simMatrix <- calculateSimMatrix(r_go_set$ID,
                                org.Zm.eg.db,
                                semdata = ZmGO,
                                ont="BP",
                                method="Wang")
levels(factor(ego_df$reg_set))
r_go_set <- ego_df %>% dplyr::filter(set=="leaf.Downregulated")


# scores <- setNames(-log10(r_go_set$p.adjust), r_go_set$ID)
scores <- setNames(r_go_set$FoldEnrich, r_go_set$ID)
# reduced_terms <- reduceSimMatrix(simMatrix,
#                                  scores = scores,
#                                  threshold= 0.5,
#                                  orgdb=org.Zm.eg.db)

runs <- rle(unname(r_go_set$semcluster))
r_go_set$run <- rep(1:length(runs$lengths), runs$lengths)


r_go_set %>% 
  dplyr::filter(ontology=="BP") %>%
  #dplyr::filter(ID %in% reduced_terms$go) %>%
  dplyr::select(-predictor,-reg_set,-Count)  %>%
  dplyr::group_by(set,semcluster) %>%
  dplyr::arrange(set,run,semcluster,-FoldEnrich) %>%
  dplyr::slice(1) %>%
  dplyr::arrange(set,p.adjust) %>% print(n=100)


data.frame()
with(ego_df,
     table(set,ontology)
)





nrow(ego_df)
ego_df$reg_set
BP_annotation_top20 <- ego_df %>%
  ungroup() %>%
  group_by(ontology,predictor,reg_set,semcluster) %>%
  arrange(-FoldEnrich, p.adjust) %>%
  dplyr::slice(1)  %>%
  ungroup() %>%
  group_by(ontology,predictor,reg_set) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice(1:20) %>%
  ungroup() %>%
  filter(reg_set!="Leaf.Unregulated" & ontology =="BP") %>%
  arrange(ontology,predictor,reg_set,p.adjust) %>% 
  dplyr::select(-ID,-predictor,-reg_set,-semcluster) %>%  print(n=100)

# go_slim <- read.table(file="~/Desktop/plant_go_slim.tab", sep ="\t", header=TRUE)

ego_df %>%
  #inner_join(go_slim %>% dplyr::select(GO_ID), by=c(ID="GO_ID"))
  #dplyr::select(-ID,-semcluster) %>% 
  group_by(ontology,predictor) %>%
  arrange(ontology,predictor,p.adjust) %>%
  dplyr::filter(Description =="aging")

r  <-  effects %>%
  dplyr::filter(predictor=="Leaf")

regulated <- r %>% 
  dplyr::filter(is_significant & is_mh_outlier) 

aging <- "GO:0007568"
regulated[regulated$gene %in% go_search(method="GO2gene",aging)[,2], ] %>%
  arrange(-logFC)


leaf_senescence <- "GO:0010150"
regulated[regulated$gene %in% go_search(method="GO2gene",leaf_senescence)[,2], ] %>%
  arrange(-logFC)



regulated[regulated$gene %in% go_search(method="GO2gene",leaf_senescence_regulation)[,2], ]

ego_df %>%dplyr::filter(Description =="lipid transport")


write.csv(BP_annotation_top20,file=("~/Desktop/BP_annotation_top20.csv"))

list_ego_results <- ego_analysis(regulated,universe)

lapply(list_ego_results, function(x) sum(x@result$p.adjust < 0.05))

quartz()
dotplot(list_ego_results$ego_BP, x="FoldEnrich",showCategory=10, title="BP")

quartz()
dotplot(list_ego_results$ego_MF,x="FoldEnrich", showCategory=10, title="MF")

quartz()
dotplot(list_ego_results$ego_CC,x="FoldEnrich",showCategory=10, title="CC")



