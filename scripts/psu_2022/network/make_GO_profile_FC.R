library(dplyr)
library(ggplot2)
library(clusterProfiler)

effect_order <-  c("Leaf","-P","inv4m")

# Imposing fold change limits

effects <- read.csv("~/Desktop/DEGs_FDR_for_table.csv") %>%
  mutate(neglogP = -log10(adj.P.Val)) %>%
  mutate(is_significant = adj.P.Val<0.05) %>%
  mutate(is_up = case_when(
    (predictor !="Leaf" ) & (logFC >= 2)    & is_significant ~ TRUE,
    (predictor =="Leaf" ) & (logFC >= 0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(is_down = case_when(
    (predictor !="Leaf" ) & (logFC <= -2)    & is_significant ~ TRUE,
    (predictor =="Leaf" ) & (logFC <= -0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(regulation = case_when(
    is_significant & is_up ~ "Upregulated",
    is_significant & is_down ~ "Downregulated",
    .default = "Unregulated"
  )) %>%
  mutate(is_DEG= is_significant & regulation != "Unregulated") %>%
  arrange(predictor,adj.P.Val)

with(effects,
table(predictor,regulation)
)

universe <- union( effects$gene[effects$predictor=="Leaf"],effects$gene[effects$predictor=="Leaf"])
length(universe)

# Load GO annotation
TERM2NAME <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2NAME.rds")

TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024_full.rds")

ego <- function(x, ontology ="BP"){ 
  n_genes <- length(x)
  print(paste("calculating enrichment for:",n_genes,"genes"))
  enricher(
    gene=x,
    pvalueCutoff = 0.05,
    pAdjustMethod = "none",
    universe=universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE= TERM2GENE[TERM2GENE$GO %in% TERM2NAME[[ontology]]$go_id,],
    TERM2NAME = TERM2NAME[[ontology]]
  )}


DEGs <- with( effects,{
  list(
    Leaf.Down = gene[predictor=="Leaf" & regulation=="Downregulated"],
    Leaf.Up = gene[predictor=="Leaf" & regulation=="Upregulated"],
   # Leaf.all = gene[predictor=="Leaf" & regulation!="Unregulated"],
    `-P.Down` = gene[predictor=="-P" & regulation=="Downregulated"],
    `-P.Up` = gene[predictor=="-P" & regulation=="Upregulated"],
  #  `-P.all` = gene[predictor=="-P" & regulation!="Unregulated"],
    inv4m.Down = gene[predictor=="inv4m" & regulation=="Downregulated"],
    inv4m.Up = gene[predictor=="inv4m" & regulation=="Upregulated"]
  #  inv4m.all = gene[predictor=="inv4m" & regulation!="Unregulated"]
  )}
)

n.genes <- lapply(DEGs, length) %>% unlist
n.genes

compGO_BP <- compareCluster(geneCluster   = DEGs,
                            fun           = ego)

compGO_BP@compareClusterResult$Cluster %>% table()

compGO_BP@compareClusterResult 

cluster_order <- names(DEGs)


levels(compGO_BP@compareClusterResult$Cluster) <- cluster_order

quartz()
dotplot(compGO_BP,size="Count",  showCategory = 10) 

curated_goa <-NULL
curated_goa <-read.csv(file="~/Desktop/inv4mP_FC_goa_Fattel_full_curated.csv")

slim <- curated_goa %>%
  dplyr::filter(to_slim2==TRUE) %>%
  #dplyr::select(ontology,ID,Description) %>%
  dplyr::select(ID,Description)  %>%
  dplyr::arrange(ID) %>%
  dplyr::distinct() 

slim <- read.csv("~/slim_to_plot.csv")
effects[effects$gene %in% c("Zm00001eb188570","Zm00001eb188030","Zm00001eb187430"),]
slim_comp <- compGO_BP
slim_comp@compareClusterResult <- slim_comp@compareClusterResult %>%
  filter(ID %in% slim$ID)
slim_comp@compareClusterResult$Cluster %>% table()

slim_comp@compareClusterResult %>% 
  filter( Cluster == "Leaf.Down" & ID =="GO:0010582") %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many")



v3_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v3_to_B73v5.tsv", sep = "\t", header= FALSE) %>%
  dplyr::rename(v3="V1",v5="V2")%>%
  tidyr::separate_longer_delim(cols = v5,delim=",")


sags <- v3_v5 %>%
  inner_join(read.csv("~/Desktop/SAG_orthologs.csv"), by=c(v3="GENE.ID..Maize."))

leaf_nat <- read.csv("~/Desktop/natural_senescence.csv")

quartz()
leaf_nat %>%
  ggplot(aes(x=log_ESL_ML,y=log_LSL_ML)) +
  geom_point()


slim_comp@compareClusterResult %>% 
   filter( Cluster == "Leaf.Up" & ID =="GO:0007568") %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many") 

# no overlap between aging and SAGs from 
aging_genes_v3 <- slim_comp@compareClusterResult %>% 
  filter( #Cluster == "Leaf.Down" &
            ID =="GO:0007568") %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID")  %>%
  inner_join(v3_v5, by = c(gene="v5")) %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many") %>% pull(v3) %>% sort() %>% unique()

aging_genes_v3



quartz()
leaf_nat %>%
  dplyr::mutate(is_DEG = leaf_nat$GeneID %in% aging_genes_v3) %>%
  arrange(is_DEG) %>%
  ggplot(aes(x=log_ESL_ML,y=log_LSL_ML, col = is_DEG)) +
  geom_point()

leaf_mh_outlier <- effects %>%
  dplyr::filter(predictor=="Leaf", is_mh_outlier==TRUE)  %>%
  inner_join(v3_v5, by = c(gene="v5")) %>%
  pull(v3) %>% sort() %>% unique()

leaf_nat %>%
  dplyr::mutate(is_DEG = leaf_nat$GeneID %in% leaf_mh_outlier) %>%
  filter(abs(log_ESL_ML)>10 | abs(log_LSL_ML)>10, is_DEG==TRUE)

quartz()
leaf_nat %>%
  dplyr::mutate(is_DEG = leaf_nat$GeneID %in% leaf_mh_outlier) %>%
  arrange(is_DEG) %>%
  ggplot(aes(x=log_ESL_ML,y=log_LSL_ML, col = is_DEG)) +
  geom_point()


slim_comp@compareClusterResult %>% 
  filter( Cluster == "-P.Up" & ID =="GO:0016036") %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many")





with_known_gene <- slim_comp@compareClusterResult %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many") %>%
  select(Cluster,ID, Description, gene, SYMBOL,locus_name) %>%
  filter(!is.na(SYMBOL) & !is.na(locus_name)) %>%
  pull(ID) %>% sort() %>% unique()

# select the minimum (p_value per term)
sorted_by_p <- slim_comp@compareClusterResult %>%  
  filter( ID %in% with_known_gene) %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  dplyr::slice(1:6) %>%
  mutate(label.order = 1:n()) %>%
  mutate(Description= forcats::fct_reorder(Description,label.order))

top_terms <- levels(sorted_by_p$Description)


# slim_to_plot <-  rbind(
#   slim %>%
#     filter(Description %in% top_terms) %>%
#     filter(ID != "GO:0010358"),
#   slim %>%
#     filter(Description =="aging"))
# write.csv(slim_to_plot,"~/slim_to_plot.csv")

# to_bold <- c("GO:0015706","GO:0016036","GO:0010073","GO:0006032","GO:0009813","GO:0006479","GO:0019374")



to_plot <- slim_comp@compareClusterResult  %>%  
  filter( ID %in% with_known_gene) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  # droplevels() %>%
  filter(Description %in% top_terms) %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  mutate(Description = factor(Description, levels = top_terms)) %>% 
  mutate(label_order = 1:n()) %>%
  # droplevels() %>%
  mutate(Description= forcats::fct_reorder(Description,label_order))

# gene_symbol <- effects %>% 
#   mutate(Cluster = case_when(
#     regulation == "Upregulated" ~ paste(predictor,"Up", sep="."),
#     regulation == "Downregulated" ~ paste(predictor,"Down",sep=".")
#   )) %>%
#   dplyr::select(Cluster,gene, SYMBOL="locus_symbol",
#                 locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
#   filter(!is.na(SYMBOL) & !is.na(locus_name) & !is.na(Cluster)) %>%
#   mutate(neglogPG=-log10(P))  %>%
#   inner_join(TERM2GENE, relationship =  "many-to-many") 

effects_symbol <- effects %>% 
  mutate(Cluster = case_when(
    regulation == "Upregulated" ~ paste(predictor,"Up", sep="."),
    regulation == "Downregulated" ~ paste(predictor,"Down",sep=".")
  )) %>%
  dplyr::select(Cluster,gene, SYMBOL="locus_symbol",
                locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
  filter(!is.na(SYMBOL) & !is.na(locus_name) & !is.na(Cluster)) %>%
  mutate(neglogPG=-log10(P)) 

# gene_symbol <- read.table("data/gene_symbol.tab",quote="",header=TRUE, sep ="\t", na.strings = "")

term2symbol <- to_plot %>%  
  filter( ID %in% with_known_gene) %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(effects_symbol, relationship = "many-to-many") %>%
  dplyr::select(Cluster,gene,SYMBOL,ID, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG) %>%
  distinct() %>%
  group_by(Cluster,ID) %>%
  arrange(Cluster,ID, -mahalanobis) %>%
  dplyr::slice(1)

levels(to_plot$Cluster) <- cluster_order
levels(term2symbol$Cluster) <- cluster_order

labels <- paste0(
  gsub("."," ",cluster_order,fixed = TRUE),
  "\n", "(",n.genes,")"
)

names(labels) <-names(n.genes)


to_plot$Cluster %>% levels()
bg <- to_plot  %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, fill= neglogP)) +
  ggtitle("GO Biological Process Annotation of DEGs ") +
  scale_fill_continuous(low="dodgerblue", high="tomato", name = "-log10(FDR)",
                        guide=guide_colorbar())+ 
  scale_y_discrete( # labels = scales::wrap_format(50),
                   limits = rev(levels(to_plot$Description))) +
  scale_x_discrete(labels=labels,drop=FALSE) +
  geom_point(shape=21)  +
  scale_size(range = c(2, 8))+
  ylab(NULL) +  xlab(NULL) +
  DOSE::theme_dose(font.size = 14) +
  theme(
    plot.title=element_text(face = "bold", size=20),
    legend.position = c(0.87,0.67),
    legend.box = "horizontal", 
    axis.text.x = element_text(vjust=0.5))



annotation_plot <- bg + # data point size
  geom_text(
    data=term2symbol,
    position=position_nudge(0.1),
    hjust=0,vjust=0.5,
    fontface="italic",
    mapping =aes(x=Cluster,y =Description, label =SYMBOL),
    inherit.aes = FALSE
  )

quartz()
annotation_plot 

kegg <- readRDS(out, file ="~/Desktop/FC_kegg_plot.RDS")


out <- ggpubr::ggarrange(kegg,annotation_plot, ncol=1,align = "hv", heights = c(0.5,0.6))
quartz(height = 6.5, width=12)
out 
ggsave(out, file ="~/Desktop/FC_annotation_plot.svg",height = 6.75, width=12)


