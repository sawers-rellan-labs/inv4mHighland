library(dplyr)
library(ggplot2)
library(clusterProfiler)

effects <- read.csv("~/Desktop/network_effects.csv") 

with(effects,
     table(predictor,regulation)
)

background <-  effects %>%
  filter(predictor=="Inv4m") %>%
  pull(gene) %>%
  sort() %>%
  unique()

length(background)

# Load GO annotation
TERM2NAME <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2NAME.rds")

TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024_full.rds")

ego <- function(x, ontology ="BP", universe=NULL,...){ 
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

# Load network from analyze MaizeNetome


DEGs <-   list(
    `Introgression`= trans_introgression_nodes,
    `Flanking` = trans_flanking_nodes,
    `Inv4m` = trans_Inv4m_nodes
  )

n.genes <- lapply(DEGs, length) %>% unlist
n.genes

# GSEA works well with large lists
# to_GSEA <- effects %>%
#   dplyr::filter(is_significant) %>%
#   arrange(-abs(logFC))
# 
# geneList <- abs(to_GSEA$logFC)
# names(geneList) <- to_GSEA$gene
# 
# y <- GSEA(geneList,
#           TERM2NAME = TERM2NAME[["BP"]],
#           TERM2GENE = TERM2GENE[TERM2GENE$GO %in% TERM2NAME[["BP"]]$go_id,])
#
# y@result$Description
# quartz()
#
# dotplot(y)


compGO_BP <- compareCluster(
  geneCluster   = DEGs,
  fun           = function(x){
    ego(x, ontology ="BP", universe=background)
    }
)


# another_bg <- DEGs %>% unlist() %>% sort() %>% unique()
# 
# compGO_BP <- compareCluster(
#   geneCluster   = DEGs,
#   fun           = function(x){
#     ego(x, ontology ="BP", universe=another_bg)
#   }
# )

compGO_BP@compareClusterResult$Cluster %>% table()

compGO_BP@compareClusterResult %>%
  filter(Cluster=="Inv4m") %>%
  filter(ID =="GO:0008283") %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  inner_join(
    effects %>%
      select(gene,SYMBOL=locus_symbol,locus_name) %>%
      distinct(),
    relationship = "many-to-many") 




cluster_order <- names(DEGs)


levels(compGO_BP@compareClusterResult$Cluster) <- cluster_order

slim_comp <- compGO_BP

slim_comp@compareClusterResult <- compGO_BP@compareClusterResult %>% 
  filter(Count >1) %>%
  #filter(ID != "GO:0051567") %>% # No protein annotated with histone demethylase molecular function
  filter(!ID %in% c("GO:0006305","GO:0006306","GO:0006304")) %>% # PCNA2 doesnt't have a DNA alkylation/methylation function 
    group_by(Cluster,geneID) %>% # this grouping reduces semantic redundancy!!! I should try it elsewhere
  arrange(p.adjust,.by_group = TRUE) %>%
  dplyr::slice(1)



slim_comp@compareClusterResult$Cluster %>% table()

quartz()
dotplot(slim_comp,size="Count",  showCategory = 10) +
  ggtitle("GO annotation of Inv4m\ntrans coexpression network DEGs") +
  scale_size(range =c(4,8),breaks = (1:4)) +
  scale_x_discrete(drop=FALSE) 

# Nice!


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

# idon't really need this bit
# I only have a few terms that can be plotted
sorted_by_p <- slim_comp@compareClusterResult %>%  
  filter( ID %in% with_known_gene) %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  dplyr::slice(1:10) %>%
  mutate(label.order = 1:n()) %>%
  mutate(Description= forcats::fct_reorder(Description,label.order))

top_terms <- levels(sorted_by_p$Description)


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

effects_symbol <- effects %>% 
  filter(is_significant) %>%
  dplyr::select(gene, SYMBOL="locus_symbol",
                locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
  filter(!is.na(SYMBOL) & !is.na(locus_name)) %>%
  mutate(neglogPG=-log10(P)) 

to_plot %>%  
  #filter( ID %in% with_known_gene) %>%
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(gene="geneID") %>%
  left_join(effects_symbol, relationship = "many-to-many") %>%
  dplyr::select(Cluster,gene,SYMBOL,ID, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG)

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

to_plot$ID
  
bg <- to_plot  %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, fill= neglogP)) +
  ggtitle("GO Biological Process Annotation of Top DEGs") +
  # scale_fill_continuous(low="dodgerblue", high="tomato", name = "-log10(FDR)",
  #                       guide=guide_colorbar()) + 
  scale_fill_gradient(low="dodgerblue", high="tomato", name = "-log10(FDR)",
                      limits = c(1,4), breaks=1:4)+
   scale_y_discrete(  labels = scales::wrap_format(60),
    limits = rev(levels(to_plot$Description))) +
  scale_x_discrete(labels=labels,drop=FALSE) +
  geom_point(shape=21)  +
  # scale_size(range =c(5,9),breaks = c(3,6,9)) +
  scale_size(range =c(6,9),breaks = 2:4) +
  ylab(NULL) +  xlab(NULL) +
  DOSE::theme_dose(font.size = 14) +
  theme(
    plot.title=element_text(face = "bold", size=20),
    legend.position = c(0.25,0.3),
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


quartz(height = 4, width=11)
annotation_plot 

ggsave(annotation_plot, file ="~/Desktop/Network_annotation_plot.svg",height = 4, width=11)

