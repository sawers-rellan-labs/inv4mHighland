library(dplyr)
library(ggplot2)
library(clusterProfiler)

effect_order <-  c("Leaf","-P","inv4m")

# Without imposing fold change limits

effects <- read.csv("~/Desktop/DEGs_FDR_for_table.csv") %>%
  mutate(neglogP = -log10(adj.P.Val)) %>%
  mutate(is_significant = adj.P.Val<0.05) %>%
  mutate(is_up = case_when(
    (logFC > 0) & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(is_down = case_when(
    (logFC < 0)  & is_significant ~ TRUE,
    .default = FALSE))  %>%
  mutate(regulation = case_when(
    is_significant & is_up ~ "Upregulated",
    is_significant & is_down ~ "Downregulated",
    .default = "Unregulated"
  )) %>%
  mutate(is_DEG= is_significant & regulation != "Unregulated") %>%
  arrange(predictor,adj.P.Val)


universe <- effects$gene[effects$predictor=="Leaf"]
length(universe)

ego <- function(x, ontology ="BP"){ 
  TERM2NAME <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2NAME.rds")
  
  # Load GO annotation
  TERM2GENE <- readRDS("/Users/fvrodriguez/Desktop/GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024_full.rds")
  
  effects$gene[effects$predictor=="Leaf"]
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
    TERM2GENE= TERM2GENE[TERM2GENE$GO %in% TERM2NAME[[ontology]]$go_id,],
    TERM2NAME = TERM2NAME[[ontology]]
  )}



DEGs <- with( effects,{
list(
  Leaf.Down = gene[predictor=="Leaf" & is_down],
  Leaf.Up = gene[predictor=="Leaf" & is_up],
  `-P.Down` = gene[predictor=="-P" & is_down],
  `-P.Up` = gene[predictor=="-P" & is_up],
  inv4m.Down = gene[predictor=="inv4m" & is_down],
  inv4m.Up = gene[predictor=="inv4m" & is_up]

)}
)

n.genes <- lapply(DEGs, length) %>% unlist

n.genes

compGO_BP <- compareCluster(geneCluster   = DEGs,
                            fun           = ego)

compGO_BP@compareClusterResult$Cluster %>% table()


# write.csv(compGO_BP,file="~/Desktop/inv4mP_FC_ALL_goa_Fattel_full.csv")

levels(compGO_BP@compareClusterResult$Cluster) <- cluster_order

quartz()
dotplot(compGO_BP) +
  scale_x_discrete(drop=FALSE) 

curated_goa <-NULL
curated_goa <-read.csv(file="~/Desktop/inv4mP_FDR_Fattel_full_goa_curated.csv")

slim <- curated_goa %>%
  dplyr::filter(to_slim2==TRUE) %>%
  #dplyr::select(ontology,ID,Description) %>%
  dplyr::select(ID,Description)  %>%
  dplyr::arrange(ID) %>%
  dplyr::distinct() 


slim_comp <- compGO_BP
slim_comp@compareClusterResult <- slim_comp@compareClusterResult %>%
  filter(ID %in% slim$ID)
slim_comp@compareClusterResult$Cluster %>% table()


cluster_order <- c("Leaf.Down","Leaf.Up","-P.Down","-P.Up","inv4m.Down","inv4m.Up")

# slect the minimum (p_value per term)
sorted_by_p <- slim_comp@compareClusterResult  %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  dplyr::slice(1:10) %>%
  mutate(label.order = 1:n()) %>%
  mutate(Description= forcats::fct_reorder(Description,label.order))

top_terms <- levels(sorted_by_p$Description)

# to_bold <- c("GO:0015706","GO:0016036","GO:0010073","GO:0006032","GO:0009813","GO:0006479","GO:0019374")

to_plot <-compGO_BP@compareClusterResult  %>%
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

gene_symbol <- effects %>% 
  mutate(Cluster = case_when(
    regulation == "Upregulated" ~ paste(predictor,"Up", sep="."),
    regulation == "Downregulated" ~ paste(predictor,"Down",sep=".")
  )) %>%
  dplyr::select(Cluster,gene, SYMBOL="locus_symbol",
                locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
  filter(!is.na(SYMBOL) & !is.na(locus_name) & !is.na(Cluster)) %>%
  mutate(neglogPG=-log10(P))  %>%
  inner_join(TERM2GENE, relationship =  "many-to-many") 


term2symbol <- to_plot %>%   
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  dplyr::rename(GO="ID") %>%
  inner_join(gene_symbol, relationship = "many-to-many") %>%
  dplyr::select(Cluster,gene,SYMBOL,GO, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG) %>%
  distinct() %>%
  group_by(Cluster,GO) %>%
  arrange(Cluster,GO, -mahalanobis) %>%
  dplyr::slice(1)

levels(to_plot$Cluster) <- cluster_order
levels(term2symbol$Cluster) <- cluster_order

labels <- paste0(
  c("Leaf Down",
    "Leaf Up", 
    "-P Down",
    "-P Up",
    "Inv4m Up",
    "Inv4m Down"),
  "\n", "(",n.genes,")" 
)


names(labels) <-c( "Leaf.Down",
                   "Leaf.Up", 
                   "-P.Down",
                   "-P.Up",
                   "inv4m.Up",
                   "inv4m.Down")
to_plot$Cluster %>% levels()
bg <- to_plot %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, fill= neglogP)) +
  ggtitle("Biological Process Annotation of DEGs ") +
  scale_fill_continuous(low="dodgerblue", high="tomato", name = "-log10(FDR)",
                        guide=guide_colorbar())+ 
  scale_y_discrete(labels = scales::wrap_format(50),
                   limits = rev(levels(to_plot$Description))) +
  scale_x_discrete(labels=labels,drop=FALSE) +
  geom_point(shape=21)  +
  scale_size(range = c(2, 8))+
  ylab(NULL) +  xlab(NULL) +
  DOSE::theme_dose() +
  theme(axis.text.x = element_text(vjust=0.5))


out <- bg + # data point size
  geom_text(
    data=term2symbol,  position=position_nudge(0.1), hjust=0,vjust=0.5,
    fontface="italic",
    mapping =aes(x=Cluster,y =Description, label =SYMBOL),
    inherit.aes = FALSE
  )

quartz(height=6, width = 12)
print(out)
ggsave(file="~/Desktop/editme.svg", plot=out, width=12, height=6)



