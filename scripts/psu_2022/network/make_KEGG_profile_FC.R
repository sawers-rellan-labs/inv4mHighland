library(dplyr)
library(ggplot2)
library(clusterProfiler)

id_map <- read.csv("/Users/fvrodriguez/Desktop/Gene_xRef.csv", header= FALSE)

colnames(id_map) <-c("gene", "v5")

effect_order <-  c("Leaf","-P","inv4m")

# Imposing fold change limits

effects <- read.csv("~/Desktop/DEGs_FDR_for_table.csv") %>%
  rename(gene="v5") %>%
  inner_join(id_map) %>%
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


DEGs <- with( effects,{
  list(
    Leaf.Down = gene[predictor=="Leaf" & regulation=="Downregulated"],
    Leaf.Up = gene[predictor=="Leaf" & regulation=="Upregulated"],
   # Leaf.all = gene[predictor=="Leaf" & regulation!="Unregulated"],
    `-P.Down` = gene[predictor=="-P" & regulation=="Downregulated"],
    `-P.Up` = gene[predictor=="-P" & regulation=="Upregulated"],
   # `-P.all` = gene[predictor=="-P" & regulation!="Unregulated"],
    inv4m.Down = gene[predictor=="inv4m" & regulation=="Downregulated"],
    inv4m.Up = gene[predictor=="inv4m" & regulation=="Upregulated"]
   # inv4m.all = gene[predictor=="inv4m" & regulation!="Unregulated"]
  )}
)



n.genes <- lapply(DEGs, length) %>% unlist()

compKEGG <-lapply(DEGs, FUN=function(x){
enrichKEGG(gene         = x,
                 organism     = 'zma',
                 pvalueCutoff = 0.05)
})


cluster_result <- lapply(names(compKEGG), FUN=function(x) {
  result <- compKEGG[[x]]@result
  result %>%
    filter(p.adjust <0.05) %>%
    mutate(Cluster=x) %>%
    select(Cluster,everything())
}) %>% dplyr::bind_rows()

cluster_order <- names(DEGs)



sorted_by_p <- cluster_result  %>%
  mutate(neglogP= -log10(p.adjust)) %>%
  mutate(Cluster=factor(Cluster, levels=cluster_order)) %>%
  droplevels() %>%
  group_by(Cluster) %>%
  arrange(Cluster,p.adjust) %>%
  dplyr::slice(1:20) %>%
  mutate(label.order = 1:n()) %>%
  mutate(Description= forcats::fct_reorder(Description,label.order))

top_terms <- levels(sorted_by_p$Description)


cluster_result




to_plot <- cluster_result  %>%
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

labels <- paste0(
  gsub("."," ",cluster_order,fixed = TRUE),
  "\n", "(",n.genes,")"
  )


effects_symbol <- effects %>% 
  mutate(Cluster = case_when(
    regulation == "Upregulated" ~ paste(predictor,"Up", sep="."),
    regulation == "Downregulated" ~ paste(predictor,"Down",sep=".")
  )) %>%
  dplyr::select(Cluster,v5,gene, SYMBOL="locus_symbol",
                locus_name, logFC, P= adj.P.Val,  mahalanobis) %>% 
 # filter(!is.na(SYMBOL) & !is.na(locus_name) & !is.na(Cluster)) %>%
  mutate(neglogPG=-log10(P)) 

# gene_symbol <- read.table("data/gene_symbol.tab",quote="",header=TRUE, sep ="\t", na.strings = "")

term2symbol <- to_plot %>%   
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  rename(geneID="gene") %>%
  mutate(gene = as.integer(gene)) %>%
  inner_join(effects_symbol) %>%
  dplyr::select(Cluster,gene,SYMBOL,ID, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG) %>%
  distinct() %>%
  group_by(Cluster,ID) %>%
  arrange(Cluster,ID, -mahalanobis) %>%
  dplyr::slice(1)



to_plot %>%   
  tidyr::separate_longer_delim(geneID,delim="/") %>%
  rename(geneID="gene") %>%
  mutate(gene = as.integer(gene)) %>%
  inner_join(effects_symbol) %>%
  dplyr::select(Cluster,v5,gene,SYMBOL,ID, Description,p.adjust,Count,mahalanobis,logFC,P,neglogPG) %>%
  distinct() %>%
  group_by(Cluster,ID) %>%
  arrange(Cluster,ID, -mahalanobis) %>%
  filter(ID=="zma00196")



bg <- to_plot %>%
  ggplot(aes(x = Cluster, y = Description, size = Count, fill= neglogP)) +
  ggtitle("KEGG Pathway Annotation of DEGs ") +
  scale_fill_continuous(low="dodgerblue", high="tomato", name = "-log10(FDR)",
                        guide=guide_colorbar())+ 
  scale_y_discrete( # labels = scales::wrap_format(50),
    limits = rev(levels(to_plot$Description))) +
  scale_x_discrete(labels=labels,drop=FALSE) +
  geom_point(shape=21)  +
  scale_size(range = c(2, 8))+
  ylab(NULL) +  xlab(NULL) +
  DOSE::theme_dose(font.size = 12) +
  theme(
    legend.position = c(0.87,0.6),
    legend.box = "horizontal", 
    axis.text.x = element_text(vjust=0.5))


out <- bg + # data point size
  geom_text(
    data=term2symbol,  position=position_nudge(0.1), hjust=0,vjust=0.5,
    fontface="italic",
    mapping =aes(x=Cluster,y =Description, label =SYMBOL),
    inherit.aes = FALSE
  )

quartz()
out

saveRDS(out, file ="~/Desktop/FC_kegg_plot.RDS")

