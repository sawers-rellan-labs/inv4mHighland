library(dplyr)
library(rtracklayer)
library(GenomicRanges)
adjust.corr <- function (x, type = c("pearson", "spearman"), use = c("complete.obs", 
                                                                     "pairwise.complete.obs"), adjust = "holm") 
{
  opt <- options(scipen = 5)
  on.exit(options(opt))
  type <- match.arg(type)
  use <- match.arg(use)
  x <- if (use == "complete.obs") 
    as.matrix(na.omit(x))
  else as.matrix(x)
  R <- Hmisc::rcorr(x, type = type)
  P <- P.unadj <- R$P
  p <- P[lower.tri(P)]
  adj.p <- p.adjust(p, method = adjust)
  P[lower.tri(P)] <- adj.p
  P[upper.tri(P)] <- 0
  P <- P + t(P)
  P <- ifelse(P < 1e-04, 0, P)
  P <- format(round(P, 4))
  diag(P) <- ""
  P[c(grep("0.0000", P), grep("^ 0$", P))] <- "<.0001"
  P[grep("0.000$", P)] <- "<.001"
  P.unadj <- ifelse(P.unadj < 1e-04, 0, P.unadj)
  P.unadj <- format(round(P.unadj, 4))
  diag(P.unadj) <- ""
  P.unadj[c(grep("0.0000$", P.unadj), grep("^ 0$", P.unadj))] <- "<.0001"
  P.unadj[grep("0.000$", P.unadj)] <- "<.001"
  result <- list(R = R, P = P, P.unadj = P.unadj, type = type, 
                 adjustr = adjust)
  class(result) <- "adjust.corr"
  result
}

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

# sum(rownames(voomR$E) %in% flanking_introgression_gene_ids)

# sum(!rownames(voomR$E) %in% no_chr & !rownames(voomR$E) %in% shared_introgression_gene_ids)

effects <- read.csv("~/Desktop/predictor_effects.csv") %>%
  dplyr::filter(predictor=="GenotypeINV4") %>%
  mutate(is_significant = adj.P.Val <0.05,
         is_topDEG = is_significant & abs(logFC) > 2,
         in_Inv4m = gene %in% inv4m_gene_ids,
         in_shared = gene %in% shared_introgression_gene_ids,
         outside = !in_shared,
         in_flanking = gene %in% flanking_introgression_gene_ids
         ) %>%
  group_by(gene) %>%
  arrange(adj.P.Val) %>%
  dplyr::slice(1) %>%  # This slice is necessary because I have a many to many relationship
  ungroup()

effects$in_inv4m <-NULL 

nrow(effects)

sum(effects$in_shared)
sum(effects$is_significant)
all(inv4m_gene_ids %in% shared_introgression_gene_ids)


# DEG enrichment by genomic region ----
## DEGs -----
# odds of being signficant in shared introgression vs outside
so <- with(effects,
           table(in_shared,is_significant)
)
so
fisher.test(so)
fisher.test(so)$p.value

# odds of being signficant inv4m vs outside
io <- with(effects %>% filter(!in_flanking),
        table(in_Inv4m, is_significant)
)

io
fisher.test(io)
fisher.test(io)$p.value

# odds of being signficant flanking vs outside

fo <- with(effects %>% filter(in_flanking|outside),
 table(in_flanking,is_significant)
)
fo
fisher.test(fo)
fisher.test(fo)$p.value


is <- with(effects %>% filter(in_shared),
           table(in_Inv4m,is_significant)
)
is
fisher.test(is)
fisher.test(is)$p.value
## Top DEGs -----

with(effects,
     table(is_topDEG)
)

so <- with(effects,
           table(in_shared,is_topDEG)
)
so
fisher.test(so)
fisher.test(so)$p.value

# odds of being signficant inv4m vs outside
io <- with(effects %>% filter(!in_flanking),
           table(in_Inv4m, is_topDEG)
)

io
fisher.test(io)
fisher.test(io)$p.value

# odds of being signficant flanking vs outside

fo <- with(effects %>% filter(in_flanking|outside),
           table(in_flanking,is_topDEG)
)
fo
fisher.test(fo)
fisher.test(fo)$p.value


is <- with(effects %>% filter(in_shared),
           table(in_Inv4m,is_topDEG)
)
is
fisher.test(is)
fisher.test(is)$p.value


effects[effects$gene=="Zm00001eb193120",] %>% t()
effects[effects$gene=="Zm00001eb384800",] %>% t()

# new_label


v4_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v4_to_B73v5.tsv", sep = "\t", header= FALSE) %>%
  dplyr::rename(v4="V1",v5="V2")%>%
  tidyr::separate_longer_delim(cols = v5,delim=",")


regulators  <- effects %>% 
  # filter( in_shared == TRUE, is_significant==TRUE) %>%
  filter( in_shared== TRUE, is_DEG==TRUE) %>%
  arrange(adj.P.Val) %>% 
  dplyr::select(gene:desc_merged, logFC,adj.P.Val) %>%
  inner_join(v4_v5, by=c(gene="v5")) %>%
  group_by(gene) %>%
  arrange(adj.P.Val) 
  

regulators$gene %>% sort() %>% unique()


# flanking_introgression_gene_ids 
targets <-  effects %>%
  # filter( outside == TRUE, is_significant==TRUE) %>%
   filter( outside == TRUE, is_DEG==TRUE) %>%
  arrange(adj.P.Val) %>% 
  dplyr::select(gene:desc_merged, logFC,adj.P.Val) %>%
  inner_join(v4_v5, by=c(gene="v5")) %>%
  group_by(gene) %>%
  arrange(adj.P.Val) 


regulators$gene %>% sort() %>% unique()
targets$gene %>% sort() %>% unique()

# oxidative stress like correlated with zcn26
regulators[regulators$gene=="Zm00001eb193120",]



# get MaizeNetome data----


v4_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v4_to_B73v5.tsv", sep = "\t", header= FALSE) %>%
  dplyr::rename(v4="V1",v5="V2")%>%
  tidyr::separate_longer_delim(cols = v5,delim=",")

grassius <- read.csv("~/Desktop/grassius_annotation_v5.csv", na.strings = "") %>%
  inner_join(v4_v5, by=c(gene.ID ="v5")) %>%
  distinct()



# I split the search into 3
#
MaizeNetome <-rbind(
   #left: 4:155195539-170747954
  read.csv("~/Desktop/inv4m/coexpression/left_MaizeNetome_coexpression_edges.csv") %>%
    tidyr::separate_longer_delim(cols = Direct.Connetation.Node,delim="\n") %>%
    mutate(network="coexpression"),
  # inv4m: 4:170747955-186145605
  read.csv("~/Desktop/inv4m/coexpression/inv4m_MaizeNetome_coexpression_edges.csv") %>%
    tidyr::separate_longer_delim(cols = Direct.Connetation.Node,delim="\n") %>%
    mutate(network="coexpression"),
  read.csv("~/Desktop/inv4m/coexpression/right_MaizeNetome_coexpression_edges.csv") %>%
  # right: 4:186145606-195900523
    tidyr::separate_longer_delim(cols = Direct.Connetation.Node,delim="\n") %>%
    mutate(network="coexpression")
) %>% dplyr::rename( Direct.Connection="Direct.Connetation.Node")

all_ids <- MaizeNetome$Element.ID %>% sort %>% unique() 

noncoding <- all_ids[!grepl("Zm00001d",all_ids)]

MaizeNetome %>% pull(Element.ID) %>% sort () %>% unique %>% length()

MaizeNetome %>% pull(Direct.Connection) %>% sort () %>% unique %>% length()

# It's in v4

myGFF <- "~/ref/zea/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3"

v4<- rtracklayer::import(myGFF)  %>%
  subset(type=="gene")

B73 <- v4

genes <- as.data.frame(B73) 

genes$ID <- gsub("gene:","",genes$ID)

# Extract these coordinates from the MCScan data

# for Zm00001eb190470
inv4m_start <- genes[genes$ID=="Zm00001d051819","start"]
inv4m_start
#170747955

# for Zm00001eb194800
inv4m_end <- genes[genes$ID=="Zm00001d052292","end"]
inv4m_end
#186145605

#4:170747955-186145605

inv4m_gene_ids <- genes %>%
  filter(
    seqnames==4,
    start >= inv4m_start,
    end <= inv4m_end
  ) %>% pull(ID)

length(inv4m_gene_ids)

# checked manually in the RNAseq genotype table

introgression_start <-  155195539 # 157012149 
introgression_end <- 193981713 # 195900523

shared_introgression_gene_ids<-  genes %>%
  filter(
    seqnames==4,
    start >= introgression_start,
    end <= introgression_end
  ) %>% pull(ID)

length(shared_introgression_gene_ids)

flanking_introgression_gene_ids <- shared_introgression_gene_ids[!(shared_introgression_gene_ids %in% inv4m_gene_ids)]




target_edges <- MaizeNetome %>%
  left_join(grassius %>% select(protein.name,v4), by=c(Element.ID="v4")) %>%
  dplyr::select(network,Element.ID, protein.name, everything()) %>%
  distinct() %>%
  filter(Direct.Connection %in% targets$v4)

nrow(target_edges)


regulator_edges <- MaizeNetome %>%
  left_join(grassius %>% select(protein.name,v4), by=c(Element.ID="v4")) %>%
  dplyr::select(network,Element.ID, protein.name, everything()) %>%
  distinct() %>%
  filter(Element.ID %in% regulators$v4) %>% distinct()

nrow(regulator_edges)

two_way <- regulator_edges%>%
  filter(network == "coexpression") %>%
  filter(Direct.Connection %in% targets$v4)  %>%
  inner_join(v4_v5, by=c(Direct.Connection="v4")) %>% 
  inner_join(effects %>% dplyr::select(gene:desc_merged, logFC), by=c(v5="gene")) %>%
  dplyr::rename( target_gene="v5", target_symbol="locus_symbol" , target_desc="desc_merged", target_logFC="logFC") %>%
  inner_join(v4_v5, by=c(Element.ID="v4")) %>%
  inner_join(effects %>% 
               dplyr::select(v5="gene", introgressed_symbol="locus_symbol", introgressed_desc="desc_merged", introgressed_logFC="logFC",
                             in_shared, in_Inv4m)) %>%
  dplyr::select(network, introgressed_gene=v5, protein.name,introgressed_symbol,
                introgressed_desc,introgressed_logFC, in_shared, in_Inv4m, target_gene,
                target_symbol,target_desc,target_logFC) %>%
  group_by(network,introgressed_gene,target_gene) %>%
  arrange(abs(introgressed_logFC)) %>%
  dplyr::slice(1) %>%  # This slice is necessary because I have a many to many relationship
                       # multiple genes in one version mapping to multiple genes in other
                       # The right solution would be to build a filter of the
                       # allowed combinations and replace this slice
  filter(abs(introgressed_logFC) & abs(target_logFC)>2) %>%
  arrange(-abs(introgressed_logFC)) %>%
  filter(introgressed_gene !="Zm00001eb187900") %>% # plant height candidate. TOM like protein, response to auxin
                                                    # is these gene on the breakpoints of introgression/inversion
                                                    # it might be in two lists at the same time?
  tibble() %>% print(n=100) 
 nrow(two_way)  



inv4m_regulators<-  two_way$introgressed_gene[two_way$in_Inv4m] %>% sort() %>% unique()
length(inv4m_regulators)

inv4m_targets<-  two_way$target_gene[two_way$in_Inv4m] %>% sort() %>% unique()
length(inv4m_targets)

flanking_regulators <-  two_way$introgressed_gene[!two_way$in_Inv4m] %>% sort() %>% unique()

intersect(flanking_regulators,inv4m_regulators)

introgressed_regulators <- union(flanking_regulators,inv4m_regulators)
regulators %>%
  filter(locus_symbol=="jmj2")

effects[effects$gene=="Zm00001eb191790",] %>% t()
two_way[(two_way$introgressed_gene=="Zm00001eb191790"),-(1:3)]

# zcn26
effects[effects$gene=="Zm00001eb384770",] %>% t()
two_way[(two_way$target_gene=="Zm00001eb384770"),]

#zap1
effects[effects$gene=="Zm00001eb118120",] %>% t()
two_way[(two_way$target_gene=="Zm00001eb118120"),]


# OXS3 like gene correlated with zcn26
effects[effects$gene=="Zm00001eb193120",] %>% t()
two_way[(two_way$introgressed_gene=="Zm00001eb193120"),] 

# OXS3 like gene shows up in plant height and flowering time candidates
# however it's in the flanking regions :/
two_way[(two_way$introgressed_gene=="Zm00001eb186670"),] 


length(flanking_regulators)

outside_targets <-  two_way$target_gene %>% sort() %>% unique()
length(outside_targets )

intersect(flanking_regulators,outside_targets)
intersect(inv4m_regulators,outside_targets)

trans_network_genes <- union(c(inv4m_regulators,flanking_regulators),outside_targets)
length(trans_network_genes)
sum(trans_network_genes %in% regulators$gene)
sum(trans_network_genes %in% targets$gene)
trans_network_genes[ !(trans_network_genes %in% regulators$gene) & !(trans_network_genes %in% targets$gene)]
trans_network_genes[ trans_network_genes %in% union(inv4m_regulators,flanking_regulators)]
trans_network_genes[ trans_network_genes %in% union(introgressed_regulators,flanking_regulators)]


# ft_xp <- psu %>% 
#   right_join(
#     cbind(
#       row=voomR$targets$row,
#       as.data.frame(t(voomR$E)[,of_interest])


# plant_csv <- "/Users/fvrodriguez/Desktop/Desktop/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"
# 
# psu <- read.csv(plant_csv) %>%
#   dplyr::select(row=P22.,Anthesis=DTA,Silking=DTS)

voomR <- readRDS("~/Desktop/normalized_expression_voom_object.rda")

of_interest <- trans_network_genes[trans_network_genes %in% rownames(voomR$E)]
outside_targets <-  outside_targets[outside_targets %in% rownames(voomR$E)]
introgressed_regulators <-introgressed_regulators[introgressed_regulators %in% rownames(voomR$E)]

xp <- t(voomR$E)[,of_interest]
order <- of_interest
length(order)
corr <- round(cor(xp, use = "pairwise.complete.obs"),2)
dim(corr)
p_adj <- adjust.corr(xp,adjust = "fdr")
p_adj$P.unadj
to_plot <- corr[order,order]
to_plot[p_adj$P[order,order]>0.05] <-0

# to_plot <- to_plot[order,order %>% rev()]

# x <- to_plot[inv4m_regulators,inv4m_targets ] !=0
outside_targets[!outside_targets %in% order]
x <- to_plot[introgressed_regulators,outside_targets] 
dim(x)

edge <- cbind(c(row(x)), c(col(x)), c(x))

idx <- rownames(x)[edge[,1]]
idy <- colnames(x)[edge[,2]]
all_edge <- data.frame(
  from = idx,
  to = idy,
  r = edge[,3]
)
n_samples <- ncol(voomR$E)
d <- function(r,n){sqrt(2 * n * (1 - abs(r)))}
d(all_edge$r,n_samples)
cost <- function(r){log2(1/(1-r^2))}
cost(all_edge$r)



edge_list <- all_edge  %>%
  filter(from %in% inv4m_regulators, to %in% outside_targets) %>% 
  mutate(e=case_when(abs(r)>0 ~1, .default = 0),
         weight=cost(r),
         distance=d(r,n_samples)) %>%
  dplyr::filter(e==1) %>%
  dplyr::filter(from!=to) %>%
  dplyr::mutate(color="lightgrey") %>%
  #dplyr::mutate(color= case_when(r>0 ~ "dodgerblue", r<0 ~ "tomato", .default="black")) %>%
  arrange(from,-abs(r)) %>% tibble()





edge_list <- all_edge  %>%
  inner_join(two_way %>% dplyr::select(from="introgressed_gene", to="target_gene")) %>% 
  mutate(e=case_when(abs(r)>0 ~1, .default = 0),
         weight=cost(r),
         distance=d(r,n_samples)) %>%
  dplyr::filter(e==1) %>%
  dplyr::filter(from!=to) %>%
  dplyr::mutate(color="lightgrey") %>%
  #dplyr::mutate(color= case_when(r>0 ~ "dodgerblue", r<0 ~ "tomato", .default="black")) %>%
  arrange(from,-abs(r)) %>% tibble()



nrow(edge_list)

length(edge_list$from %>% sort %>% unique())

length(edge_list$to%>% sort %>% unique())

library(igraph)
library(visNetwork)
edge_list

g <- graph_from_edgelist(edge_list[,1:2] %>% as.matrix(),directed = TRUE) %>%
  set_edge_attr("color", value = edge_list$color) %>%
  set_edge_attr("arrows.to.type", value = c("bar","arrow")[(edge_list$r>0)+1]) %>%
  set_edge_attr("width", value = 5) 
  # set_edge_attr("width", value = ceiling(3*edge_list$r)+2) 

V(g)
e <-E(g)
in_inv4m <- names(V(g)) %in% inv4m_regulators
in_flanking <- names(V(g)) %in% flanking_regulators

upregulated <- names(V(g)) %in% effects$gene[effects$logFC>0]
region <- case_when(
  names(V(g)) %in% inv4m_regulators ~ "inv4m",
  names(V(g)) %in% flanking_regulators ~ "flank",
  .default ="out"
) %>% as.factor()



g <- g %>%
  set_vertex_attr("region", value = region) %>%
  set_vertex_attr("color.border", value = c("#663399","gold")[region]) %>%
  set_vertex_attr("borderWidth", value = 3) %>%
  set_vertex_attr("color.background", 
                  value = case_when(
                    upregulated & in_inv4m ~ "#663399",
                    upregulated & !in_inv4m ~ "gold",
                    !upregulated ~ "white",
                    .default = "white")
  )

g <- g %>%
  set_vertex_attr("region", value = region) %>%
  set_vertex_attr("color.border", value = c("deeppink","#663399","gold")[region]) %>%
  set_vertex_attr("borderWidth", value = 3) %>%
  set_vertex_attr("color.background", 
                  value = case_when(
                    upregulated & in_inv4m ~ "#663399",
                    upregulated & in_flanking ~ "deeppink",
                    upregulated & !in_inv4m ~ "gold",
                    !upregulated ~ "white",
                    .default = "white")
                  )

#c("triangleDown","triangle")[upregulated +1]) 


gene_locus <- data.frame(gene=names(V(g))) %>%
  left_join(effects  %>%
              group_by(gene) %>%
              dplyr::slice(1)
            ) %>%
  ungroup()
nrow(gene_locus)
node_gene <-gene_locus$gene
node_symbol <- gene_locus$locus_symbol
node_label <- node_gene
length(node_label)
node_label[!is.na(node_symbol)] <-node_symbol[!is.na(node_symbol)]
names(node_label) <- node_gene
length(node_label)
node_label
new_label <-c()
# I need to make a table with these names.
# maybe the same I need to add to the manuscript
new_label["Zm00001eb384770"] <-"zcn26"
new_label["Zm00001eb192330"] <-"roc4-1"
new_label["Zm00001eb193260"] <- "trxl2"
new_label["Zm00001eb194300"] <- "paspa" # nana-like
new_label["Zm00001eb213980"] <- "dywl"
new_label["Zm00001eb194080"] <- "ugt85a2"
new_label["Zm00001eb112470"] <- "engd2"
new_label["Zm00001eb332450"] <- "roc4-2"
new_label["Zm00001eb193270"] <- "trxl5"
new_label["Zm00001eb066620"] <- "tipin"
new_label["Zm00001eb241410"] <- "etfa"
new_label["Zm00001eb013800"] <- "mrlk"
new_label["Zm00001eb249540"] <- "prpl29"
new_label["Zm00001eb000540"] <- "prd3"
new_label["Zm00001eb060540"] <- "actin-1"
new_label["Zm00001eb192580"] <- "map4k8"
new_label["Zm00001eb193120"] <- "o3l1"

replace_label <- new_label[names(new_label) %in% gene_locus$gene]
node_label[names(replace_label)] <- replace_label

g <- g %>%
  set_vertex_attr("label", value = node_label) 

trans_introgression_nodes <- names(V(g))

library(visNetwork)

visNetwork::visIgraph(g,idToLabel = FALSE,layout = "layout.sphere")

visNetwork::visIgraph(mst(g),idToLabel = FALSE)



data <- toVisNetworkData(g)
data$edges$arrows.to.type
data$nodes$color.border <- V(g)$color.border
data$nodes$color.background <- V(g)$color.background
data$nodes$borderWidth <- 4


data$nodes$label <- paste("<i>",
  node_label[data$nodes$id],
  "</i>"
)


edges <- data$edges%>% filter(from %in% flanking_regulators,!to %in% flanking_regulators)

nrow(edges)

edges$from %>% sort() %>% unique() %>% length()
edges$to %>% sort() %>% unique() %>% length()

flanking_edges <- edges

nodes <- data$nodes %>% filter( id %in% union(edges$from, edges$to))

trans_flanking_nodes <- nodes$id

visNetwork(nodes = nodes, edges = edges, width = "100%", height = 700)  %>%
  visEdges(arrows ="to") %>%
  visNodes(font = list(size = 20, multi = "html"))

edges <- data$edges%>% filter(!from %in% flanking_regulators,!to %in% flanking_regulators)
nodes <- data$nodes %>% filter(!id %in% flanking_regulators,
                               id %in% edges$to | id %in% edges$from)

nrow(edges)

edges$from %>% sort() %>% unique() %>% length()
edges$to %>% sort() %>% unique() %>% length()

visNetwork(nodes = nodes,
           edges =  edges,
           width = "100%", height = 700)  %>%
  visEdges(arrows ="to") %>%
  visNodes(font = list(size = 25, multi = "html"))

trans_Inv4m_nodes <- nodes$id


edges <- data$edges%>% filter(!from %in% flanking_regulators,
                              !to %in% flanking_regulators,
                              !to %in% flanking_edges$to)

edges$from %>% sort() %>% unique() %>% length()
edges$to %>% sort() %>% unique() %>% length()

flanking_edges <- edges

nodes <- data$nodes %>% filter(!id %in% flanking_regulators,
                               id %in% edges$to | id %in% edges$from)


visNetwork(nodes = nodes,
           edges =  edges,
           width = "100%", height = 700)  %>%
  visEdges(arrows ="to") %>%
  visNodes(font = list(size = 25, multi = "html"))


network_effects  <- effects  %>%
  dplyr::filter(predictor=="GenotypeINV4") %>%
  mutate(predictor="Inv4m") %>%
  mutate(is_significant = adj.P.Val <0.05,
         is_topDEG = is_significant & abs(logFC) > 2,
         in_Inv4m = gene %in% inv4m_gene_ids,
         in_shared = gene %in% shared_introgression_gene_ids,
         outside = !in_shared,
         in_flanking = gene %in% flanking_introgression_gene_ids
  ) %>%
  group_by(gene) %>%
  arrange(adj.P.Val) %>%
  dplyr::slice(1) %>%  # This slice is necessary because I have a many to many relationship
  ungroup()

write.csv(network_effects,file="~/Desktop/network_effects.csv")

