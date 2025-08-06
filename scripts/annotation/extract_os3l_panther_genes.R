library(dplyr)
library(rbioapi)

v4_v5 <- read.table(file="/Users/fvrodriguez/Desktop/B73_gene_xref/B73v4_to_B73v5.tsv", sep = "\t", header= FALSE) %>%
  dplyr::rename(v4="V1",v5="V2")%>%
  tidyr::separate_longer_delim(cols = v5,delim=",")

panther <- read.table("~/Desktop/PTHR33172.tab", sep="\t", header=TRUE)
nrow(panther)

target_organisms <- c(
  `Marchantia polymorpha` = 3197,
  `Arabidopsis thaliana` =3702,
  `Oryza sativa` = 4530,
  `Zea mays` = 4577
)

rba_panther_family(
  "PTHR33172",
  what="ortholog", 
  # target_organisms = target_organisms
)
rba_options(timeout = 10)

tree <- rba_panther_family(
  "PTHR33172",
  what="tree", 
  # target_organisms = target_organisms
)

# I don't know how to deal with this tree structure
# It's a nested list of NAMED dataframes from the json output of the panther API.
# however I'm here for the uniprot identifiers:

unlisted <- unlist(tree)

unlisted[grepl(":SF",unlisted,fixed=TRUE)] %>% unname()

unlisted[grepl("MARPO|EnsemblGenome",unlisted, fixed=TRUE)] 
long_name <- "^annotation_node.children.annotation_node.children.annotation_node.children.annotation_node.children.annotation_node.node_name1"
unlisted[grepl(long_name,names(unlisted), perl = TRUE)]

cat(protein$uniprot)
cat(protein$locus)
protein <- data.frame(
  Gene_ID = c(
  unlisted[grepl("MARPO|EnsemblGenome",unlisted, fixed=TRUE)] %>%  unname(),
  unlisted[grepl("ARATH",unlisted, fixed=TRUE)] %>%  unname(),
  unlisted[grepl("MAIZE|EnsemblGenome",unlisted, fixed=TRUE)] %>%  unname(),
  unlisted[grepl("ORYSJ|EnsemblGenome",unlisted, fixed=TRUE)] %>%  unname()
)
) %>%
  mutate(Gene_ID = gsub("=locus=","=",Gene_ID)) %>%
  mutate(Gene_ID = gsub("|UniProtKB=","=",Gene_ID,fixed = TRUE)) %>%
  mutate(Gene_ID = gsub("|","=",Gene_ID,fixed = TRUE))  %>%
  tidyr::separate(Gene_ID, sep="=", into=c("orgcode","db", "locus","uniprot"))

#Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Interacts with	Gene Ontology (biological process)	Subcellular location [CC]	PubMed ID	Domain [CC]	Protein families	Domain [FT]	Sequence

nrow(protein)

library(httr)

res = GET("https://rest.uniprot.org/uniprotkb/search",
    query = list(
      query=paste0("accession:",
                   paste(protein$uniprot,collapse = " OR accession:")), 
      fields=paste("organism_name, organism_id","gene_names","lit_pubmed_id","accession","sequence", sep=","),
      format="tsv",
      size =nrow(protein)), verbose()
)


res$request
res$status_code
res$content

selected <-protein %>%
  left_join(
  read.table(
  text= memDecompress(res$content,
                     type = "gzip",
                     asChar = TRUE), 
           header=TRUE, sep ="\t"
  ) %>% arrange(Organism..ID.,Entry),
  by=c(uniprot="Entry") 
) %>%
  dplyr::select(orgcode,taxid="Organism..ID." ,db:uniprot, ,
                Gene.Names, Sequence,PMID="PubMed.ID")

write.csv(selected, file="~/Desktop/PTHR33172_selected_protein.csv")



# I'd like to know where did I get the
panther <- read.table("~/Desktop/PTHR33172.tab", sep="\t", header=TRUE)

panther %>% select(Gene:Organism,Publications,Gene.ID,Protein.name) %>%
  arrange(-Publications) %>%
  filter(Organism%in% c(
    "Arabidopsis thaliana",
    "Oryza sativa",
    "Zea mays",
    "Marchantia polymorpha")
  ) %>% inner_join(v4_v5, by= c(Gene.ID="v4"))

