library(dplyr)

fulldata <- read.table("~/ref/zea/Zm00001eb.1.fulldata.txt", 
                       na.strings = c(""," "),sep= "\t", quote='"',fill = TRUE, header=TRUE)

fulldata$locus_symbol <- gsub(" ","",fulldata$locus_symbol)
fulldata$locus_symbol[fulldata$locus_symbol==""] <- NA
fulldata$locus_name <- gsub("^ +","",fulldata$locus_name, perl =TRUE)
fulldata$locus_name <- gsub(" +$","",fulldata$locus_name, perl =TRUE)
fulldata$locus_name[fulldata$locus_name==""] <- NA
fulldata$canonical_transcript_name <- gsub(" ","",fulldata$canonical_transcript_name)
fulldata$canonical_transcript_name[fulldata$canonical_transcript_name==""] <- NA

head(fulldata)

colnames(fulldata)
model_names <- fulldata[,c("gene_model","canonical_transcript_name","locus_symbol","locus_name")]

# GRAMENE PANNZER annotation 2022 I don't know Why the track data for gbrowse, which contains 
# a more recent PANNZER gene descripion is not aavaliable somewhere

pannzer_desc <-read.table("~/ref/zea/B73_DE.out", 
                          na.strings = c(""," "),sep= "\t", quote='"',fill = TRUE, header=TRUE)
  
# select best annotation by GSZ
# https://academic.oup.com/bioinformatics/article/31/10/1544/176441

# Check that it doesn't match any non-number
numbers_only <- function(x) !grepl("\\D", x)


gene_desc <- model_names %>%
  mutate(qpid = gsub("_T","_P",canonical_transcript_name)) %>%
  left_join(pannzer_desc, relationship = "many-to-many") %>%
  group_by(qpid) %>%
  arrange(-cluster_GSZ) %>%
  dplyr::slice(1) %>% 
  ungroup() %>%
  dplyr::select(gene_model,locus_symbol,locus_name,desc,genename) %>%  
  mutate(label = case_when(!is.na(locus_name) ~ locus_symbol,
                          is.na(locus_symbol) & !is.na(genename) & !numbers_only(genename)   ~ genename,
                           .default = NA)) %>%
  mutate(description= coalesce(locus_name,desc)) %>%
  dplyr::select(gene_model, label, description) %>%
  left_join(grassius %>% select(protein.name, gene_model="gene.ID")) %>%
  mutate(label= coalesce(label,protein.name)) %>%
  select(-protein.name) %>% 
  mutate(label=tolower(label)) %>%
  distinct() %>%
  print(n=100)

gene_desc[gene_desc$gene_model %in% effects$gene[effects$is_DEG],] %>% print(n=500)


