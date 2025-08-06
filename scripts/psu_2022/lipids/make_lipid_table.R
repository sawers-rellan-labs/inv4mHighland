library(dplyr)
library(janitor)
library(ggplot2)


pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")


# filter internal standards and odd chains
to_exclude<- c("DG_12_0", "LPC_17_0", "LPE_17_1",
               "PC_25_0", "PG_17_0","SM_35_1",
               "Sphingosine_17_1","TG_17_0","TG_57_6" 
)

#lipid_csv <- "data/PSU_Normalized_MSDial_240416.csv"

lipid_csv <- "data/PSU_RawData_MSDial_240416.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) %>%
  dplyr::select(-all_of(to_exclude)) %>%
  filter(!grepl("Methanol",sampleID)) %>%
  filter(!grepl("Pool",sampleID))
nrow(raw)

colnames(raw)
m <- raw %>% select(DGDG_34_0:TG_58_5) %>% as.matrix()
row.names(m) <- raw$sampleID

plant_csv <- "data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"

psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", rowid = "P22." ) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

sampleInfo <- read.csv('data/inv4mRNAseq_metadata.csv') %>%
  rename(Genotype=genotype) %>% 
  rename(rowid=row)



head_group <- c(
  DG ="glycolipid", 
  DGDG="glycolipid",
  DGGA="glycolipid",
  LPC="phospholipid",
  LPE="phospholipid",
  MGDG="glycolipid",
  PC="phospholipid",
  PE="phospholipid",
  PG="phospholipid",
  PI="phospholipid",
  SQDG="glycolipid",
  TG="glycolipid")

side_chains <- c(
  DG = 2, 
  DGDG=2,
  DGGA=2,
  LPC=1,
  LPE=1,
  MGDG=2,
  PC=2,
  PE=2,
  PG=2,
  PI=2,
  SQDG=2,
  TG=3)

# - col_name: string for identifying the lipid in row and columns in R, NO COLONS
# - class: type of lipid (colname no number)
# - IUB_name: official IUB name with colons
# - head_group: type of head group glycolipid(inlcuding glycerol) or phospholipid
# - chain_n: number of side chains, how many faty acid side chains 1 2 3
# - total_C: total fatty acid  carbon atoms
# - total_FAU: total fatty acid unstaturations

lipid_meta <- data.frame(
  # col_name: string for identifying the lipid in row and columns in R, NO COLONS
  col_name = colnames(m)
) %>% 
  tidyr::separate(col_name, c("class", "n_C","n_desat"), remove=FALSE) %>%
  # IUB_name: official IUB name with colons
  mutate(IUB_name = paste0(class," ",n_C,":",n_desat)) %>%
  mutate(head_group = head_group[class]) %>%
  mutate(n_chains = side_chains[class]) %>%
  select(head_group,IUB_name,everything())

lipid_meta$n_C <- as.numeric(lipid_meta$n_C)
lipid_meta$n_desat<- as.numeric(lipid_meta$n_desat)

lipid_meta$desat_idx = lipid_meta$n_desat/lipid_meta$n_C

glyco <- lipid_meta$col_name[lipid_meta$head_group == "glycolipid"]
phospho <- lipid_meta$col_name[lipid_meta$head_group == "phospholipid"]


# IUB_name: official IUB name with colons
# head_group: type of head group glycolipid (inlcuding glycerol) or phospholipid
# chain_n: number of side chains, how many faty acid side chains 1 2 3
# n_C: total fatty acid  carbon atoms
# n_desat: total fatty acid unstaturations

# table of possible reactions
rx <-  read.csv("data/rx.csv", na.strings = c("","#N/A","NA","Inf"))



r1 <- lipid_meta %>%
  filter(class=="DGDG") %>%
  select(col_name,n_C) %>%
  full_join(
    lipid_meta %>%
    filter(class=="MGDG") %>%
      select(col_name,n_C),
    by = character()
  ) %>%
  filter(n_C.x== n_C.y) %>%
  select(x=col_name.x, y=col_name.y) %>%
  arrange(x,y) %>%
  distinct()
  
r2 <- lipid_meta %>%
  filter(class=="DG") %>%
  select(col_name,n_C) %>%
  full_join(
    lipid_meta %>%
      filter(class %in% c("DGGA", "SQDG", "MGDG", "PC", "PE")) %>%
      select(col_name,n_C),
    by = character()
  ) %>%
  filter(n_C.x== n_C.y) %>%
  select(x=col_name.x, y=col_name.y) %>%
  arrange(x,y) %>%
  distinct()


# destaturatiions

r3 <- lapply(lipid_meta$class, FUN = function(x){
  lipid_meta %>%
  filter(class==x) %>%
  select(col_name,n_C,n_desat) %>%
  full_join(
    lipid_meta %>%
      filter(class==x) %>%
      select(col_name,n_C,n_desat),
    by = character()
  ) %>%
  filter(n_C.x==n_C.y) %>%
  filter(abs(n_desat.x-n_desat.y)==1) %>%
  select(x=col_name.x, y=col_name.y)
}) %>% dplyr::bind_rows() %>%
  arrange() %>%
  distinct(x,y)



rx2 <- rbind(rx,r1,r2,r3) %>% as.matrix() %>% t() %>%
  apply(2, sort) %>% t() %>%
  as.data.frame() %>%
  rename(x="V1",y="V2") %>%
  arrange() %>%
  distinct() %>%
  filter(x %in% lipid_meta$col_name) %>%
  filter(y %in% lipid_meta$col_name)


write.csv(lipid_meta,"data/lipid_metadata.csv", quote=FALSE, row.names = FALSE)
write.csv(rx2,"data/rx2.csv", quote=FALSE, row.names = FALSE)
