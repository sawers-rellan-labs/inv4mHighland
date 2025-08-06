#' Check and Process SRA Metadata for Crow 2020 Experiment
#'
#' This script processes SRA metadata to identify samples suitable for
#' differential expression analysis, focusing on proliferating tissues
#' (SAM and root) from NIL lines.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load common configuration
source("common_config.R")

# Load additional required libraries
library(tidyr)

# Load SRA metadata
meta1 <- read.csv(file = file.path(DATA_DIR, "SraRunInfo.csv")) %>%
  tidyr::separate_wider_delim(
    cols = "LibraryName",
    names = c("base", "parent", "ID"), 
    delim = "_",
    cols_remove = FALSE
  )


with(meta1,
     table(parent)) %>%sum()
with(meta1,
     table(parent,base))

meta2 <- read.csv(file = file.path(DATA_DIR, "gc7_sample_submission.csv"))

with(meta2,
     table(old_line,full_geno))
meta2$old_line %>% factor %>% levels()

table(meta1$LibraryName)

(factor(meta2$samp) %>% levels() )%in% (factor(meta1$LibraryName) %>% levels())

nrow(meta1)

nrow(meta2)

meta2$old_line


meta2 %>%
  rename(LibraryName="samp") %>%
  left_join(
    meta1 %>% select(Run, LibraryName)) %>%
  select(old_line, LibraryName,Run)



with(meta2,
     table(old_line)) 


meta2$
with(meta2 %>% filter(old_line %in% c("PT_NIL","Mi21_NIL") & grepl("SAM|root",tissue)),
     table(old_line,tissue)) 

apical_tissues <- meta2 %>% 
  filter(old_line %in% c("PT_NIL","Mi21_NIL") & grepl("SAM|root",tissue)) %>%
  select(old_line,inv_temp,tissue,samp) %>%
  tidyr::separate_wider_delim(inv_temp, "_", names = c("inv_gt", "temp")) %>%
inner_join(                 
meta1 %>%
  select(samp=LibraryName,Run)
)

with( apical_tissues %>%
  filter(old_line =="PT_NIL"),
  table(inv_gt,temp,tissue)
)

# Save processed metadata to data directory
output_file <- file.path(DATA_DIR, "crow2020_apical_tissue_samples.tab")
write.table(
  file = output_file, 
  apical_tissues,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Processed sample metadata saved to:", output_file, "\n")
cat("Total samples for analysis:", nrow(apical_tissues), "\n")
cat("Genotype distribution:\n")
print(table(apical_tissues$old_line))
cat("\nTissue distribution:\n")
print(table(apical_tissues$tissue))
cat("\nTemperature distribution:\n")
print(table(apical_tissues$inv_gt))



with(meta2 %>% filter(old_line %in% c("PT_NIL","Mi21_NIL") & grepl("SAM|root",tissue)),
     table()) 
with(meta2 %>% filter(old_line=="F1" & grepl("Leaf|leaf",tissue)),
     table(old_line,tissue))

with(meta2 %>% filter(old_line=="F1" & grepl("Leaf|leaf",tissue)),
     table(old_line,tissue))

with(meta2 %>%
     table(old_line,tissue))

with(meta2,
     table(parent,base))
