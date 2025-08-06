#' Validate Crow 2020 Genes Against Current Analysis
#'
#' This script validates genes identified in the Crow 2020 paper against
#' the current inv4m analysis results, focusing on conservation of effects
#' and identification of target genes for temperature response analysis.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load common configuration
source("common_config.R")

# Load Crow 2020 supplementary data
crow_file <- file.path(DATA_DIR, "Crow2020_table_S2.tab")
if (!file.exists(crow_file)) {
  stop("Crow 2020 supplementary table not found: ", crow_file)
}

crow <- read.table(file = crow_file, header = TRUE, sep = "\t")

# Load gene ID mapping if available
if (exists("v4_v5")) {
  crow <- crow %>%
    left_join(v4_v5, by = c(Gene = "v4"))
}

# Load functional annotation if available
if (exists("pannzer")) {
  crow <- crow %>%
    left_join(pannzer, by = c(v5 = "gene_model"))
}

# Analyze Crow 2020 gene list
crow_genes <- crow$Gene %>% sort() %>% unique()
cat("Total Crow 2020 genes:", length(crow_genes), "\n")

# Check conservation with current results (if available)
if (exists("effects")) {
  conserved <- crow %>%
    tibble() %>%
    filter(
      if ("v5" %in% names(.)) {
        v5 %in% effects$gene[effects$is_significant]
      } else {
        Gene %in% effects$gene[effects$is_significant]
      }
    )
  
  cat("Genes conserved in current analysis:", nrow(conserved), "\n")
  if (nrow(conserved) > 0) {
    print(conserved)
  }
} else {
  cat("Current analysis results not available for comparison\n")
}

# Identify target genes for temperature response analysis
target_analysis_genes <- crow_genes[crow_genes %in% unlist(TARGET_GENES)]

cat("\nTarget genes found in Crow 2020 data:", length(target_analysis_genes), "\n")
if (length(target_analysis_genes) > 0) {
  for (gene_name in names(TARGET_GENES)) {
    gene_ids <- TARGET_GENES[[gene_name]]
    found_genes <- intersect(gene_ids, crow_genes)
    if (length(found_genes) > 0) {
      cat(paste0(gene_name, ": "), paste(found_genes, collapse = ", "), "\n")
    }
  }
}

# Check genomic location of Crow genes (if genomic data available)
if (exists("genes") && "start" %in% names(genes)) {
  # Map Crow genes to genomic coordinates
  crow_with_coords <- crow
  if ("v5" %in% names(crow)) {
    crow_with_coords <- crow %>%
      left_join(
        genes %>% select(ID, seqnames, start, end),
        by = c("v5" = "ID")
      )
  }
  
  # Check if any are in inv4m region
  if (exists("inv4m_start") && exists("inv4m_end")) {
    inv4m_genes <- crow_with_coords %>%
      filter(
        seqnames == 4,
        start >= inv4m_start,
        end <= inv4m_end
      )
    
    cat("\nCrow genes in inv4m region:", nrow(inv4m_genes), "\n")
    if (nrow(inv4m_genes) > 0) {
      print(inv4m_genes %>% select(Gene, v5, start, end))
    }
  }
}

cat("\nValidation complete. Use these target genes in detect_crow2020_consistency_degs.R\n")
