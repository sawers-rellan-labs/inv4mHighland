#' Common Configuration for PSU 2022 Analysis Scripts
#' 
#' This file contains shared configuration settings and file paths
#' used across multiple PSU 2022 analysis scripts. Source this file
#' to ensure consistent paths and settings.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Directory paths --------------------------------------------------------
DATA_DIR <- "../../data"
OUTPUT_DIR <- "results_psu_2022" 
REF_DIR <- "~/ref/zea"

# Common data file paths -------------------------------------------------
PHENOTYPE_FILES <- list(
  plant = file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"),
  ear = file.path(DATA_DIR, "22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv"),
  expression_counts = file.path(DATA_DIR, "inv4mRNAseq_gene_sample_exp.csv"),
  expression_metadata = file.path(DATA_DIR, "PSU-PHO22_Metadata.csv"),
  ionome = file.path(DATA_DIR, "PSU_inv4m_ionome_all.csv"),
  lipids = file.path(DATA_DIR, "PSU_inv4m_lipids.csv"),
  gene_symbol = file.path(DATA_DIR, "gene_symbol.tab")
)

# Reference files --------------------------------------------------------
REFERENCE_FILES <- list(
  gff = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"),
  pannzer = "~/Desktop/PANNZER/PANNZER_DESC.tab"  # Optional external annotation
)

# Analysis parameters ----------------------------------------------------
ANALYSIS_PARAMS <- list(
  min_lib_size = 2e7,
  fdr_threshold = 0.05,
  logfc_threshold_general = 2,
  logfc_threshold_tissue = 0.7,
  consistency_threshold = 0.5
)

# Color palettes ---------------------------------------------------------
COLOR_PALETTES <- list(
  genotype = c("CTRL" = "gold", "Inv4m" = "#4a0f82"),
  treatment = c("+P" = "#2E8B57", "-P" = "#CD5C5C")
)

# Helper functions -------------------------------------------------------

#' Validate file paths and warn about missing files
#' @param file_list Named list of file paths
#' @param required Logical indicating if files are required (stops on missing)
validate_files <- function(file_list, required = FALSE) {
  missing_files <- file_list[!file.exists(unlist(file_list))]
  
  if (length(missing_files) > 0) {
    missing_names <- names(missing_files)
    message_text <- paste("Missing files:", paste(missing_names, collapse = ", "))
    
    if (required) {
      stop(message_text)
    } else {
      warning(message_text)
    }
  }
  
  return(file_list[file.exists(unlist(file_list))])
}

#' Create output directory if it doesn't exist
#' @param output_dir Output directory path
ensure_output_dir <- function(output_dir = OUTPUT_DIR) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  return(output_dir)
}

# Initialize output directory
ensure_output_dir()

cat("PSU 2022 common configuration loaded\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")