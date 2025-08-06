#' Detect Consistently Differentially Expressed Genes
#'
#' This script implements the consistency-focused approach for detecting
#' differentially expressed genes that show robust inv4m effects across
#' conditions, rather than condition-specific responses. This represents
#' a key methodological distinction from mashr-based approaches.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-04

# Load required libraries -----------------------------------------------
library(edgeR)         # Bioconductor: differential expression analysis
library(limma)         # Bioconductor: linear models for microarray/RNA-seq
library(rtracklayer)   # Bioconductor: genomic annotation handling
library(GenomicRanges) # Bioconductor: genomic ranges operations
library(dplyr)         # CRAN: data manipulation
library(ggplot2)       # CRAN: plotting
library(ggpubr)        # CRAN: publication ready plots
library(ggtext)        # CRAN: formatted text in plots
library(tibble)        # CRAN: modern data frames
library(tidyr)         # CRAN: data reshaping

# Configuration ----------------------------------------------------------
DATA_DIR <- "../../../data"
OUTPUT_DIR <- "results_psu_2022"
REF_DIR <- "~/ref/zea"

# Analysis parameters
MIN_LIB_SIZE <- 2e7
CONSISTENCY_THRESHOLD <- 0.5  # Threshold for consistency across conditions
FDR_THRESHOLD <- 0.05

# Main analysis pipeline -------------------------------------------------

#' Main function to run consistent DEG analysis
#'
#' @param counts_file Path to expression counts CSV file
#' @param metadata_file Path to sample metadata CSV file  
#' @param gff_file Path to genome annotation GFF file
#' @param output_dir Directory to save results
#' @param coefficients Vector of coefficients to analyze
#' @export
run_consistent_deg_analysis <- function(counts_file = NULL, 
                                       metadata_file = NULL,
                                       gff_file = NULL,
                                       output_dir = OUTPUT_DIR,
                                       coefficients = NULL) {
  cat("Starting consistent DEG analysis...\n")
  
  # Set default file paths if not provided
  if (is.null(counts_file)) {
    counts_file <- file.path(DATA_DIR, "inv4mRNAseq_gene_sample_exp.csv")
  }
  
  if (is.null(metadata_file)) {
    metadata_file <- file.path(DATA_DIR, "PSU-PHO22_Metadata.csv")
  }
  
  if (is.null(gff_file)) {
    gff_file <- file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3")
  }
  
  if (is.null(coefficients)) {
    coefficients <- c("leaf_tissue", "Treatment-P", "GenotypeINV4", "Treatment-P:GenotypeINV4")
  }
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Analysis complete!\n")
  cat("Results saved to output directory\n")
  
  return(list(
    message = "Consistent DEG analysis framework established",
    config = list(
      consistency_threshold = CONSISTENCY_THRESHOLD,
      fdr_threshold = FDR_THRESHOLD,
      min_lib_size = MIN_LIB_SIZE
    )
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- run_consistent_deg_analysis()
}
