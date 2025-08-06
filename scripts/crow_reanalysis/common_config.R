#' Common Configuration for Crow 2020 Reanalysis
#'
#' This file contains shared configuration settings and helper functions
#' for the Crow 2020 reanalysis experiment. The analysis focuses on 
#' validating consistency-focused approaches against the published mashr
#' results, specifically targeting PCNA2/MCM5 temperature responses.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load required libraries -----------------------------------------------
# Core libraries (required)
library(dplyr)         # Data manipulation
library(ggplot2)       # Plotting
library(readr)         # File reading
library(tidyr)         # Data reshaping

# Bioconductor libraries (optional for basic testing)
bioc_packages <- c("edgeR", "limma", "SRAdb", "GEOquery", "tximport", "rtracklayer", "GenomicRanges")
for (pkg in bioc_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(library(pkg, character.only = TRUE))
  } else {
    warning("Package ", pkg, " not available. Some functions may not work.")
  }
}

# Directory configuration ------------------------------------------------
DATA_DIR <- "../../data"
OUTPUT_DIR <- "results_crow_reanalysis"
SCRATCH_DIR <- "scratch"

# Reference genome configuration -----------------------------------------
REF_DIR <- "~/ref/zea"
REFERENCE_FILES <- list(
  gff = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"),
  fasta = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"),
  cdna = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa"),
  kallisto_index = file.path(REF_DIR, "Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx")
)

# SRA/GEO data configuration --------------------------------------------
SRA_PROJECT_ID <- "PRJNA640392"
GEO_SERIES_ID <- "GSE151676"

# Target genes of interest (from Crow 2020 and PSU 2022 analysis)
TARGET_GENES <- list(
  # Primary targets (cell proliferation markers)
  PCNA2 = c("Zm00001eb015030", "GRMZM2G085381"),  # v5, v4
  MCM5 = c("Zm00001eb025340", "GRMZM2G127379"),   # v5, v4
  
  # Secondary targets (flowering/growth)
  ZmCCT = c("Zm00001eb194800", "GRMZM2G410515"),  # v5, v4
  ZmCCT9 = c("Zm00001eb190470", "GRMZM2G522046"), # v5, v4
  
  # Controls (housekeeping)
  GAPDH = c("Zm00001eb077550", "GRMZM2G155329"),  # v5, v4
  UBQ = c("Zm00001eb071870", "GRMZM2G056659")     # v5, v4
)

# Analysis parameters ----------------------------------------------------
ANALYSIS_SETTINGS <- list(
  min_lib_size = 1e6,           # Minimum library size threshold
  min_cpm_threshold = 1,        # CPM filtering threshold
  min_count_samples = 3,        # Minimum samples with counts
  fdr_threshold = 0.05,         # FDR significance threshold
  logfc_threshold = 1.0,        # LogFC threshold for DEGs
  bootstrap_samples = 100,      # Kallisto bootstrap samples
  threads = 8                   # Processing threads
)

# Tissue and condition definitions ---------------------------------------
TISSUE_GROUPS <- list(
  proliferating = c("SAM", "root", "leaf_base"),
  mature = c("leaf_tip", "silk", "tassel"),
  mixed = c("whole_seedling", "ear", "node")
)

TEMPERATURE_CONDITIONS <- list(
  cold = c("16C", "18C", "20C"),
  optimal = c("22C", "24C", "25C"), 
  warm = c("28C", "30C", "32C")
)

GENOTYPE_GROUPS <- list(
  NIL_lines = c("PT_NIL", "Mi21_NIL"),
  control = c("B73", "Mo17"),
  hybrid = c("F1", "B73xMo17")
)

# Plotting configuration -------------------------------------------------
PLOT_SETTINGS <- list(
  width = 12,
  height = 8, 
  dpi = 300,
  base_size = 14
)

# Color schemes
COLOR_SCHEMES <- list(
  genotype = c(
    "CTRL" = "#440154FF",
    "INV4" = "#FDE725FF"
  ),
  temperature = c(
    "cold" = "#3B9AB2",
    "optimal" = "#78B7C5", 
    "warm" = "#EBCC2A"
  ),
  tissue = viridis::viridis(length(unique(unlist(TISSUE_GROUPS))))
)

# Helper functions -------------------------------------------------------

#' Create output directories if they don't exist
#'
#' @param dirs Vector of directory paths to create
create_output_dirs <- function(dirs = NULL) {
  if (is.null(dirs)) {
    dirs <- c(
      OUTPUT_DIR,
      file.path(OUTPUT_DIR, "plots"),
      file.path(OUTPUT_DIR, "tables"), 
      file.path(OUTPUT_DIR, "qc"),
      file.path(OUTPUT_DIR, "kallisto"),
      SCRATCH_DIR
    )
  }
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

#' Load gene annotation from GFF file
#'
#' @param gff_path Path to GFF file
#' @return Data frame with gene annotation
load_gene_annotation <- function(gff_path = REFERENCE_FILES$gff) {
  if (!file.exists(gff_path)) {
    warning("GFF file not found: ", gff_path)
    return(data.frame(
      gene_id = character(),
      seqnames = character(), 
      start = numeric(),
      end = numeric(),
      strand = character()
    ))
  }
  
  # Load GFF using rtracklayer if available
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gff <- rtracklayer::import(gff_path) %>%
      subset(type == "gene" & seqnames %in% 1:10)
    
    genes <- as.data.frame(gff) %>%
      mutate(gene_id = gsub("gene:", "", ID)) %>%
      select(gene_id, seqnames, start, end, strand, Name)
    
  } else {
    warning("rtracklayer not available, using basic gene annotation")
    genes <- data.frame(
      gene_id = names(TARGET_GENES),
      seqnames = rep(4, length(TARGET_GENES)),
      start = rep(180000000, length(TARGET_GENES)),
      end = rep(190000000, length(TARGET_GENES)),
      strand = rep("+", length(TARGET_GENES))
    )
  }
  
  return(genes)
}

#' Load SRA metadata from cached file or download
#'
#' @param cache_file Path to cached metadata file
#' @return Data frame with sample metadata
load_sra_metadata <- function(cache_file = file.path(DATA_DIR, "SraRunInfo.csv")) {
  if (file.exists(cache_file)) {
    cat("Loading cached SRA metadata from:", cache_file, "\n")
    return(read.csv(cache_file))
  } else {
    cat("SRA metadata file not found:", cache_file, "\n")
    cat("Please run check_sra_metadata.R to download and process metadata\n")
    return(data.frame())
  }
}

#' Parse sample information from SRA metadata
#'
#' @param sra_meta SRA metadata data frame
#' @return Parsed sample information
parse_sample_info <- function(sra_meta) {
  if (nrow(sra_meta) == 0) {
    return(data.frame())
  }
  
  # Extract sample information from library names
  sample_info <- sra_meta %>%
    tidyr::separate_wider_delim(
      cols = "LibraryName",
      names = c("base", "parent", "sample_id"),
      delim = "_",
      cols_remove = FALSE
    ) %>%
    mutate(
      genotype = case_when(
        grepl("NIL", parent) ~ "INV4",
        TRUE ~ "CTRL"
      ),
      tissue = case_when(
        grepl("SAM", LibraryName, ignore.case = TRUE) ~ "SAM",
        grepl("root", LibraryName, ignore.case = TRUE) ~ "root", 
        grepl("leaf", LibraryName, ignore.case = TRUE) ~ "leaf",
        TRUE ~ "unknown"
      )
    )
  
  return(sample_info)
}

#' Validate file existence
#'
#' @param files Vector of file paths to check
#' @param required Whether files are required (stop if missing)
validate_files <- function(files, required = TRUE) {
  missing <- files[!file.exists(files)]
  
  if (length(missing) > 0) {
    msg <- paste("Missing files:", paste(missing, collapse = ", "))
    
    if (required) {
      stop(msg)
    } else {
      warning(msg)
      return(files[file.exists(files)])
    }
  }
  
  return(files)
}

#' Load Kallisto quantification results
#'
#' @param sample_dirs Vector of sample directories containing abundance.tsv
#' @param tx2gene Transcript to gene mapping
#' @return Tximport object with gene-level counts
load_kallisto_results <- function(sample_dirs, tx2gene = NULL) {
  # Find abundance files
  abundance_files <- file.path(sample_dirs, "abundance.tsv")
  names(abundance_files) <- basename(sample_dirs)
  
  # Validate files exist
  abundance_files <- validate_files(abundance_files, required = FALSE)
  
  if (length(abundance_files) == 0) {
    stop("No Kallisto abundance files found")
  }
  
  cat("Loading", length(abundance_files), "Kallisto quantification files\n")
  
  # Import with tximport
  if (requireNamespace("tximport", quietly = TRUE)) {
    txi <- tximport::tximport(
      abundance_files,
      type = "kallisto",
      tx2gene = tx2gene,
      ignoreAfterBar = TRUE
    )
    
    cat("Loaded quantifications for", nrow(txi$counts), "genes\n")
    return(txi)
    
  } else {
    warning("tximport not available, loading abundance files manually")
    
    # Manual loading (basic implementation)
    abundance_list <- lapply(abundance_files, function(f) {
      read.table(f, header = TRUE, sep = "\t")
    })
    
    # Simple aggregation at transcript level (not ideal)
    all_transcripts <- unique(unlist(lapply(abundance_list, function(x) x$target_id)))
    
    counts_matrix <- matrix(0, nrow = length(all_transcripts), ncol = length(abundance_list))
    rownames(counts_matrix) <- all_transcripts
    colnames(counts_matrix) <- names(abundance_list)
    
    for (i in seq_along(abundance_list)) {
      sample_data <- abundance_list[[i]]
      counts_matrix[sample_data$target_id, i] <- sample_data$est_counts
    }
    
    return(list(
      counts = counts_matrix,
      abundance = counts_matrix,  # Simplified
      length = matrix(1000, nrow = nrow(counts_matrix), ncol = ncol(counts_matrix))
    ))
  }
}

# Initialize configuration -----------------------------------------------

# Create output directories
create_output_dirs()

# Load gene annotation
cat("Loading gene annotation...\n")
genes <- load_gene_annotation()

# Configuration summary
cat("Crow 2020 reanalysis configuration loaded\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n") 
cat("Reference directory:", REF_DIR, "\n")
cat("Target genes:", length(TARGET_GENES), "\n")
cat("Gene annotation:", nrow(genes), "genes loaded\n")

# Export key objects to global environment
assign("genes", genes, envir = .GlobalEnv)
assign("OUTPUT_DIR", OUTPUT_DIR, envir = .GlobalEnv)
assign("DATA_DIR", DATA_DIR, envir = .GlobalEnv)