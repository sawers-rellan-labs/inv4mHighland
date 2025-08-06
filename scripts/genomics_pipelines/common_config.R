#' Common Configuration for Genomics Pipelines
#' 
#' This file contains shared configuration settings and file paths
#' used across genomics pipeline R scripts. Source this file
#' to ensure consistent paths and settings.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Directory paths --------------------------------------------------------
DATA_DIR <- "../../data"
OUTPUT_DIR <- "results_genomics_pipelines"
REF_DIR <- "~/ref/zea"
SCRATCH_DIR <- "~/scratch/genomics"

# Reference genome files -------------------------------------------------
REFERENCE_FILES <- list(
  gff = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"),
  fasta = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"),
  annotation = file.path(REF_DIR, "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.gtf")
)

# Common data files ------------------------------------------------------
DATA_FILES <- list(
  metadata = file.path(DATA_DIR, "inv4mRNAseq_metadata.csv"),
  sample_info = file.path(DATA_DIR, "PSU-PHO22_Metadata.csv"),
  plots = file.path(DATA_DIR, "PSU_PHO22_PLOTS.csv"),
  tube_info = file.path(DATA_DIR, "tuberna_ctn.tab")
)

# Genomic coordinates (from MCScan analysis) ----------------------------
GENOMIC_REGIONS <- list(
  inv4m = list(
    chr = "4",
    start_gene = "Zm00001eb190470",  # Will be resolved to coordinate
    end_gene = "Zm00001eb194800",    # Will be resolved to coordinate
    start = 172883675,               # Known coordinates
    end = 188132113
  ),
  introgression = list(
    chr = "4", 
    start = 157012149,
    end = 195900523
  )
)

# Analysis parameters ----------------------------------------------------
ANALYSIS_PARAMS <- list(
  min_mapping_quality = 20,
  min_coverage = 10,
  max_missing = 0.2,
  min_allele_freq = 0.05
)

# Visualization settings -------------------------------------------------
PLOT_SETTINGS <- list(
  width = 10,
  height = 8,
  dpi = 300,
  base_size = 12
)

# Color palettes ---------------------------------------------------------
COLOR_PALETTES <- list(
  genotype = c("B73" = "#1f77b4", "TIL18" = "#ff7f0e", "PT" = "#2ca02c"),
  quality = c("High" = "#2ca02c", "Medium" = "#ff7f0e", "Low" = "#d62728"),
  chromosome = RColorBrewer::brewer.pal(10, "Set3")
)

# Helper functions -------------------------------------------------------

#' Validate file paths and create missing directories
#' @param file_list Named list of file paths
#' @param create_dirs Logical, create missing directories for output files
validate_files <- function(file_list, create_dirs = TRUE) {
  missing_files <- file_list[!file.exists(unlist(file_list))]
  
  if (length(missing_files) > 0) {
    warning("Missing files: ", paste(names(missing_files), collapse = ", "))
  }
  
  if (create_dirs) {
    # Create output directories
    output_dirs <- unique(dirname(unlist(file_list)))
    for (dir in output_dirs) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        cat("Created directory:", dir, "\n")
      }
    }
  }
  
  return(file_list[file.exists(unlist(file_list))])
}

#' Load gene annotation and resolve genomic coordinates
#' @param gff_path Path to GFF3 file
#' @return List with gene coordinates and annotation
load_gene_annotation <- function(gff_path = REFERENCE_FILES$gff) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer package required for gene annotation")
  }
  
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package required for gene annotation")
  }
  
  if (!file.exists(gff_path)) {
    stop("GFF file not found: ", gff_path)
  }
  
  # Load gene annotation
  gff <- rtracklayer::import(gff_path) %>%
    subset(type == "gene" & seqnames %in% 1:10)
  
  genes <- as.data.frame(gff)
  genes$ID <- gsub("gene:", "", genes$ID)
  
  # Resolve genomic region coordinates
  inv4m_start <- genes[genes$ID == GENOMIC_REGIONS$inv4m$start_gene, "start"]
  inv4m_end <- genes[genes$ID == GENOMIC_REGIONS$inv4m$end_gene, "end"]
  
  # Create GRanges objects
  inv4m_region <- GenomicRanges::GRanges(
    seqnames = GENOMIC_REGIONS$inv4m$chr,
    ranges = GenomicRanges::IRanges(start = inv4m_start, end = inv4m_end),
    strand = "+"
  )
  
  introgression_region <- GenomicRanges::GRanges(
    seqnames = GENOMIC_REGIONS$introgression$chr,
    ranges = GenomicRanges::IRanges(
      start = GENOMIC_REGIONS$introgression$start, 
      end = GENOMIC_REGIONS$introgression$end
    ),
    strand = "+"
  )
  
  # Find overlapping genes
  inv4m_overlap <- GenomicRanges::findOverlaps(inv4m_region, gff, ignore.strand = TRUE)
  inv4m_genes <- genes$ID[S4Vectors::subjectHits(inv4m_overlap)]
  
  introgression_overlap <- GenomicRanges::findOverlaps(introgression_region, gff, ignore.strand = TRUE)
  introgression_genes <- genes$ID[S4Vectors::subjectHits(introgression_overlap)]
  
  flanking_genes <- setdiff(introgression_genes, inv4m_genes)
  
  return(list(
    genes = genes,
    gff = gff,
    regions = list(
      inv4m = inv4m_region,
      introgression = introgression_region
    ),
    gene_sets = list(
      inv4m = inv4m_genes,
      introgression = introgression_genes,
      flanking = flanking_genes
    ),
    coordinates = list(
      inv4m_start = inv4m_start,
      inv4m_end = inv4m_end
    )
  ))
}

#' Create standard output directory structure
ensure_output_dirs <- function() {
  dirs <- c(
    OUTPUT_DIR,
    file.path(OUTPUT_DIR, "plots"),
    file.path(OUTPUT_DIR, "tables"),
    file.path(OUTPUT_DIR, "qc")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
  
  return(dirs)
}

# Initialize directories
ensure_output_dirs()

cat("Genomics pipelines configuration loaded\n")
cat("Data directory:", DATA_DIR, "\n") 
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Reference directory:", REF_DIR, "\n")