#' Analyze Kallisto RNA-seq Quantification Statistics
#'
#' This script analyzes RNA-seq quantification statistics from Kallisto,
#' integrates sample metadata, and generates quality control reports.
#' It processes tube information, sample metadata, and plot assignments
#' to create comprehensive sample tracking.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load required libraries -----------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readr)

# Source common configuration
source("../common_config.R")

# Helper functions -------------------------------------------------------

#' Load and validate RNA-seq metadata files
#'
#' @param tube_file Path to tube container file
#' @param sample_file Path to sample information file  
#' @param plots_file Path to plots information file
#' @return List containing loaded data frames
load_metadata_files <- function(tube_file, sample_file, plots_file) {
  # Validate files exist
  files <- c(tube = tube_file, sample = sample_file, plots = plots_file)
  missing <- files[!file.exists(files)]
  
  if (length(missing) > 0) {
    stop("Missing metadata files: ", paste(names(missing), collapse = ", "))
  }
  
  # Load data with error handling
  tube_data <- tryCatch({
    read.table(tube_file, sep = "\t", quote = "", header = TRUE)
  }, error = function(e) {
    stop("Failed to load tube data: ", e$message)
  })
  
  sample_data <- tryCatch({
    read.csv(sample_file, quote = "", header = TRUE)
  }, error = function(e) {
    stop("Failed to load sample data: ", e$message)
  })
  
  plots_data <- tryCatch({
    read.csv(plots_file)
  }, error = function(e) {
    stop("Failed to load plots data: ", e$message)
  })
  
  return(list(
    tube = tube_data,
    sample = sample_data, 
    plots = plots_data
  ))
}

#' Integrate metadata from multiple sources
#'
#' @param metadata_list List containing tube, sample, and plots data
#' @return Integrated metadata data frame
integrate_metadata <- function(metadata_list) {
  # Join metadata sources
  integrated <- metadata_list$tube %>%
    inner_join(metadata_list$sample, by = c(tube = "top_tag")) %>%
    inner_join(metadata_list$plots, by = c(row = "PHO22"))

  
  # Apply known corrections for swapped samples
  integrated <- apply_sample_corrections(integrated)
  
  # Process time information
  integrated <- process_time_data(integrated)
  
  return(integrated)
}

#' Apply known sample corrections
#'
#' @param data Metadata data frame
#' @return Corrected metadata data frame
apply_sample_corrections <- function(data) {
  # Correct R33 and R34 - they have the leaves swapped
  data$leaf_tissue[data$tube == 'R33'] <- 2
  data$leaf_tissue[data$tube == 'R34'] <- 1
  
  cat("Applied sample corrections for R33 and R34\n")
  return(data)
}

#' Process time data and create decimal time
#'
#' @param data Metadata data frame
#' @return Data frame with processed time information
process_time_data <- function(data) {
  data <- within(data, {
    # Fix time format (convert 1: to 13: for PM times)
    TIME <- sub("1:", "13:", TIME)
    TIME <- hm(TIME)
    decimal_time <- hour(TIME) + minute(TIME)/60 + second(TIME)/3600
  })
  
  return(data)
}


#' Load and process Kallisto QC statistics
#'
#' @param qc_file Path to Kallisto QC statistics file
#' @return Processed QC statistics data frame
load_kallisto_qc <- function(qc_file) {
  if (!file.exists(qc_file)) {
    stop("Kallisto QC file not found: ", qc_file)
  }
  
  # Load QC data
  qc_data <- read.table(
    qc_file, 
    sep = "\t", 
    quote = "", 
    header = FALSE, 
    skip = 1
  )
  
  # Set column names
  colnames(qc_data) <- c(
    "gsl_sample", "pseudoaligned_pct", 
    "pseudoaligned", "processed"
  )
  
  # Extract tube information
  qc_data$tube <- substr(qc_data$gsl_sample, 1, 3)
  
  return(qc_data)
}

#' Generate QC summary statistics
#'
#' @param qc_data Kallisto QC data frame
#' @return Summary statistics data frame
generate_qc_summary <- function(qc_data) {
  summary_stats <- qc_data %>%
    summarise(
      total_samples = n(),
      mean_pseudoaligned_pct = mean(pseudoaligned_pct, na.rm = TRUE),
      median_pseudoaligned_pct = median(pseudoaligned_pct, na.rm = TRUE),
      min_pseudoaligned_pct = min(pseudoaligned_pct, na.rm = TRUE),
      max_pseudoaligned_pct = max(pseudoaligned_pct, na.rm = TRUE),
      mean_processed = mean(processed, na.rm = TRUE),
      low_quality_samples = sum(pseudoaligned_pct < 70, na.rm = TRUE)
    )
  
  return(summary_stats)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to analyze Kallisto statistics
#'
#' @param tube_file Path to tube container file
#' @param sample_file Path to sample information file
#' @param plots_file Path to plots file
#' @param qc_file Path to Kallisto QC file
#' @export
analyze_kallisto_stats <- function(tube_file = NULL, sample_file = NULL, 
                                  plots_file = NULL, qc_file = NULL) {
  
  cat("Starting Kallisto statistics analysis...\n")
  
  # Set default file paths
  if (is.null(tube_file)) {
    tube_file <- file.path(DATA_DIR, "tuberna_ctn.tab")
  }
  if (is.null(sample_file)) {
    sample_file <- file.path(DATA_DIR, "PENN-PHO22-TUBES.csv")
  }
  if (is.null(plots_file)) {
    plots_file <- file.path(DATA_DIR, "PENN_PHO22_PLOTS.csv")
  }
  if (is.null(qc_file)) {
    qc_file <- file.path(OUTPUT_DIR, "kallisto_qc.tab")
  }
  
  # Load and integrate metadata
  cat("Loading metadata files...\n")
  metadata_list <- load_metadata_files(tube_file, sample_file, plots_file)
  
  cat("Integrating metadata...\n")
  metadata <- integrate_metadata(metadata_list)
  
  # Save integrated metadata
  output_metadata <- metadata %>%
    select(
      tube:NCSU_RNA_plant,
      Treatment, genotype, leaf_tissue, side_tag,
      TIME, decimal_time, everything()
    )
  
  write.csv(
    output_metadata,
    file = file.path(DATA_DIR, "inv4mRNAseq_metadata.csv"),
    row.names = FALSE
  )
  
  cat("Saved integrated metadata to inv4mRNAseq_metadata.csv\n")
  
  # Load and process Kallisto QC data
  if (file.exists(qc_file)) {
    cat("Loading Kallisto QC statistics...\n")
    qc_data <- load_kallisto_qc(qc_file)
    
    # Generate summary statistics
    qc_summary <- generate_qc_summary(qc_data)
    
    cat("QC Summary:\n")
    print(qc_summary)
    
    # Save QC results
    write.csv(
      qc_data,
      file = file.path(OUTPUT_DIR, "kallisto_qc_processed.csv"),
      row.names = FALSE
    )
    
    write.csv(
      qc_summary,
      file = file.path(OUTPUT_DIR, "kallisto_qc_summary.csv"),
      row.names = FALSE
    )
  } else {
    warning("Kallisto QC file not found: ", qc_file)
    qc_data <- NULL
    qc_summary <- NULL
  }
  
  cat("Analysis complete!\n")
  
  return(list(
    metadata = metadata,
    qc_data = qc_data,
    qc_summary = qc_summary
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- analyze_kallisto_stats()
}

# check leaf number
ggplot(sampleInfo, aes(y=leaf_number, x =genotype)) + 
  ggbeeswarm::geom_quasirandom() + facet_wrap(~Treatment)
anova(lm(data=sampleInfo,leaf_number ~ Treatment*genotype))


#sum over genes
# make gene _sample table

# if you are thinikg in adding readcounts from poor qc libraries together
# sampleInfo[sampleInfo$genotype=="CTRL" & sampleInfo$Treatment=="Low_P" & sampleInfo$leaf_tissue ==1, c("tube","row","Rep")]

# this is for runing in the server

samples <- dir("./quant_out")

all_exp <-lapply(samples, function(x){
  print(x)
  sample_exp <- read.table(file.path("quant_out", x,"abundance.tsv"), sep = "\t", quote ="", header = TRUE)
  sample_exp$est_counts <-as.integer(sample_exp$est_counts)
  sample_exp$gene = factor(sub("_T.+$","",sample_exp$target_id, perl =TRUE))

  out <- sample_exp %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(counts = sum(est_counts))
  out$sample <- x
  out
}) %>% dplyr::bind_rows()



gene_sample_exp <- all_exp %>% pivot_wider( names_from = "sample", values_from ="counts" )


write.csv(gene_sample_exp,"gene_sample_exp.csv")











kalqc %>% filter(pseudoaligned_pct < 60) 
  
quartz()
kalqc %>%
  tidyr::pivot_longer(pseudoaligned_pct:processed, names_to = "var", values_to = "value" ) %>%
  ggplot2::ggplot(aes(x=value))+
  ggtitle("kallisto pseudoalignment") +
  ylab("# tissue libraries") +
  xlab("read stat") +
  geom_histogram()+
  facet_wrap( ~ var, scales = "free")+
  ggpubr::theme_classic2()
