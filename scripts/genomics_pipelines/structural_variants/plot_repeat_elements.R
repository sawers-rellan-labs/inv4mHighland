#' Plot Repeat Elements Across Genomes
#'
#' This script analyzes BLAST results for repeat elements (knobs and TR-1)
#' across different maize genomes and creates visualizations showing their
#' distribution and similarity scores along chromosomes.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load required libraries -----------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
library(readr)

# Source common configuration
source("../common_config.R")

# Analysis configuration -------------------------------------------------
GENOMES <- c("TIL18", "PT", "B73")
BLAST_COLUMNS <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

# Default blast file locations
DEFAULT_BLAST_FILES <- list(
  knob = c(
    TIL18 = file.path(SCRATCH_DIR, "knob_til18.blast"),
    PT = file.path(SCRATCH_DIR, "knob_pt.blast"),
    B73 = file.path(SCRATCH_DIR, "knob_b73.blast")
  ),
  tr1 = c(
    TIL18 = file.path(SCRATCH_DIR, "TR-1_til18.blast"),
    PT = file.path(SCRATCH_DIR, "TR-1_pt.blast"),
    B73 = file.path(SCRATCH_DIR, "TR-1_b73.blast")
  )
)

# Helper functions -------------------------------------------------------

#' Load BLAST results for a repeat type
#'
#' @param blast_files Named vector of BLAST file paths
#' @param repeat_name Name of the repeat element
#' @return Data frame with combined BLAST results
load_blast_results <- function(blast_files, repeat_name) {
  
  # Validate files exist
  missing_files <- blast_files[!file.exists(blast_files)]
  if (length(missing_files) > 0) {
    warning("Missing BLAST files for ", repeat_name, ": ", 
            paste(names(missing_files), collapse = ", "))
    blast_files <- blast_files[file.exists(blast_files)]
  }
  
  if (length(blast_files) == 0) {
    stop("No BLAST files found for ", repeat_name)
  }
  
  # Load and combine results
  results <- lapply(names(blast_files), function(genome) {
    cat("Loading", repeat_name, "results for", genome, "...\n")
    
    blast_data <- tryCatch({
      read.table(
        blast_files[genome], 
        sep = "\t", 
        header = FALSE,
        col.names = BLAST_COLUMNS,
        colClasses = c(qseqid = "character", sseqid = "character")
      )
    }, error = function(e) {
      warning("Failed to load BLAST file for ", genome, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(blast_data)) {
      blast_data$queryg <- genome
      blast_data$repeat_type <- repeat_name
      return(blast_data)
    }
    return(NULL)
  })
  
  # Remove NULL results and combine
  results <- results[!sapply(results, is.null)]
  
  if (length(results) > 0) {
    combined <- bind_rows(results) %>%
      mutate(
        bitscaled = scales::rescale(bitscore, to = c(0, 100)),
        queryg = factor(queryg, levels = GENOMES)
      )
    
    cat("Loaded", nrow(combined), "BLAST hits for", repeat_name, "\n")
    return(combined)
  } else {
    return(NULL)
  }
}

#' Create repeat element visualization
#'
#' @param repeat_data Combined repeat element data
#' @param inversion_boundaries Data frame with inversion boundaries
#' @return ggplot object
create_repeat_plot <- function(repeat_data, inversion_boundaries = NULL) {
  
  # Create base plot
  p <- repeat_data %>%
    ggplot(aes(x = qstart, y = bitscaled, color = repeat_type)) +
    geom_point(alpha = 0.7) +
    facet_wrap(
      ~queryg, 
      nrow = length(GENOMES), 
      ncol = 1,
      scales = "free_x"
    ) +
    scale_color_manual(
      values = c("knob180" = "#1d7f7a", "TL-R1" = "gold"),
      name = "Repeat Type"
    ) +
    scale_x_continuous(
      labels = unit_format(unit = "Mb", scale = 1e-6),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(
      x = "Genomic Position",
      y = "Normalized Repeat Match Score",
      title = "Repeat Element Distribution Across Genomes"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
  
  # Add inversion boundaries if provided
  if (!is.null(inversion_boundaries)) {
    p <- p + 
      geom_vline(
        data = inversion_boundaries,
        aes(xintercept = position),
        color = "#3232a0",
        linewidth = 1.2,
        alpha = 0.7
      )
  }
  
  return(p)
}

# Main analysis pipeline -------------------------------------------------

#' Main function to analyze and plot repeat elements
#'
#' @param knob_files Named vector of knob BLAST files
#' @param tr1_files Named vector of TR-1 BLAST files
#' @param output_prefix Prefix for output files
#' @export
plot_repeat_elements <- function(knob_files = NULL, tr1_files = NULL,
                                output_prefix = "repeat_elements") {
  
  cat("Starting repeat elements analysis...\n")
  
  # Use default files if not provided
  if (is.null(knob_files)) {
    knob_files <- DEFAULT_BLAST_FILES$knob
  }
  if (is.null(tr1_files)) {
    tr1_files <- DEFAULT_BLAST_FILES$tr1
  }
  
  # Load BLAST results
  knob_data <- load_blast_results(knob_files, "knob180")
  tr1_data <- load_blast_results(tr1_files, "TL-R1")
  
  # Combine data
  combined_data <- bind_rows(knob_data, tr1_data)
  
  if (nrow(combined_data) == 0) {
    stop("No repeat element data loaded")
  }
  
  # Create inversion boundary annotations
  inversion_boundaries <- data.frame(
    queryg = rep(GENOMES, each = 2),
    position = c(
      # TIL18
      180365316, 193570651,
      # PT  
      173486186, 186925483,
      # B73
      172882309, 188131461
    ),
    boundary = rep(c("start", "end"), length(GENOMES))
  )
  
  # Create visualization
  cat("Creating repeat element plot...\n")
  repeat_plot <- create_repeat_plot(combined_data, inversion_boundaries)
  
  # Save plot
  output_file <- file.path(OUTPUT_DIR, "plots", paste0(output_prefix, ".pdf"))
  ggsave(
    output_file,
    repeat_plot,
    width = PLOT_SETTINGS$width,
    height = PLOT_SETTINGS$height * 1.5,  # Taller for multiple panels
    dpi = PLOT_SETTINGS$dpi
  )
  
  cat("Plot saved to:", output_file, "\n")
  
  # Save data
  write.csv(
    combined_data,
    file = file.path(OUTPUT_DIR, "tables", paste0(output_prefix, "_data.csv")),
    row.names = FALSE
  )
  
  write.csv(
    inversion_boundaries,
    file = file.path(OUTPUT_DIR, "tables", paste0(output_prefix, "_boundaries.csv")),
    row.names = FALSE
  )
  
  # Generate summary statistics
  summary_stats <- combined_data %>%
    group_by(queryg, repeat_type) %>%
    summarise(
      total_hits = n(),
      mean_bitscore = mean(bitscore, na.rm = TRUE),
      mean_identity = mean(pident, na.rm = TRUE),
      mean_length = mean(length, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(
    summary_stats,
    file = file.path(OUTPUT_DIR, "tables", paste0(output_prefix, "_summary.csv")),
    row.names = FALSE
  )
  
  cat("Analysis complete!\n")
  cat("Summary statistics:\n")
  print(summary_stats)
  
  return(list(
    data = combined_data,
    plot = repeat_plot,
    summary = summary_stats,
    boundaries = inversion_boundaries
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- plot_repeat_elements()
}

