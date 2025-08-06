#' Clayton 2025 Spatial Analysis
#'
#' Main function to run the complete spatial analysis pipeline for Clayton field experiment.

#' Run Clayton 2025 comprehensive spatial analysis
#' 
#' Executes the complete spatial analysis pipeline reproducing the results from 
#' docs/notebooks/inv4m_field_modelling.Rmd as a standalone analysis function.
#' 
#' @param data_file Path to the CLY25_Inv4m.csv data file (default: "data/CLY25_Inv4m.csv")
#' @param output_dir Directory for saving results (default: "results_clayton_2025")
#' @param phenotypes Vector of phenotype names to analyze (default: all available)
#' @param create_plots Logical, whether to create diagnostic plots (default: TRUE)
#' @param verbose Logical, whether to print detailed progress messages (default: TRUE)
#' 
#' @return List containing all analysis results
#' @export
#' @import dplyr
#' @import ggplot2
#' @import tidyverse
#' @import purrr
run_clayton_spatial_analysis <- function(data_file = "data/CLY25_Inv4m.csv",
                                        output_dir = "results_clayton_2025",
                                        phenotypes = NULL,
                                        create_plots = TRUE,
                                        verbose = TRUE) {
  
  if (verbose) {
    cat("=== CLAYTON 2025 FIELD SPATIAL ANALYSIS ===\n")
    cat("Multi-trait spatial analysis pipeline\n")
    cat("Date:", as.character(Sys.Date()), "\n\n")
  }
  
  # Validate setup
  validate_analysis_setup()
  
  # Load and clean data
  if (!file.exists(data_file)) {
    stop("Error: Data file not found at: ", data_file)
  }
  
  field_data_raw <- read.csv(data_file, na.strings = c("","#N/A","NA"))
  
  # Clean and prepare data following Rmd approach
  field_data <- field_data_raw %>%
    # Calculate Estimated Blade Area (EBA) = 0.75 * BL * BW FIRST
    mutate(EBA = 0.75 * BL * BW) %>%
    rename(
      plant_id = plant,
      block = rep,
      x = X_pos,
      y = Y_pos,
      inv4m = inv4m_gt
    ) %>%
    # Convert relevant columns to factors
    mutate(
      plot_id = as.factor(plot_id),
      block = as.factor(block),
      donor = as.factor(donor),
      inv4m = as.factor(inv4m)
    ) %>%
    # Add small amount of noise to x coordinates to avoid identical positions
    mutate(x = x + runif(n(), min = 0.0, max = 0.01)) %>%
    # Create centered coordinates and other required variables
    mutate(
      x_c = x - mean(x, na.rm = TRUE),
      y_c = y - mean(y, na.rm = TRUE),
      plot_row = as.numeric(as.character(plot_row)),
      plot_col = as.numeric(as.character(plot_col))
    )
  
  # Define comprehensive phenotypes list - all traits from DTA to BW + EBA
  if (is.null(phenotypes)) {
    phenotypes <- c("DTA", "DTS", "LAE", "PH", "EN", "SL", "BL", "BW", "EBA")
  }
  
  # Verify all phenotypes exist in the dataset
  available_phenotypes <- intersect(phenotypes, names(field_data))
  phenotypes <- available_phenotypes
  
  if (verbose) {
    cat("Data loaded:", nrow(field_data), "observations\n")
    cat("Phenotypes for analysis:", paste(phenotypes, collapse = ", "), "\n")
    cat("Number of phenotypes:", length(phenotypes), "\n\n")
  }
  
  # Missing Data Assessment
  if (verbose) cat("=== MISSING DATA ASSESSMENT ===\n")
  
  missing_summary <- field_data %>%
    select(all_of(phenotypes)) %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "phenotype", values_to = "missing_count") %>%
    mutate(
      total_obs = nrow(field_data),
      missing_pct = round(missing_count / total_obs * 100, 1),
      available_n = total_obs - missing_count
    )
  
  if (verbose) print(missing_summary)
  
  # Scaled Variogram Analysis
  if (verbose) cat("\n=== SCALED VARIOGRAM ANALYSIS ===\n")
  
  variogram_results <- list()
  
  for (trait in phenotypes) {
    if (verbose) cat("Processing variogram for", trait, "...\n")
    
    trait_data <- field_data %>%
      filter(!is.na(.data[[trait]]) & 
             !is.na(x) & !is.na(y) & 
             !is.na(donor) & !is.na(inv4m) & !is.na(block))
    
    if (verbose) cat("  Sample size:", nrow(trait_data), "\n")
    
    if (nrow(trait_data) > 10) {
      tryCatch({
        variogram_results[[trait]] <- calculate_scaled_variogram(trait_data, trait)
        if (verbose) cat("  Successfully calculated variogram\n")
      }, error = function(e) {
        if (verbose) cat("  ERROR calculating variogram:", e$message, "\n")
      })
    } else {
      if (verbose) cat("  Warning: Insufficient data for", trait, "\n")
    }
  }
  
  # Model Fitting
  if (verbose) cat("\n=== COMPREHENSIVE MODEL COMPARISON ===\n")
  
  all_models <- list()
  all_model_stats <- tibble()
  
  for (trait in phenotypes) {
    if (verbose) cat("\n=== FITTING MODELS FOR", trait, "===\n")
    
    tryCatch({
      trait_data <- field_data %>%
        filter(!is.na(.data[[trait]]) & 
               !is.na(x) & !is.na(y) & 
               !is.na(donor) & !is.na(inv4m) & 
               !is.na(block) & !is.na(plot_id) &
               !is.na(plot_row) & !is.na(plot_col))
      
      if (verbose) cat("  Sample size:", nrow(trait_data), "\n")
      
      if (nrow(trait_data) > 50) {
        if (verbose) cat("  Fitting models...\n")
        trait_models <- fit_all_models(trait_data, trait)
        all_models[[trait]] <- trait_models
        
        successful_models <- names(trait_models)[!sapply(trait_models, is.null)]
        
        if (length(successful_models) > 0) {
          if (verbose) cat("  Extracting model statistics...\n")
          trait_stats <- extract_model_stats(trait_models, trait)
          all_model_stats <- bind_rows(all_model_stats, trait_stats)
        }
      }
    }, error = function(e) {
      if (verbose) cat("  ERROR processing", trait, ":", e$message, "\n")
    })
  }
  
  # Export Results
  if (verbose) cat("\n=== EXPORTING RESULTS ===\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  if (nrow(all_model_stats) > 0) {
    write.csv(all_model_stats, file.path(output_dir, "all_model_statistics.csv"), row.names = FALSE)
  }
  
  write.csv(missing_summary, file.path(output_dir, "missing_data_summary.csv"), row.names = FALSE)
  
  if (length(variogram_results) > 0) {
    vgm_summary <- map_dfr(variogram_results, function(x) {
      tibble(
        phenotype = x$phenotype,
        n_obs = x$n_obs,
        max_semivariance = round(x$max_gamma, 3)
      )
    })
    write.csv(vgm_summary, file.path(output_dir, "variogram_summary.csv"), row.names = FALSE)
  }
  
  if (verbose) {
    cat("Results exported to:", output_dir, "\n")
    cat("\n=== ANALYSIS COMPLETE ===\n")
  }
  
  return(list(
    field_data = field_data,
    missing_summary = missing_summary,
    variogram_results = variogram_results,
    all_models = all_models,
    model_stats = all_model_stats,
    phenotypes = phenotypes
  ))
}