#' Clayton 2025 Spatial Analysis - Multi-trait Pipeline
#' 
#' This script performs comprehensive spatial analysis of the Clayton 2025 field experiment
#' using the inv4mHighland package functions. It contains the Clayton-specific data loading
#' and processing logic while leveraging the general spatial analysis functions.
#' 
#' @author Fausto Rodriguez
#' @date 2025-08-07

# Load required packages ----
library(inv4mHighland)
library(tidyverse)
library(ggpubr)

# Clayton-specific data loading function ----
load_clayton_data <- function(data_file) {
  if (!file.exists(data_file)) {
    stop("Error: Data file not found at: ", data_file)
  }
  
  field_data_raw <- read.csv(data_file, na.strings = c("","#N/A","NA"))
  
  # Clean and prepare data following Clayton 2025 approach
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
  
  return(field_data)
}

# Setup and Data Preparation ----
cat("=== CLAYTON 2025 FIELD SPATIAL ANALYSIS ===\n")
cat("Multi-trait spatial analysis pipeline\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Analysis Parameters ----
file_path <- "../../../data/CLY25_Inv4m.csv"
output_dir <- "results_clayton_2025"
phenotypes_to_analyze <- c("DTA", "DTS", "LAE", "PH", "EN", "SL", "BL", "BW", "EBA")


# Load and clean Clayton 2025 data
cat("Loading and preprocessing Clayton 2025 data...\n")
field_data <- load_clayton_data(file_path)

# Verify all phenotypes exist in the dataset
available_phenotypes <- intersect(phenotypes_to_analyze, names(field_data))
phenotypes_to_analyze <- available_phenotypes

cat("Data loaded:", nrow(field_data), "observations\n")
cat("Phenotypes for analysis:", paste(phenotypes_to_analyze, collapse = ", "), "\n")
cat("Number of phenotypes:", length(phenotypes_to_analyze), "\n\n")

# Missing Data Assessment ----
cat("=== MISSING DATA ASSESSMENT ===\n")
missing_results <- assess_missing_data(field_data, phenotypes_to_analyze, verbose = TRUE)

# Scaled Variogram Analysis ----
cat("\n=== SCALED VARIOGRAM ANALYSIS ===\n")
variogram_results <- calculate_scaled_variogram(field_data, phenotypes_to_analyze)

if (length(variogram_results) > 0) {
  cat("Variograms calculated for", length(variogram_results), "traits\n")
  
  # Create variogram summary plot
  all_vgm_data <- map_dfr(variogram_results, function(x) {
    x$variogram %>%
      mutate(phenotype = x$phenotype,
             n_obs = x$n_obs)
  })
  
  variogram_plot <- ggplot(all_vgm_data, aes(x = dist, y = gamma_scaled, color = phenotype)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_color_brewer(type = "qual", palette = "Set3") +
    labs(
      title = "Scaled Empirical Variograms - All Phenotypes",
      subtitle = "Scaled to 0-100 for direct comparison of spatial autocorrelation patterns",
      x = "Distance", y = "Scaled Semivariance (0-100)", color = "Phenotype"
    ) +
    theme_classic() +
    theme(legend.position = "right")
  
  print(variogram_plot)
  
  # Save variogram plot
  ggsave(file.path(output_dir, "scaled_variograms_all_traits.png"), 
         variogram_plot, width = 12, height = 8, dpi = 300)
}

# Comprehensive Model Comparison ----
cat("\n=== COMPREHENSIVE MODEL COMPARISON ===\n")
all_models <- analyze_multiple_traits(field_data, phenotypes_to_analyze, verbose = TRUE)

# Extract model statistics
all_model_stats <- tibble()
for (trait in names(all_models)) {
  if (!is.null(all_models[[trait]])) {
    trait_stats <- extract_model_stats(all_models[[trait]], trait)
    all_model_stats <- bind_rows(all_model_stats, trait_stats)
  }
}

# Model comparison and best model selection
best_models <- tibble()
if (nrow(all_model_stats) > 0) {
  # Find best model for each trait (excluding model_1 - fixed effects only)
  best_models <- all_model_stats %>%
    filter(!is.na(BIC) & model != "model_1") %>%
    group_by(phenotype) %>%
    slice_min(BIC, n = 1) %>%
    select(phenotype, best_model = model, best_BIC = BIC)
  
  cat("\n=== BEST MODEL SELECTION ===\n")
  print(best_models)
  
  # Model selection frequency
  if (nrow(best_models) > 0) {
    model_frequency <- best_models %>%
      count(best_model, sort = TRUE) %>%
      mutate(frequency = n / nrow(best_models) * 100)
    
    cat("\n=== MODEL SELECTION FREQUENCY ===\n")
    print(model_frequency)
  }
}

# Treatment Effects Analysis ----
cat("\n=== TREATMENT EFFECTS ANALYSIS ===\n")
treatment_effects <- list()
significant_effects <- tibble()

if (nrow(best_models) > 0) {
  for (trait in names(all_models)) {
    best_model_info <- best_models %>% filter(phenotype == trait)
    
    if (nrow(best_model_info) > 0) {
      best_model_name <- best_model_info$best_model
      
      effects <- extract_treatment_effects_emmeans(
        all_models[[trait]], 
        best_model_name, 
        trait
      )
      
      if (!is.null(effects)) {
        treatment_effects[[trait]] <- effects
      }
    }
  }
  
  # Summary of significant effects
  if (length(treatment_effects) > 0) {
    significant_effects <- map_dfr(treatment_effects, function(x) {
      if (!is.null(x$contrast_summary)) {
        x$contrast_summary %>%
          as_tibble() %>%
          mutate(phenotype = x$phenotype,
                 model_used = x$model_used) %>%
          filter(p.value < 0.05) %>%
          select(phenotype, model_used, donor, estimate, SE, p.value)
      }
    })
    
    if (nrow(significant_effects) > 0) {
      cat("\n=== SIGNIFICANT TREATMENT EFFECTS (p < 0.05) ===\n")
      print(significant_effects)
    } else {
      cat("\n=== NO SIGNIFICANT TREATMENT EFFECTS DETECTED ===\n")
    }
  }
}

# Spatial Distribution Visualization ----
cat("\n=== SPATIAL DISTRIBUTION VISUALIZATION ===\n")
show_spatial_distribution(field_data, phenotypes_to_analyze)

# Generate Diagnostic Plots for Best Models ----
cat("\n=== GENERATING DIAGNOSTIC PLOTS ===\n")

if (nrow(best_models) > 0) {
  for (trait in names(all_models)) {
    if (trait %in% best_models$phenotype) {
      best_model_name <- best_models$best_model[best_models$phenotype == trait]
      best_model <- all_models[[trait]][[best_model_name]]
      
      if (!is.null(best_model)) {
        cat("Generating diagnostics for", trait, "using", best_model_name, "...\n")
        
        # Create comprehensive residual diagnostics
        trait_data <- field_data[!is.na(field_data[[trait]]), ]
        diagnostic_plots <- residual_diagnostics(best_model, trait_data, trait)
        
        # Display plots
        print(diagnostic_plots$resid_fitted)
        if (!is.null(diagnostic_plots$spatial)) {
          print(diagnostic_plots$spatial)
        }
        print(diagnostic_plots$qq)
      }
    }
  }
}

# Export Results ----
cat("\n=== EXPORTING RESULTS ===\n")

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare comprehensive results list
results <- list(
  field_data = field_data,
  missing_summary = missing_results$summary,
  missing_by_treatment = missing_results$by_treatment,
  variogram_results = variogram_results,
  all_models = all_models,
  model_stats = all_model_stats,
  best_models = best_models,
  treatment_effects = treatment_effects,
  significant_effects = significant_effects,
  phenotypes = phenotypes_to_analyze
)

# Export using package function
export_spatial_results(results, output_dir, verbose = TRUE)

# Save additional Clayton-specific results
if (nrow(best_models) > 0) {
  write.csv(best_models, file.path(output_dir, "best_models_selection.csv"), row.names = FALSE)
}

if (nrow(significant_effects) > 0) {
  write.csv(significant_effects, file.path(output_dir, "significant_treatment_effects.csv"), row.names = FALSE)
}

if (exists("model_frequency")) {
  write.csv(model_frequency, file.path(output_dir, "model_selection_frequency.csv"), row.names = FALSE)
}

# Generate comprehensive summary report
summary_text <- paste0(
  "CLAYTON 2025 SPATIAL ANALYSIS SUMMARY\n",
  "====================================\n\n",
  "Analysis date: ", Sys.Date(), "\n",
  "Total observations: ", nrow(field_data), "\n",
  "Phenotypes analyzed: ", paste(phenotypes_to_analyze, collapse = ", "), "\n",
  "Number of phenotypes: ", length(phenotypes_to_analyze), "\n\n"
)

# Add missing data summary
summary_text <- paste0(summary_text, "MISSING DATA PATTERNS\n", "--------------------\n")
for (i in 1:nrow(missing_results$summary)) {
  row <- missing_results$summary[i, ]
  summary_text <- paste0(summary_text,
    row$phenotype, ": ", row$missing_count, " missing (", row$missing_pct, "%), ",
    row$available_n, " available\n"
  )
}

# Add variogram summary
if (length(variogram_results) > 0) {
  summary_text <- paste0(summary_text, "\nVARIOGRAM ANALYSIS\n", "-----------------\n")
  summary_text <- paste0(summary_text, "Variograms calculated for ", length(variogram_results), " traits\n")
  summary_text <- paste0(summary_text, "All traits showed evidence of spatial structure\n")
}

# Add model selection results
if (nrow(best_models) > 0) {
  summary_text <- paste0(summary_text, "\nBEST MODELS BY PHENOTYPE\n", "------------------------\n")
  for (i in 1:nrow(best_models)) {
    row <- best_models[i, ]
    summary_text <- paste0(summary_text,
      row$phenotype, ": ", row$best_model, " (BIC: ", round(row$best_BIC, 2), ")\n"
    )
  }
}

# Add treatment effects summary
if (nrow(significant_effects) > 0) {
  summary_text <- paste0(summary_text, "\nSIGNIFICANT TREATMENT EFFECTS (p < 0.05)\n", 
                        "---------------------------------------\n")
  for (i in 1:nrow(significant_effects)) {
    row <- significant_effects[i, ]
    summary_text <- paste0(summary_text,
      row$phenotype, " (", row$donor, "): ", round(row$estimate, 3),
      " ± ", round(row$SE, 3), " (p = ", format.pval(row$p.value, digits = 3), ")\n"
    )
  }
} else {
  summary_text <- paste0(summary_text, "\nTREATMENT EFFECTS\n", "-----------------\n",
                        "No significant inv4m effects detected at p < 0.05\n")
}

writeLines(summary_text, file.path(output_dir, "clayton_2025_analysis_summary.txt"))

# Final Summary ----
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results exported to:", output_dir, "\n")
cat("\nKey outputs:\n")
cat("- all_model_statistics.csv: Complete model comparison\n")
cat("- best_models_selection.csv: Optimal model for each trait\n")
cat("- significant_treatment_effects.csv: Treatment effects (if any)\n")
cat("- missing_data_summary.csv: Missing data patterns\n")
cat("- variogram_summary.csv: Spatial autocorrelation summary\n")
cat("- clayton_2025_analysis_summary.txt: Human-readable report\n")
cat("- scaled_variograms_all_traits.png: Variogram visualization\n")
cat("- Multiple spatial distribution and diagnostic plots\n\n")

cat("This analysis successfully demonstrates the power of the inv4mHighland package\n")
cat("for comprehensive spatial analysis while keeping dataset-specific logic\n")
cat("in the analysis script where it belongs.\n")