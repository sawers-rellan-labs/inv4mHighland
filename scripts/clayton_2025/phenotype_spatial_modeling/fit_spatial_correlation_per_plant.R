#' Multi-Trait Spatial Analysis of Replicated Latin Square Field Experiment
#' 
#' This script reproduces the comprehensive spatial analysis from 
#' docs/inv4m_field_modelling.Rmd as a standalone R script.
#' Analyzes all measured phenotypes in the inv4m Clayton field experiment using:
#' 1. Proper missing data handling for valid model comparisons
#' 2. Scaled variogram analysis across traits
#' 3. Six hierarchical model structures for each phenotype 
#' 4. Treatment effect quantification using optimal models
#'
#' @author Fausto Rodriguez
#' @date 2025-08-06

# Load required libraries ----
library(tidyverse)
library(dplyr)
library(ggplot2)
library(nlme)
library(gstat)
library(emmeans)
library(ggpubr)
library(ggtext)
library(VIM)
library(mgcv)
library(ape)

# Source helper functions ----
source("spatial_correlation_helpers.R")

# Setup and Data Preparation ----
cat("=== CLAYTON 2025 FIELD SPATIAL ANALYSIS ===\n")
cat("Multi-trait spatial analysis pipeline\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# Load and clean data
file_path <- "../../../data/CLY25_Inv4m.csv"
if (!file.exists(file_path)) {
  stop("Error: CLY25_Inv4m.csv not found at: ", file_path)
}

field_data_raw <- read.csv(file_path, na.strings = c("","#N/A","NA"))

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
phenotypes <- c("DTA", "DTS", "LAE", "PH", "EN", "SL", "BL", "BW", "EBA")

# Verify all phenotypes exist in the dataset
missing_cols <- setdiff(phenotypes, names(field_data))
if (length(missing_cols) > 0) {
  cat("Warning: These phenotype columns are missing from the dataset:\n")
  print(missing_cols)
}

available_phenotypes <- intersect(phenotypes, names(field_data))
phenotypes <- available_phenotypes

cat("Data loaded:", nrow(field_data), "observations\n")
cat("Phenotypes for analysis:", paste(phenotypes, collapse = ", "), "\n")
cat("Number of phenotypes:", length(phenotypes), "\n\n")

# Missing Data Assessment ----
cat("=== MISSING DATA ASSESSMENT ===\n")

# Create missing data summary
missing_summary <- field_data %>%
  select(all_of(phenotypes)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "phenotype", values_to = "missing_count") %>%
  mutate(
    total_obs = nrow(field_data),
    missing_pct = round(missing_count / total_obs * 100, 1),
    available_n = total_obs - missing_count
  )

print(missing_summary)

# Check for systematic missing data
missing_by_treatment <- field_data %>%
  select(donor, inv4m, all_of(phenotypes)) %>%
  group_by(donor, inv4m) %>%
  summarise(across(all_of(phenotypes), ~sum(is.na(.))), .groups = 'drop')

cat("\nMissing data counts by treatment combination:\n")
print(missing_by_treatment)

# Spatial distribution of missing data
missing_spatial_plot <- field_data %>%
  select(x, y, all_of(phenotypes)) %>%
  mutate(any_missing = rowSums(is.na(select(., all_of(phenotypes)))) > 0) %>%
  ggplot(aes(x = x, y = y, color = any_missing)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  labs(title = "Spatial distribution of missing observations",
       color = "Has missing\ndata") +
  theme_classic() +
  coord_equal()

print(missing_spatial_plot)

# Scaled Variogram Analysis ----
cat("\n=== SCALED VARIOGRAM ANALYSIS ===\n")

# Calculate scaled variograms for each phenotype
variogram_results <- list()

for (trait in phenotypes) {
  cat("Processing variogram for", trait, "...\n")
  
  # Check if trait exists in data
  if (!trait %in% names(field_data)) {
    cat("  ERROR: Column", trait, "not found in field_data\n")
    next
  }
  
  # Create complete-case dataset for this trait
  trait_data <- field_data %>%
    filter(!is.na(.data[[trait]]) & 
           !is.na(x) & !is.na(y) & 
           !is.na(donor) & !is.na(inv4m) & !is.na(block))
  
  cat("  Sample size:", nrow(trait_data), "\n")
  
  if (nrow(trait_data) > 10) {
    tryCatch({
      variogram_results[[trait]] <- calculate_scaled_variogram(trait_data, trait)
      cat("  Successfully calculated variogram\n")
    }, error = function(e) {
      cat("  ERROR calculating variogram:", e$message, "\n")
    })
  } else {
    cat("  Warning: Insufficient data for", trait, "\n")
  }
}

cat("\nVariograms successfully calculated for:", length(variogram_results), "traits\n")

# Create variogram summary
if (length(variogram_results) > 0) {
  vgm_summary <- map_dfr(variogram_results, function(x) {
    tibble(
      phenotype = x$phenotype,
      n_obs = x$n_obs,
      max_semivariance = round(x$max_gamma, 3)
    )
  })
  
  cat("\n=== VARIOGRAM SUMMARY STATISTICS ===\n")
  print(vgm_summary)
  
  # Plot scaled variograms
  all_vgm_data <- map_dfr(variogram_results, function(x) {
    x$variogram %>%
      mutate(phenotype = x$phenotype,
             n_obs = x$n_obs)
  })
  
  # Create combined variogram plot
  variogram_plot <- ggplot(all_vgm_data, aes(x = dist, y = gamma_scaled, color = phenotype)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_color_brewer(type = "qual", palette = "Set3") +
    labs(
      title = "Scaled Empirical Variograms - All Phenotypes",
      subtitle = "Scaled to 0-100 for direct comparison of spatial autocorrelation patterns",
      x = "Distance",
      y = "Scaled Semivariance (0-100)",
      color = "Phenotype"
    ) +
    theme_classic2(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
  
  print(variogram_plot)
} else {
  cat("No variograms could be calculated!\n")
}

# Comprehensive Model Comparison ----
cat("\n=== COMPREHENSIVE MODEL COMPARISON ===\n")

# Initialize results storage
all_models <- list()
all_model_stats <- tibble()
model_fitting_diagnostics <- tibble()

# Fit models for each phenotype
for (trait in phenotypes) {
  cat("\n=== FITTING MODELS FOR", trait, "===\n")
  
  tryCatch({
    # Create complete-case dataset
    trait_data <- field_data %>%
      filter(!is.na(.data[[trait]]) & 
             !is.na(x) & !is.na(y) & 
             !is.na(donor) & !is.na(inv4m) & 
             !is.na(block) & !is.na(plot_id) &
             !is.na(plot_row) & !is.na(plot_col))
    
    cat("  Sample size:", nrow(trait_data), "\n")
    
    # Record diagnostics for this trait
    trait_diagnostics <- tibble(
      phenotype = trait,
      sample_size = nrow(trait_data),
      n_donors = length(unique(trait_data$donor)),
      n_inv4m_levels = length(unique(trait_data$inv4m)),
      n_blocks = length(unique(trait_data$block)),
      n_plots = length(unique(trait_data$plot_id)),
      trait_range = paste(round(range(trait_data[[trait]], na.rm = TRUE), 2), collapse = " - "),
      models_attempted = 0,
      models_successful = 0,
      failed_models = "",
      error_messages = ""
    )
    
    if (nrow(trait_data) > 50) {
      # Fit all models using comprehensive hierarchy
      cat("  Fitting models...\n")
      trait_models <- fit_all_models(trait_data, trait)
      all_models[[trait]] <- trait_models
      
      # Count successful models
      successful_models <- names(trait_models)[!sapply(trait_models, is.null)]
      failed_models <- names(trait_models)[sapply(trait_models, is.null)]
      
      trait_diagnostics$models_attempted <- 6
      trait_diagnostics$models_successful <- length(successful_models)
      trait_diagnostics$failed_models <- paste(failed_models, collapse = ", ")
      
      cat("  Successful models:", length(successful_models), "out of 6\n")
      if (length(failed_models) > 0) {
        cat("  Failed models:", paste(failed_models, collapse = ", "), "\n")
      }
      
      # Extract statistics
      if (length(successful_models) > 0) {
        cat("  Extracting model statistics...\n")
        trait_stats <- extract_model_stats(trait_models, trait)
        cat("    Extracted", nrow(trait_stats), "rows of statistics\n")
        all_model_stats <- bind_rows(all_model_stats, trait_stats)
        cat("    Total model stats now:", nrow(all_model_stats), "rows\n")
      } else {
        cat("  No successful models for", trait, "\n")
      }
    } else {
      cat("  Warning: Insufficient data for", trait, "\n")
      trait_diagnostics$error_messages <- "Insufficient sample size (< 50 observations)"
    }
    
    # Add diagnostics
    model_fitting_diagnostics <- bind_rows(model_fitting_diagnostics, trait_diagnostics)
    
  }, error = function(e) {
    cat("  MAJOR ERROR processing", trait, ":", e$message, "\n")
    # Still add basic diagnostics even if everything failed
    error_diagnostics <- tibble(
      phenotype = trait,
      sample_size = NA,
      n_donors = NA,
      n_inv4m_levels = NA,
      n_blocks = NA,
      n_plots = NA,
      trait_range = NA,
      models_attempted = 0,
      models_successful = 0,
      failed_models = "",
      error_messages = paste("Major error:", e$message)
    )
    model_fitting_diagnostics <<- bind_rows(model_fitting_diagnostics, error_diagnostics)
  })
}

cat("\n=== FINAL SUMMARY ===\n")
cat("Processed", length(all_models), "phenotypes out of", length(phenotypes), "requested\n")
cat("Phenotypes with models:", names(all_models), "\n")
cat("Missing phenotypes:", setdiff(phenotypes, names(all_models)), "\n")

# Display model fitting diagnostics
cat("\n=== MODEL FITTING DIAGNOSTICS ===\n")
print(model_fitting_diagnostics)

# Model Comparison Results
if (nrow(all_model_stats) > 0) {
  # BIC comparison table - exclude model_1 from comparison
  bic_matrix <- all_model_stats %>%
    filter(model != "model_1") %>%
    select(phenotype, model, BIC) %>%
    complete(phenotype, model, fill = list(BIC = NA)) %>%
    pivot_wider(names_from = model, values_from = BIC)
  
  # Convert to matrix for display
  if (nrow(bic_matrix) > 0 && ncol(bic_matrix) > 1) {
    if ("phenotype" %in% names(bic_matrix)) {
      bic_table <- bic_matrix %>%
        select(-phenotype) %>%
        as.data.frame()
      
      rownames(bic_table) <- bic_matrix$phenotype
      
      cat("\n=== BIC VALUES FOR MIXED-EFFECTS MODELS (LOWER IS BETTER) ===\n")
      cat("Note: Model_1 (fixed effects only) excluded from comparison\n")
      print(round(bic_table, 1))
    }
  }
  
  # Best model for each phenotype (excluding model_1)
  best_models <- all_model_stats %>%
    filter(!is.na(BIC) & model != "model_1") %>%
    group_by(phenotype) %>%
    slice_min(BIC, n = 1) %>%
    select(phenotype, best_model = model, best_BIC = BIC)
  
  if (nrow(best_models) > 0) {
    cat("\n=== BEST MIXED-EFFECTS MODEL (LOWEST BIC) FOR EACH PHENOTYPE ===\n")
    print(best_models)
    
    # Model selection frequency
    model_frequency <- best_models %>%
      count(best_model, sort = TRUE) %>%
      mutate(frequency = n / nrow(best_models) * 100)
    
    cat("\n=== FREQUENCY OF EACH MODEL BEING SELECTED AS BEST ===\n")
    print(model_frequency)
  }
} else {
  cat("No model statistics available\n")
  best_models <- tibble()
}

# Treatment Effects Analysis ----
cat("\n=== TREATMENT EFFECTS ANALYSIS ===\n")

# Extract treatment effects using best models
treatment_effects <- list()
significant_effects <- tibble()

if (nrow(best_models) > 0) {
  phenotypes_for_effects <- names(all_models)
  
  for (trait in phenotypes_for_effects) {
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
  
  # Summary of significant contrasts
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
}

if (nrow(significant_effects) > 0) {
  cat("\n=== SIGNIFICANT INV4M EFFECTS (P < 0.05) BY DONOR BACKGROUND ===\n")
  print(significant_effects)
} else {
  cat("No significant inv4m effects detected at p < 0.05\n")
}

# Create forest plot if we have treatment effects
if (length(treatment_effects) > 0) {
  forest_data <- map_dfr(treatment_effects, function(x) {
    if (!is.null(x$contrast_summary)) {
      x$contrast_summary %>%
        as_tibble() %>%
        mutate(
          phenotype = x$phenotype,
          model_used = x$model_used,
          std_effect = estimate / SE,
          ci_lower = (estimate - 1.96 * SE) / SE,
          ci_upper = (estimate + 1.96 * SE) / SE
        ) %>%
        select(phenotype, donor, std_effect, ci_lower, ci_upper, p.value)
    }
  })
  
  if (nrow(forest_data) > 0) {
    forest_data <- forest_data %>%
      mutate(p.adj = p.adjust(p.value, method = "bonferroni")) %>%
      filter(p.value < 0.05) %>%
      mutate(bonferroni_sig = p.adj < 0.05)
    
    if (nrow(forest_data) > 0) {
      forest_data <- forest_data %>%
        mutate(
          trait_label = case_when(
            phenotype == "DTA" ~ "Days to Anthesis",
            phenotype == "DTS" ~ "Days to Silking", 
            phenotype == "LAE" ~ "Leaves Above Ear",
            phenotype == "PH" ~ "Plant Height",
            phenotype == "EN" ~ "Nodes with Ears",
            phenotype == "SL" ~ "Sheath Length",
            phenotype == "BL" ~ "Blade Length",
            phenotype == "BW" ~ "Blade Width",
            phenotype == "EBA" ~ "Estimated Blade Area",
            TRUE ~ phenotype
          ),
          predictor = paste0("Inv4m from ", donor)
        ) %>%
        arrange(desc(abs(std_effect))) %>%
        mutate(trait_label = factor(trait_label, levels = unique(trait_label)))
      
      # Create forest plot
      forest_plot <- forest_data %>%
        ggplot(aes(x = std_effect, y = trait_label)) +
        geom_vline(xintercept = 0, lty = 2) +
        geom_point(size = 4, aes(color = bonferroni_sig)) +
        geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, size = 0.7) +
        facet_wrap(. ~ predictor, ncol = 2) +
        scale_color_manual(
          values = c("TRUE" = "red", "FALSE" = "black"),
          labels = c("TRUE" = "Bonferroni significant", "FALSE" = "Nominally significant"),
          name = "Significance"
        ) +
        labs(
          title = "Forest Plot of Standardized inv4m Effects",
          x = "Standardized Effect Size",
          y = "Trait",
          caption = "All nominally significant effects shown (p < 0.05). Red points survive Bonferroni correction."
        ) +
        theme_classic2(base_size = 12) +
        theme(
          strip.background = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom"
        )
      
      print(forest_plot)
      
      cat("\n=== ALL NOMINALLY SIGNIFICANT EFFECTS (p < 0.05) ===\n")
      print(forest_data %>% select(phenotype, donor, std_effect, p.value, p.adj, bonferroni_sig) %>% arrange(p.value))
      
      bonferroni_sig_effects <- forest_data %>% filter(bonferroni_sig)
      if (nrow(bonferroni_sig_effects) > 0) {
        cat("\n=== BONFERRONI-CORRECTED SIGNIFICANT EFFECTS (p.adj < 0.05) ===\n")
        print(bonferroni_sig_effects %>% select(phenotype, donor, std_effect, p.value, p.adj) %>% arrange(p.adj))
      } else {
        cat("\n=== NO EFFECTS SURVIVE BONFERRONI CORRECTION ===\n")
      }
    }
  }
}

# Spatial Distribution Visualization ----
cat("\n=== SPATIAL DISTRIBUTION VISUALIZATION ===\n")

# Create spatial plots for key traits
if ("DTA" %in% names(field_data)) {
  plot_dta <- create_spatial_plot(field_data, DTA, "Days to Anthesis", "days")
  print(plot_dta)
}

if ("DTS" %in% names(field_data)) {
  plot_dts <- create_spatial_plot(field_data, DTS, "Days to Silking", "days")
  print(plot_dts)
}

if ("PH" %in% names(field_data)) {
  plot_ph <- create_spatial_plot(field_data, PH, "Plant Height", "cm")
  print(plot_ph)
}

if ("EBA" %in% names(field_data)) {
  plot_eba <- create_spatial_plot(field_data, EBA, "Estimated Blade Area", "cm²")
  print(plot_eba)
}

# Generate diagnostic plots for best models
cat("\n=== GENERATING DIAGNOSTIC PLOTS ===\n")

for (trait in names(all_models)) {
  if (trait %in% best_models$phenotype) {
    best_model_name <- best_models$best_model[best_models$phenotype == trait]
    best_model <- all_models[[trait]][[best_model_name]]
    
    if (!is.null(best_model)) {
      cat("\nGenerating diagnostics for", trait, "...\n")
      
      # Create comprehensive residual diagnostics
      trait_data <- field_data[!is.na(field_data[[trait]]), ]
      
      diagnostic_plots <- residual_diagnostics(best_model, trait_data, trait)
      
      # Display plots
      print(diagnostic_plots$resid_fitted)
      if (!is.null(diagnostic_plots$spatial)) {
        print(diagnostic_plots$spatial)
      }
      print(diagnostic_plots$qq)
      
      # Test for remaining spatial autocorrelation
      if (nrow(trait_data) > 10) {
        trait_data$residuals <- residuals(best_model, type = "normalized")
        
        coords <- cbind(trait_data$x, trait_data$y)
        coords_complete <- coords[complete.cases(coords), ]
        residuals_complete <- trait_data$residuals[complete.cases(trait_data$residuals)]
        
        if (length(residuals_complete) > 10 && nrow(coords_complete) > 10) {
          dists <- as.matrix(dist(coords_complete))
          diag(dists) <- 1
          
          moran_test <- tryCatch({
            Moran.I(residuals_complete, weight = 1/dists)
          }, error = function(e) {
            cat("  Moran's I test failed:", e$message, "\n")
            NULL
          })
          
          if (!is.null(moran_test)) {
            cat("  Moran's I:", round(moran_test$observed, 4), 
                "(p-value:", round(moran_test$p.value, 4), ")\n")
            
            if (moran_test$p.value < 0.05) {
              cat("  WARNING: Significant spatial autocorrelation remains\n")
            } else {
              cat("  Good: No significant spatial autocorrelation detected\n")
            }
          }
        }
      }
    }
  }
}

# Summary and Results Export ----
cat("\n=== ANALYSIS SUMMARY ===\n")

cat("1. MISSING DATA:\n")
print(missing_summary)

cat("\n2. MODEL SELECTION FREQUENCY:\n")
if (exists("model_frequency") && nrow(model_frequency) > 0) {
  print(model_frequency)
} else {
  cat("No successful model fits available\n")
}

cat("\n3. SIGNIFICANT TREATMENT EFFECTS:\n")
if (nrow(significant_effects) > 0) {
  print(significant_effects)
} else {
  cat("No significant effects detected at p < 0.05\n")
}

cat("\n4. SPATIAL AUTOCORRELATION:\n")
if (length(variogram_results) > 0) {
  cat("Variograms calculated for", length(variogram_results), "phenotypes\n")
  cat("All phenotypes showed evidence of spatial structure in their variograms\n")
} else {
  cat("No variograms could be calculated\n")
}

# Export Results ----
cat("\n=== EXPORTING RESULTS ===\n")

# Create output directory
output_dir <- "results_clayton_2025"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save model comparison results
if (nrow(all_model_stats) > 0) {
  write.csv(all_model_stats, file.path(output_dir, "all_model_statistics.csv"), row.names = FALSE)
}

# Save best models
if (nrow(best_models) > 0) {
  write.csv(best_models, file.path(output_dir, "best_models_selection.csv"), row.names = FALSE)
}

# Save treatment effects
if (nrow(significant_effects) > 0) {
  write.csv(significant_effects, file.path(output_dir, "significant_treatment_effects.csv"), row.names = FALSE)
}

# Save missing data summary
write.csv(missing_summary, file.path(output_dir, "missing_data_summary.csv"), row.names = FALSE)

# Save variogram results summary
if (exists("vgm_summary")) {
  write.csv(vgm_summary, file.path(output_dir, "variogram_summary.csv"), row.names = FALSE)
}

# Create comprehensive summary report
summary_text <- paste0(
  "CLAYTON 2025 MULTI-TRAIT SPATIAL ANALYSIS SUMMARY\n",
  "================================================\n\n",
  "Analysis date: ", Sys.Date(), "\n",
  "Total observations: ", nrow(field_data), "\n",
  "Phenotypes analyzed: ", paste(phenotypes, collapse = ", "), "\n",
  "Number of phenotypes: ", length(phenotypes), "\n\n",
  
  "MISSING DATA PATTERNS\n",
  "--------------------\n"
)

for (i in 1:nrow(missing_summary)) {
  row <- missing_summary[i, ]
  summary_text <- paste0(summary_text,
    row$phenotype, ": ", row$missing_count, " missing (", row$missing_pct, "%), ",
    row$available_n, " available\n"
  )
}

summary_text <- paste0(summary_text, "\n")

if (nrow(best_models) > 0) {
  summary_text <- paste0(summary_text,
    "BEST MODELS BY PHENOTYPE\n",
    "------------------------\n"
  )
  
  for (i in 1:nrow(best_models)) {
    row <- best_models[i, ]
    summary_text <- paste0(summary_text,
      row$phenotype, ": ", row$best_model, " (BIC: ", round(row$best_BIC, 2), ")\n"
    )
  }
  
  summary_text <- paste0(summary_text, "\n")
}

if (nrow(significant_effects) > 0) {
  summary_text <- paste0(summary_text,
    "SIGNIFICANT TREATMENT EFFECTS (p < 0.05)\n",
    "---------------------------------------\n"
  )
  
  for (i in 1:nrow(significant_effects)) {
    row <- significant_effects[i, ]
    summary_text <- paste0(summary_text,
      row$phenotype, " (", row$donor, "): ", round(row$estimate, 3),
      " ± ", round(row$SE, 3), " (p = ", format.pval(row$p.value, digits = 3), ")\n"
    )
  }
} else {
  summary_text <- paste0(summary_text,
    "TREATMENT EFFECTS\n",
    "-----------------\n",
    "No significant inv4m effects detected at p < 0.05\n"
  )
}

writeLines(summary_text, file.path(output_dir, "comprehensive_analysis_summary.txt"))

cat("Results exported to:", output_dir, "\n")
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("This script has reproduced the comprehensive spatial analysis\n")
cat("from docs/inv4m_field_modelling.Rmd as a standalone R script.\n")