#' Applying Univariate Response Spatial Models to Multiple Phenotypes.
#' Functions for analyzing multiple traits simultaneously.

#' 
#' @description This function applies unviariate response spatial mixed models to multiple traits
#' @param data Dataset containing trait data
#' @param traits Vector of trait names to analyze
#' 
#' @return List of analysis results for each trait
#' @export
#' @import dplyr
#' @import nlme
#' @import emmeans
fit_spatial_models_by_trait() <- function(data, traits) {
  results <- list()
  
  for (trait in traits) {
    # Skip if insufficient data
    if (sum(!is.na(data[[trait]])) < 100) {
      cat("Skipping", trait, "- insufficient data\n")
      next
    }
    
    cat("\nAnalyzing", trait, "...\n")
    
    # Fit model (plot_re)
    formula <- as.formula(paste(trait, "~ donor * inv4m + block + poly(x_c, 2) + poly(y_c, 2)"))
    
    mod <- tryCatch({
      lme(formula,
          random = ~ 1 | plot_id,
          data = data,
          method = "REML",
          na.action = na.omit)
    }, error = function(e) {
      cat("  Model failed for", trait, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(mod)) {
      results[[trait]] <- list(
        model = mod,
        effects = extract_treatment_effects(mod, trait),
        emmeans = emmeans(mod, ~ donor * inv4m),
        contrasts = pairs(emmeans(mod, ~ donor * inv4m), adjust = "tukey")
      )
    }
  }
  
  results
}



#' Create summary report
#' 
#' Creates a comprehensive summary report from analysis results
#' 
#' @param models List of fitted models
#' @param var_comp Variance components data frame
#' @param trait_results Multi-trait analysis results
#' @param output_dir Output directory for results
#' 
#' @return Summary text (also saves to file)
#' @export
create_summary_report <- function(models, var_comp, trait_results, output_dir) {
  
  # Best model identification
  model_comp <- create_model_comparison(models)
  best_model_name <- model_comp$Model[1]
  best_model <- models[[best_model_name]]
  
  # Create summary text
  summary_text <- paste0(
    "SPATIAL MIXED MODEL ANALYSIS SUMMARY\n",
    "====================================\n\n",
    "Date: ", Sys.Date(), "\n",
    "Best Model (by AIC): ", best_model_name, "\n",
    "Model AIC: ", round(model_comp$AIC[1], 2), "\n\n",
    
    "VARIANCE COMPONENTS\n",
    "-------------------\n"
  )
  
  if (!is.null(var_comp) && best_model_name %in% var_comp$Model) {
    vc <- var_comp[var_comp$Model == best_model_name, ]
    summary_text <- paste0(summary_text,
      "Plot-level variance: ", round(vc$Plot_Variance, 3), "\n",
      "Residual variance: ", round(vc$Residual_Variance, 3), "\n",
      "ICC: ", round(vc$ICC, 3), " (", round(vc$Percent_Plot, 1), "% at plot level)\n\n"
    )
  }
  
  summary_text <- paste0(summary_text,
    "TREATMENT EFFECTS\n",
    "-----------------\n"
  )
  
  # Add treatment effects
  effects <- extract_treatment_effects(best_model)
  if (!is.null(effects)) {
    for (i in 1:nrow(effects)) {
      summary_text <- paste0(summary_text,
        effects$effect[i], ": ", round(effects$estimate[i], 3), 
        " (SE = ", round(effects$se[i], 3), ", p = ", 
        format.pval(effects$p_value[i], digits = 3), ")\n"
      )
    }
  }
  
  # Save summary
  writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  
  return(summary_text)
}