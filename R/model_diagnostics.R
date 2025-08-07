#' Model Diagnostics Functions
#'
#' Functions for creating diagnostic plots, testing model assumptions, and assessing missing data.

#' Comprehensive residual diagnostics
#' 
#' Creates a comprehensive set of diagnostic plots for model validation
#' 
#' @param model Fitted model object
#' @param data Dataset used for model fitting
#' @param trait Name of the trait being analyzed
#' 
#' @return List of diagnostic plots
#' @export
#' @import ggplot2
#' @import dplyr
residual_diagnostics <- function(model, data, trait) {
  # Add fitted values and residuals to data
  data$fitted <- fitted(model)
  data$residuals <- residuals(model, type = "normalized")
  
  # Create diagnostic plots
  plots <- list()
  
  # Residuals vs Fitted
  plots$resid_fitted <- ggplot(data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = paste("Residuals vs Fitted Values -", trait),
         x = "Fitted Values", y = "Normalized Residuals") +
    theme_minimal()
  
  # Q-Q plot
  plots$qq <- ggplot(data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = paste("Q-Q Plot of Residuals -", trait)) +
    theme_minimal()
  
  # Spatial residuals
  plots$spatial <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = residuals), size = 2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    facet_wrap(~block) +
    labs(title = paste("Spatial Distribution of Residuals -", trait),
         color = "Residual") +
    coord_equal() +
    theme_minimal()
  
  # Residuals by treatment (if treatment column exists)
  if ("treatment" %in% names(data)) {
    plots$treatment <- ggplot(data, aes(x = treatment, y = residuals)) +
      geom_boxplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Residuals by Treatment -", trait),
           x = "Treatment", y = "Normalized Residuals") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    plots$treatment <- NULL
  }
  
  return(plots)
}

#' Extract variance components
#' 
#' Extracts variance components from mixed-effects models
#' 
#' @param models List of fitted mixed-effects models
#' 
#' @return Data frame with variance components
#' @export
#' @import nlme
#' @import dplyr
extract_variance_components <- function(models) {
  var_comp_list <- list()
  
  for (mod_name in names(models)) {
    mod <- models[[mod_name]]
    if (inherits(mod, "lme")) {
      vc <- VarCorr(mod)
      
      # Extract variance values
      if ("plot_id" %in% rownames(vc)) {
        var_plot <- as.numeric(vc["plot_id", "Variance"])
        var_resid <- as.numeric(vc["Residual", "Variance"])
        
        var_comp_list[[mod_name]] <- data.frame(
          Model = mod_name,
          Plot_Variance = var_plot,
          Residual_Variance = var_resid,
          Total_Variance = var_plot + var_resid,
          ICC = var_plot / (var_plot + var_resid),
          Percent_Plot = 100 * var_plot / (var_plot + var_resid)
        )
      }
    }
  }
  
  if (length(var_comp_list) > 0) {
    var_comp_df <- do.call(rbind, var_comp_list)
    row.names(var_comp_df) <- NULL
    return(var_comp_df)
  } else {
    return(NULL)
  }
}

#' Assess missing data patterns
#' 
#' Creates comprehensive missing data assessment for multiple phenotypes
#' 
#' @param data Dataset to assess
#' @param phenotypes Vector of phenotype column names
#' @param verbose Whether to print progress messages
#' 
#' @return List containing missing data summary and diagnostic plots
#' @export
#' @import dplyr
#' @import ggplot2
assess_missing_data <- function(data, phenotypes, verbose = TRUE) {
  if (verbose) cat("=== MISSING DATA ASSESSMENT ===\n")
  
  # Create missing data summary
  missing_summary <- data %>%
    select(all_of(phenotypes)) %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "phenotype", values_to = "missing_count") %>%
    mutate(
      total_obs = nrow(data),
      missing_pct = round(missing_count / total_obs * 100, 1),
      available_n = total_obs - missing_count
    )
  
  if (verbose) print(missing_summary)
  
  # Check for systematic missing data by treatment
  if (all(c("donor", "inv4m") %in% names(data))) {
    missing_by_treatment <- data %>%
      select(donor, inv4m, all_of(phenotypes)) %>%
      group_by(donor, inv4m) %>%
      summarise(across(all_of(phenotypes), ~sum(is.na(.))), .groups = 'drop')
    
    if (verbose) {
      cat("\nMissing data counts by treatment combination:\n")
      print(missing_by_treatment)
    }
  } else {
    missing_by_treatment <- NULL
  }
  
  # Spatial distribution of missing data
  missing_spatial_plot <- NULL
  if (all(c("x", "y") %in% names(data))) {
    missing_spatial_plot <- data %>%
      select(x, y, all_of(phenotypes)) %>%
      mutate(any_missing = rowSums(is.na(select(., all_of(phenotypes)))) > 0) %>%
      ggplot(aes(x = x, y = y, color = any_missing)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
      labs(title = "Spatial distribution of missing observations",
           color = "Has missing\ndata") +
      theme_classic() +
      coord_equal()
    
    if (verbose) print(missing_spatial_plot)
  }
  
  return(list(
    summary = missing_summary,
    by_treatment = missing_by_treatment,
    spatial_plot = missing_spatial_plot
  ))
}