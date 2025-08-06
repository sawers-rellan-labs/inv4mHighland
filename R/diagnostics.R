#' Model Diagnostics Functions
#'
#' Functions for creating diagnostic plots and testing model assumptions.

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