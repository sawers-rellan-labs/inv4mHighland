#' Treatment Effects Analysis Functions
#'
#' Functions for extracting and analyzing treatment effects from spatial models.

#' Extract treatment effects from any model
#' 
#' Extracts fixed effects related to treatments from model objects
#' 
#' @param model Fitted model object
#' @param model_name Optional name for the model
#' 
#' @return Data frame with treatment effects
#' @export
#' @import nlme
extract_treatment_effects <- function(model, model_name = NULL) {
  fe <- safe_fixef(model)
  if (is.null(fe)) return(NULL)
  
  vcov_mat <- safe_vcov(model)
  if (is.null(vcov_mat)) return(NULL)
  
  se <- sqrt(diag(vcov_mat))
  
  # Get indices for treatment effects
  idx <- grep("donor|inv4m", names(fe))
  
  if (length(idx) > 0) {
    df <- data.frame(
      effect = names(fe)[idx],
      estimate = fe[idx],
      se = se[idx],
      t_value = fe[idx] / se[idx],
      p_value = 2 * pnorm(-abs(fe[idx] / se[idx])),
      row.names = NULL
    )
    
    if (!is.null(model_name)) {
      df$model <- model_name
    }
    
    return(df)
  } else {
    return(NULL)
  }
}

#' Extract treatment effects using emmeans
#' 
#' Uses emmeans to extract marginal means and contrasts for treatment effects
#' 
#' @param models List of fitted models
#' @param best_model_name Name of the best model to use
#' @param phenotype_name Name of the phenotype
#' 
#' @return List containing emmeans results and contrasts
#' @export
#' @import emmeans
#' @import dplyr
extract_treatment_effects_emmeans <- function(models, best_model_name, phenotype_name) {
  
  model <- models[[best_model_name]]
  
  if (is.null(model)) {
    return(NULL)
  }
  
  tryCatch({
    # Calculate estimated marginal means
    emm <- emmeans(model, specs = ~ inv4m | donor)
    
    # Calculate contrasts with CTRL as reference (reverse = TRUE)
    contrasts <- pairs(emm, by = "donor", reverse = TRUE)
    
    return(list(
      phenotype = phenotype_name,
      model_used = best_model_name,
      emmeans = emm,
      contrasts = contrasts,
      emm_summary = summary(emm),
      contrast_summary = summary(contrasts)
    ))
    
  }, error = function(e) {
    cat("Error extracting effects for", phenotype_name, ":", e$message, "\n")
    return(NULL)
  })
}

#' Analyze BLUPs from random effects models
#' 
#' Extracts and analyzes Best Linear Unbiased Predictors
#' 
#' @param model Mixed-effects model with random effects
#' @param plot_data Plot-level data
#' 
#' @return List with BLUP analysis results
#' @export
#' @import nlme
#' @import dplyr
analyze_blups <- function(model, plot_data) {
  # Extract BLUPs
  blups <- ranef(model)
  
  if ("plot_id" %in% names(blups)) {
    plot_blups <- blups$plot_id
    
    # Merge with plot information
    plot_info <- plot_data %>%
      select(plot_id, donor, inv4m, treatment)
    
    plot_info$blup <- plot_blups[match(plot_info$plot_id, rownames(plot_blups)), 1]
    
    # Test for correlation with treatment
    blup_aov <- aov(blup ~ donor * inv4m, data = plot_info)
    
    return(list(
      plot_info = plot_info,
      anova = blup_aov,
      summary = summary(blup_aov)
    ))
  } else {
    return(NULL)
  }
}