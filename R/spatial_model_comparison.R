#' Spatial Model Fitting Functions
#'
#' Functions for fitting hierarchical spatial models and extracting statistics.

#' Fit all six model structures for a given phenotype
#' 
#' Fits a hierarchy of 6 spatial models from simple fixed effects to 
#' complex mixed-effects with spatial correlation
#' 
#' @param data Complete-case dataset for the phenotype
#' @param phenotype_name Name of the phenotype column
#' 
#' @return List of fitted model objects
#' @export
#' @import nlme
#' @import dplyr
fit_all_models <- function(data, phenotype_name) {
  
  formula_str <- paste(phenotype_name, "~ donor * inv4m")
  models <- list()
  
  # Create plot-level data for Model 1
  plot_data <- data %>%
    group_by(plot_id, plot_row, plot_col, block, donor, inv4m) %>%
    summarise(trait_mean = mean(.data[[phenotype_name]], na.rm = TRUE), 
              .groups = 'drop')
  
  # Model 1: Plot means baseline
  cat("    Attempting Model 1...\n")
  model_1 <- tryCatch({
    lm(trait_mean ~ poly(plot_row, 2) + poly(plot_col, 2) + block + donor * inv4m, 
       data = plot_data)
  }, error = function(e) {
    cat("    Model 1 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_1"]] <- model_1
  if (!is.null(model_1)) cat("    Model 1: SUCCESS\n")
  
  # Model 2: Spatial correlation only  
  cat("    Attempting Model 2...\n")
  model_2 <- tryCatch({
    gls(as.formula(formula_str),
        correlation = corSpher(form = ~ x + y | block/plot_id, nugget = TRUE),
        data = data, method = "REML")
  }, error = function(e) {
    cat("    Model 2 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_2"]] <- model_2
  if (!is.null(model_2)) cat("    Model 2: SUCCESS\n")
  
  # Model 3: Plot random effects
  cat("    Attempting Model 3...\n")
  model_3 <- tryCatch({
    lme(as.formula(paste(phenotype_name, "~ donor * inv4m + block")),
        random = ~ 1 | plot_id, data = data, method = "REML")
  }, error = function(e) {
    cat("    Model 3 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_3"]] <- model_3
  if (!is.null(model_3)) cat("    Model 3: SUCCESS\n")
  
  # Model 4: Plot random + gradients
  cat("    Attempting Model 4...\n")
  model_4 <- tryCatch({
    lme(as.formula(paste(phenotype_name, "~ donor * inv4m + block + plot_row + plot_col")),
        random = ~ 1 | plot_id, data = data, method = "REML")
  }, error = function(e) {
    cat("    Model 4 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_4"]] <- model_4
  if (!is.null(model_4)) cat("    Model 4: SUCCESS\n")
  
  # Model 5: Block random + spatial correlation
  cat("    Attempting Model 5...\n")
  model_5 <- tryCatch({
    lme(as.formula(formula_str),
        random = ~ 1 | block,
        correlation = corSpher(form = ~ x + y | block),
        data = data, method = "REML")
  }, error = function(e) {
    cat("    Model 5 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_5"]] <- model_5
  if (!is.null(model_5)) cat("    Model 5: SUCCESS\n")
  
  # Model 6: Full hierarchical + spatial polynomials
  cat("    Attempting Model 6...\n")
  model_6 <- tryCatch({
    lme(as.formula(paste(phenotype_name, "~ donor * inv4m + poly(x, 2) + poly(y, 2)")),
        random = ~ 1 | block/plot_id, data = data, method = "REML")
  }, error = function(e) {
    cat("    Model 6 FAILED:", e$message, "\n")
    return(NULL)
  })
  models[["model_6"]] <- model_6
  if (!is.null(model_6)) cat("    Model 6: SUCCESS\n")
  
  return(models)
}

#' Extract model comparison statistics
#' 
#' Extracts AIC, BIC, and other model comparison statistics
#' 
#' @param models List of fitted models
#' @param phenotype_name Name of the phenotype
#' 
#' @return Data frame with model comparison metrics
#' @export
#' @import dplyr
extract_model_stats <- function(models, phenotype_name) {
  
  cat("    Extracting stats for", phenotype_name, "\n")
  cat("    Models available:", names(models), "\n")
  cat("    Non-null models:", sum(!sapply(models, is.null)), "\n")
  
  # Initialize tibble with proper column structure
  model_stats <- tibble(
    phenotype = character(),
    model = character(),
    AIC = numeric(),
    BIC = numeric(),
    logLik = numeric(),
    converged = logical()
  )
  
  for (i in 1:6) {
    model_name <- paste0("model_", i)
    model <- models[[model_name]]
    
    cat("      Processing", model_name, "...")
    
    # Check if model exists and is a proper model object
    if (!is.null(model)) {
      cat(" exists, class:", class(model)[1])
      tryCatch({
        # Extract stats based on model class
        if (inherits(model, "lm")) {
          # For lm objects
          aic_val <- AIC(model)
          bic_val <- BIC(model)
          ll_val <- as.numeric(logLik(model))
          
          model_stats <- bind_rows(model_stats, tibble(
            phenotype = phenotype_name,
            model = model_name,
            AIC = aic_val,
            BIC = bic_val,
            logLik = ll_val,
            converged = TRUE
          ))
          cat(" - SUCCESS (AIC:", round(aic_val, 1), ", BIC:", round(bic_val, 1), ")\n")
          
        } else if (inherits(model, c("lme", "gls"))) {
          # For lme/gls objects - check if they have proper logLik
          aic_val <- AIC(model)
          bic_val <- BIC(model)
          ll_val <- as.numeric(logLik(model))
          
          model_stats <- bind_rows(model_stats, tibble(
            phenotype = phenotype_name,
            model = model_name,
            AIC = aic_val,
            BIC = bic_val,
            logLik = ll_val,
            converged = TRUE
          ))
          cat(" - SUCCESS (AIC:", round(aic_val, 1), ", BIC:", round(bic_val, 1), ")\n")
          
        } else {
          cat(" - Unknown class:", class(model), "\n")
        }
      }, error = function(e) {
        cat(" - ERROR:", e$message, "\n")
        # Add row with NAs for failed extraction
        model_stats <<- bind_rows(model_stats, tibble(
          phenotype = phenotype_name,
          model = model_name,
          AIC = NA_real_,
          BIC = NA_real_,
          logLik = NA_real_,
          converged = FALSE
        ))
      })
    } else {
      cat(" - NULL (failed to fit)\n")
      # Model is NULL (failed to fit)
      model_stats <- bind_rows(model_stats, tibble(
        phenotype = phenotype_name,
        model = model_name,
        AIC = NA_real_,
        BIC = NA_real_,
        logLik = NA_real_,
        converged = FALSE
      ))
    }
  }
  
  cat("    Final stats for", phenotype_name, ":", nrow(model_stats), "rows\n")
  return(model_stats)
}

#' Create model comparison table
#' 
#' Creates a formatted comparison table with delta AIC/BIC
#' 
#' @param models List of fitted models
#' 
#' @return Data frame with model comparison results
#' @export
create_model_comparison <- function(models) {
  comparison <- data.frame(
    Model = names(models),
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    logLik = sapply(models, logLik),
    df = sapply(models, function(x) attr(logLik(x), "df"))
  )
  
  # Add delta AIC and BIC
  comparison$deltaAIC <- comparison$AIC - min(comparison$AIC)
  comparison$deltaBIC <- comparison$BIC - min(comparison$BIC)
  
  # Sort by AIC
  comparison <- comparison[order(comparison$AIC), ]
  
  return(comparison)
}

#' Safely extract fixed effects
#' 
#' Safely extracts fixed effects from model objects with error handling
#' 
#' @param model Model object
#' 
#' @return Vector of fixed effects or NULL
#' @export
safe_fixef <- function(model) {
  tryCatch({
    if (inherits(model, c("lme", "gls"))) {
      return(fixef(model))
    } else if (inherits(model, "lm")) {
      return(coef(model))
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

#' Safely extract variance-covariance matrix
#' 
#' Safely extracts vcov matrix with error handling
#' 
#' @param model Model object
#' 
#' @return Variance-covariance matrix or NULL
#' @export
safe_vcov <- function(model) {
  tryCatch({
    return(vcov(model))
  }, error = function(e) {
    return(NULL)
  })
}

#' Export spatial model results
#' 
#' Exports spatial analysis results to specified directory
#' 
#' @param results List containing spatial analysis results
#' @param output_dir Directory for saving results
#' @param verbose Whether to print progress messages
#' 
#' @export
export_spatial_results <- function(results, output_dir, verbose = TRUE) {
  if (verbose) cat("=== EXPORTING SPATIAL RESULTS ===\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save model statistics
  if (!is.null(results$model_stats) && nrow(results$model_stats) > 0) {
    write.csv(results$model_stats, file.path(output_dir, "all_model_statistics.csv"), row.names = FALSE)
  }
  
  # Save missing data summary  
  if (!is.null(results$missing_summary)) {
    write.csv(results$missing_summary, file.path(output_dir, "missing_data_summary.csv"), row.names = FALSE)
  }
  
  # Save variogram summary
  if (!is.null(results$variogram_results) && length(results$variogram_results) > 0) {
    vgm_summary <- purrr::map_dfr(results$variogram_results, function(x) {
      tibble(
        phenotype = x$phenotype,
        n_obs = x$n_obs,
        max_semivariance = round(x$max_gamma, 3)
      )
    })
    write.csv(vgm_summary, file.path(output_dir, "variogram_summary.csv"), row.names = FALSE)
  }
  
  if (verbose) {
    cat("Spatial analysis results exported to:", output_dir, "\n")
  }
}