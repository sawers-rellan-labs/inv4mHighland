# ============================================================================
# COMPREHENSIVE SPATIAL ANALYSIS FUNCTIONS FOR INV4M EXPERIMENT
# ============================================================================
# This script contains all helper functions needed to reproduce the comprehensive 
# spatial analysis from docs/inv4m_field_modelling.Rmd as a standalone R script.
# Includes functions for variogram analysis, model fitting, emmeans analysis, and visualization.

# Function to create spatial heat maps
plot_trait_spatial <- function(data, trait, title, midpoint = NULL) {
  if (is.null(midpoint)) {
    midpoint <- median(data[[trait]], na.rm = TRUE)
  }
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(color = get(trait)), size = 2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = midpoint) +
    facet_wrap(~block, labeller = labeller(block = function(x) paste("Block", x))) +
    labs(title = title, color = trait, x = "X Position", y = "Y Position") +
    coord_equal() +
    theme_minimal()
}

# Function to extract treatment effects from any model
extract_treatment_effects <- function(model, model_name = NULL) {
  fe <- fixef(model)
  se <- sqrt(diag(vcov(model)))
  
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

# Function to create directional variograms
create_directional_variograms <- function(data, response, angles = c(0, 45, 90, 135)) {
  sp_data <- data
  coordinates(sp_data) <- ~x+y
  
  vario_dir <- variogram(as.formula(paste(response, "~ 1")), 
                         data = sp_data,
                         alpha = angles,
                         cutoff = 30)
  
  return(vario_dir)
}

# Function to fit model hierarchy
fit_model_hierarchy <- function(data, response = "DTA", verbose = TRUE) {
  
  models <- list()
  
  # Base formula components
  base_formula <- paste(response, "~ donor * inv4m")
  
  # Model 0: Null
  if (verbose) cat("Fitting Null model...\n")
  models$Null <- gls(as.formula(paste(response, "~ 1")), 
                     data = data, method = "REML")
  
  # Model 1: Treatment only
  if (verbose) cat("Fitting Treatment model...\n")
  models$Treatment <- gls(as.formula(base_formula), 
                          data = data, method = "REML")
  
  # Model 2: Treatment + Block
  if (verbose) cat("Fitting Treatment + Block model...\n")
  models$`Treatment + Block` <- gls(as.formula(paste(base_formula, "+ block")), 
                                    data = data, method = "REML")
  
  # Model 3: Treatment + Block + Trends
  if (verbose) cat("Fitting Treatment + Block + Trends model...\n")
  models$`Treatment + Block + Trends` <- gls(
    as.formula(paste(base_formula, "+ block + poly(x_c, 2) + poly(y_c, 2)")), 
    data = data, method = "REML"
  )
  
  # Model 4: Add plot random effect
  if (verbose) cat("Fitting model with Plot random effects...\n")
  models$`T + B + Trends + Plot RE` <- lme(
    as.formula(paste(base_formula, "+ block + poly(x_c, 2) + poly(y_c, 2)")),
    random = ~ 1 | plot_id,
    data = data, method = "REML"
  )
  
  # Model 5: Spatial only
  if (verbose) cat("Fitting Spatial correlation model...\n")
  models$`T + B + Trends + Spatial` <- tryCatch({
    gls(as.formula(paste(base_formula, "+ block + poly(x_c, 2) + poly(y_c, 2)")),
        correlation = corSpher(form = ~ x + y | block, nugget = TRUE),
        data = data, method = "REML")
  }, error = function(e) {
    if (verbose) cat("  Spatial model failed:", e$message, "\n")
    NULL
  })
  
  # Model 6: Plot RE + Spatial
  if (verbose) cat("Fitting combined Plot RE + Spatial model...\n")
  models$`T + B + Trends + Plot RE + Spatial` <- tryCatch({
    lme(as.formula(paste(base_formula, "+ block + poly(x_c, 2) + poly(y_c, 2)")),
        random = ~ 1 | plot_id,
        correlation = corSpher(form = ~ x + y | block, nugget = TRUE),
        data = data, method = "REML",
        control = lmeControl(opt = "optim", maxIter = 200))
  }, error = function(e) {
    if (verbose) cat("  Combined model failed:", e$message, "\n")
    NULL
  })
  
  # Remove NULL models
  models <- models[!sapply(models, is.null)]
  
  return(models)
}

# Function to create model comparison table
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

# Function to extract variance components
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

# Function to perform BLUP analysis
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

# Function for comprehensive residual diagnostics
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

# Function to perform multi-trait analysis
analyze_multiple_traits <- function(data, traits, model_type = "plot_re") {
  results <- list()
  
  for (trait in traits) {
    # Skip if insufficient data
    if (sum(!is.na(data[[trait]])) < 100) {
      cat("Skipping", trait, "- insufficient data\n")
      next
    }
    
    cat("\nAnalyzing", trait, "...\n")
    
    # Fit model based on type
    if (model_type == "plot_re") {
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
    }
    
    if (!is.null(mod)) {
      results[[trait]] <- list(
        model = mod,
        effects = extract_treatment_effects(mod, trait),
        emmeans = emmeans(mod, ~ donor * inv4m),
        contrasts = pairs(emmeans(mod, ~ donor * inv4m), adjust = "tukey")
      )
    }
  }
  
  return(results)
}

# Function to create summary report
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
  for (i in 1:nrow(effects)) {
    summary_text <- paste0(summary_text,
      effects$effect[i], ": ", round(effects$estimate[i], 3), 
      " (SE = ", round(effects$se[i], 3), ", p = ", 
      format.pval(effects$p_value[i], digits = 3), ")\n"
    )
  }
  
  # Save summary
  writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  
  return(summary_text)
}

# Utility function for creating publication-ready tables
create_publication_table <- function(model, traits = NULL) {
  if (is.null(traits)) {
    # Single model table
    effects <- extract_treatment_effects(model)
    if (!is.null(effects)) {
      effects$CI_lower <- effects$estimate - 1.96 * effects$se
      effects$CI_upper <- effects$estimate + 1.96 * effects$se
      effects$p_value <- format.pval(effects$p_value, digits = 3)
      
      return(effects)
    } else {
      return(NULL)
    }
  } else {
    # Multi-trait table
    # Implementation for multiple traits
    # ... (additional code as needed)
    return(NULL)
  }
}

# ============================================================================
# NEW FUNCTIONS FROM RMD
# ============================================================================

# Function to create spatial plots (updated for Rmd compatibility)
create_spatial_plot <- function(data, col_var, plot_title, legend_name, x_lab = "") {
  col_sym <- ensym(col_var)
  
  data %>%
    filter(!is.na(!!col_sym)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = !!col_sym), size = 1.2, alpha = 0.8) +
    scale_color_distiller(palette = "RdYlGn", direction = 1, name = legend_name) +
    labs(title = plot_title, x = x_lab, y = "Field Y position") +
    theme_classic2(base_size = 10) +
    theme(legend.position = "right") +
    coord_equal()
}

# Function to calculate scaled variogram for a phenotype (from Rmd)
calculate_scaled_variogram <- function(data, phenotype_name) {
  
  # Fit simple model to get residuals
  formula_str <- paste(phenotype_name, "~ donor * inv4m")
  m_simple <- lm(as.formula(formula_str), data = data)
  
  # Add residuals to data
  data$resids <- NA
  data$resids[as.integer(rownames(m_simple$model))] <- residuals(m_simple)
  
  # Create gstat object and calculate variogram
  g_obj <- gstat(id = "resids", 
                 formula = resids ~ 1, 
                 locations = ~x+y,
                 data = data %>% filter(!is.na(resids)))
  
  vgm_data <- variogram(g_obj)
  
  # Scale to 0-100
  max_gamma <- max(vgm_data$gamma)
  vgm_data$gamma_scaled <- (vgm_data$gamma / max_gamma) * 100
  
  return(list(
    variogram = vgm_data,
    max_gamma = max_gamma,
    n_obs = nrow(data),
    phenotype = phenotype_name
  ))
}

# Function to fit all six model structures for a given phenotype (from Rmd)
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

# Function to extract model comparison statistics (from Rmd)
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

# Function to extract treatment effects using emmeans (from Rmd)
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

# Additional utility functions

# Function to check and load required packages
check_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it first."))
    }
  }
}

# Function to safely extract fixed effects
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

# Function to safely extract variance-covariance matrix
safe_vcov <- function(model) {
  tryCatch({
    return(vcov(model))
  }, error = function(e) {
    return(NULL)
  })
}

# Final validation function
validate_analysis_setup <- function() {
  required_packages <- c("tidyverse", "dplyr", "ggplot2", "nlme", 
                         "gstat", "emmeans", "ggpubr", "ggtext", 
                         "VIM", "mgcv", "ape")
  
  cat("Checking required packages...\n")
  
  missing_packages <- character(0)
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    stop(paste("The following packages are missing:", paste(missing_packages, collapse = ", "), 
               "\nPlease install them using: install.packages(c('", 
               paste(missing_packages, collapse = "', '"), "'))"))
  }
  
  cat("All required packages are available.\n")
  return(TRUE)
}
