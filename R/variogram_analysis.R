#' Variogram Analysis Functions
#'
#' Functions for spatial variogram calculation and analysis.

#' Calculate scaled variogram for one or multiple phenotypes
#' 
#' Calculates empirical variogram of model residuals and scales to 0-100
#' for comparison across traits. Can handle single trait or multiple traits.
#' 
#' @param data Complete-case dataset for the phenotype(s)
#' @param phenotype_names Name of the phenotype column (character) or vector of phenotype names
#' 
#' @return List containing variogram data and scaling info. If multiple phenotypes 
#'         provided, returns named list with results for each successful calculation.
#' @export
#' @import gstat
#' @import dplyr
calculate_scaled_variogram <- function(data, phenotype_names) {
  
  # Handle single phenotype (backward compatibility)
  if (length(phenotype_names) == 1) {
    return(calculate_single_variogram(data, phenotype_names))
  }
  
  # Handle multiple phenotypes
  variogram_results <- list()
  
  cat("=== CALCULATING VARIOGRAMS FOR", length(phenotype_names), "TRAITS ===\n")
  
  for (trait in phenotype_names) {
    cat("Processing variogram for", trait, "...\n")
    
    # Create complete-case dataset for this trait
    trait_data <- data %>%
      filter(!is.na(.data[[trait]]) & 
             !is.na(x) & !is.na(y) & 
             !is.na(donor) & !is.na(inv4m) & 
             !is.na(block))
    
    cat("  Sample size:", nrow(trait_data), "\n")
    
    if (nrow(trait_data) > 10) {  # Minimum observations for variogram
      tryCatch({
        variogram_results[[trait]] <- calculate_single_variogram(trait_data, trait)
        cat("  Successfully calculated variogram\n")
      }, error = function(e) {
        cat("  ERROR calculating variogram:", e$message, "\n")
      })
    } else {
      cat("  Warning: Insufficient data for", trait, "(need > 10 observations)\n")
    }
  }
  
  cat("\nVariograms calculated for:", length(variogram_results), "out of", length(phenotype_names), "traits\n")
  return(variogram_results)
}

#' Internal function to calculate variogram for a single phenotype
#' 
#' @param data Complete-case dataset for the phenotype
#' @param phenotype_name Name of the phenotype column
#' @return List containing variogram data and scaling info
#' @keywords internal
calculate_single_variogram <- function(data, phenotype_name) {
  
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

#' Create directional variograms
#' 
#' Creates directional variograms to assess anisotropy
#' 
#' @param data Spatial data frame
#' @param response Response variable name
#' @param angles Vector of angles for directional analysis
#' 
#' @return Variogram object
#' @export
#' @import gstat
#' @import sp
create_directional_variograms <- function(data, response, angles = c(0, 45, 90, 135)) {
  sp_data <- data
  coordinates(sp_data) <- ~x+y
  
  vario_dir <- variogram(as.formula(paste(response, "~ 1")), 
                         data = sp_data,
                         alpha = angles,
                         cutoff = 30)
  
  return(vario_dir)
}