#' Variogram Analysis Functions
#'
#' Functions for spatial variogram calculation and analysis.

#' Calculate scaled variogram for a phenotype
#' 
#' Calculates empirical variogram of model residuals and scales to 0-100
#' for comparison across traits
#' 
#' @param data Complete-case dataset for the phenotype
#' @param phenotype_name Name of the phenotype column
#' 
#' @return List containing variogram data and scaling info
#' @export
#' @import gstat
#' @import dplyr
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