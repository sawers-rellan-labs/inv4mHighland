#' Utility Functions
#'
#' Helper functions for package operations and validation.

#' Check and load required packages
#' 
#' Validates that all required packages are available
#' 
#' @param packages Character vector of package names
#' 
#' @return Invisibly returns TRUE if all packages available
#' @export
check_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it first."))
    }
  }
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

#' Validate analysis setup
#' 
#' Validates that all required packages and dependencies are available
#' 
#' @return Invisibly returns TRUE if setup is valid
#' @export
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

#' Create publication-ready tables
#' 
#' Creates formatted tables suitable for publication
#' 
#' @param model Single model object
#' @param traits Optional vector of trait names for multi-trait tables
#' 
#' @return Formatted data frame
#' @export
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