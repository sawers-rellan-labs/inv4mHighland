#!/usr/bin/env Rscript

# ============================================================================
# SETUP SCRIPT FOR SPATIAL ANALYSIS
# ============================================================================
# This script ensures all required packages are installed

cat("Checking and installing required packages for spatial analysis...\n\n")

# List of required packages
required_packages <- c(
  # Core analysis packages
  "nlme",      # Mixed models with spatial correlation
  "lme4",      # Alternative mixed models
  "gstat",     # Variogram analysis
  "sp",        # Spatial data classes and methods
  "ape",       # Phylogenetic analysis (for Moran's I)
  
  # Data manipulation
  "dplyr",     # Data wrangling
  "tidyr",     # Data reshaping
  
  # Statistical tools
  "emmeans",   # Estimated marginal means
  "DHARMa",    # Residual diagnostics for mixed models
  "car",       # Regression diagnostics
  "rstatix",   # Statistical tests for ggplot2
  
  # Visualization
  "ggplot2",   # Graphics
  "ggpubr",    # Publication ready plots
  "gridExtra", # Arrange multiple plots
  
  # Reporting
  "knitr",     # R Markdown processing
  "rmarkdown", # R Markdown rendering
  "DT"         # Interactive tables
)

# Check which packages are not installed
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("The following packages need to be installed:\n")
  cat(paste(" -", missing_packages), sep = "\n")
  cat("\nInstalling packages...\n")
  
  # Install missing packages
  install.packages(missing_packages, repos = "https://cran.rstudio.com/")
  
  cat("\nPackage installation complete!\n")
} else {
  cat("All required packages are already installed.\n")
}

# Verify installation
all_installed <- all(required_packages %in% installed.packages()[,"Package"])

if (all_installed) {
  cat("\n✓ Setup complete! All packages are ready.\n")
  cat("\nYou can now run the analysis by:\n")
  cat("1. Opening 'spatial_analysis_inv4m.Rmd' in RStudio and clicking 'Knit'\n")
  cat("2. Or running: source('render_notebook.R')\n")
} else {
  cat("\n✗ Some packages failed to install. Please check the error messages above.\n")
}

# Test data file existence
data_file <- "../../../data/CLY25_Inv4m.csv"
if (file.exists(data_file)) {
  cat("\n✓ Data file found at:", data_file, "\n")
} else {
  cat("\n✗ Warning: Data file not found at:", data_file, "\n")
  cat("  Please ensure the data file is in the correct location.\n")
  cat("  Expected path: data/CLY25_Inv4m.csv in project root\n")
}
