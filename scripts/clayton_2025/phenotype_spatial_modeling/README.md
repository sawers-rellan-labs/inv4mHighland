# Clayton 2025 Phenotype Spatial Modeling

This directory contains refactored scripts for spatial analysis of the Clayton 2025 field experiment. The refactoring implements the standardized spatial analysis framework described in `docs/analysis.md`.

## Scripts Overview

### Core Analysis Scripts

- **`fit_cly25_spatial_correlation_models.R`** - Main spatial modeling script
  - Fits hierarchical model comparison for multiple traits
  - Uses standardized spatial analysis framework from PSU 2022
  - Generates comprehensive diagnostic plots and model comparisons
  - Exports results to `results_clayton_2025/` directory

- **`analyze_cly25_phenotypes.R`** - Phenotype visualization and basic statistics
  - Creates publication-ready trait plots with statistical annotations
  - Generates summary statistics across genotype and donor combinations
  - Exports individual plots and summary tables

### Helper Functions and Setup

- **`spatial_correlation_helpers.R`** - Reusable function library
  - Contains standardized functions for spatial analysis
  - Implements hierarchical model fitting with multiple correlation structures
  - Provides comprehensive diagnostic tools
  - Uses consistent API across experiments

- **`setup_package_dependencies.R`** - Package installation and validation
  - Ensures all required packages are installed
  - Validates data file availability
  - Provides setup instructions

## Key Refactoring Changes

### 1. Removed Redundancy
- Eliminated duplicate `_FORMATTED` version of helper functions
- Consolidated repetitive model fitting code into reusable functions

### 2. Standardized Framework
- Implemented consistent spatial analysis approach from PSU 2022 experiment
- Uses hierarchical model comparison (Null → Treatment → Block → Trends → Random Effects → Spatial)
- Standardized diagnostic plots and model evaluation metrics

### 3. Improved Data Handling
- Fixed hardcoded file paths to use relative paths from project root
- Added data validation and error handling
- Consistent factor level ordering

### 4. Enhanced Documentation
- Added comprehensive function documentation with parameters and return values
- Improved code comments and structure
- Created standardized analysis workflow

### 5. Error Fixes
- Fixed typo in residual plot code (`alPHa` → `alpha`)
- Added input validation to prevent common errors
- Improved error messages and debugging information

## Usage

### Prerequisites
```bash
# From the phenotype_spatial_modeling directory
Rscript setup_package_dependencies.R
```

### Run Complete Analysis
```r
# Load and run spatial correlation analysis
source("fit_cly25_spatial_correlation_models.R")

# Generate phenotype plots and summaries  
source("analyze_cly25_phenotypes.R")
```

### Use Helper Functions
```r
# Load helper functions
source("spatial_correlation_helpers.R")

# Example usage
models <- fit_model_hierarchy(data, response = "PH")
comparison <- create_model_comparison(models)
diagnostics <- residual_diagnostics(models[[1]], data, "PH")
```

## Output Structure

Results are saved to `results_clayton_2025/`:
- `*_model_comparison.csv` - AIC/BIC comparison tables for each trait
- `*_treatment_effects.csv` - Treatment effect estimates and statistics
- `*_plot.pdf` - Individual trait visualization plots
- `trait_summary_statistics.csv` - Summary statistics by treatment
- `analysis_summary.txt` - Comprehensive analysis report

## Integration with Project Framework

This refactored code integrates with the overall inv4m project framework:

1. **Consistent with PSU 2022 Analysis**: Uses the same statistical approach and helper functions
2. **Follows Naming Conventions**: Verb-based script names, noun-based helper functions  
3. **Standardized Outputs**: Compatible with cross-experiment comparison workflows
4. **Documentation Standards**: Follows project documentation and commenting standards

## Next Steps

1. **Cross-Experiment Integration**: Compare Clayton 2025 results with PSU 2022 using standardized framework
2. **Method Validation**: Validate spatial correlation approach effectiveness
3. **Results Integration**: Incorporate findings into comprehensive inv4m analysis

This refactoring establishes Clayton 2025 as a fully integrated component of the inv4m research program's spatial analysis framework.