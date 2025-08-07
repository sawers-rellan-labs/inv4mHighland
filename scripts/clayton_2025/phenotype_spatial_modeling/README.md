# Clayton 2025 Phenotype Spatial Modeling

This directory demonstrates the modernized spatial analysis workflow using the **inv4mHighland R package**. The previous scattered script approach has been replaced with a streamlined, package-based methodology.

## Current Status: Package-Based Workflow

### Active Script

- **`fit_spatial_model_per_plant.R`** - **UPDATED** to use inv4mHighland package
  - Single function call: `run_clayton_spatial_analysis()`
  - Reduced from 700+ lines to ~100 lines
  - Automatic error handling and progress reporting
  - Comprehensive analysis pipeline with all original functionality
  - Exports results to `results_clayton_2025_package/` directory

### Deprecated Scripts (Functionality Now in Package)

- **`plot_pairwise_comparisons.R.deprecated`** - Replaced by package visualization functions
  - Functionality now available through `run_clayton_spatial_analysis()`
  - Treatment effect plots generated automatically
  - More robust statistical testing and visualization

- **`plot_spatial_layout_CLY25.R.deprecated`** - Replaced by `show_spatial_distribution()`
  - Multi-trait spatial visualization with organized layout
  - Custom trait labeling and consistent formatting
  - Automatic trait selection and display optimization

## Package-Based Analysis Workflow

### Prerequisites

```r
# Install and load the inv4mHighland package
# (Package should be installed from the project root)
library(inv4mHighland)
library(tidyverse)
```

### Run Complete Analysis

```r
# Single function call replaces all previous scripts
results <- run_clayton_spatial_analysis(
  data_file = "../../../data/CLY25_Inv4m.csv",
  output_dir = "results_clayton_2025_package",
  trait_names = c("DTA", "DTS", "LAE", "PH", "EN", "SL", "BL", "BW", "EBA"),
  use_variogram = TRUE,
  create_plots = TRUE,
  export_results = TRUE
)
```

### Advanced Usage (Individual Functions)

```r
# Use package functions independently if needed
library(inv4mHighland)

# Load and clean data
field_data <- load_clayton_data("../../../data/CLY25_Inv4m.csv")

# Calculate variograms for multiple traits
vgm_results <- calculate_scaled_variogram(field_data, c("DTA", "DTS", "PH"))

# Fit spatial models with model comparison
models <- fit_all_models(field_data, "PH")
model_stats <- extract_model_stats(models, "PH")

# Create visualizations
show_spatial_distribution(field_data, c("DTA", "DTS", "PH"))

# Extract treatment effects
effects <- extract_treatment_effects_emmeans(models, "model_6", "PH")
```

## Package Functions Overview

### Core Analysis Functions
- `run_clayton_spatial_analysis()` - Master workflow function
- `fit_all_models()` - Hierarchical model comparison  
- `calculate_scaled_variogram()` - Multi-trait variogram analysis
- `extract_model_stats()` - Model comparison statistics
- `extract_treatment_effects()` - Treatment effect quantification

### Visualization Functions
- `show_spatial_distribution()` - Multi-trait spatial plots (3 per row)
- `create_spatial_plot()` - Individual trait spatial visualization
- `residual_diagnostics()` - Model diagnostic plots

### Utility Functions
- `load_clayton_data()` - Data loading with validation
- `validate_analysis_setup()` - Dependency checking
- `create_model_comparison()` - Model comparison tables
- `create_publication_table()` - Publication-ready results

## Analysis Pipeline

The package implements a comprehensive analysis pipeline:

1. **Data Loading & Cleaning**
   - Automatic factor conversion and coordinate centering
   - EBA calculation (Estimated Blade Area = 0.75 × BL × BW)
   - Missing data assessment and validation

2. **Spatial Analysis**
   - Scaled variogram calculation for spatial autocorrelation assessment
   - Six hierarchical models: Null → Treatment → Block → Trends → Random Effects → Spatial
   - Automatic model selection using BIC

3. **Treatment Effects**
   - Treatment effect extraction using optimal models
   - Emmeans-based contrasts for inv4m vs. control comparisons
   - Multiple testing correction and significance assessment

4. **Visualization & Diagnostics**
   - Spatial distribution plots for all traits
   - Residual diagnostics for model validation
   - Forest plots for treatment effects
   - Comprehensive diagnostic assessments

5. **Results Export**
   - Model comparison tables
   - Treatment effect summaries
   - Missing data reports
   - Comprehensive analysis reports

## Output Structure

Results are saved to `results_clayton_2025_package/`:
- `all_model_statistics.csv` - Complete model comparison results
- `best_models_selection.csv` - Best model for each trait
- `significant_treatment_effects.csv` - Significant treatment effects
- `missing_data_summary.csv` - Missing data patterns
- `variogram_summary.csv` - Spatial autocorrelation patterns
- `comprehensive_analysis_summary.txt` - Human-readable report
- Multiple diagnostic and visualization plots

## Advantages of Package-Based Approach

✅ **Code Reduction**: 700+ lines reduced to ~20 functional lines
✅ **Error Handling**: Automatic validation and robust error management
✅ **Consistency**: Standardized parameters and outputs across analyses
✅ **Documentation**: Complete roxygen2 documentation for all functions
✅ **Modularity**: Functions can be used independently or as complete pipeline
✅ **Reproducibility**: Identical results across different users and systems
✅ **Maintainability**: Updates to methodology automatically available to all scripts
✅ **Extensibility**: Easy adaptation to new experiments and traits

## Integration with Project Framework

This package-based approach provides a template for modernizing other experimental analyses:

1. **Consistent Methodology**: Same spatial analysis approach across all experiments
2. **Standardized Interface**: Common function signatures and parameters
3. **Cross-Experiment Compatibility**: Results can be easily compared across studies
4. **Future Development**: New features automatically available to all analyses

## Migration Notes

**For users familiar with the old scripts:**
- All original functionality is preserved and enhanced
- Results format remains compatible with downstream analyses
- Individual functions can still be called for custom workflows
- Deprecated scripts remain available with `.deprecated` extension for reference

**For new users:**
- Start with `run_clayton_spatial_analysis()` for complete analysis
- Refer to function documentation: `?run_clayton_spatial_analysis`
- Use `show_spatial_distribution()` for quick trait visualization
- Check package vignettes for detailed examples

This modernization establishes Clayton 2025 as the reference implementation for the inv4mHighland spatial analysis framework.