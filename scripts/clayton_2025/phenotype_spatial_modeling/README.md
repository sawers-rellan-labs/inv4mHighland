# Clayton 2025 Phenotype Spatial Modeling

This directory demonstrates the modernized spatial analysis workflow using the **inv4mHighland R package**. The previous scattered script approach has been replaced with a streamlined, package-based methodology.

## Current Status: Package-Based Workflow

### Active Script

- **`fit_spatial_model_per_plant.R`** - **UPDATED** to use inv4mHighland package
  - Single function call: `run_clayton_spatial_analysis()`
  - Reduced from 700+ lines to ~100 lines
  - Automatic error handling and progress reporting
  - Comprehensive analysis pipeline with all original functionality
  - Exports results to a directory specified in the script (e.g., `results_clayton_2025_package/`)

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

To run the analysis, execute the `fit_spatial_model_per_plant.R` script. This script will use the `inv4mHighland` R package to perform the spatial analysis. Make sure the package is installed and the input data is available in the expected location.
