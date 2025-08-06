# inv4mHighland

[![R Package](https://img.shields.io/badge/R%20Package-v0.1.0-blue)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/Documentation-roxygen2-green)](man/)

An R package for comprehensive spatial analysis of the **Inv4m** chromosomal inversion effects in maize varieties adapted to Mexican highlands. The package provides tools for multi-omics spatial field experiments including RNA-seq, lipidomics, and ionomics data analysis with proper spatial correlation modeling.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("username/inv4m-highland")
```

## Quick Start

```r
library(inv4mHighland)

# Run the complete Clayton 2025 spatial analysis
results <- run_clayton_spatial_analysis(
  data_file = "data/CLY25_Inv4m.csv",
  output_dir = "results_clayton_2025"
)

# Or analyze specific phenotypes
results <- run_clayton_spatial_analysis(
  data_file = "data/CLY25_Inv4m.csv",
  phenotypes = c("DTA", "DTS", "PH", "EBA"),
  create_plots = TRUE
)

# Individual function usage
# Calculate scaled variograms
variogram_result <- calculate_scaled_variogram(trait_data, "PH")

# Fit spatial model hierarchy
models <- fit_all_models(trait_data, "DTA")

# Extract treatment effects with emmeans
effects <- extract_treatment_effects_emmeans(models, "model_3", "DTA")

# Create spatial visualizations
spatial_plot <- create_spatial_plot(field_data, PH, "Plant Height", "cm")
```

## Package Structure

The package is organized into several functional modules:

### Core Analysis Functions
- **`run_clayton_spatial_analysis()`** - Main function to execute complete analysis pipeline
- **`fit_all_models()`** - Fits 6 hierarchical spatial model structures  
- **`extract_model_stats()`** - Extracts model comparison statistics
- **`extract_treatment_effects_emmeans()`** - Treatment effects via emmeans

### Spatial Analysis
- **`calculate_scaled_variogram()`** - Scaled variogram analysis across traits
- **`create_directional_variograms()`** - Directional variogram analysis
- **`create_spatial_plot()`** - Spatial distribution visualization

### Model Diagnostics
- **`residual_diagnostics()`** - Comprehensive model diagnostics
- **`extract_variance_components()`** - Variance component analysis
- **`analyze_blups()`** - BLUP analysis for random effects

### Multi-trait Analysis
- **`analyze_multiple_traits()`** - Standardized multi-trait analysis
- **`create_summary_report()`** - Comprehensive analysis reports

## Research Context

This package implements the spatial analysis framework for investigating **Inv4m**, a large chromosomal inversion found in maize varieties adapted to Mexican highlands. Through multi-omics field experiments, we explore how this inversion contributes to enhanced plant development.

### Key Research Findings

- **Inv4m accelerates flowering and increases plant height** regardless of phosphorus availability
- **Transcriptional effects are highly localized** to genes within the inversion  
- **No evidence for phosphorus starvation advantage** - responses independent of inversion
- **Modulates cell proliferation network** including flowering genes `pcna2`, `mcm5`, `jmj2`, `zcn26`

### Methodological Innovations

- **Spatial-aware field genomics**: Proper spatial correlation modeling for agricultural experiments
- **Consistency-focused analysis**: Robust FDR-based approach over mashr/lsfr methods
- **6-model hierarchy**: From simple fixed effects to complex mixed-effects with spatial correlation
- **Multi-environment validation**: Replication across climatic conditions

## File Organization

```
├── R/                          # Package functions
│   ├── clayton_spatial_analysis.R
│   ├── spatial_models.R
│   ├── variogram_analysis.R
│   ├── spatial_visualization.R
│   ├── treatment_effects.R
│   ├── diagnostics.R
│   ├── multi_trait_analysis.R
│   └── utils.R
├── man/                        # Documentation (generated)
├── docs/
│   └── notebooks/             # Analysis notebooks
│       └── inv4m_field_modelling.Rmd
├── data/                      # Experimental datasets  
├── scripts/                   # Legacy analysis scripts
├── inst/                      # Package installation files
├── tests/                     # Unit tests
└── vignettes/                 # Package vignettes
```

## Datasets

The package works with several experimental datasets:

- **CLY25_Inv4m.csv** - Clayton 2025 field phenotype data
- **inv4mRNAseq_gene_sample_exp.csv** - RNA-seq expression matrix
- **PSU_inv4m_ionome_all.csv** - Ionomics data
- **PSU_inv4m_lipids.csv** - Lipidomics data

## Spatial Model Hierarchy

The package implements a systematic 6-model hierarchy:

1. **Model 1**: Plot means baseline (fixed effects on plot averages)
2. **Model 2**: Spatial correlation only (spherical correlation structure)
3. **Model 3**: Plot random effects (plot-to-plot variation)
4. **Model 4**: Plot random + field gradients (row/column trends)
5. **Model 5**: Block random + spatial correlation
6. **Model 6**: Full hierarchical + spatial polynomials

## Dependencies

Key R packages required:
- `nlme` (≥3.1.0) - Mixed-effects models
- `gstat` (≥2.0.0) - Geostatistical analysis  
- `emmeans` (≥1.6.0) - Marginal means and contrasts
- `tidyverse` (≥1.3.0) - Data manipulation and visualization
- `ggplot2` (≥3.3.0) - Plotting
- `dplyr` (≥1.0.0) - Data manipulation

See `DESCRIPTION` file for complete dependency list.

## Usage Examples

### Basic Spatial Analysis
```r
library(inv4mHighland)

# Load data
data <- read.csv("data/CLY25_Inv4m.csv")

# Run analysis for specific traits
results <- run_clayton_spatial_analysis(
  data_file = "data/CLY25_Inv4m.csv", 
  phenotypes = c("DTA", "PH", "BW"),
  verbose = TRUE
)

# Access results
missing_data <- results$missing_summary
variograms <- results$variogram_results
models <- results$all_models
```

### Custom Spatial Modeling
```r
# Fit model hierarchy for a specific trait
trait_data <- prepare_trait_data(data, "PH")
models <- fit_all_models(trait_data, "PH")

# Compare models
model_stats <- extract_model_stats(models, "PH")
best_model <- model_stats[which.min(model_stats$BIC), ]

# Extract treatment effects
effects <- extract_treatment_effects_emmeans(
  models, best_model$model, "PH"
)
```

## Documentation

- Function documentation: `?function_name` or `help(function_name)`
- Package vignettes: `browseVignettes("inv4mHighland")`
- Analysis notebook: `docs/notebooks/inv4m_field_modelling.Rmd`

## Citation

If you use this package in your research, please cite:

```bibtex
@Manual{inv4mhighland,
  title = {inv4mHighland: Spatial Analysis of Inv4m Chromosomal Inversion Effects in Highland Maize},
  author = {Fausto Rodriguez and Highland Maize Research Team},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/username/inv4m-highland}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Please read our contributing guidelines and submit issues or pull requests through GitHub.

## Contact

For questions about the package or research, please open an issue on GitHub or contact the maintainers.