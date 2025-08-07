# inv4mHighland

[![R Package](https://img.shields.io/badge/R%20Package-v0.1.0-blue)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/Documentation-roxygen2-green)](man/)

An R package for comprehensive spatial analysis of the **Inv4m** chromosomal inversion effects in maize varieties adapted to Mexican highlands. The package provides tools for multi-omics spatial field experiments including RNA-seq, lipidomics, and ionomics data analysis with proper spatial correlation modeling.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("your-username/inv4mHighland")
```

## Quick Start

To demonstrate the package's functionality, here is a basic example of how to use the `run_clayton_spatial_analysis` function. Please note that you will need to provide your own data file.

```r
library(inv4mHighland)

# Create a dummy data frame for demonstration purposes
field_data <- data.frame(
  Genotype = rep(c("B73", "Inv4m"), each = 50),
  Rep = rep(1:10, each = 10),
  Row = rep(1:10, 10),
  Col = rep(rep(1:10, each = 10), 1),
  DTA = rnorm(100, 60, 5),
  PH = rnorm(100, 200, 20)
)

# Run the complete Clayton 2025 spatial analysis
# Replace 'field_data' with your actual data
results <- run_clayton_spatial_analysis(
  data = field_data,
  output_dir = "results_clayton_2025"
)

# Or analyze specific phenotypes
results <- run_clayton_spatial_analysis(
  data = field_data,
  phenotypes = c("DTA", "PH"),
  create_plots = TRUE
)
```

## Package Structure

The package is organized into several functional modules, with R scripts located in the `R/` directory:

- **`diagnostics.R`**: Functions for model diagnostics, including `residual_diagnostics()` and `analyze_blups()`.
- **`multi_trait_analysis.R`**: Tools for analyzing multiple traits, such as `analyze_multiple_traits()` and `create_summary_report()`.
- **`spatial_models.R`**: Core functions for fitting spatial models, including `fit_all_models()`, `extract_model_stats()`, and `extract_variance_components()`.
- **`spatial_visualization.R`**: Functions for creating spatial plots, such as `create_spatial_plot()` and `show_spatial_distribution()`.
- **`treatment_effects.R`**: Tools for extracting treatment effects, including `extract_treatment_effects_emmeans()`.
- **`variogram_analysis.R`**: Functions for variogram analysis, such as `calculate_scaled_variogram()` and `create_directional_variograms()`.

## File Organization

```
├── R/                          # Package functions
│   ├── diagnostics.R
│   ├── multi_trait_analysis.R
│   ├── spatial_models.R
│   ├── spatial_visualization.R
│   ├── treatment_effects.R
│   └── variogram_analysis.R
├── man/                        # Documentation (generated)
├── docs/
│   └── notebooks/             # Analysis notebooks
│       └── inv4m_field_modelling.Rmd
├── data/                      # Experimental datasets (add your own)
├── scripts/                   # Legacy analysis scripts
├── inst/                      # Package installation files
├── tests/                     # Unit tests
└── vignettes/                 # Package vignettes
```

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
  url = {https://github.com/your-username/inv4mHighland}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Please read our contributing guidelines and submit issues or pull requests through GitHub.

## Contact

For questions about the package or research, please open an issue on GitHub or contact the maintainers.
