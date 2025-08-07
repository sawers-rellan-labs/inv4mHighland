# PSU 2022 Field Experiment Analysis

This directory contains scripts for the comprehensive multi-omics analysis of the PSU 2022 field experiment. The experiment was designed to investigate the effects of the Inv4m chromosomal inversion on maize phenotypes under different environmental conditions.

## Scripts Overview

This directory is organized by data type.

### Expression Analysis (`expression/`)
- **`add_locus_label_descriptions.R`**: Adds descriptions to locus labels.
- **`detect_rnaseq_consistent_degs.R`**: Detects consistently differentially expressed genes.
- **`detect_rnaseq_fdr_degs.R`**: Detects differentially expressed genes using an FDR approach.
- **`estimate_phenotype_mashr_effects.R`**: Estimates phenotype effects using mashr.
- **`fastman_plots.R`**: Generates Manhattan plots.
- **`make_gene_set_GO_analysis.R`**: Performs Gene Ontology analysis on gene sets.
- **`make_manhattan_plots.R`**: Creates Manhattan plots.
- **`make_volcano_plot.R`**: Creates volcano plots.
- **`plot_deg_manhattan.R`**: Plots Manhattan plots for DEGs.
- **`run_deg_linear_models.R`**: Runs linear models for differential expression.
- **`run_rnaseq_mashr_leaf_analysis.R`**: Runs a mashr analysis on leaf RNA-seq data.

### Ionomics (`ionomics/`)
- **`analyze_psu_ionome_multivariate.R`**: Performs multivariate analysis of ionome data.
- **`correct_ionome_rack_effects.R`**: Corrects for rack effects in ionome data.
- **`plot_p31_traditional_varieties.R`**: Plots P31 data for traditional varieties.
- **`plot_p31_traditional_varieties_v2.R`**: An updated version of the P31 plotting script.

### Lipids (`lipids/`)
- **`check_internal_standards.R`**: Checks internal standards in lipidomics data.
- **`clean_lipid_data_ms_dial_normalized.R`**: Cleans and normalizes lipidomics data from MS-DIAL.

### Other
- **`common_config.R`**: Contains common configuration settings for the PSU 2022 analysis.

## Configuration and Setup

### Common Configuration (`common_config.R`)

This file is essential for all scripts in this directory. It sets up paths, loads required libraries, and defines common variables.

**Key settings:**
- `BASE_DIR`: The base directory for the project.
- `DATA_DIR`: The directory where data is stored.
- `RESULTS_DIR`: The directory where results are saved.
- `load_libraries()`: A function to load all necessary R packages.

### Required Data

To run these analyses, you will need the following data files (not included in this repository):

- **Expression Data**: Raw RNA-seq data and corresponding metadata.
- **Ionomics Data**: Ionomics measurements and sample information.
- **Lipidomics Data**: Lipidomics profiles and associated metadata.

These files should be placed in the appropriate subdirectories within the `data/` directory.

### Running Analyses

To run any of the analysis scripts, first ensure that the `common_config.R` file is correctly configured and that the required data is in place. Then, you can source the desired script from within an R session.

**Example:**
```R
# Set the working directory to the psu_2022 directory
setwd("scripts/psu_2022")

# Source the common configuration
source("common_config.R")

# Run an analysis
source("expression/detect_rnaseq_consistent_degs.R")
```

## Expected Outputs

The analysis scripts will generate various outputs, including:

- **Tables**: CSV files containing differential expression results, ionomics data, etc.
- **Plots**: PDF and PNG files of Manhattan plots, volcano plots, and other visualizations.
- **R Objects**: RDS files containing processed data and model fits.

These outputs will be saved in the `results/` directory, organized by data type.

## Notes

- The `_OLD` scripts are deprecated and should not be used for new analyses.
- The `mashr` related scripts represent an alternative analysis approach that was explored but is not the primary method used in the final analysis.
- The `ionomics` and `lipids` directories contain specialized scripts for those data types.
- All scripts are designed to be run from the `psu_2022` directory.
