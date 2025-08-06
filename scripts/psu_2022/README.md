# PSU 2022 Field Experiment Analysis

This directory contains refactored scripts for the comprehensive multi-omics analysis of the PSU 2022 field experiment. The refactoring implements the standardized analytical framework described in `docs/analysis.md`, emphasizing consistency-focused approaches over condition-specific methods.

## Scripts Overview

### Core Spatial Analysis (`spatial_modeling/`)

- **`analyze_psu_spatial_phenotypes.R`** - Main spatial phenotype analysis
  - Implements spherical spatial correlation modeling for field experiments
  - Uses comprehensive helper functions for multi-trait analysis
  - Generates publication-ready visualizations and diagnostic plots
  - **Status**: ✅ Refactored with FORMATTED version

- **`analyze_psu_ionome_spatial.R`** - Spatial ionome analysis  
  - Spatial mixed-effects models for mineral concentration data
  - Handles seed and stover tissue analysis with ratio calculations
  - Comprehensive statistical validation and outlier detection
  - **Status**: ✅ Refactored with FORMATTED version

### RNA-seq Analysis (`expression/`)

- **`detect_rnaseq_fdr_degs.R`** - Main differential expression analysis
  - Standard FDR-based approach (NOT mashr) for consistent inv4m effects
  - Comprehensive quality control with MDS analysis
  - Genomic coordinate integration and cis/trans classification
  - Mahalanobis outlier detection for robust gene identification
  - **Status**: ✅ Refactored with FORMATTED version

- **`detect_rnaseq_consistent_degs.R`** - Consistency-focused DEG detection  
  - Framework for detecting robust effects across conditions
  - Methodological distinction from mashr-based approaches
  - **Status**: ✅ Refactored (framework established)

### Supporting Analysis Scripts

#### Ionomics (`ionomics/`)
- **`analyze_psu_ionome_multivariate.R`** - Multivariate ionome analysis with PCA
- **`correct_ionome_rack_effects.R`** - Technical batch effect correction
- **`plot_p31_traditional_varieties*.R`** - Traditional variety comparisons
- **Status**: ✅ Paths updated for consistency

#### Lipidomics (`lipids/`)
- **`run_lipidomics_mashr_analysis.R`** - mashr application to lipid profiles
- **`clean_lipid_data*.R`** - Data processing and normalization pipelines  
- **`lipid_ratio_analysis.R`** - Ratio-based lipid analysis
- **Status**: ✅ Paths updated for consistency

#### Network Analysis (`network/`)
- **`build_netome_coexpression_networks.R`** - B73 co-expression network construction
- **`make_GO_profile_*.R`** - Gene Ontology enrichment analysis
- **`extract_gwas_network_hits.R`** - GWAS-network integration
- **Status**: ✅ Ready for integration

## Key Refactoring Improvements

### 1. Consolidated Redundant Versions
- Replaced original scripts with superior FORMATTED versions
- Removed duplicate code and inconsistent implementations
- Preserved original versions as `*_OLD.R` for reference

### 2. Standardized File Paths
- Updated hardcoded paths to use relative paths from project root
- Created `common_config.R` for shared configuration settings
- Fixed issues with `/Users/fvrodriguez/Desktop/Desktop/` patterns

### 3. Enhanced Documentation  
- Added comprehensive script headers with author and date
- Improved function documentation with parameters and return values
- Consistent commenting and code organization

### 4. Methodological Consistency
- Implemented consistent spatial correlation modeling across experiments
- Standardized statistical approaches (FDR over mashr for robustness)
- Unified diagnostic and visualization approaches

### 5. Improved Error Handling
- Added file validation and error checking
- Graceful handling of missing optional files (e.g., PANNZER annotations)
- Better progress reporting and logging

## Configuration and Setup

### Common Configuration (`common_config.R`)
```r
# Shared settings for all PSU 2022 scripts
source("common_config.R")

# Key variables available:
# - DATA_DIR: "../../data"  
# - OUTPUT_DIR: "results_psu_2022"
# - PHENOTYPE_FILES: List of common data files
# - ANALYSIS_PARAMS: Standard thresholds and parameters
```

### Required Data Files
All scripts expect data files in `../../data/` relative to script location:
- `inv4mRNAseq_gene_sample_exp.csv` - RNA-seq count matrix
- `PSU-PHO22_Metadata.csv` - RNA-seq sample metadata
- `22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv` - Field phenotype data
- `22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv` - Ear phenotype data
- `PSU_inv4m_ionome_all.csv` - Ionome measurements
- `PSU_inv4m_lipids.csv` - Lipidomics data
- `gene_symbol.tab` - Gene annotation

### Output Organization
All results are saved to `results_psu_2022/` with standardized naming:
- Model comparison tables: `*_model_comparison.csv`
- Treatment effects: `*_treatment_effects.csv`  
- Visualization plots: `*_plot.pdf`
- Summary statistics: `*_summary_statistics.csv`

## Analysis Workflow

### 1. Spatial Phenotype Analysis
```bash
cd spatial_modeling/
Rscript analyze_psu_spatial_phenotypes.R
```

### 2. RNA-seq Differential Expression  
```bash
cd expression/
Rscript detect_rnaseq_fdr_degs.R
```

### 3. Ionome Analysis
```bash
cd spatial_modeling/
Rscript analyze_psu_ionome_spatial.R
```

### 4. Multi-omics Integration
```bash
cd network/
Rscript build_netome_coexpression_networks.R
```

## Integration with Project Framework

This refactored PSU 2022 analysis integrates with the broader inv4m project:

1. **Consistent with Clayton 2025**: Uses same statistical framework and helper functions
2. **Cross-Experiment Comparison**: Enables systematic comparison of inv4m effects across environments  
3. **Standardized Outputs**: Compatible with meta-analysis and integration workflows
4. **Documentation Standards**: Follows project-wide documentation conventions

## Methodological Innovations

The PSU 2022 analysis implements several key methodological advances:

1. **Spatial-Aware Field Genomics**: Proper spherical correlation modeling for agricultural settings
2. **Consistency-Focused Statistics**: Robust effects across conditions rather than condition-specific responses  
3. **Integrated Multi-Omics**: Consistent analytical framework across RNA-seq, ionomics, and lipidomics
4. **Network-Validated Results**: B73 co-expression integration for biological validation

## Next Steps

1. **Execute Complete Analysis Pipeline**: Run all refactored scripts with current data
2. **Cross-Experiment Integration**: Compare PSU 2022 and Clayton 2025 results using consistent framework
3. **Method Validation**: Validate spatial correlation approach effectiveness
4. **Results Integration**: Incorporate into comprehensive inv4m analysis

This refactored PSU 2022 analysis establishes the foundation for robust, reproducible multi-omics analysis of inv4m adaptive effects in field conditions.