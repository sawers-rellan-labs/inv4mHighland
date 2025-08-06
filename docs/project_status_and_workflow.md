# Project Status and Workflow Analysis (v2)

## 1. Executive Summary

This document provides a detailed analysis of the project's analytical workflows. Unlike a simple file inventory, this report synthesizes the logic from the source code to map the flow of data and the sequence of operations. The analysis reveals a set of robust but disconnected workflows with significant opportunities for improvement.

**Key Findings:**
- **Fragmented Workflows:** The analysis is composed of many standalone scripts. Data processing logic, especially for the main `22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv` file, is repeated across numerous scripts, leading to redundancy and potential for inconsistency.
- **Hardcoded Paths:** A critical issue is the widespread use of hardcoded, absolute local paths (e.g., `/Users/fvrodriguez/Desktop/Desktop/...`), which makes the project non-portable and non-reproducible.
- **Implicit Dependencies:** The execution order of scripts is not explicitly defined, relying on the user to run them in the correct sequence.
- **Mixed Methodologies:** The project uses both standard FDR-based DEG analysis and `mashr`-based analysis, but the relationship and rationale for using both are not clearly documented in the code structure.

This detailed analysis forms the basis for a refactoring plan that will centralize data processing, fix pathing issues, and create a more coherent and reproducible analysis pipeline.

## 2. Core Analytical Workflows

### 2.1. PSU 2022 - Spatial Phenotype Analysis

This workflow corrects for spatial variation in the field to get accurate estimates of genotype and treatment effects on plant phenotypes.

**Data Flow:**
1.  **Input:** `data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv`
2.  **Script:** `scripts/psu_2022/spatial_modeling/analyze_psu_spatial_phenotypes.R`
    - **Action:**
        - Loads the raw phenotype data directly.
        - Performs data cleaning and type conversion.
        - Fits a series of linear mixed-effects models using `nlme`.
        - Critically, it compares different spatial correlation structures (`corLin`, `corSpher`, etc.) to find the best fit for the data based on AIC.
        - Estimates the fixed effects for genotype, phosphorus treatment, and their interaction from the best model.
    - **Issue:** This script contains its own data loading and cleaning logic, which is duplicated in many other scripts. The path to the input CSV is hardcoded.
3.  **Output:** The script prints model summaries and effect estimates to the console but does not appear to save a structured output file (e.g., a results table).

### 2.2. PSU 2022 - RNA-seq Differential Expression Analysis

This workflow identifies genes that are differentially expressed in response to genotype and phosphorus treatment. It appears to have two parallel methodologies.

**Workflow A: Standard FDR-based Analysis**

1.  **Input:**
    - Kallisto quantification results (not directly referenced, but implied).
    - `data/inv4mRNAseq_metadata.csv` (Metadata for RNA-seq samples).
2.  **Script:** `scripts/psu_2022/expression/run_deg_linear_models.R`
    - **Action:** Uses the `limma-voom` pipeline to fit linear models for each gene, testing for the effects of genotype, treatment, and their interaction.
    - **Output:** Saves DEG results to `results/limma_voom_results.csv`.
3.  **Script:** `scripts/psu_2022/expression/detect_rnaseq_fdr_degs.R`
    - **Action:** Takes the output from the linear models and applies a False Discovery Rate (FDR) correction to adjust p-values.
    - **Output:** Saves FDR-corrected results.
4.  **Visualization:**
    - `scripts/psu_2022/expression/make_volcano_plot.R` and `plot_deg_manhattan.R` take the corrected DEG tables and generate visualizations.

**Workflow B: `mashr` Analysis**

1.  **Script:** `scripts/psu_2022/expression/run_rnaseq_mashr_leaf_analysis.R`
    - **Action:** Implements a more complex analysis using `mashr` (Multivariate Adaptive Shrinkage) to analyze effects across different leaf samples, likely to identify condition-specific effects.
    - **Output:** Saves `mashr` results objects.

### 2.3. PSU 2022 - Ionomics and Lipidomics Analysis

These workflows follow a similar pattern to the phenotype analysis.

**Ionomics Workflow:**
1.  **Input:** `data/PSU_inv4m_ionome_all.csv` and the main `22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv` for metadata.
2.  **Script:** `scripts/psu_2022/ionomics/correct_ionome_rack_effects.R`
    - **Action:** Performs initial data cleaning and corrects for experimental batch effects (rack effects).
3.  **Script:** `scripts/psu_2022/ionomics/analyze_psu_ionome_multivariate.R`
    - **Action:** Performs multivariate analysis on the corrected ionome data.
    - **Issue:** Both scripts reload and process the main phenotype CSV, duplicating effort. Paths are hardcoded.

**Lipidomics Workflow:**
1.  **Input:** `data/PSU_inv4m_lipids.csv` and `data/inv4mRNAseq_metadata.csv`.
2.  **Scripts:** A series of scripts in `scripts/psu_2022/lipids/` perform cleaning (`clean_lipid_data.R`), dimensionality reduction (`lipid_dimensionality_reduction.R`), and differential abundance analysis (`run_lipidomics_mashr_analysis.R`).
    - **Issue:** Significant code duplication for data loading and cleaning across these scripts.

### 2.4. Clayton 2025 - Spatial Phenotype Analysis

This workflow applies the spatial modeling methodology to a different experiment.

1.  **Input:** `CLY25_Inv4m.csv`
    - **CRITICAL ISSUE:** This file is not in the project. All related scripts reference it from a hardcoded local path: `~/Desktop/CLY25_Inv4m.csv`.
2.  **Script:** `scripts/clayton_2025/phenotype_spatial_modeling/fit_cly25_spatial_correlation_models.R`
    - **Action:** Similar to the PSU spatial analysis, it fits mixed-effects models with spatial correlation structures.
3.  **Helper Script:** `scripts/clayton_2025/phenotype_spatial_modeling/spatial_correlation_helpers.R`
    - **Action:** Contains helper functions used by the main analysis script. This is a good example of separating library code.

## 3. Synthesis and Refactoring Implications

The current structure is a collection of powerful but isolated analyses. The most significant barrier to usability and reproducibility is the decentralized and duplicated data processing and the use of hardcoded paths.

A successful refactoring must prioritize:
1.  **Centralizing Data Processing:** Create a single script (`scripts/0_process_data.R`) that loads all raw data, cleans it, merges it, and saves processed versions to `data/processed/`. All downstream analysis scripts should **only** read from this processed data directory.
2.  **Eliminating Hardcoded Paths:** All file paths must be relative to the project root to ensure portability.
3.  **Defining a Clear Execution Order:** The refactored directory structure and naming convention should make the workflow sequence obvious (e.g., `1.1_...`, `1.2_...`).
4.  **Consolidating Helper Functions:** Move reusable functions (like the spatial correlation helpers) into a central `scripts/lib/` directory.
5.  **Clarifying Methodologies:** The documentation and structure should clarify when and why different statistical approaches (`limma` vs. `mashr`) are used.