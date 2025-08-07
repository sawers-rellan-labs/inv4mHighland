# Analysis Workflow for Inv4m Adaptive Effects in Maize

## 1. Project Summary

This project investigates the adaptive effects of the *Inv4m* chromosomal inversion in maize using a **multi-experiment research program** with Near Isogenic Lines (NILs) in a B73 background. The core methodological innovation is the application of **consistency-focused statistical approaches** with **spatial correlation modeling** rather than condition-specific mashr/lsfr methods. This approach detects robust, replicable inv4m effects across different environments and conditions.

## 2. Experimental Program Overview

The project follows an **experiment-centric organization** with four main components:

### 2.1. Core Experiments

1. **PSU Field Experiment 2022 (Central/Core Experiment)**
   - **Status**: ✅ Complete analysis pipeline, results published in `docs/paper.md/.tex`
   - **Design**: Complete Randomized Block Design (CRBD) with spatial correlation
   - **Data**: Multi-omics (RNA-seq, ionomics, lipidomics, phenotypics)
   - **Innovation**: Consistency-based DEG detection vs. mashr approaches

2. **Clayton Field Experiment 2025 (Environmental Comparison)**
   - **Status**: 🔄 Active development, spatial analysis implemented
   - **Design**: Replicated Latin Square in North Carolina climate
   - **Input**: `data/CLY25_Inv4m.csv` (✅ Now available in data directory)
   - **Purpose**: Test environmental robustness of inv4m effects

3. **Chamber Experiments (Developmental Analysis)**
   - **Status**: ⏳ Scripts developed, ready for integration
   - **Components**: Seedling analysis, SAM morphometry
   - **Purpose**: Link molecular findings to cellular/developmental phenotypes

4. **Crow 2020 Reanalysis (Framework Validation)**
   - **Status**: ⏳ Framework prepared for temperature-response analysis
   - **Target**: PCNA2, MCM5 temperature responses in proliferating tissues
   - **Purpose**: Validate consistency approach against published mashr results
   - **Notes for HPC Execution**:
     - **Repository Setup**: Clone the entire repository to your dedicated working space on the HPC cluster.
     - **Dependency Verification**: On the HPC, confirm the presence of Kallisto reference files (`Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx` or `Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa`) in `~/ref/zea/`. Also, ensure SRA Toolkit (`prefetch`, `fastq-dump`) and `kallisto` are installed and in your PATH.
     - **Conda Environment**: Activate the appropriate Conda environment containing the necessary tools before running scripts.
     - **Execution Command**: Navigate to `scripts/crow_reanalysis/` and run `./process_crow2020_rnaseq.sh --threads <num_threads> --bootstrap <num_bootstrap_samples>`. Consider submitting this as an LSF job for resource management.

## 3. Core Methodological Framework

### 3.1. Statistical Innovation: Consistency over Condition-Specificity

**Traditional Approach (Crow 2020):**
```
mashr + lsfr → Detect ANY condition-specific effects
```

**Our Approach:**
```
Linear Mixed Models + Spatial Correlation → Detect CONSISTENT effects across conditions
```

### 3.2. Key Statistical Implementation

#### Spatial Field Analysis (Primary Innovation)
```r
# Template from scripts/psu_2022/spatial_modeling/analyze_psu_spatial_phenotypes.R
lme(trait ~ Genotype * Treatment,
    random = ~ 1 | Rep,
    correlation = corSpher(form = ~ Plot_Row + Plot_Column),
    data = field_data)
```

#### RNA-seq Consistency Analysis
```r
# From scripts/psu_2022/expression/detect_rnaseq_consistent_degs.R
# Standard FDR-based approach (NOT mashr/lsfr)
design <- model.matrix(~ Plot_Column + Plot_Row + leaf_tissue + Treatment*Genotype)
voom_result <- voom(y_filtered, design=design)
fit <- lmFit(voom_result)
eb_fit <- eBayes(fit, robust=TRUE)
```

## 4. Current Project Structure and Implementation Status

### 4.1. Experiment-Centric Directory Organization

```
scripts/
├── psu_2022/           # ✅ Complete multi-omics pipeline
│   ├── expression/     # RNA-seq consistency analysis (11 files)
│   ├── spatial_modeling/  # Core spatial innovation (4 files)
│   ├── network/        # B73 co-expression integration (7 files)
│   ├── ionomics/       # Elemental analysis (4 files)
│   └── lipids/         # Metabolomic profiling (11 files)
├── clayton_2025/      # 🔄 Active spatial analysis (5 files + subdir)
│   └── phenotype_spatial_modeling/  # Spatial framework application (6 files)
├── chamber_experiments/  # ⏳ Developmental analysis (3 files)
├── crow_reanalysis/    # ⏳ Framework validation (3 files + subdir)
├── genomics_pipelines/ # ✅ Complete infrastructure (3 subdirs, ~35 files)
│   ├── rnaseq_processing/   # STAR/Kallisto quantification
│   ├── variant_calling/     # Complete GATK pipeline
│   └── structural_variants/ # Inversion characterization
└── annotation/         # ✅ Functional annotation (5 files)
```

### 4.2. Analysis-Plot Integration Implementation

**Successfully Achieved:**
- **Expression analysis** + Manhattan plots in `psu_2022/expression/`
- **Network analysis** + GO/KEGG plots in `psu_2022/network/`
- **Spatial models** + diagnostic plots in spatial_modeling directories
- **Variant calling** + QC plots in `genomics_pipelines/variant_calling/`

### 4.3. Script Naming Convention Status

**✅ Correctly Applied:**
- **Verb-based (Actions)**: `analyze_*`, `detect_*`, `fit_*`, `build_*`, `plot_*`, `validate_*`
- **Noun-based (Libraries)**: `*_helpers.R`, `fastman_plots.R`, `*_template.R`

**Examples:**
```
detect_rnaseq_consistent_degs.R     # Action: DEG detection
fit_cly25_spatial_correlation_models.R  # Action: Model fitting
spatial_correlation_helpers.R       # Library: Helper functions
fastman_plots.R                     # Library: Plotting functions
```

## 5. Core Analysis Workflows

### 5.1. PSU 2022 Multi-Omics Analysis (✅ Complete)

#### 5.1.1. Spatial Phenotype Analysis
- **Script**: `scripts/psu_2022/spatial_modeling/analyze_psu_spatial_phenotypes.R`
- **Input**: `data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv`
- **Method**: nlme-based spatial correlation with spherical correlation structure
- **Output**: Environment-corrected treatment effects, spatial diagnostic plots
- **Status**: Complete implementation, serves as template for other experiments

#### 5.1.2. RNA-seq Consistency Analysis  
- **Script**: `scripts/psu_2022/expression/detect_rnaseq_consistent_degs.R`
- **Input**: `data/inv4mRNAseq_gene_sample_exp.csv`, `data/inv4mRNAseq_metadata.csv`
- **Method**: limma/edgeR with spatial covariates, FDR correction
- **Innovation**: Consistent effects across phosphorus treatments (not condition-specific)
- **Output**: DEG lists, MDS plots, Manhattan plots
- **Status**: Complete implementation with visualization pipeline

#### 5.1.3. Multi-Omics Integration
- **Network Analysis**: `scripts/psu_2022/network/build_netome_coexpression_networks.R`
- **Ionomics**: `scripts/psu_2022/ionomics/analyze_psu_ionome_multivariate.R`
- **Lipidomics**: `scripts/psu_2022/lipids/run_lipidomics_mashr_analysis.R`
- **Status**: All components functional, integrated through consistent statistical framework

### 5.2. Clayton 2025 Environmental Comparison (🔄 Active)

#### 5.2.1. Active Analysis Implementation
- **Primary Script**: `docs/inv4m_field_modelling.Rmd` (R Markdown notebook)
- **Core Model Script**: `scripts/clayton_2025/phenotype_spatial_modeling/fit_cly25_spatial_correlation_models.R`
- **Input**: `data/CLY25_Inv4m.csv` (✅ Now available in data directory)
- **Design**: Replicated Latin Square with spatial correlation
- **Method**: Model comparison (corSpher, corExp, corGaus) with AIC/BIC selection

#### 5.2.2. Spatial Framework Application
```r
# From fit_cly25_spatial_correlation_models.R
lme(trait ~ donor * genotype,
    random = ~ 1 | rep,
    correlation = corSpher(form = ~ X_pos + Y_pos | rep),
    data = flowering_data)
```

#### 5.2.3. Supporting Analysis
- **Layout Visualization**: `visualize_cly25_spatial_layout.R`
- **Phenotype Analysis**: `analyze_cly25_phenotypes.R`
- **Helper Functions**: `spatial_correlation_helpers.R`
- **Status**: Framework implemented, active data analysis in progress

### 5.3. Chamber Experiments Integration (⏳ Ready for Implementation)

#### 5.3.1. Seedling Analysis
- **Script**: `scripts/chamber_experiments/analyze_seedling_root_shoot.R`
- **Purpose**: 3-day seedling growth measurements linking to PCNA2/MCM5 expression
- **Measurements**: Root length, coleoptile length, cell proliferation indices
- **Status**: Analysis script developed, needs integration with field results

#### 5.3.2. SAM Morphometry
- **Script**: `scripts/chamber_experiments/model_sam_morphometry.R`
- **Purpose**: Shoot apical meristem size/elongation from microscopy
- **Connection**: Links transcriptomic cell proliferation signature to developmental phenotype
- **Status**: Morphometric analysis implemented, ready for correlation analysis

### 5.4. Crow 2020 Reanalysis (⏳ Framework Prepared)

#### 5.4.1. Target Gene Validation
- **Script**: `scripts/crow_reanalysis/validate_crow2020_genes.R`
- **Approach**: Apply consistency framework to published temperature-response data
- **Target Genes**: PCNA2, MCM5 (cell cycle/proliferation markers identified in PSU 2022)
- **Tissues**: SAM and primary root (proliferating tissues)
- **Status**: Methodological framework defined, ready for execution

#### 5.4.2. Method Comparison
- **Reference Scripts**: `scripts/crow_reanalysis/legacy_analysis/mashr_tutorial.R`
- **Purpose**: Direct comparison of consistency approach vs. original mashr results
- **Expected Outcome**: Validate temperature-dependent inv4m effects using robust statistics

## 6. Data Management and File Organization

### 6.1. Core Datasets (Complete)

#### PSU 2022 Experiment
- ✅ `data/inv4mRNAseq_gene_sample_exp.csv` - Expression matrix
- ✅ `data/inv4mRNAseq_metadata.csv` - Sample metadata with spatial coordinates
- ✅ `data/PSU_inv4m_ionome_all.csv` - Ionomics measurements
- ✅ `data/PSU_inv4m_lipids.csv` - Lipidomics profiles
- ✅ `data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv` - Field phenotypes

#### Clayton 2025 Experiment  
- ✅ `data/CLY25_Inv4m.csv` - Multi-trait phenotypes with spatial coordinates

#### Supporting Data
- ✅ Reference annotations, gene lists, coexpression networks (16 additional files)

### 6.2. Results Organization
```
results/
├── mds_genotype.png              # Expression analysis QC
├── mds_treatment_leaf.png         # Sample quality visualization  
└── mds_variables.pdf              # Technical factor assessment
```

## 7. Implementation Priorities and Next Steps

### 7.1. Immediate Actions (Complete Active Workflows)

#### High Priority
1. **✅ COMPLETED**: Move CLY25_Inv4m.csv to data directory
2. **✅ COMPLETED**: Execute Clayton 2025 Analysis (Primary: `docs/inv4m_field_modelling.Rmd`, Supporting: `scripts/clayton_2025/phenotype_spatial_modeling/`)
3. **Generate Cross-Environment Results**: Compare PSU vs Clayton inv4m effects
4. **Get Crow 2020 Gene Expression Matrix**: Obtain the gene expression matrix for Crow 2020 data.

### 7.2. Short-Term Goals (Cross-Experiment Integration)

#### Analytical Integration
1. **Execute Crow Reanalysis**: Apply consistency framework to temperature data
   - Script: `scripts/crow_reanalysis/validate_crow2020_genes.R`
   - Focus: PCNA2, MCM5 temperature responses
   
2. **Chamber Experiment Integration**: Connect developmental phenotypes with field results
   - Link SAM morphometry to field plant height phenotypes
   - Correlate seedling growth with field flowering time acceleration

3. **Network Integration**: Systematic co-expression network integration
   - Apply B73 network prioritization across all experiments
   - Cross-validate DEG findings through network modules

#### Documentation and Reproducibility
1. **Update Analysis Documentation**: Align with current implementation
2. **Consolidate Script Versions**: Remove redundant `_FORMATTED` versions  
3. **Create Cross-Experiment Vignettes**: Document integrated analysis workflows
4. **Document Conda Environments**: Specify Conda environments required for HPC cluster scripts.
   To ensure reproducibility and avoid version clashes, we use a minimal Conda environment strategy, separating R packages from command-line bioinformatics tools.

   **Conda Environment Strategy:**

   *   **`inv4m_r_env` (for R and R packages)**:
       - **Purpose**: Houses R and all R packages for analysis scripts (e.g., `clayton_spatial_analysis.R`, `detect_crow2020_consistency_degs.R`, `inv4m_field_modelling.Rmd`).
       - **Rationale**: R packages have specific R version requirements and can conflict if mixed with other language environments or compiled tools.
       - **`inv4m_r_env.yml`** (create this file in your project root):
         ```yaml
         name: inv4m_r_env
         channels:
           - conda-forge
           - bioconda
           - defaults
         dependencies:
           - r-base=4.3.2 # Specify R version to match your development environment
           - r-essentials # Includes many common R packages
           - r-dplyr
           - r-ggplot2
           - r-ape
           - r-emmeans
           - r-gstat
           - r-nlme
           - r-rlang
           - r-tidyr
           - r-readr
           - bioconductor-edger # For edgeR
           - bioconductor-limma # For limma
           - bioconductor-tximport # For tximport
           - bioconductor-rtracklayer # For rtracklayer
           - bioconductor-genomicranges # For GenomicRanges
           - r-vim # For VIM package
           - r-knitr
           - r-rmarkdown
           - r-ggpubr
           - r-ggtext
           - r-viridis # For viridis color palettes
           - r-corrplot # For correlation plots
           - r-pheatmap # For heatmaps
           # Add any other R packages you might use that are not covered by r-essentials
         ```

   *   **`inv4m_bio_tools_env` (for Bioinformatics Command-Line Tools)**:
       - **Purpose**: Contains all command-line bioinformatics tools used in shell scripts (e.g., `process_crow2020_rnaseq.sh`, GATK pipelines).
       - **Rationale**: Tools like SRA Toolkit, Kallisto, GATK, and STAR are often compiled binaries with specific system library dependencies. Isolating them prevents conflicts.
       - **`inv4m_bio_tools_env.yml`** (create this file in your project root):
         ```yaml
         name: inv4m_bio_tools_env
         channels:
           - conda-forge
           - bioconda
           - defaults
         dependencies:
           - sra-tools # Includes prefetch, fastq-dump
           - kallisto
           - gatk # For GATK tools (variant calling, etc.)
           - bwa # If used for alignment
           - samtools # Essential for BAM/SAM manipulation
           - picard # Often used with GATK for BAM processing
           - star # For STAR aligner (RNA-seq alignment)
           # Add any other command-line tools you use (e.g., bedtools, vcftools)
         ```

   **Implementation Steps:**
   1.  **Create the `environment.yml` files**: Copy the content above into `inv4m_r_env.yml` and `inv4m_bio_tools_env.yml` in your project's root directory.
   2.  **Install the environments**: On your HPC cluster (after loading the `conda` module, if necessary), navigate to your project directory and run:
       ```bash
       conda env create -f inv4m_r_env.yml
       conda env create -f inv4m_bio_tools_env.yml
       ```
   3.  **Activate the environments**: Before running any R script, activate the R environment: `conda activate inv4m_r_env`. Before running any shell script that uses bioinformatics tools, activate the bioinformatics tools environment: `conda activate inv4m_bio_tools_env`.

### 7.3. Medium-Term Objectives (Advanced Framework Development)

#### Methodological Advances
1. **Cross-Environment Comparison Framework**: 
   - Systematic PSU ↔ Clayton environment × genotype interactions
   - Climate adaptation signatures
   
2. **Developmental Trajectory Analysis**:
   - Seedling → SAM → Field developmental progression
   - Proliferation-to-maturation molecular trajectories
   
3. **Temperature-Response Modeling**:
   - Crow data integration with field temperature variation
   - Seasonal adaptation signatures

#### Integration and Validation
1. **Multi-Experiment Meta-Analysis**: Effect size consistency across studies
2. **Biological Validation Pipeline**: Systematic gene function validation
3. **Cross-Population Validation**: Apply framework to additional genetic backgrounds

## 8. Methodological Contributions and Impact

### 8.1. Statistical Innovations Implemented
1. **Spatial-Aware Field Genomics**: Proper spatial correlation in agricultural genomics
2. **Consistency over Condition-Specificity**: More robust biological inference
3. **Multi-Environment Validation**: True environmental replication
4. **Integrated Multi-Omics**: Consistent statistical framework across data types

### 8.2. Biological Discoveries Enabled
1. **Consistent Inv4m Effects**: Accelerated flowering and growth across environments
2. **Cell Proliferation Network**: PCNA2, MCM5, and florigen regulation
3. **Developmental Coordination**: SAM morphology linked to field phenotypes
4. **Environmental Robustness**: Adaptive effects maintained across climates

### 8.3. Reproducible Research Framework
- **Experiment-centric organization** facilitates replication
- **Analysis-plot integration** ensures visualization reproducibility
- **Consistent statistical approaches** enable cross-study comparison
- **Complete data processing pipelines** from raw data to results

## 9. Quality Control and Validation

### 9.1. Statistical Validation
- **Spatial autocorrelation testing**: Moran's I for residual analysis
- **Model comparison**: AIC/BIC selection for correlation structures  
- **Effect size estimation**: Confidence intervals and significance testing
- **Multiple testing correction**: FDR control across experiments

### 9.2. Biological Validation
- **Cross-experiment consistency**: Effect replication across studies
- **Network integration**: Co-expression network validation of DEGs
- **Literature validation**: Cross-reference with published findings
- **Independent dataset validation**: Crow reanalysis framework

This analysis workflow provides a comprehensive framework for investigating inv4m adaptive effects through robust, spatially-aware statistical methods with systematic cross-experiment validation. The experiment-centric organization and consistency-focused approach represent significant methodological advances over previous condition-specific approaches.