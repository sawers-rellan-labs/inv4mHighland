# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based genomics research project analyzing **Inv4m**, a chromosomal inversion in maize varieties adapted to Mexican highlands. The project uses multi-omics approaches (RNA-seq, lipidomics, ionomics) across multiple field and controlled experiments to understand how this inversion contributes to enhanced plant development.

## Common Development Commands

### R Environment
```bash
# Open R project in RStudio
open inv4mRNA.Rproj

# Run R scripts directly
Rscript scripts/path/to/script.R

# Install required packages (check script headers for dependencies)
R -e "install.packages(c('dplyr', 'ggplot2', 'edgeR', 'limma'))"
```

### Batch Job Management (HPC) with LSF/quiyt

```bash
# Submit jobs to cluster (from batch/ directory)
bsub < q_STAR_align.sh
bsub < q_kallisto_idx.sh
bsub < q_kal_quant_sample.sh

# Check job status
bjobs
```

### Data Processing Pipelines
```bash
# RNA-seq quantification
bash scripts/genomics_pipelines/rnaseq_processing/kallisto_quantify_sample.sh
bash scripts/genomics_pipelines/rnaseq_processing/star_quantify_sample.sh

# Variant calling pipeline (GATK workflow)
bash scripts/genomics_pipelines/variant_calling/align_sample.sh
bash scripts/genomics_pipelines/variant_calling/add_read_group.sh
bash scripts/genomics_pipelines/variant_calling/mark_duplicates.sh
bash scripts/genomics_pipelines/variant_calling/call_variants.sh
```

## Code Architecture

### Core Statistical Framework
The project implements a **consistency-focused approach** that differs from standard mashr/lsfr methods:

1. **Spatial Field Analysis** (`scripts/psu_2022/spatial_modeling/`)
   - Uses spherical spatial correlation modeling for field experiments
   - Accounts for environmental gradients in agricultural settings
   - Key script: `analyze_psu_spatial_phenotypes.R`

2. **RNA-seq Consistency Analysis** (`scripts/psu_2022/expression/`)
   - Standard FDR-based DEG detection (NOT mashr/lsfr)
   - Focuses on genes with consistent Inv4m effects across conditions
   - Key script: `detect_rnaseq_fdr_degs.R`

3. **Multi-Omics Integration** (`scripts/psu_2022/network/`)
   - Integrates results with B73 co-expression networks
   - Uses 72 tissue conditions for biological validation
   - Prioritizes DEGs based on co-expression patterns

### Experimental Structure

#### Primary Experiments
- **`psu_2022/`** - Core field experiment (Pennsylvania, multi-omics)
- **`clayton_2025/`** - Environmental comparison experiment (North Carolina)
- **`chamber_experiments/`** - Controlled environment studies (SAM morphology)
- **`crow_reanalysis/`** - Reanalysis of Crow 2020 data using project's framework

#### Supporting Infrastructure
- **`genomics_pipelines/`** - Bioinformatics workflows (variant calling, RNA-seq, structural variants)
- **`annotation/`** - Functional annotation scripts
- **`batch/`** - HPC cluster job submission scripts

### Script Naming Conventions

**Verb-based naming** (executable analysis scripts):
- `analyze_*` - Statistical analysis
- `detect_*` - Pattern/gene identification  
- `fit_*` - Model fitting
- `plot_*`, `visualize_*` - Create visualizations
- `run_*` - Execute pipelines
- `validate_*` - Verification analysis

**Noun-based naming** (function libraries):
- `*_helpers.R` - Function libraries
- `*_template.R` - Code templates
- `fastman_plots.R` - Plotting function definitions

## Key Data Files

### Input Data
- `data/inv4mRNAseq_gene_sample_exp.csv` - RNA-seq expression matrix
- `data/PSU-PHO22_Metadata.csv` - Sample metadata for RNA-seq
- `data/PSU_inv4m_ionome_all.csv` - Ionomics data
- `data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv` - Field phenotype data

### Configuration
- `config.yaml` - Reference genome and file paths configuration
- `inv4mRNA.Rproj` - R project settings

### Reference Data
- AGPv5 (NAM-5.0) maize reference genome
- Ensembl Plants annotation (release 56)

## Statistical Methodology

### Key Differences from Standard Approaches
1. **Spatial Correction**: Uses spherical spatial correlation for field experiments instead of ignoring spatial structure
2. **Consistency Focus**: Prioritizes genes with consistent effects across treatments rather than condition-specific responses
3. **FDR over mashr**: Uses standard FDR correction for more robust, interpretable results
4. **Network Integration**: Leverages extensive B73 co-expression data for biological validation

### Analysis Workflow
```
Raw field data → Spatial correlation modeling → Mixed-effects models → Consistent treatment effects
RNA-seq counts → Standard normalization → Linear models → FDR correction → Consistent DEGs
Consistent DEGs → B73 co-expression network → Functional modules → Biological interpretation
```

## Important Dependencies

### R Packages
Core analysis packages commonly used across scripts:
- **Bioconductor**: `edgeR`, `limma`, `rtracklayer`, `GenomicRanges`
- **CRAN**: `dplyr`, `ggplot2`, `ggpubr`, `ggtext`

### External Tools
- **STAR**: RNA-seq alignment
- **Kallisto**: RNA-seq quantification  
- **GATK**: Variant calling pipeline
- **AnchorWave**: Structural variant detection

## Development Notes

### File Paths
- Scripts contain hardcoded paths to `/rsstu/users/r/rrellan/sara/` (HPC environment)
- Local development may require path adjustments in scripts
- Reference data stored in `config.yaml`

### HPC Integration
- Most compute-intensive scripts designed for SLURM job submission
- Batch scripts in `batch/` directory use `sbatch` commands
- Queue scripts (`queue_*.sh`) handle job dependencies

### Data Processing Pipeline Order
1. Genome indexing and alignment (`genomics_pipelines/`)
2. Quantification and quality control
3. Statistical analysis (`psu_2022/expression/`, `psu_2022/spatial_modeling/`)
4. Network integration and visualization (`psu_2022/network/`)

This repository represents a comprehensive genomics analysis framework emphasizing spatial-aware statistics and consistency-focused differential expression analysis for agricultural field experiments.