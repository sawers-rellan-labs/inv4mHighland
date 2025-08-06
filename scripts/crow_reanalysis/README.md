# Crow 2020 Reanalysis - Consistency-Focused Differential Expression

This directory contains scripts for reanalyzing the Crow 2020 temperature response dataset using a consistency-focused approach rather than mashr/lsfr methods. The analysis focuses on identifying robust inv4m effects across temperature treatments, particularly targeting PCNA2 and MCM5 temperature responses in proliferating tissues.

## Overview

The Crow 2020 reanalysis validates the consistency-focused statistical framework developed in this project against published mashr results. Rather than detecting condition-specific effects, we focus on genes with consistent inv4m effects across temperature treatments.

## Directory Structure

```
crow_reanalysis/
├── common_config.R                        # Shared configuration and helper functions
├── process_crow2020_rnaseq.sh             # Kallisto-based RNA-seq processing pipeline
├── detect_crow2020_consistency_degs.R     # Main DEG analysis (FDR-based, consistency-focused)
├── check_sra_metadata.R                   # SRA metadata processing and sample selection
├── validate_crow2020_genes.R              # Gene validation against current results
├── legacy_analysis/                       # Original mashr-based approaches
│   ├── mashr_cross_validation_template.R  # mashr cross-validation template
│   └── mashr_tutorial.R                   # mashr methodology tutorial
├── results_crow_reanalysis/               # Analysis outputs
│   ├── plots/                             # Visualizations
│   ├── tables/                            # Result tables
│   ├── qc/                                # Quality control outputs
│   └── kallisto/                          # Quantification results
└── scratch/                               # Temporary files and downloads
```

## Analysis Pipeline

### 1. Data Acquisition and Processing

#### Download and Process SRA Metadata
```r
# Process SRA metadata to identify target samples
source("check_sra_metadata.R")
```

**Output:** `../../data/crow2020_apical_tissue_samples.tab` - Processed sample metadata focusing on proliferating tissues (SAM and root) from NIL lines.

#### RNA-seq Quantification
```bash
# Download SRA data and perform Kallisto quantification
./process_crow2020_rnaseq.sh

# Options:
./process_crow2020_rnaseq.sh --threads 16 --bootstrap 200  # More resources
./process_crow2020_rnaseq.sh --samples-only               # List samples only
./process_crow2020_rnaseq.sh my_samples.txt              # Process specific samples
```

**Requirements:**
- SRA toolkit (fastq-dump, prefetch)
- Kallisto
- Reference transcriptome index

**Output:** Kallisto quantification results in `results_crow_reanalysis/kallisto/`

### 2. Differential Expression Analysis

#### Main Consistency-Focused Analysis
```r
# Run complete DEG analysis
source("detect_crow2020_consistency_degs.R")
```

**Analysis Features:**
- **Consistency-focused approach:** Identifies genes with robust effects across temperature treatments
- **Target gene focus:** Special attention to PCNA2, MCM5, ZmCCT, ZmCCT9
- **Temperature × Genotype interactions:** Tests for inv4m effects that vary with temperature
- **Proliferating tissue focus:** Primary analysis on SAM and root tissues

**Statistical Model:**
```r
# Linear mixed effects model
~ old_line + genotype * temperature
```

**Key Predictors:**
- `genotypeINV4`: Main inv4m effect
- `temperatureoptimal`, `temperaturewarm`: Temperature effects
- `genotypeINV4:temperatureoptimal`, `genotypeINV4:temperaturewarm`: Inv4m × Temperature interactions

### 3. Gene Validation and Comparison

#### Validate Against Published Results
```r
# Compare with Crow 2020 findings
source("validate_crow2020_genes.R")
```

**Validation Features:**
- Cross-reference with published Crow 2020 gene lists
- Identify conservation of effects between studies
- Map target genes to genomic coordinates
- Check for inv4m region enrichment

## Key Configuration

### Target Genes (defined in `common_config.R`)
- **PCNA2**: Cell cycle progression, proliferating cell marker
- **MCM5**: DNA replication, cell cycle control
- **ZmCCT**: Flowering time regulation
- **ZmCCT9**: Circadian clock and flowering

### Analysis Parameters
- **FDR threshold**: 0.05
- **LogFC threshold**: 1.0 (general), 0.5 (temperature effects)
- **Minimum CPM**: 1.0
- **Minimum samples**: 3

### Tissue Classification
- **Proliferating**: SAM, root, leaf base
- **Mature**: Leaf tip, silk, tassel
- **Mixed**: Whole seedling, ear, node

## Expected Outputs

### Tables
- `crow2020_all_de_results.csv`: Complete differential expression results
- `crow2020_significant_degs.csv`: Significant DEGs with classifications
- `crow2020_target_gene_results.csv`: Results for target genes (PCNA2, MCM5, etc.)

### Visualizations
- `crow2020_mds_genotype.pdf`: Sample clustering by genotype
- `crow2020_mds_temperature.pdf`: Sample clustering by temperature
- `crow2020_target_genes_heatmap.pdf`: Expression patterns of target genes

### QC and Intermediate Files
- `crow2020_voom_object.rds`: Normalized expression data
- `crow2020_ebayes_fit.rds`: Fitted statistical model

## Methodological Innovation

### Consistency vs. Condition-Specificity

**Traditional Approach (Crow 2020):**
```
mashr + lsfr → Detect ANY temperature-specific effects
```

**Our Approach:**
```
Linear Mixed Models + FDR → Detect CONSISTENT effects across temperatures
```

### Key Advantages
1. **Robust biological inference**: Focus on reproducible effects
2. **Cross-environment validation**: Effects that replicate across conditions
3. **Statistical clarity**: Interpretable linear model coefficients
4. **Integration-ready**: Compatible with field experiment results

## Usage Examples

### Complete Analysis Pipeline
```bash
# 1. Process metadata (if SRA data available)
R -e "source('check_sra_metadata.R')"

# 2. Download and quantify samples
./process_crow2020_rnaseq.sh --threads 8 --bootstrap 100

# 3. Run differential expression analysis
R -e "source('detect_crow2020_consistency_degs.R')"

# 4. Validate against published results
R -e "source('validate_crow2020_genes.R')"
```

### Testing Without SRA Data
```r
# Test configuration and functions
source("common_config.R")

# The main analysis script includes simulated data for testing
source("detect_crow2020_consistency_degs.R")
```

## Integration with Other Experiments

The crow_reanalysis results integrate with:

1. **PSU 2022 Field Experiment**: Cross-validate inv4m effects
2. **Clayton 2025 Experiment**: Environmental consistency
3. **Chamber Experiments**: Link temperature response to developmental phenotypes

## Requirements

### R Packages
- **Core**: dplyr, ggplot2, readr, tidyr
- **Bioconductor**: edgeR, limma, tximport, rtracklayer, GenomicRanges
- **Optional**: ggtext, ggpubr, corrplot, pheatmap

### External Tools
- **SRA Toolkit**: For downloading sequencing data
- **Kallisto**: For transcript quantification
- **Reference Files**: Maize B73 transcriptome and annotation

### Data Files
- `SraRunInfo.csv`: SRA metadata (downloadable)
- `gc7_sample_submission.csv`: Sample submission metadata
- `Crow2020_table_S2.tab`: Published gene lists (optional)

## Notes

- The analysis gracefully handles missing packages and data files for testing
- Simulated data is generated when real data is not available
- All scripts follow the established naming convention and configuration patterns
- Results are compatible with the broader inv4m analysis framework

This reanalysis provides a methodological validation of the consistency-focused approach and identifies key temperature-responsive genes that may mediate inv4m adaptive effects across different thermal environments.