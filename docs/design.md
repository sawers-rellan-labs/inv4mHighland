# Comprehensive Multi-Experiment Design

## Executive Summary

This repository contains analysis code for a **comprehensive multi-experiment research program** investigating the **inv4m chromosomal inversion** across different environments, genetic backgrounds, and developmental stages. The central focus is identifying **consistent inv4m responses** across conditions rather than condition-specific effects, which fundamentally differs from previous mashr-based approaches like Crow 2020.

## Experimental Design Overview

### Current Completed/Ongoing Experiments

#### 1. **PSU Field Experiment 2022** (Central/Core Experiment)
- **Location**: Rock Springs, Pennsylvania
- **Design**: Complete Randomized Block Design (CRBD) with 16 replicates
- **Treatments**: 
  - 2 adjacent fields: +P (35 ppm) vs -P (5 ppm) phosphorus
  - 2 genotypes: CTRL vs Inv4m (Near Isogenic Lines)
- **Sampling**: 
  - Field phenotypes
  - RNA-seq: Every other leaf from top leaf with fully developed collar
  - Lipidomics: Same tissue sampling strategy
- **Rationale**: This experiment serves as the core test for the effects of *Inv4m* in a temperate field environment. The `paper.md` findings confirmed that *Inv4m* introgression into a B73 background significantly accelerates flowering time and increases plant height, providing a strong basis for this field validation.
- **Key Feature**: Spatial variation requiring spherical correlation modeling

#### 2. **Growth Chamber Experiment 2024**
- **Conditions**: 3 weeks controlled growth chamber
- **Tissue**: Shoot apical meristem (SAM) dissection
- **Analysis**: Microscopy measurements
- **Purpose**: Developmental/cellular level analysis
- **Rationale**: The  PSU 2022 experiment revealed that *Inv4m* plants have significantly taller and more elongated shoot apical meristems (SAMs). This experiment is designed to investigate the cellular underpinnings of this developmental phenotype by performing detailed microscopy on dissected SAMs.

#### 3. **Clayton Field Experiment 2025** (NC State)
- **Design**: Replicated Latin Square
- **Factors**: 
  - 2 genotypes: CTRL vs Inv4m
  - 2 donors: MI21 vs TMEx (teosinte mexicana backgrounds)
- **Sampling**: Field phenotypes
- **Purpose**: Environmental and genetic background effects
- **Rationale**: To understand the environmental dependency of *Inv4m*'s effects, this experiment compares its performance in a hotter, more humid climate (Clayton) against the milder temperate climate of Pennsylvania. This directly tests whether the accelerated flowering and growth phenotypes observed in `paper.md` are robust across different environments.

#### 4. **Seedling Growth Chamber 2025**
- **Conditions**: 3-day seedlings, growth chamber
- **Measurements**: 
  - Cell proliferation index (root apical meristem)
  - Primary root length
  - Coleoptile length
- **Factors**: 2 genotypes × 2 donors (CTRL/Inv4m × MI21/TMEx)
- **Rationale**: The  PSU 2022 experiment analysis identified a trans-coexpression network modulated by *Inv4m* that is enriched for cell proliferation genes (*pcna2*, *mcm5*). This experiment aims to capture the earliest manifestation of this effect by measuring cell proliferation indices and growth metrics in young seedlings, linking the molecular findings to early developmental phenotypes.

### Planned Future Experiments

#### 5. **Pennsylvania Field Experiment 2025**
- **Purpose**: Environmental comparison with Clayton 2025
- **Rationale**: This experiment serves as a direct environmental comparison to Clayton 2025 and a temporal replicate of PSU 2022. Since in PSU 202 established the baseline effects of *Inv4m* on flowering and height in a temperate setting, this experiment will test the consistency of these effects and their interaction with a different climate.
- **Hypothesis**: Stronger inv4m effects in Pennsylvania climate

#### 6. **Crow 2020 Data Reanalysis**

- **Objective**: Apply your statistical framework to Crow's proliferating tissue data
- **Rationale**: The findings from the  PSU 2022 experiment —that *Inv4m* accelerates flowering, promotes growth, and modulates a cell proliferation network (*pcna2*, *mcm5*) and a key florigen (*zcn26*)—provide a strong, targeted hypothesis for re-examining the Crow 2020 data. The original study collected transcriptomes from highly proliferative tissues (SAM, primary root) under different temperatures, making it an ideal dataset to test if the expression of these specific developmental genes is temperature-dependent.
- **Target genes**: PCNA2, MCM5 (cell cycle/proliferation markers)
- **Tissues**: SAM and primary root
- **Key difference**: Temperature-dependent inv4m effects using your approach

##### Crow et al. 2020 Experimental Design Summary for Reanalysis

The original study aimed to dissect the functional importance of the Inv4m inversion using transcriptomic data. Key aspects of the experimental design are summarized below for reanalysis purposes.

*   **Genetic Material and Crosses:**
    *   **Parental Lines:**
        *   Two highland maize landraces carrying the highland (inverted) Inv4m haplotype:
            1.  **Palomero Toluqueno (PT)** from accession mexi5.
            2.  **Cónico (Mi21)** from accession Michoacán 21.
        *   The B73 maize inbred line, which carries the lowland (standard) Inv4m haplotype.
    *   **Crossing Scheme:**
        *   Each landrace (PT and Mi21) was crossed with B73.
        *   The resulting F1 progeny were backcrossed to B73 for five generations (BC5), selecting for the Inv4m inversion at each stage using a diagnostic CAPS marker.
        *   Two heterozygous BC5 individuals from each landrace background were self-pollinated to generate four BC5S1 families segregating for the Inv4m locus.

*   **Growth Conditions and Treatments:**
    *   **Environments:** Plants were grown in controlled growth chambers under two temperature regimes to simulate different altitudes:
        *   **Warm (Lowland):** 32°C day / 22°C night with a 12-hour light cycle.
        *   **Cold (Highland):** 22°C day / 11°C night with a 12-hour light cycle.
    *   **Replication:** The entire experiment was replicated twice, with growth chambers switched between replicates to account for instrument variation.

*   **Tissue Sampling and RNA Sequencing:**
    *   **Developmental Stages & Tissues for Reanalysis:** The reanalysis will focus on proliferating tissues:
        *   **V3 stage:** Stem apical meristem (SAM)
        *   **V1 & V3 stages:** Primary root
    *   **Original Sampling Strategy:**
        *   A total of nine distinct tissues were collected across the two developmental stages.
        *   Three biological replicates were sampled for each unique combination of `Inv4m donor` x `Inv4m genotype` x `temperature` x `tissue`.
    *   **Sequencing:**
        *   mRNA was isolated and prepared into strand-specific libraries (BRaD-seq protocol).
        *   Libraries were sequenced on an Illumina HiSeq X platform.

*   **Original Data Analysis Pipeline:**
    *   **Alignment:** RNA-seq reads were aligned to the B73 v4 reference genome using HISAT2.
    *   **Quantification:** Gene expression was quantified using Kallisto.
    *   **Differential Expression:** The `limma` R package (with `voom`) was used for differential expression analysis, and `mashr` was used to combine results across conditions to identify any condition-specific effects.
    *   **Genotyping:** Residual landrace alleles in the BC5S1 families were identified by calling SNPs from the RNA-seq data using GATK.

## Critical Statistical Framework Distinction

### Your Approach vs. Crow 2020
- **Crow 2020**: mashr + lsfr to detect **any** condition-specific differential expression
- **Your Approach**: Linear mixed models with spatial correction to detect **consistent** inv4m responses across conditions
- **Core Scripts**: 
  - `PSU_pheno.R`: Spherical spatial correlation modeling
  - `get_DEG_FDR.R`: Standard FDR-based DEG detection (NO mashr, NO lsfr)

### Statistical Rationale
1. **Field experiments** require spatial correlation modeling due to environmental gradients
2. **Consistent effects** across phosphorus conditions are more biologically meaningful than condition-specific responses
3. **Standard FDR approach** more appropriate for detecting robust, replicable effects
4. **Gene co-expression network integration** (72 tissue conditions in B73) for biological validation

## Revised Repository Structure Assessment

### Core Analysis Workflows Identified

#### 1. **Spatial Field Analysis Pipeline** (Primary)
```
Raw field data → Spatial correlation modeling → Mixed-effects models → Consistent treatment effects
```
- **Key script**: `PSU_pheno.R`
- **Method**: Spherical spatial correlation structure
- **Output**: Environment-corrected treatment effects

#### 2. **RNA-seq Consistency Analysis** (Primary)
```
RNA-seq counts → Standard normalization → Linear models → FDR correction → Consistent DEGs
```
- **Key script**: `get_DEG_FDR.R`
- **Method**: Standard limma/edgeR workflow with FDR
- **NOT using**: mashr or lsfr approaches
- **Focus**: Genes with consistent inv4m effects across phosphorus treatments

#### 3. **Multi-Omics Integration with Network Validation**
```
Consistent DEGs → B73 co-expression network → Functional modules → Biological interpretation
```
- **Network data**: 72 tissue conditions in B73 (recurrent parent)
- **Purpose**: Prioritize DEGs based on co-expression patterns
- **Validation**: Cross-reference with known developmental pathways

#### 4. **Cross-Environment Comparison Pipeline** (Planned)
```
Pennsylvania 2025 ↔ Clayton 2025 → Environment × Genotype interactions → Climate adaptation signatures
```

#### 5. **Developmental Axis Analysis** (Integration of experiments 2, 4, 6)
```
Seedling → SAM → Field developmental stages → Proliferation-to-maturation trajectory
```

## Key Methodological Priorities

### High Priority (Core to Your Approach)
1. **Standardize spatial correlation modeling** across all field experiments
2. **Implement consistent FDR-based DEG detection** (avoiding mashr complexity)
3. **Establish cross-condition consistency metrics** for effect size comparison
4. **Integrate B73 co-expression network** for DEG prioritization

### Medium Priority (Cross-Experiment Integration)
1. **Develop environment comparison framework** (Pennsylvania vs Clayton)
2. **Create developmental trajectory analysis** (seedling → SAM → field)
3. **Implement genetic background comparison** (MI21 vs TMEx donors)
4. **Establish temperature-response modeling** for Crow reanalysis

### Unique Analytical Contributions

Your approach makes several important methodological advances:

1. **Spatial-aware field genomics**: Proper spatial correlation in agricultural genomics
2. **Consistency over condition-specificity**: More robust biological inference
3. **Multi-environment validation**: True replication across climates
4. **Network-integrated prioritization**: Leveraging extensive co-expression data
5. **Proliferation-focused reanalysis**: Targeted hypothesis testing with better statistical framework

## Implementation Strategy

### Phase 1: Core Framework (PSU 2022 analysis refinement)
- Standardize `PSU_pheno.R` spatial modeling
- Optimize `get_DEG_FDR.R` consistency detection
- Integrate B73 co-expression network

### Phase 2: Cross-Experiment Framework
- Develop Clayton vs PSU comparison tools
- Create developmental trajectory analysis
- Implement genetic background effect modeling

### Phase 3: Crow Reanalysis
- Apply your framework to proliferating tissues
- Focus on PCNA2/MCM5 temperature responses
- Compare with original mashr results

This summary provides the foundational information needed to apply the new statistical framework to the relevant subsets of the Crow et al. 2020 data.

This is a much more sophisticated and biologically meaningful research program than initially apparent. Your focus on consistent effects and proper spatial modeling represents a significant methodological advance over previous approaches.