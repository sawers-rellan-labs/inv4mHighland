#!/bin/bash

#' Common Configuration for Genomics Pipeline Shell Scripts
#'
#' This file contains shared configuration settings and paths
#' used across genomics pipeline shell scripts. Source this file
#' to ensure consistent settings.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Base directories -------------------------------------------------------
export BASE_DIR="/rsstu/users/r/rrellan/sara"
export REF_DIR="${BASE_DIR}/ref"
export DATA_DIR="${BASE_DIR}/RNAseq/RellanAlvarez"
export OUTPUT_BASE_DIR="$(pwd)/results_genomics_pipelines"

# Reference genome files -------------------------------------------------
export GENOME_DIR="${REF_DIR}/NAM5_CHR"
export GTF_FILE="${REF_DIR}/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf"
export GFF_FILE="${REF_DIR}/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
export TRANSCRIPTS_FA="${REF_DIR}/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa"
export GENOME_FA="${REF_DIR}/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"

# Kallisto configuration -------------------------------------------------
export KALLISTO_INDEX="Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx"
export KALLISTO_THREADS=8
export KALLISTO_BOOTSTRAP=100

# STAR configuration -----------------------------------------------------
export STAR_THREADS=16
export STAR_GENOME_DIR="${GENOME_DIR}"
export STAR_ANNOTATION="${GTF_FILE}"

# GATK configuration -----------------------------------------------------
export GATK_JAVA_OPTIONS="-Xmx8g"
export PICARD_JAVA_OPTIONS="-Xmx8g"

# File naming patterns ---------------------------------------------------
export READ1_SUFFIX="_R1_001.fastq.gz"
export READ2_SUFFIX="_R2_001.fastq.gz"

# Tool paths (adjust based on system) -----------------------------------
export STAR_PATH="/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR"
export GATK_PATH="gatk"
export PICARD_PATH="picard"
export KALLISTO_PATH="kallisto"

# Environment setup ------------------------------------------------------
# Uncomment and modify based on your system
# conda activate /usr/local/usrapps/maize/sorghum/conda/envs/rnaseq
# module load STAR/2.7.7a-foss-2020b
# module load GATK/4.2.0.0-Java-11

# Helper functions -------------------------------------------------------

#' Create output directory if it doesn't exist
#' @param $1 Directory path to create
create_output_dir() {
    local dir_path="$1"
    if [[ ! -d "$dir_path" ]]; then
        mkdir -p "$dir_path"
        echo "Created directory: $dir_path"
    fi
}

#' Validate required files exist
#' @param $@ List of file paths to validate
validate_files() {
    local missing_files=()
    for file in "$@"; do
        if [[ ! -f "$file" ]]; then
            missing_files+=("$file")
        fi
    done
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        echo "ERROR: Missing required files:"
        printf ' - %s\n' "${missing_files[@]}"
        exit 1
    fi
}

#' Log script execution details
#' @param $1 Script name
#' @param $2 Sample name (optional)
log_execution() {
    local script_name="$1"
    local sample_name="${2:-}"
    
    echo "============================================="
    echo "Script: $script_name"
    echo "Date: $(date)"
    echo "Host: $(hostname)"
    if [[ -n "$sample_name" ]]; then
        echo "Sample: $sample_name"
    fi
    echo "============================================="
}

# Create base output directories
create_output_dir "$OUTPUT_BASE_DIR"
create_output_dir "$OUTPUT_BASE_DIR/logs"

echo "Genomics pipeline configuration loaded"
echo "Base directory: $BASE_DIR"
echo "Reference directory: $REF_DIR"
echo "Output directory: $OUTPUT_BASE_DIR"