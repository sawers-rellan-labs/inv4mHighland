#!/bin/bash

#' Kallisto RNA-seq Quantification for Single Sample
#'
#' This script performs RNA-seq quantification using Kallisto for a single sample.
#' It uses the common configuration for consistent paths and settings.
#'
#' Usage: ./kallisto_quantify_sample.sh SAMPLE_NAME
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load common configuration
source ../common_config.sh

# Load required modules (uncomment based on your system)
# conda activate /usr/local/usrapps/maize/sorghum/conda/envs/rnaseq

# Validate arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 S001"
    exit 1
fi

# Get sample name from command line
SAMPLE="$1"

# Log execution details
log_execution "kallisto_quantify_sample.sh" "$SAMPLE"

# Set up file paths
OUTPUT_DIR="${OUTPUT_BASE_DIR}/kallisto/${SAMPLE}"
READ1="${DATA_DIR}/${SAMPLE}${READ1_SUFFIX}"
READ2="${DATA_DIR}/${SAMPLE}${READ2_SUFFIX}"
INDEX_PATH="${REF_DIR}/${KALLISTO_INDEX}"

# Validate required files
validate_files "$READ1" "$READ2" "$INDEX_PATH"

# Create output directory
create_output_dir "$OUTPUT_DIR"

# Check if index exists, create if needed
if [[ ! -f "$INDEX_PATH" ]]; then
    echo "Index file does not exist. Running kallisto index..."
    $KALLISTO_PATH index -i "$INDEX_PATH" "$TRANSCRIPTS_FA"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to create kallisto index"
        exit 1
    fi
else
    echo "Index file exists: $INDEX_PATH"
fi

# Log quantification details
echo "Running kallisto quant for sample: $SAMPLE"
echo "Read 1: $READ1"
echo "Read 2: $READ2" 
echo "Output directory: $OUTPUT_DIR"
echo "Index: $INDEX_PATH"
echo "Threads: $KALLISTO_THREADS"
echo "Bootstrap samples: $KALLISTO_BOOTSTRAP"

# Run kallisto quantification
echo "Starting kallisto quantification..."
$KALLISTO_PATH quant \
    -i "$INDEX_PATH" \
    -o "$OUTPUT_DIR" \
    -t $KALLISTO_THREADS \
    --bootstrap-samples=$KALLISTO_BOOTSTRAP \
    "$READ1" "$READ2"

# Check exit status
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Kallisto quantification completed for sample: $SAMPLE"
    echo "Results saved to: $OUTPUT_DIR"
else
    echo "ERROR: Kallisto quantification failed for sample: $SAMPLE"
    exit 1
fi

echo "Kallisto quantification script completed at $(date)"
echo "============================================="