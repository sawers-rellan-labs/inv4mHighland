#!/bin/bash

#' STAR RNA-seq Alignment for Single Sample
#'
#' This script performs RNA-seq alignment using STAR for a single sample.
#' It uses the common configuration for consistent paths and settings.
#'
#' Usage: ./align_sample.sh SAMPLE_NAME
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load common configuration
source ../common_config.sh

# Load required modules (uncomment based on your system)
# conda activate /usr/local/usrapps/maize/sorghum/conda/envs/rnaseq
# module load STAR/2.7.7a-foss-2020b

# Validate arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 S001"
    exit 1
fi

# Get sample name from command line
SAMPLE="$1"

# Log execution details
log_execution "align_sample.sh" "$SAMPLE"

# Set up file paths
OUTPUT_DIR="${OUTPUT_BASE_DIR}/star/${SAMPLE}"
READ1="${DATA_DIR}/${SAMPLE}${READ1_SUFFIX}"
READ2="${DATA_DIR}/${SAMPLE}${READ2_SUFFIX}"

# Validate required files
validate_files "$READ1" "$READ2" "$STAR_GENOME_DIR" "$STAR_ANNOTATION"

# Create output directory
create_output_dir "$OUTPUT_DIR"

# Log alignment details
echo "Running STAR alignment for sample: $SAMPLE"
echo "Read 1: $READ1"
echo "Read 2: $READ2"
echo "Output directory: $OUTPUT_DIR"
echo "Genome directory: $STAR_GENOME_DIR"
echo "Annotation: $STAR_ANNOTATION"
echo "Threads: $STAR_THREADS"

# Run STAR alignment
echo "Starting STAR alignment..."
$STAR_PATH \
    --runThreadN $STAR_THREADS \
    --genomeDir "$STAR_GENOME_DIR" \
    --sjdbGTFfile "$STAR_ANNOTATION" \
    --readFilesCommand zcat \
    --readFilesIn "$READ1" "$READ2" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUTPUT_DIR}/" \
    --outReadsUnmapped Fastx

# Check exit status
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: STAR alignment completed for sample: $SAMPLE"
    echo "Results saved to: $OUTPUT_DIR"
else
    echo "ERROR: STAR alignment failed for sample: $SAMPLE"
    exit 1
fi

echo "STAR alignment script completed at $(date)"
echo "============================================="