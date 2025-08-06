#!/bin/bash

#' Crow 2020 RNA-seq Processing Pipeline
#'
#' This script processes RNA-seq data from the Crow 2020 experiment using Kallisto
#' quantification. It downloads SRA data for the specified samples and performs
#' transcript quantification for downstream differential expression analysis.
#'
#' Usage: ./process_crow2020_rnaseq.sh [SAMPLE_LIST_FILE]
#'
#' @author Francisco Rodriguez  
#' @date 2025-08-06

# Configuration ----------------------------------------------------------

# Base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="../../data"
OUTPUT_DIR="results_crow_reanalysis"
SCRATCH_DIR="scratch"

# Reference files (adjust paths as needed)
REF_DIR="~/ref/zea"
KALLISTO_INDEX="${REF_DIR}/Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx"
TRANSCRIPTS_FA="${REF_DIR}/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa"

# SRA project information
SRA_PROJECT="PRJNA640392"
SAMPLE_METADATA="${DATA_DIR}/crow2020_apical_tissue_samples.tab"

# Tool paths (adjust based on system)
KALLISTO="kallisto"
FASTQ_DUMP="fastq-dump"
PREFETCH="prefetch"

# Processing parameters
KALLISTO_THREADS=8
KALLISTO_BOOTSTRAP=100

# Helper functions -------------------------------------------------------

#' Print usage information
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS] [SAMPLE_LIST_FILE]

Process Crow 2020 RNA-seq data using Kallisto quantification.

OPTIONS:
    -h, --help              Show this help message
    -t, --threads N         Number of threads for Kallisto (default: 8)
    -b, --bootstrap N       Bootstrap samples for Kallisto (default: 100)
    -o, --output DIR        Output directory (default: results_crow_reanalysis)
    --skip-download         Skip SRA download (use existing FASTQ files)
    --samples-only          Only list samples, don't process

ARGUMENTS:
    SAMPLE_LIST_FILE        Optional file with specific sample IDs to process
                           (default: use all samples from metadata file)

EXAMPLES:
    $0                                      # Process all samples
    $0 --threads 16 --bootstrap 200        # Use more resources  
    $0 my_samples.txt                      # Process specific samples
    $0 --samples-only                      # List available samples

EOF
}

#' Log message with timestamp
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

#' Create directory if it doesn't exist
create_dir() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        log_message "Created directory: $dir"
    fi
}

#' Validate required tools
validate_tools() {
    local missing_tools=()
    
    # Check SRA toolkit
    if ! command -v "$FASTQ_DUMP" &> /dev/null; then
        missing_tools+=("fastq-dump (SRA toolkit)")
    fi
    
    if ! command -v "$PREFETCH" &> /dev/null; then
        missing_tools+=("prefetch (SRA toolkit)")
    fi
    
    # Check Kallisto
    if ! command -v "$KALLISTO" &> /dev/null; then
        missing_tools+=("kallisto")
    fi
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_message "ERROR: Missing required tools:"
        printf '  - %s\n' "${missing_tools[@]}"
        log_message "Please install missing tools and ensure they are in PATH"
        exit 1
    fi
    
    log_message "All required tools found"
}

#' Validate required files
validate_files() {
    local missing_files=()
    
    # Check sample metadata
    if [[ ! -f "$SAMPLE_METADATA" ]]; then
        missing_files+=("$SAMPLE_METADATA")
    fi
    
    # Check reference files
    if [[ ! -f "$KALLISTO_INDEX" ]]; then
        if [[ ! -f "$TRANSCRIPTS_FA" ]]; then
            missing_files+=("$TRANSCRIPTS_FA (needed to build index)")
        else
            log_message "Kallisto index not found, will build from: $TRANSCRIPTS_FA"
        fi
    fi
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log_message "ERROR: Missing required files:"
        printf '  - %s\n' "${missing_files[@]}"
        exit 1
    fi
    
    log_message "Required files validated"
}

#' Build Kallisto index if needed
build_kallisto_index() {
    if [[ ! -f "$KALLISTO_INDEX" ]]; then
        log_message "Building Kallisto index..."
        log_message "Input: $TRANSCRIPTS_FA"
        log_message "Output: $KALLISTO_INDEX"
        
        $KALLISTO index \
            -i "$KALLISTO_INDEX" \
            "$TRANSCRIPTS_FA"
        
        if [[ $? -eq 0 ]]; then
            log_message "Kallisto index built successfully"
        else
            log_message "ERROR: Failed to build Kallisto index"
            exit 1
        fi
    else
        log_message "Kallisto index found: $KALLISTO_INDEX"
    fi
}

#' Download SRA data for a sample
download_sra_sample() {
    local sample_id="$1"
    local fastq_dir="$2"
    
    log_message "Processing sample: $sample_id"
    
    # Check if FASTQ files already exist
    local r1_file="${fastq_dir}/${sample_id}_1.fastq"
    local r2_file="${fastq_dir}/${sample_id}_2.fastq"
    
    if [[ -f "$r1_file" && -f "$r2_file" ]]; then
        log_message "FASTQ files already exist for $sample_id, skipping download"
        return 0
    fi
    
    # Download SRA file
    log_message "Downloading SRA data for $sample_id"
    if ! $PREFETCH "$sample_id"; then
        log_message "WARNING: Failed to prefetch $sample_id, skipping"
        return 1
    fi
    
    # Convert to FASTQ
    log_message "Converting to FASTQ: $sample_id"
    cd "$fastq_dir" || exit 1
    
    if ! $FASTQ_DUMP --split-files "$sample_id"; then
        log_message "WARNING: Failed to convert $sample_id to FASTQ"
        cd "$SCRIPT_DIR" || exit 1
        return 1
    fi
    
    cd "$SCRIPT_DIR" || exit 1
    
    # Verify FASTQ files were created
    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        log_message "WARNING: FASTQ files not found after conversion for $sample_id"
        return 1
    fi
    
    log_message "Successfully downloaded and converted: $sample_id"
    return 0
}

#' Run Kallisto quantification for a sample
quantify_sample() {
    local sample_id="$1"
    local fastq_dir="$2"
    local output_base_dir="$3"
    
    local r1_file="${fastq_dir}/${sample_id}_1.fastq"
    local r2_file="${fastq_dir}/${sample_id}_2.fastq"
    local output_dir="${output_base_dir}/kallisto/${sample_id}"
    
    # Check if files exist
    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        log_message "WARNING: FASTQ files not found for $sample_id, skipping quantification"
        return 1
    fi
    
    # Check if quantification already done
    if [[ -f "${output_dir}/abundance.tsv" ]]; then
        log_message "Quantification already exists for $sample_id, skipping"
        return 0
    fi
    
    # Create output directory
    create_dir "$output_dir"
    
    # Run Kallisto quantification
    log_message "Running Kallisto quantification for $sample_id"
    log_message "Input: $r1_file, $r2_file"
    log_message "Output: $output_dir"
    
    $KALLISTO quant \
        -i "$KALLISTO_INDEX" \
        -o "$output_dir" \
        -t "$KALLISTO_THREADS" \
        --bootstrap-samples="$KALLISTO_BOOTSTRAP" \
        "$r1_file" "$r2_file"
    
    if [[ $? -eq 0 ]]; then
        log_message "Quantification completed for $sample_id"
        
        # Compress FASTQ files to save space (optional)
        if command -v gzip &> /dev/null; then
            log_message "Compressing FASTQ files for $sample_id"
            gzip "$r1_file" "$r2_file" &
        fi
        
        return 0
    else
        log_message "ERROR: Quantification failed for $sample_id"
        return 1
    fi
}

#' Get list of samples to process
get_sample_list() {
    local sample_file="$1"
    
    if [[ -n "$sample_file" && -f "$sample_file" ]]; then
        # Use provided sample list
        log_message "Using sample list from: $sample_file"
        cat "$sample_file"
    elif [[ -f "$SAMPLE_METADATA" ]]; then
        # Extract from metadata file (skip header, get Run column)
        log_message "Extracting sample list from metadata: $SAMPLE_METADATA"
        tail -n +2 "$SAMPLE_METADATA" | cut -f4
    else
        log_message "ERROR: No sample list available"
        exit 1
    fi
}

#' List available samples
list_samples() {
    log_message "Available samples in $SAMPLE_METADATA:"
    
    if [[ -f "$SAMPLE_METADATA" ]]; then
        echo "Run_ID      Sample_Name     Genotype    Tissue"
        echo "------      -----------     --------    ------"
        tail -n +2 "$SAMPLE_METADATA" | while IFS=$'\t' read -r line_name inv_temp tissue sample run_id; do
            genotype=$(echo "$line_name" | grep -o "NIL" && echo "INV4" || echo "CTRL")
            printf "%-10s  %-12s    %-8s    %-10s\n" "$run_id" "$sample" "$genotype" "$tissue"
        done
        echo ""
        echo "Total samples: $(tail -n +2 "$SAMPLE_METADATA" | wc -l)"
    else
        log_message "ERROR: Sample metadata file not found: $SAMPLE_METADATA"
        exit 1
    fi
}

# Main pipeline ----------------------------------------------------------

main() {
    log_message "Starting Crow 2020 RNA-seq processing pipeline"
    log_message "Script directory: $SCRIPT_DIR"
    log_message "Data directory: $DATA_DIR"
    log_message "Output directory: $OUTPUT_DIR"
    
    # Parse command line arguments
    SKIP_DOWNLOAD=false
    SAMPLES_ONLY=false
    SAMPLE_LIST_FILE=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -t|--threads)
                KALLISTO_THREADS="$2"
                shift 2
                ;;
            -b|--bootstrap)
                KALLISTO_BOOTSTRAP="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --skip-download)
                SKIP_DOWNLOAD=true
                shift
                ;;
            --samples-only)
                SAMPLES_ONLY=true
                shift
                ;;
            -*)
                log_message "ERROR: Unknown option $1"
                show_usage
                exit 1
                ;;
            *)
                SAMPLE_LIST_FILE="$1"
                shift
                ;;
        esac
    done
    
    # List samples only if requested
    if [[ "$SAMPLES_ONLY" == "true" ]]; then
        list_samples
        exit 0
    fi
    
    log_message "Configuration:"
    log_message "  Threads: $KALLISTO_THREADS"
    log_message "  Bootstrap samples: $KALLISTO_BOOTSTRAP"
    log_message "  Skip download: $SKIP_DOWNLOAD"
    
    # Validate environment
    validate_tools
    validate_files
    
    # Create output directories
    create_dir "$OUTPUT_DIR"
    create_dir "$OUTPUT_DIR/kallisto"
    create_dir "$OUTPUT_DIR/logs"
    create_dir "$SCRATCH_DIR"
    
    # Build Kallisto index if needed
    build_kallisto_index
    
    # Get sample list
    samples=($(get_sample_list "$SAMPLE_LIST_FILE"))
    log_message "Found ${#samples[@]} samples to process"
    
    if [[ ${#samples[@]} -eq 0 ]]; then
        log_message "ERROR: No samples found to process"
        exit 1
    fi
    
    # Set up directories
    FASTQ_DIR="$SCRATCH_DIR/fastq"
    create_dir "$FASTQ_DIR"
    
    # Process samples
    successful=0
    failed=0
    
    for sample in "${samples[@]}"; do
        log_message "Processing sample $((successful + failed + 1))/${#samples[@]}: $sample"
        
        # Download SRA data
        if [[ "$SKIP_DOWNLOAD" == "false" ]]; then
            if ! download_sra_sample "$sample" "$FASTQ_DIR"; then
                ((failed++))
                continue
            fi
        fi
        
        # Quantify with Kallisto
        if quantify_sample "$sample" "$FASTQ_DIR" "$OUTPUT_DIR"; then
            ((successful++))
        else
            ((failed++))
        fi
    done
    
    # Summary
    log_message "Pipeline completed"
    log_message "Successful: $successful"
    log_message "Failed: $failed"
    log_message "Total: ${#samples[@]}"
    
    if [[ $successful -gt 0 ]]; then
        log_message "Results saved to: $OUTPUT_DIR/kallisto/"
        log_message "Next step: Run detect_crow2020_consistency_degs.R for differential expression analysis"
    fi
    
    if [[ $failed -gt 0 ]]; then
        log_message "WARNING: $failed samples failed processing"
        exit 1
    fi
}

# Execute main function with all arguments
main "$@"