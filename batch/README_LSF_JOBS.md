# LSF Batch Job Submission Guide

This directory contains LSF batch submission scripts for the inv4mRNA genomics pipeline. All scripts use the LSF `bsub` command with appropriate resource allocations for the `sara` queue.

## Available Batch Scripts

### RNA-seq Processing

#### Primary RNA-seq Pipelines
- **`q_star_align_batch.sh`** - Modern STAR alignment pipeline
  - Usage: `./q_star_align_batch.sh SAMPLE_LIST_FILE`
  - Resources: 16 CPUs, 32GB RAM, 4h runtime per sample
  - Uses: `../scripts/genomics_pipelines/rnaseq_processing/star_quantify_sample.sh`

- **`q_kallisto_quant_batch.sh`** - Modern Kallisto quantification pipeline  
  - Usage: `./q_kallisto_quant_batch.sh SAMPLE_LIST_FILE`
  - Resources: 8 CPUs, 8GB RAM, 2h runtime per sample
  - Uses: `../scripts/genomics_pipelines/rnaseq_processing/kallisto_quantify_sample.sh`

#### Legacy RNA-seq Scripts
- **`q_STAR_align.sh`** - Legacy STAR alignment
  - Usage: `./q_STAR_align.sh SAMPLE_LIST_FILE`
  - Resources: 16 CPUs, 32GB RAM, 4h runtime per sample
  - Uses: `./align_sample.sh` (should be in same directory)

- **`q_kal_quant_sample.sh`** - Legacy Kallisto quantification
  - Usage: `./q_kal_quant_sample.sh SAMPLE_LIST_FILE` 
  - Resources: 8 CPUs, 4GB RAM, 2h runtime per sample
  - Uses: `./quant_sample.sh` (should be in same directory)

- **`q_kallisto_idx.sh`** - Kallisto index building
  - Usage: `bsub < q_kallisto_idx.sh`
  - Resources: 1 CPU, 2GB RAM, 20min runtime
  - Direct LSF submission script

### Variant Calling Pipeline

- **`q_variant_calling_pipeline.sh`** - Complete GATK variant calling pipeline
  - Usage: `./q_variant_calling_pipeline.sh SAMPLE_LIST_FILE`
  - Submits dependent jobs: align → add_read_group → mark_duplicates → call_variants
  - Resources: Variable by step (see script for details)

#### Individual Variant Calling Steps
- **`q_add_read_group.sh`** - Add read groups to BAM files
  - Usage: `./q_add_read_group.sh INPUT_FILE`
  - Resources: 1 CPU, 1GB RAM, 1h runtime per sample

- **`q_mark_duplicates.sh`** - Mark PCR duplicates  
  - Usage: `./q_mark_duplicates.sh INPUT_FILE`
  - Resources: 1 CPU, 8GB RAM, 2h runtime per sample

### Structural Variants

- **`q_structural_variants.sh`** - AnchorWave structural variant detection
  - Usage: `./q_structural_variants.sh SAMPLE_LIST_FILE`
  - Submits dependent jobs: anchorwave_part1 → anchorwave_part2 → plotting
  - Resources: Variable by step (see script for details)

### Specialized Analysis

- **`q_crow_reanalysis.sh`** - Crow 2020 data reanalysis  
  - Usage: `bsub < q_crow_reanalysis.sh`
  - Resources: 8 CPUs, 16GB RAM, 12h runtime
  - Direct LSF submission script

## LSF Resource Specifications

All scripts use standardized LSF parameters:

- **Queue**: `sara` (specified in CLAUDE.md)
- **Output**: `%J.stdout` and `%J.stderr` files  
- **Memory**: Specified with `-R 'rusage[mem=XGB]'`
- **CPUs**: Specified with `-n X`
- **Runtime**: Specified with `-W H:MM`
- **Host locality**: `-R 'span[hosts=1]'` for multi-CPU jobs

## Usage Patterns

### Single Job Submission
```bash
# Direct submission (for scripts with #BSUB directives)
bsub < q_kallisto_idx.sh
bsub < q_crow_reanalysis.sh

# Interactive submission
./q_star_align_batch.sh my_samples.txt
./q_variant_calling_pipeline.sh my_samples.txt
```

### Sample List Format
Most scripts expect a sample list file with one sample ID per line:
```
Sample_001
Sample_002  
Sample_003
```

### Checking Job Status
```bash
# Check all your jobs
bjobs

# Check specific job
bjobs JOBID

# Check job output
cat JOBID.stdout
cat JOBID.stderr
```

## Best Practices

1. **Test with small sample sets** before submitting large batches
2. **Check resource requirements** match your data size
3. **Monitor job outputs** for errors or resource issues
4. **Use appropriate queues** - all scripts default to `sara` queue
5. **Validate input files** exist before submission

## Troubleshooting

- **"Command not found"**: Ensure scripts are executable (`chmod +x *.sh`)
- **"File not found"**: Check that referenced scripts exist in correct locations
- **Resource limits**: Adjust memory/CPU requirements in scripts if jobs fail
- **Path issues**: Scripts assume specific directory structure - check paths

## Script Dependencies

Scripts in this directory reference other scripts in the repository:
- `../scripts/genomics_pipelines/` - Modern pipeline scripts
- `./` - Legacy scripts (should be in same directory as batch scripts)

Ensure all referenced scripts exist and are executable before job submission.