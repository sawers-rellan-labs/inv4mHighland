# Step 1: Download all the necessary run files from the project.
# The SRA toolkit will automatically find the correct files in its cache
# when you run fastq-dump later.
prefetch PRJNA640392

# Step 2: Create a list of the specific run IDs to process from your file.
# This command skips the header row and extracts the 4th column ('Run').
RUN_IDS=$(tail -n +2 crow2020_apical_tissue_samples.tab | cut -f 4)

# Step 3: Loop through each specific run ID and convert it to FASTQ format.
# The --split-files option creates separate files for paired-end reads.
for run_id in $RUN_IDS; do
  echo "Processing $run_id..."
  fastq-dump --split-files "$run_id"
done

echo "✅ All specified runs have been converted to FASTQ."
