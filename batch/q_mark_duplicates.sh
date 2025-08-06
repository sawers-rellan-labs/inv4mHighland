#!/bin/tcsh

# LSF Batch Submission Script for Mark Duplicates
# Usage: ./q_mark_duplicates.sh INPUT_FILE
# 
# This script submits mark duplicates jobs using GATK/Picard

set input_file=$1

# Validate input file
if ( ! -f "$input_file" ) then
    echo "ERROR: Input file not found: $input_file"
    echo "Usage: $0 INPUT_FILE"
    exit 1
endif

# Picard JAR file location
set picard="/usr/local/usrapps/maize/gatk-4.5.0.0/picard.jar"

echo "Submitting mark duplicates jobs from: $input_file"

# Loop through each line in the file
foreach line ("`cat $input_file`")
    echo "Processing line: $line"
    # Extract the fields 2 and 3
    set field2 = `echo $line | awk '{print $2}'| sed 's/_L002//g'`
    echo "Sample ID: $field2"
    set field3 = `echo $line | awk '{print $3}'`
    echo "Field 3: $field3"

    # LSF parameters for mark duplicates
    set par="-q sara -n 1 -W 2:00 -o md_${field2}_%J.stdout -e md_${field2}_%J.stderr"

    # Submit the bsub command to the job scheduler
    bsub $par -R 'rusage[mem=8GB]' ./mark_duplicates.sh ${field2} ${field3}
end

echo "All mark duplicates jobs submitted"
