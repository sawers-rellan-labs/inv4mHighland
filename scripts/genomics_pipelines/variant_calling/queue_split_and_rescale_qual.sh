#!/bin/tcsh

set input_file=$1

# Loop through each line in the file
# (I feel dirty)
foreach sampleid ("`cut -f3 input_file|sort -u`")
    echo $sampleid
# Construct the bsub command to run the alignment job for the current sample
    set  par="-n 1 -W 2:00 -o %J.stdout -e %J.stderr"
    echo $par

    # Submit the bsub command to the job scheduler
    bsub $par  -R 'rusage[mem=1GB]' $par ./apply_BQSR.sh $sampleid
end
