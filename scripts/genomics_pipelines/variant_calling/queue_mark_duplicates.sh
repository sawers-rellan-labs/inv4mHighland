#!/bin/tcsh

set input_file=$1


# Loop through each line in the file
# (I feel dirty)
foreach line ("`cat $input_file`")
    echo $line
    # Extract the fields 2 and 3
    set field2 = `echo $line | awk '{print $2}'| sed 's/_L002//g'`
    echo "field2" $field2
    set field3 = `echo $line | awk '{print $3}'`
    echo "field3" $field3

# Construct the bsub command to run the alignment job for the current sample
    set  par="-n 1 -W 1:00 -o %J.stdout -e %J.stderr"
    # echo $par

    # Submit the bsub command to the job scheduler
    bsub  -R 'rusage[mem=1GB]' $par  ./mark_duplicates.sh ${field2} ${field3}
end
