#!/bin/bash
#
# Goes through all data_ folders and parses subsets of OBJ files for interface
# and knots at given iterations/timepoints.
#
# Usage: sh parseobjs.sh
# - Run at the root folder with all the data_* folders, with parseobj binary
#   present at the root as well.
#

# Timepoints to process
TIMEPOINTS=(4 20 36 68)

for d in data_*; 
do
    cd $d
    for t in "${TIMEPOINTS[@]}";
    do
        # Output files names with the folder name instead of the run ID.
        SOURCE=$t'_*.obj'
        TARGET=$t'_'$d'.obj'
        cp $SOURCE $TARGET
        ../parseobj $TARGET interface knots
        rm $TARGET
    done
    cd ..
done
