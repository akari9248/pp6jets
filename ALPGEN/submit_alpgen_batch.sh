#!/bin/bash

# Usage: ./submit_alpgen_batch.sh <start_part> <end_part>
# Example: ./submit_alpgen_batch.sh 6 10

if [ $# -ne 2 ]; then
    echo "Usage: $0 <start_part> <end_part>"
    exit 1
fi

START=$1
END=$2
JDL_FILE="condor_algen.jdl"

for PART in $(seq $START $END); do
    JDL_TEMP="condor_alpgen_part${PART}.jdl"
    sed "s|/Part[0-9]*$|/Part${PART}|" ${JDL_FILE} > ${JDL_TEMP}
    
    echo "Submitting Part${PART}..."
    condor_submit ${JDL_TEMP}
    rm ${JDL_TEMP}

    sleep 10
done

echo "All jobs submitted!"
