#!/bin/bash

# Usage: ./submit_amplicol_batch.sh <start_part> <end_part> <njet> [queue_num]
# Example: ./submit_amplicol_batch.sh 1 5 6        # submit Part1-5, 1000 jobs each for 6j
#          ./submit_amplicol_batch.sh 1 1 4 500    # submit Part1, 500 jobs for 4j

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <start_part> <end_part> <njet> [queue_num]"
    exit 1
fi

START=$1
END=$2
NJET=$3
QUEUE_NUM=${4:-1000}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONDOR_JDL="${SCRIPT_DIR}/condor_amplicol.jdl"

for PART_NUM in $(seq $START $END); do
    echo "[Part${PART_NUM}] Submitting ${QUEUE_NUM} jobs for ${NJET}j..."
    condor_submit "${CONDOR_JDL}" PART_NUM="${PART_NUM}" QUEUE_NUM="${QUEUE_NUM}" NJET="${NJET}"
    sleep 10
done

echo "All parts submitted!"
