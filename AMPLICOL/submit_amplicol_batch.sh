#!/bin/bash

# Usage: ./submit_amplicol_batch.sh <start_part> <end_part> [queue_num]
# Example: ./submit_amplicol_batch.sh 1 5        # submit Part1-5, 1000 jobs each
#          ./submit_amplicol_batch.sh 1 1 500     # submit Part1, 500 jobs

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <start_part> <end_part> [queue_num]"
    exit 1
fi

START=$1
END=$2
QUEUE_NUM=${3:-1000}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONDOR_JDL="${SCRIPT_DIR}/condor_amplicol.jdl"

mkdir -p "${SCRIPT_DIR}/log"

for PART_NUM in $(seq $START $END); do
    echo "[Part${PART_NUM}] Submitting ${QUEUE_NUM} jobs..."
    condor_submit "${CONDOR_JDL}" PART_NUM="${PART_NUM}" QUEUE_NUM="${QUEUE_NUM}"
    sleep 10
done

echo "All parts submitted!"
