#!/bin/bash
# ./submit_condor_auto.sh <Part> [quene]
# ./submit_condor_auto.sh 2 500
# for i in {1..5}; do ./submit_condor_auto.sh $i; done

set -e

PART_NUM=$1
QUEUE_NUM=${2:-1000}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDOR_JDL="${SCRIPT_DIR}/condor.jdl"

cd "${SCRIPT_DIR}"
mkdir -p log
echo "[Part${PART_NUM}] Submit tasks (queue: ${QUEUE_NUM})..."
condor_submit "${CONDOR_JDL}" PART_NUM="${PART_NUM}" QUEUE_NUM="${QUEUE_NUM}"

echo "[Part${PART_NUM}] Done!"

