#!/bin/bash

# All operations run inside the AmpliCol directory.
# Run directly on el9 (no container needed).
# Manually set LHAPDF runtime paths so the CMS GCC 7 libstdc++ is NOT loaded,
# allowing the binary compiled with el9 system GCC to run correctly.

LHAPDF_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3"
export LD_LIBRARY_PATH="${LHAPDF_BASE}/lib"
export LHAPDF_DATA_PATH="${LHAPDF_BASE}/share/LHAPDF"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}" || { echo "Cannot cd to ${SCRIPT_DIR}"; exit 1; }

# On Condor worker nodes, the executable is placed at the IWD root
# while Transfer_Input_Files puts AmpliCol/ as a subdirectory.
if [ ! -f "amplicol_generate" ] && [ -f "AmpliCol/amplicol_generate" ]; then
    cd AmpliCol
fi

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Please provide part number and process number"
    echo "Usage: $0 <part_num> <process>"
    exit 1
fi

PART_NUM=$1
PROCESS=$2
NEVENTS=1000

# Unique seed and tag across all parts and processes
SEED=$((PART_NUM * 10000 + PROCESS * 10 + RANDOM % 10))
TAG="part${PART_NUM}_chunk${PROCESS}"

echo "========================================"
echo "  AmpliCol pp -> 6j Event Generation"
echo "========================================"
echo "Part:     ${PART_NUM}"
echo "Process:  ${PROCESS}"
echo "Events:   ${NEVENTS}"
echo "Seed:     ${SEED}"
echo "Tag:      ${TAG}"
echo "Work dir: ${SCRIPT_DIR}"
echo "========================================"

# Step 1: Process generation (creates processes.txt, only once)
if [ ! -f "processes.txt" ]; then
    echo "Running process generation..."
    python3 ./process_list.py 'p p > 6j'
    if [ $? -ne 0 ]; then
        echo "Error: process_list.py failed"
        exit 1
    fi
    echo "processes.txt generated."
else
    echo "processes.txt already exists, skipping process generation."
fi

mkdir -p Outputs

# Step 2: Generate leading-colour events
echo "Running event generation (${TAG})..."
./amplicol_generate --nevents=${NEVENTS} --seed=${SEED} --tag=${TAG} --process=processes.txt

LC_EVENT_FILE="Outputs/${TAG}_events.lhe"
if [ ! -f "${LC_EVENT_FILE}" ]; then
    echo "Error: Leading-colour event file '${LC_EVENT_FILE}' not generated"
    exit 1
fi
echo "Leading-colour events generated: ${LC_EVENT_FILE}"

# Step 3: Reweight to full-colour accuracy
if [ ! -f "amplicol_reweight" ]; then
    echo "Compiling amplicol_reweight..."
    make amplicol_reweight
    if [ $? -ne 0 ]; then
        echo "Error: Compilation of amplicol_reweight failed"
        exit 1
    fi
fi

echo "Running reweighting to full-colour accuracy..."
./amplicol_reweight "${LC_EVENT_FILE}"
if [ $? -ne 0 ]; then
    echo "Error: Reweighting failed"
    exit 1
fi

FC_EVENT_FILE="${LC_EVENT_FILE}.rwgt"
if [ ! -f "${FC_EVENT_FILE}" ]; then
    echo "Error: Full-colour event file '${FC_EVENT_FILE}' not generated"
    exit 1
fi

cp Outputs/${TAG}_events.lhe.rwgt ${SCRIPT_DIR}/

echo "========================================"
echo "  Completed!"
echo "  LC events: Outputs/${TAG}_events.lhe"
echo "  FC events: Outputs/${TAG}_events.lhe.rwgt"
echo "========================================"
