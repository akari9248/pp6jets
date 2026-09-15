#!/bin/bash
set -euo pipefail

FIRST_PART=1
LAST_PART=5

cd "$(dirname "$(realpath "$0")")"
SUBMITTER=$(id -un)
PROXY="$HOME/private/x509up_u$(id -u)"
install -m 600 "/tmp/x509up_u$(id -u)" "$PROXY"
mkdir -p log
BASE=/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV
OUTPUT_BASE=$BASE
if [[ "$SUBMITTER" != zhye ]]; then
  OUTPUT_BASE="$BASE/$SUBMITTER"
fi
for ((part=FIRST_PART; part<=LAST_PART; part++)); do
  mkdir -p "${OUTPUT_BASE}/LHE/Part${part}"
  echo "提交 LHE Part${part}：1000 个任务"
  condor_submit condor.jdl PART="$part" PROXY="$PROXY" OUTPUT_BASE="root://eoscms.cern.ch/${OUTPUT_BASE}"
done
