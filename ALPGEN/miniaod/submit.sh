#!/bin/bash
set -euo pipefail
# 只改这里：每个 part、每个 tune 固定 1000 个任务。
FIRST_PART=1
LAST_PART=10
# on: apply loose gen-jet veto; off: record decisions but simulate all events.
JET_FILTER=on
case "$JET_FILTER" in
  on) SAMPLE=er3p0_GenJetFilter_v1 ;;
  off) SAMPLE=er3p0_GenJetMonitor_v1 ;;
  *) echo 'JET_FILTER must be on or off' >&2; exit 2 ;;
esac

cd "$(dirname "$(realpath "$0")")"
((FIRST_PART >= 1 && LAST_PART >= FIRST_PART && LAST_PART <= 79)) || {
  echo 'MiniAOD part 范围应为 1..79（现有随机种子范围）。' >&2; exit 2;
}
# Each submitter uses their own proxy; bigbird reads the shared AFS copy.
SUBMITTER=$(id -un)
PROXY="$HOME/private/x509up_u$(id -u)"
install -m 600 "/tmp/x509up_u$(id -u)" "$PROXY"
mkdir -p log
BASE=/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV
OUTPUT_BASE=$BASE
if [[ "$SUBMITTER" != zhye ]]; then
  OUTPUT_BASE="$BASE/$SUBMITTER"
fi
for tune in CP5; do
  for ((part=FIRST_PART; part<=LAST_PART; part++)); do
    lhe_base="root://eoscms.cern.ch/${BASE}/LHE_er3p0/Part${part}/chunk"
    mkdir -p "${OUTPUT_BASE}/MiniAOD_${tune}_${SAMPLE}/Part${part}"
    echo "提交 ${tune} ${SAMPLE} Part${part}：1000 个任务"
    condor_submit condor.jdl PART="$part" TUNE="$tune" FILTER_MODE="$JET_FILTER" SAMPLE="$SAMPLE" LHE_BASE="$lhe_base" PROXY="$PROXY" OUTPUT_BASE="root://eoscms.cern.ch/${OUTPUT_BASE}"
  done
done
