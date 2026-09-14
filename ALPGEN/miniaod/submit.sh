#!/bin/bash
set -euo pipefail
# 只改这里：每个 part、每个 tune 固定 1000 个任务。
FIRST_PART=1
LAST_PART=5

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
for tune in CP2 CP5; do
  for ((part=FIRST_PART; part<=LAST_PART; part++)); do
    # 当前 5000 份旧命名输入直接对应 Part1–5，无需搬动文件。
    if ((part <= 5)); then
      lhe_base="root://eoscms.cern.ch/${BASE}/LHE_5M/chunk16823503_"
    else
      lhe_base="root://eoscms.cern.ch/${BASE}/LHE/Part${part}/chunk"
    fi
    mkdir -p "${OUTPUT_BASE}/MiniAOD_${tune}_AQCDUP_v1/Part${part}"
    echo "提交 ${tune} Part${part}：1000 个任务"
    condor_submit condor.jdl PART="$part" TUNE="$tune" LHE_BASE="$lhe_base" PROXY="$PROXY" OUTPUT_BASE="root://eoscms.cern.ch/${OUTPUT_BASE}"
  done
done
