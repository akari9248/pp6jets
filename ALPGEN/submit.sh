#!/bin/bash
set -euo pipefail
# 只改这里：每个 part 固定 1000 个任务。
# Part1–5 已由 cluster 16823503 生产，新 LHE 从 Part6 开始。
FIRST_PART=6
LAST_PART=10

cd "$(dirname "$(realpath "$0")")"
((FIRST_PART >= 6 && LAST_PART >= FIRST_PART && LAST_PART <= 30)) || {
  echo 'LHE part 范围应为 6..30（现有随机种子范围）；Part1–5 已保留。' >&2; exit 2;
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
for ((part=FIRST_PART; part<=LAST_PART; part++)); do
  mkdir -p "${OUTPUT_BASE}/LHE/Part${part}"
  echo "提交 LHE Part${part}：1000 个任务"
  condor_submit condor.jdl PART="$part" PROXY="$PROXY" OUTPUT_BASE="root://eoscms.cern.ch/${OUTPUT_BASE}"
done
