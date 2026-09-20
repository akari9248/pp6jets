#!/bin/bash
set -euo pipefail

# LHE_er3p0 独立生产；续产使用该目录下未用过的 part。
FIRST_PART=16
LAST_PART=30

cd "$(dirname "$(realpath "$0")")"
# Complete 1000-job parts with max seed 110000 + INDEX < 2147483399.
((FIRST_PART >= 1 && LAST_PART >= FIRST_PART && LAST_PART <= 2147373)) || {
  echo 'LHE part 范围应为 1..2147373（ALPGEN 随机种子取值上限）。' >&2; exit 2;
}
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
  mkdir -p "${OUTPUT_BASE}/LHE_er3p0/Part${part}"
  echo "提交 LHE_er3p0 Part${part}：1000 个任务"
  condor_submit condor.jdl PART="$part" PROXY="$PROXY" OUTPUT_BASE="root://eoscms.cern.ch/${OUTPUT_BASE}"
done
