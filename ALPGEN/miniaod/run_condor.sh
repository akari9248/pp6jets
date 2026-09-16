#!/bin/bash
# Run inside Condor's EL8 container. Inputs are transferred by Condor.
# Arguments: job_id LHE_URL CP2_or_CP5 (set by condor.jdl)
set -euo pipefail
JOB_ID=${1:?}; LHE_URL=${2:?}; TUNE=${3:?}
[[ "$JOB_ID" =~ ^[1-9][0-9]*$ ]] && ((JOB_ID < 80000)) || exit 2
INPUT=${LHE_URL##*/}
STEM=chunk$((JOB_ID - 1))
[[ "$TUNE" == CP2 || "$TUNE" == CP5 ]] || exit 2
[[ "$INPUT" == "${STEM}.lhe" ]] || { echo "Expected ${STEM}.lhe, got: $INPUT" >&2; exit 2; }
START=$PWD
WORK=$START/work
mkdir -p "$WORK"
# Always package available diagnostics, including when a stage fails.
archive_logs() {
  local status=$?
  trap - EXIT
  echo "$status" > "$WORK/exit_status.txt"
  cp fullsim.sh ALPGEN6j_Run2026C_cfi.py PU_Run2026C_cff.py minbias_files.txt "$WORK/"
  # Preserve the PU configuration; only production ROOT/LHE intermediates are excluded.
  cp Run2026C_PU.root "$WORK/Run2026C_PU.root.input"
  tar --exclude='*.root' --exclude='*.lhe' -czf "fullsim_${TUNE}_${STEM}.tgz" -C "$WORK" .
  exit "$status"
}
trap archive_logs EXIT
[[ -s "$INPUT" ]] || { echo "Missing transferred LHE: $INPUT" >&2; exit 1; }
[[ -r "${X509_USER_PROXY:-}" ]] || { echo 'Missing CMS proxy for pileup input' >&2; exit 1; }
export ALPGEN_EL8=1
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
printf 'Input=%s\nJobID=%s\nMaxLHEEvents=-1\nTune=%s\n' "$LHE_URL" "$JOB_ID" "$TUNE" > "$WORK/input.txt"
/usr/bin/time -v -o "$WORK/resources.txt" bash ./fullsim.sh "$START/$INPUT" "$WORK" -1 all "$JOB_ID" "$TUNE"
python3 - "$WORK" <<'PY'
import sys, pathlib, xml.etree.ElementTree as ET
work=pathlib.Path(sys.argv[1])
lines=[]
for stage in ['lhe','gensim','digihlt','reco','mini']:
    report=ET.parse(work/(stage+'.xml')).getroot()
    assert not report.findall('FrameworkError'), stage+' has framework errors'
    counts=[int(f.findtext('TotalEvents')) for f in report.findall('File')]
    assert len(counts)==1, stage+' must have one output'
    lines.append(f'{stage}: {counts[0]} events')
(work/'event_counts.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
PY
mv "$WORK/miniaod.root" "MiniAOD_${TUNE}_${STEM}.root"
echo "MiniAOD ready: MiniAOD_${TUNE}_${STEM}.root"
