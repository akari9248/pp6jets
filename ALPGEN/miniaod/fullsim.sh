#!/bin/bash
# Usage: fullsim.sh input.lhe work_directory [events=10] [all|lhe|gensim|digihlt|reco|mini|configs] [job_id=0] [CP2|CP5=CP5] [jet_filter=on|off]
set -euo pipefail
SCRIPT=$(realpath "${BASH_SOURCE[0]}")
BASE=$(dirname "$SCRIPT")
if (($# < 2)); then
  echo "Usage: $0 input.lhe work_directory [events=10] [all|lhe|gensim|digihlt|reco|mini|configs] [job_id=0] [CP2|CP5=CP5] [jet_filter=on|off]" >&2
  exit 2
fi
INPUT=$(realpath "$1")
WORK=$(realpath -m "$2")
EVENTS=${3:-10}
STAGE=${4:-all}
JOB_ID=${5:-0}
TUNE=${6:-CP5}
JET_FILTER=${7:-on}
[[ "$JET_FILTER" == on || "$JET_FILTER" == off ]] || { echo 'jet_filter must be on or off' >&2; exit 2; }
[[ "$TUNE" == CP2 || "$TUNE" == CP5 ]] || { echo "Tune must be CP2 or CP5" >&2; exit 2; }
[[ "$EVENTS" == -1 || "$EVENTS" =~ ^[1-9][0-9]*$ ]] || { echo "events must be positive or -1 (all)" >&2; exit 2; }
[[ "$JOB_ID" =~ ^[0-9]+$ ]] && ((JOB_ID < 80000)) || { echo "job_id must be 0..79999" >&2; exit 2; }
case "$STAGE" in all|lhe|gensim|digihlt|reco|mini|configs) ;; *) echo "Unknown stage: $STAGE" >&2; exit 2;; esac

if [[ ${ALPGEN_EL8:-0} != 1 ]]; then
  PROXY=${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}
  exec apptainer exec --cleanenv --bind /afs,/eos,/cvmfs,/tmp \
    /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
    env ALPGEN_EL8=1 X509_USER_PROXY="$PROXY" \
    X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates /bin/bash "$SCRIPT" "$INPUT" "$WORK" "$EVENTS" "$STAGE" "$JOB_ID" "$TUNE" "$JET_FILTER"
fi
mkdir -p "$WORK"
cd "$WORK"
if [[ -f tune.txt && $(cat tune.txt) != "$TUNE" ]]; then
  echo 'Use a separate work directory for each tune' >&2; exit 2
fi
printf '%s\n' "$TUNE" > tune.txt
export ALPGEN_TUNE="$TUNE"
if [[ -f jet_filter_mode.txt && $(cat jet_filter_mode.txt) != "$JET_FILTER" ]]; then
  echo 'Use a separate work directory for each jet filter mode' >&2; exit 2
fi
if [[ ! -f jet_filter_mode.txt && -f gensim.py ]]; then
  echo 'Use a new work directory for the central jet filter configuration' >&2; exit 2
fi
printf '%s\n' "$JET_FILTER" > jet_filter_mode.txt
export ALPGEN_JET_FILTER="$JET_FILTER"

# Reject accidentally reused Run2 input. The LHE header holds the actual energy.
python3 - "$INPUT" <<'PY'
import sys, xml.etree.ElementTree as ET
for _, node in ET.iterparse(sys.argv[1], events=('end',)):
    if node.tag == 'init':
        fields=node.text.split()
        assert all(abs(float(x.replace('D','E'))-6800.) < .01 for x in fields[2:4]), 'LHE must have 6800 GeV per beam'
        break
else:
    raise RuntimeError('Missing LHE init block')
PY

GT=160X_mcRun3_2026_lowPU_v3
GEN_RELEASE=/cvmfs/cms.cern.ch/el8_amd64_gcc13/cms/cmssw/CMSSW_16_0_8
RECO_RELEASE=/cvmfs/cms.cern.ch/el8_amd64_gcc13/cms/cmssw/CMSSW_16_0_6
PU_FILES=$BASE/minbias_files.txt
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
# Carry the fragment with the job; no private CMSSW area is needed on workers.
mkdir -p localgen
touch localgen/__init__.py
cp "$BASE/ALPGEN6j_Run2026C_cfi.py" "$BASE/PU_Run2026C_cff.py" "$BASE/CentralGenJetFilter_cff.py" localgen/
cp "$BASE/Run2026C_PU.root" "$WORK/"

COMMON=(--conditions "$GT" --era Run3_2026 --geometry DB:Extended --mc -n "$EVENTS" --no_exec)

# Each stage has its own release, Python config, log and FrameworkJobReport.
# All intermediate ROOT files are kept. Select one stage to rerun just that step.
step() (
  NAME=$1; RELEASE=$2; shift 2
  [[ "$STAGE" == all || "$STAGE" == configs || "$STAGE" == "$NAME" ]] || exit 0
  set +u
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  export SCRAM_ARCH=el8_amd64_gcc13
  cd "$RELEASE/src"
  eval "$(scram runtime -sh)"
  set -u
  cd "$WORK"
  export PYTHONPATH="$WORK:${PYTHONPATH:-}"
  echo "=== $NAME : $CMSSW_VERSION ==="
  cmsDriver.py "$@" "${COMMON[@]}" --python_filename "$NAME.py" > "${NAME}_config.log" 2>&1 || { tail -50 "${NAME}_config.log"; exit 1; }
  if [[ "$NAME" != lhe ]]; then
    cat >> "$NAME.py" <<'PYFILTER'

# Retain prefilter decisions and lumi-level weighted counters through MiniAOD.
from localgen.CentralGenJetFilter_cff import keep_filter_products
process = keep_filter_products(process)
PYFILTER
  fi
  python3 - "$NAME.py" <<'PY'
import ast, sys
ast.parse(open(sys.argv[1]).read())
PY
  # job_id=0 preserves the original local debugging seeds.
  # A batch job has one unique lumi and a separate seed range for each stage.
  if ((JOB_ID > 0)); then
    python3 - "$NAME.py" "$JOB_ID" "$NAME" <<'PYSEED'
import sys
path, job, stage = sys.argv[1], int(sys.argv[2]), sys.argv[3]
index = ['lhe', 'gensim', 'digihlt', 'reco', 'mini'].index(stage)
seed = 100000 + job * 10000 + index * 1000
with open(path, 'a') as out:
    out.write("\n# Batch identity and reproducible random seeds.\n")
    if stage == 'lhe':
        out.write(f"process.source.firstLuminosityBlock = cms.untracked.uint32({job})\n")
        out.write("process.source.numberEventsInLuminosityBlock = cms.untracked.uint32(1000000000)\n")
    out.write("if hasattr(process, 'RandomNumberGeneratorService'):\n")
    out.write("    from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper\n")
    out.write("    _rng = RandomNumberServiceHelper(process.RandomNumberGeneratorService)\n")
    out.write("    assert _rng.countSeeds() < 1000\n")
    out.write(f"    _rng.insertSeeds(*range({seed}, {seed} + _rng.countSeeds()))\n")
PYSEED
  fi
  if [[ "$STAGE" != configs ]]; then
    cmsRun -j "$NAME.xml" "$NAME.py" > "$NAME.log" 2>&1 || { tail -50 "$NAME.log"; exit 1; }
    if [[ "$NAME" == gensim ]]; then
      python3 "$BASE/filter_summary.py" sim.root central_jet_filter_summary.json --mode "$JET_FILTER" > central_jet_filter_summary.log 2>&1 || { cat central_jet_filter_summary.log; exit 1; }
    fi
    echo "$NAME done; log: $WORK/$NAME.log"
  fi
)

step lhe "$GEN_RELEASE" MCDBtoEDM -s NONE --eventcontent LHE --datatier LHE \
  --filein "file:$INPUT" --fileout file:lhe.root

step gensim "$GEN_RELEASE" localgen/ALPGEN6j_Run2026C_cfi.py \
  --step GEN,SIM --beamspot DBrealistic --eventcontent RAWSIM --datatier GEN-SIM --nThreads 1 \
  --customise localgen/CentralGenJetFilter_cff.customise \
  --filein file:lhe.root --fileout file:sim.root

# The base scenario supplies 25 ns and BX -5..3; customise replaces mean=5 with the 0<=mu<10 histogram.
step digihlt "$RECO_RELEASE" --step DIGI,L1,DIGI2RAW,HLT:2026v11 \
  --pileup E7TeV_AVE_5_BX2808 --pileup_input "filelist:$PU_FILES" \
  --customise localgen/PU_Run2026C_cff.customise \
  --eventcontent RAWSIM --datatier GEN-SIM-RAW --nThreads 1 \
  --filein file:sim.root --fileout file:hlt.root

step reco "$RECO_RELEASE" --step RAW2DIGI,L1Reco,RECO,RECOSIM \
  --eventcontent AODSIM --datatier AODSIM --nThreads 1 \
  --filein file:hlt.root --fileout file:reco.root

step mini "$GEN_RELEASE" --step PAT --eventcontent MINIAODSIM --datatier MINIAODSIM \
  --filein file:reco.root --fileout file:miniaod.root

echo "Finished requested stage: $STAGE. Files and logs: $WORK"
