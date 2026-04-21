#!/bin/bash
set -euo pipefail

case "$USER" in
  zhye)
    home="/afs/cern.ch/user/z/zhye/"
    ;;
  shuangyu)
    home="/afs/cern.ch/user/s/shuangyu/"
    ;;
  *)
    echo "Error: unknown USER '$USER', cannot set home path"
    exit 1
    ;;
esac

CMSSW_BASE="${CMSSW_BASE_OVERRIDE:-${home}CMSSW_10_6_28_patch1}"
SCRAM_ARCH_VALUE="${SCRAM_ARCH_OVERRIDE:-slc7_amd64_gcc700}"

if ! command -v cmssw-el7 >/dev/null 2>&1; then
  echo "Error: cmssw-el7 is not available on this machine."
  exit 1
fi

if [[ ! -d "${CMSSW_BASE}/src/GeneratorInterface/Pythia8Interface" ]]; then
  echo "Error: cannot find GeneratorInterface/Pythia8Interface under ${CMSSW_BASE}/src"
  exit 1
fi

source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "Checking the patched Pythia8 plugin inside the el7 runtime..."
cmssw-el7 -- bash -lc "source /cvmfs/cms.cern.ch/cmsset_default.sh && export SCRAM_ARCH=${SCRAM_ARCH_VALUE} && cd ${CMSSW_BASE}/src && eval \$(scramv1 runtime -sh) && edmPluginHelp -p Pythia8HadronizerFilter >/dev/null"

echo "Success: the patched Pythia8 plugin loads correctly in the el7 runtime."
