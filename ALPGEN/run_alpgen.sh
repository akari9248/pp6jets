#!/bin/bash
# Usage: ./ALPGEN/run_alpgen.sh [work_directory]
set -euo pipefail
SCRIPT=$(realpath "${BASH_SOURCE[0]}")
BASE=$(dirname "$SCRIPT")
WORK=$(realpath -m "${1:-$BASE/work/first}")

# The existing ALPGEN binary uses the old Fortran/LHAPDF libraries.
if [[ ${ALPGEN_EL7:-0} != 1 ]]; then
  exec apptainer exec --cleanenv --bind /afs,/eos,/cvmfs,/tmp \
    /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el7:x86_64 \
    env ALPGEN_EL7=1 /bin/bash "$SCRIPT" "$WORK"
fi
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-patch/CMSSW_10_6_28_patch1/src
eval "$(scram runtime -sh)"
set -u

mkdir -p "$WORK"
[[ ! -e "$WORK/alpgen" ]] || { echo "Use a new work directory: $WORK/alpgen already exists" >&2; exit 1; }
cp -a "$BASE/alpgen" "$WORK/alpgen"
cp "$BASE/input_step1.dat" "$WORK/input_step1.dat"
cp "$BASE/input_step2.dat" "$WORK/input_step2.dat"
cd "$WORK/alpgen/Njetwork"

echo "ALPGEN integration and weighted generation; log: $WORK/alpgen_step1.log"
./Njetgen < "$WORK/input_step1.dat" > "$WORK/alpgen_step1.log" 2>&1
echo "ALPGEN unweighting; log: $WORK/alpgen_step2.log"
./Njetgen < "$WORK/input_step2.dat" > "$WORK/alpgen_step2.log" 2>&1
[[ -s sixjets.lhe ]] || { tail -30 "$WORK/alpgen_step2.log"; exit 1; }
cp sixjets.lhe "$WORK/sixjets.lhe"
# This ALPGEN version can omit the final XML closing tag.
if ! grep -q '</LesHouchesEvents>' "$WORK/sixjets.lhe"; then
  printf '\n</LesHouchesEvents>\n' >> "$WORK/sixjets.lhe"
fi
echo "LHE ready: $WORK/sixjets.lhe"
echo "Unweighted events: $(grep -c '<event>' "$WORK/sixjets.lhe")"
