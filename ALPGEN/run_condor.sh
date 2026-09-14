#!/bin/bash
# Condor runs this in the EL7 image, in its private scratch directory.
# Arguments: global_chunk_index cluster_id (set by condor.jdl)
set -euo pipefail
INDEX=${1:?}; CLUSTER=${2:?}
[[ "$INDEX" =~ ^[0-9]+$ && "$CLUSTER" =~ ^[0-9]+$ ]] || exit 2
((INDEX < 30000)) || exit 2
ATTEMPTS=20000000
WARMUP=20000
ITERATIONS=1000
START=$PWD
LABEL=chunk${INDEX}
mkdir -p "$LABEL"
# Part and process determine one global index; different parts use different seeds.
S1=$((10000 + INDEX)); S2=$((30000 + INDEX))
S3=$((50000 + INDEX)); S4=$((70000 + INDEX))
awk -v warmup="$WARMUP" -v iterations="$ITERATIONS" -v attempts="$ATTEMPTS" -v s1="$S1" -v s2="$S2" '
  NR==4 {$0=warmup ", " iterations} NR==5 {$0=attempts}
  $1=="iseed1" {$0="iseed1 " s1} $1=="iseed2" {$0="iseed2 " s2} {print}
' input_step1.dat > "$LABEL/input_step1.dat"
awk -v s3="$S3" -v s4="$S4" '
  $1=="iseed3" {$0="iseed3 " s3} $1=="iseed4" {$0="iseed4 " s4} {print}
' input_step2.dat > "$LABEL/input_step2.dat"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw-patch/CMSSW_10_6_28_patch1/src
eval "$(scram runtime -sh)"
set -u
cd "$START/alpgen/Njetwork"
echo "Cluster $CLUSTER chunk $INDEX: attempts=$ATTEMPTS warmup=$WARMUP,$ITERATIONS seeds=$S1,$S2,$S3,$S4"
./Njetgen < "$START/$LABEL/input_step1.dat" > "$START/$LABEL/step1.log" 2>&1 || { tail -40 "$START/$LABEL/step1.log"; exit 1; }
./Njetgen < "$START/$LABEL/input_step2.dat" > "$START/$LABEL/step2.log" 2>&1 || { tail -40 "$START/$LABEL/step2.log"; exit 1; }
[[ -s sixjets.lhe ]] || { tail -40 "$START/$LABEL/step2.log"; exit 1; }
if ! grep -q '</LesHouchesEvents>' sixjets.lhe; then
  printf '\n</LesHouchesEvents>\n' >> sixjets.lhe
fi
EVENTS=$(grep -c '<event>' sixjets.lhe || true)
((EVENTS > 0)) || { echo 'No unweighted events generated' >&2; exit 1; }
cp sixjets.stat sixjets.par sixjets.grid1 sixjets.grid2 "$START/$LABEL/"
printf 'cluster=%s\nchunk=%s\nlhe_events=%s\nattempts=%s\n' "$CLUSTER" "$INDEX" "$EVENTS" "$ATTEMPTS" > "$START/$LABEL/summary.txt"
mv sixjets.lhe "$START/$LABEL.lhe"
cd "$START"
tar -czf "${LABEL}_logs.tgz" "$LABEL"
echo "LHE_EVENTS=$EVENTS"
tail -5 "$LABEL/step2.log"
echo "Ready for Condor transfer: $LABEL.lhe ${LABEL}_logs.tgz"
