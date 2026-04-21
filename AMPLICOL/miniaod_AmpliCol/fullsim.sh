#!/bin/bash
# 94X_mc2017_realistic_v15
set -e

# Derive paths from script location to avoid USER-based hardcoding on condor nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PP6JETS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
CMSSW_10_6_28_ROOT="$(cd "${PP6JETS_DIR}/../.." && pwd -P)"
home="$(cd "${CMSSW_10_6_28_ROOT}/.." && pwd -P)/"

FRAG_SRC="${SCRIPT_DIR}/JME-RunIISummer20UL17GEN-00006-fragment.py"
FRAG_LINK="${CMSSW_10_6_28_ROOT}/src/Configuration/GenProduction/python/JME-RunIISummer20UL17GEN-00006-fragment-AmpliCol.py"
if [[ ! -f "$FRAG_SRC" ]]; then
  echo "Error: fragment source not found: $FRAG_SRC"
  exit 1
fi
FRAG_SRC_REALPATH=$(realpath "$FRAG_SRC")
FRAG_LINK_REALPATH=""
if [[ -e "$FRAG_LINK" || -L "$FRAG_LINK" ]]; then
  FRAG_LINK_REALPATH=$(realpath "$FRAG_LINK" 2>/dev/null || true)
fi
if [[ "$FRAG_LINK_REALPATH" != "$FRAG_SRC_REALPATH" ]]; then
  ln -sf "$FRAG_SRC" "$FRAG_LINK"
fi

GT="106X_mc2017_realistic_v9For2017H_v1"
PART_NUM="${1}"
CHUNK_NUM="${2}"
inputfile="/eos/cms/store/group/phys_smp/ec/shuangyu/Amplicol/pp6j_20GeV_SC2/Part${PART_NUM}/part${PART_NUM}_chunk${CHUNK_NUM}_events.lhe.rwgt"
DEBUG_MODE=0

if [[ "$*" == *"--DEBUG"* ]]; then
  DEBUG_MODE=1
  timestamp=$(date +%Y%m%d_%H%M%S)
  logfile="fullsim_${timestamp}.log"
  exec > >(tee -a "${logfile}") 2>&1
  echo "DEBUG mode enabled, logging to ${logfile}"
fi


eventsnum=$(grep -c '<event>' "${inputfile}")
echo "Found ${eventsnum} events in ${inputfile}"
if [[ "${DEBUG_MODE}" -eq 1 ]]; then
  eventsnum=100
fi


ln -sf "${inputfile}" input.lhe

# Use 10_6_28 for LHE->GEN setup and fragment handling (UL-era GEN configs).
cd ${home}"CMSSW_10_6_28_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py MCDBtoEDM --conditions ${GT} -s NONE --eventcontent LHE --datatier LHE --filein file:input.lhe --fileout file:MCDBtoEDM_NONE.root -n ${eventsnum} --mc
cmsDriver.py Configuration/GenProduction/python/JME-RunIISummer20UL17GEN-00006-fragment-AmpliCol.py --python_filename GEN.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:gen.root --conditions ${GT} --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --era Run2_2017 --mc -n ${eventsnum} --filein file:MCDBtoEDM_NONE.root
rm MCDBtoEDM_NONE.root

# Use 10_6_17 for Run2 UL SIM and DIGI chain compatibility.
cd ${home}"CMSSW_10_6_17_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename SIM.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:sim.root --conditions ${GT} --beamspot Realistic25ns13TeVEarly2017Collision --step SIM --nThreads 4 --geometry DB:Extended --filein file:gen.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm gen.root
cmsDriver.py --python_filename DIGI2RAW.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI --fileout file:digi.root --pileup_input filelist:existing_PUfiles_RunIISummer20UL17.txt --conditions 106X_mc2017_realistic_v9For2017H_v1 --customise_commands "process.mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(0,1,2,3,4,5,6,7,8,9,10)\nprocess.mix.input.nbPileupEvents.probValue = cms.vdouble(0.00151109, 0.01743738, 0.4441798, 0.4967324, 0.0300071, 0.00585244, 0.001825343, 0.000985512, 0.000730211, 0.000545979, 0.000192059)" --step DIGI,L1,DIGI2RAW --nThreads 4 --geometry DB:Extended --filein file:sim.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum} --pileup 2017_25ns_UltraLegacy_PoissonOOTPU
rm sim.root

# HLT menu 2e34v40 is tied to 94X, so run HLT in 9_4_14_UL_patch1.
cd ${home}"CMSSW_9_4_14_UL_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename HLT.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --fileout file:hlt.root --conditions 94X_mc2017_realistic_v15 --customise_commands "process.source.bypassVersionCheck = cms.untracked.bool(True)" --step HLT:2e34v40 --nThreads 4 --geometry DB:Extended --filein file:digi.root --era Run2_2017 --mc -n ${eventsnum}
rm digi.root

# Switch back to 10_6_17 for UL RECO on 2017 conditions/geometry.
cd ${home}"CMSSW_10_6_17_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename RECO.py --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --fileout file:reco.root --conditions ${GT} --step RAW2DIGI,L1Reco,RECO,RECOSIM --nThreads 4 --geometry DB:Extended --filein file:hlt.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm hlt.root

# Final MiniAOD step uses 10_6_20 with run2_miniAOD_UL modifier.
cd ${home}"CMSSW_10_6_20"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename MiniAOD.py --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --fileout file:JME-RunIISummer20UL17MiniAODv2-${CHUNK_NUM}.root --conditions ${GT} --step PAT --procModifiers run2_miniAOD_UL --nThreads 4 --geometry DB:Extended --filein file:reco.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm reco.root
rm input.lhe
rm MCDBtoEDM*.py
rm GEN.py
rm SIM.py
rm DIGI2RAW.py
rm HLT.py
rm RECO.py
rm MiniAOD.py
