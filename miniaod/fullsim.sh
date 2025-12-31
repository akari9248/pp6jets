#!/bin/bash

home="/afs/cern.ch/user/s/shuangyu/"
eventsnum=1
cd ${home}"CMSSW_10_6_28_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py MCDBtoEDM --conditions 106X_mc2017_realistic_v6 -s NONE --eventcontent RAWSIM --datatier GEN --filein file:/eos/cms/store/group/phys_smp/ec/shuangyu/pp6j_15GeV/Part1/chunk${1}.lhe -n ${eventsnum}
cmsDriver.py Configuration/GenProduction/python/JME-RunIISummer20UL17GEN-00006-fragment.py --python_filename GEN.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:gen.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --era Run2_2017 --mc -n ${eventsnum} --filein file:MCDBtoEDM_NONE.root
rm MCDBtoEDM_NONE.root

cd ${home}"CMSSW_10_6_17_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename SIM.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:sim.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step SIM --nThreads 4 --geometry DB:Extended --filein file:gen.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm gen.root
cmsDriver.py --python_filename DIGI2RAW.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI --fileout file:digi.root --pileup_input filelist:existing_PUfiles_RunIISummer20UL17.txt --conditions 106X_mc2017_realistic_v6 --customise_commands "process.mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(0,1,2,3,4,5,6,7,8,9,10)\nprocess.mix.input.nbPileupEvents.probValue = cms.vdouble(0.00151109, 0.01743738, 0.4441798, 0.4967324, 0.0300071, 0.00585244, 0.001825343, 0.000985512, 0.000730211, 0.000545979, 0.000192059)" --step DIGI,L1,DIGI2RAW --nThreads 4 --geometry DB:Extended --filein file:sim.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum} --pileup 2017_25ns_UltraLegacy_PoissonOOTPU
rm sim.root

cd ${home}"CMSSW_9_4_14_UL_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename HLT.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --fileout file:hlt.root --conditions 94X_mc2017_realistic_v15 --customise_commands "process.source.bypassVersionCheck = cms.untracked.bool(True)" --step HLT:2e34v40 --nThreads 4 --geometry DB:Extended --filein file:digi.root --era Run2_2017 --mc -n ${eventsnum}
rm digi.root

cd ${home}"CMSSW_10_6_17_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename RECO.py --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --fileout file:reco.root --conditions 106X_mc2017_realistic_v6 --step RAW2DIGI,L1Reco,RECO,RECOSIM --nThreads 4 --geometry DB:Extended --filein file:hlt.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm hlt.root

cd ${home}"CMSSW_10_6_20"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
cmsDriver.py --python_filename MiniAOD.py --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --fileout file:JME-RunIISummer20UL17MiniAODv2-${1}.root --conditions 106X_mc2017_realistic_v9 --step PAT --procModifiers run2_miniAOD_UL --nThreads 4 --geometry DB:Extended --filein file:reco.root --era Run2_2017 --runUnscheduled --mc -n ${eventsnum}
rm reco.root
rm MCDBtoEDM.py
rm GEN.py
rm SIM.py
rm DIGI2RAW.py
rm HLT.py
rm RECO.py
rm MiniAOD.py
