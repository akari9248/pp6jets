# pp6jets

This repository contains two **parallel event-generator pipelines**:

- `ALPGEN` generator chain
- `AMPLICOL` generator chain

Both produce LHE events independently, and then each branch runs the **same MiniAOD production logic** (GEN/SIM/DIGI/HLT/RECO/MiniAOD with equivalent CMSSW-step structure).

## Workflow Overview

1. Generate LHE events with **either** ALPGEN or AMPLICOL.
2. Run full simulation to MiniAOD for the selected branch:
   - `ALPGEN/miniaod/fullsim.sh`
   - `AMPLICOL/miniaod_AmpliCol/fullsim.sh`
3. Submit Condor jobs with the matching `condor.jdl` in that branch.

## Start From Scratch (lxplus)

### 1) Initialize the environment

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmssw-el7
```

### 2) Create all required CMSSW releases

The MiniAOD scripts switch across these releases:

- `CMSSW_10_6_28_patch1` (LHE->GEN setup + fragment handling)
- `CMSSW_10_6_17_patch1` (SIM/DIGI2RAW/RECO)
- `CMSSW_9_4_14_UL_patch1` (HLT menu in 94X)
- `CMSSW_10_6_20` (final MiniAOD step)

Create all of them:

```bash
cmsrel CMSSW_10_6_28_patch1
cmsrel CMSSW_10_6_17_patch1
cmsrel CMSSW_9_4_14_UL_patch1
cmsrel CMSSW_10_6_20
```

### 3) Clone this repository

```bash
cd CMSSW_10_6_28_patch1/src
git clone git@github.com:akari9248/pp6jets.git
cd pp6jets
```

## Branch A: ALPGEN

### A1) Build ALPGEN

```bash
cmsenv
scram tool info lhapdf
```

Set `CONFIG_FILE_DIR=...` in `alpgen/alplib/alpgen.mk` to the LHAPDF `PATH` directory from the command output, then:

```bash
make cleanall
make gen
./Njetgen
```

### A2) Run MiniAOD chain for ALPGEN LHE

Working directory: `ALPGEN/miniaod`

Before submission:

- set ALPGEN LHE input path in `ALPGEN/miniaod/fullsim.sh` (`inputfile=...`)
- set EOS output path in `ALPGEN/miniaod/condor.jdl` (`output_destination`)
- ensure `x509userproxy` points to a valid proxy file

Submit with the batch helper script:

```bash
cd ALPGEN/miniaod
./submit_condor_auto.sh <PART_NUM> [QUEUE_NUM]
```

Examples:

```bash
./submit_condor_auto.sh 1
./submit_condor_auto.sh 2 500
for i in {1..5}; do ./submit_condor_auto.sh $i; done
```

## Branch B: AMPLICOL

### B1) Produce AMPLICOL LHE

Use the scripts/config under `AMPLICOL` to produce your LHE samples.

### B2) Run MiniAOD chain for AMPLICOL LHE

Working directory: `AMPLICOL/miniaod_AmpliCol`

Before submission:

- set AMPLICOL LHE input path in `AMPLICOL/miniaod_AmpliCol/fullsim.sh` (`inputfile=...`)
- set EOS output path in `AMPLICOL/miniaod_AmpliCol/condor.jdl` (`output_destination`)
- ensure `x509userproxy` points to a valid proxy file

Submit with the batch helper script:

```bash
cd AMPLICOL/miniaod_AmpliCol
./submit_condor_auto.sh <PART_NUM> [QUEUE_NUM]
```

Examples:

```bash
./submit_condor_auto.sh 1
./submit_condor_auto.sh 2 500
for i in {1..5}; do ./submit_condor_auto.sh $i; done
```

## VOMS Proxy

Create/update proxy before Condor submission:

```bash
voms-proxy-init -voms cms -rfc -out x509up
```

## Useful Checks

- detect schedd nodes where you currently have jobs:
  ```bash
  ./detect_active_schedds.sh
  ```
- check specific nodes only:
  ```bash
  ./detect_active_schedds.sh bigbird13.cern.ch bigbird19.cern.ch
  ```
