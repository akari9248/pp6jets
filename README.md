# pp6jets

A workflow for generating pp → 6 jets events using ALPGEN + PYTHIA within the CMS software framework.

## Installation

This repository must be placed inside a CMSSW release. For 2017 simulation, use `CMSSW_10_6_28_patch1`:

```bash
cmssw-el7
cmsrel CMSSW_10_6_28_patch1
cd CMSSW_10_6_28_patch1/src
git clone git@github.com:akari9248/pp6jets.git
```


## ALPGEN Setup

### 1. Configure LHAPDF Path

After entering the `pp6jets` directory, run:

```bash
cmsenv
scram tool info lhapdf
```

You will see output similar to:

```
Tool info as configured in location /afs/cern.ch/user/s/shuangyu/CMSSW_10_6_28_patch1
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Name : lhapdf
Version : 6.2.1-pafccj3
++++++++++++++++++++
SCRAM_PROJECT=no
LHAPDF_BASE=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3
LIB=LHAPDF
LIBDIR=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/lib
INCLUDE=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/include
USE=yaml-cpp root_cxxdefaults
LHAPDF_DATA_PATH=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/share/LHAPDF
ROOT_INCLUDE_PATH=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/include
PATH=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/bin
```

### 2. Update Makefile

Open `alpgen/alplib/alpgen.mk` and find the line `CONFIG_FILE_DIR=xxx`. Replace it with the `PATH` value from the output above.

### 3. Compile ALPGEN

```bash
make cleanall
make gen
```

After successful compilation, test by running `./Njetgen`.

### 4. Configure and Run

- The main executable is `run_alpgen.sh`
- ALPGEN input parameters are defined in two sections:
  - `cat > config_step1.txt << EOF`
  - `cat > config_step2.txt << EOF`

---

## PYTHIA Integration

### 1. Install the Fragment File

Copy `JME-RunIISummer20UL17GEN-00006-fragment.py` to:

```
CMSSW_10_6_28_patch1/src/Configuration/GenProduction/python/
```

> **Note:** Create this directory if it doesn't exist.

If you have custom configuration files, place them in the same directory and update the corresponding path in `miniaod/fullsim.sh`:

```bash
cmsDriver.py Configuration/GenProduction/python/JME-RunIISummer20UL17GEN-00006-fragment.py \
  --python_filename GEN.py \
  --eventcontent RAWSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier GEN \
  --fileout file:gen.root \
  --conditions 106X_mc2017_realistic_v6 \
  --beamspot Realistic25ns13TeVEarly2017Collision \
  --step GEN \
  --geometry DB:Extended \
  --era Run2_2017 \
  --mc \
  -n ${eventsnum} \
  --filein file:MCDBtoEDM_NONE.root
```

### 2. Compile CMSSW

After placing the fragment file, compile:

```bash
scram b -j 8
```

---

## Condor Job Submission

### ALPGEN Jobs

1. Exit the Singularity container using `Ctrl + A + D`
2. Submit with: `condor_submit condor_algen.jdl`
3. **Important:** Modify `output_destination` before each submission to avoid overwriting previous results

### Full Simulation (MiniAOD) Jobs

Navigate to the `miniaod` folder and complete the following steps before submission:

1. **Update home path** in `fullsim.sh` to your own directory, and install the required CMSSW versions listed there

2. **Generate VOMS proxy:**
   ```bash
   voms-proxy-init -voms cms -rfc -out x509up
   ```

3. **Set output path:** Modify `output_destination` in `condor.jdl`

4. **Configure LHE input:** Update the `--filein` parameter in `fullsim.sh`:
   ```bash
   cmsDriver.py MCDBtoEDM \
     --conditions 106X_mc2017_realistic_v6 \
     -s NONE \
     --eventcontent RAWSIM \
     --datatier GEN \
     --filein file:/eos/cms/store/group/phys_smp/ec/shuangyu/pp6j_15GeV/Part1/chunk${1}.lhe \
     -n ${eventsnum}
   ```
   > **Note:** Update the path before `/chunk${1}.lhe` to point to your LHE file location. Keep the `/chunk${1}.lhe` suffix unchanged.
