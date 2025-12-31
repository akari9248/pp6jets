# pp6jets
## 安装说明
这个git仓库要放在CMSSW下面，目前对于2017年的模拟，选择的是CMSSW_10_6_28_patch1:
```
cmssw-el7
cmsrel CMSSW_10_6_28_patch1
cd CMSSW_10_6_28_patch1/src
git clone git@github.com:akari9248/pp6jets.git
```

## ALPGEN
在进入到pp6jets文件夹之后，运行如下命令：
```
cmsenv
scram tool info lhapdf
```
会得到类似于如下的输出：
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
然后打开`alpgen/alplib/alpgen.mk`，找到`CONFIG_FILE_DIR=xxx`这一行，把它修改成你上面得到的`PATH`，随后运行：
```
make cleanall
make gen
```
编译成功，可以运行`./Njetgen`测试一下。
可执行文件为`run_alpgen.sh`，ALPGEN的输入在其中修改，具体在`cat > config_step1.txt << EOF` 和 `cat > config_step2.txt << EOF` 这两处。
在提交condor任务的时候记得 Ctrl + A + D 退出singularity，运行`condor_submit condor_algen.jdl`即可。记得在每次提交任务的修改`output_destination`防止覆盖。

# PYTHIA
随后是LHE文件对接到PYTHIA，首先把`JME-RunIISummer20UL17GEN-00006-fragment.py`放在`CMSSW_10_6_28_patch1/src/Configuration/GenProduction/python/`下面（没有这个文件夹的话自己创建一个）。如果你有其他的配置文件，也记得放在这个文件夹下面，并且修改`miniaod`文件夹中fullsim.sh对应的
```
cmsDriver.py Configuration/GenProduction/python/JME-RunIISummer20UL17GEN-00006-fragment.py --python_filename GEN.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN --fileout file:gen.root --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN --geometry DB:Extended --era Run2_2017 --mc -n ${eventsnum} --filein file:MCDBtoEDM_NONE.root
```
在放完fragment之后需要编译：`scram b - j 8`。
然后进入`miniaod`文件夹，在提交condor之前有几个注意点：
1. 记得把`fullsim.sh`中home路径改成你自己的，并且安装下面的那些CMSSW
2. 运行`voms-proxy-init -voms cms -rfc -out x509up `
3. 修改`condor.jdl`中的`output_destination`，这是文件输出路径
4. 修改`fullsim.sh`中
```
cmsDriver.py MCDBtoEDM --conditions 106X_mc2017_realistic_v6 -s NONE --eventcontent RAWSIM --datatier GEN --filein file:/eos/cms/store/group/phys_smp/ec/shuangyu/pp6j_15GeV/Part1/chunk${1}.lhe -n ${eventsnum}
```
这一行命令，主要是`--filein`后面的参数，这是LHE文件的input，`/chunk${1}.lhe`不用管，前面记得改。
