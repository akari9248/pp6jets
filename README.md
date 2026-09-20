# ALPGEN 六喷注 → Run2026C MiniAOD

只使用 **HTCondor** 生产。每个 Part 固定1000个任务，LHE、MiniAOD及对应日志均使用 Part 内的 `chunk0`–`chunk999`；缺失任务保留缺号。内部种子索引 `(part-1)*1000+chunk` 不作为文件名。

## 提交

在 lxplus 上使用自己的工作区和 proxy：

```bash
voms-proxy-init -voms cms -rfc -valid 192:00
```

修改对应 `submit.sh` 顶部的 `FIRST_PART`、`LAST_PART`，再运行：

| 阶段 | 提交入口 | 设置 |
|---|---|---|
| ALPGEN → LHE | `cd ALPGEN && ./submit.sh` | 输出 LHE_er3p0；当前范围 Part16–30 |
| LHE → MiniAOD | `cd ALPGEN/miniaod && ./submit.sh` | 当前 Part1–5、CP2、`JET_FILTER=off`；按需修改 |

续产使用未占用的 Part；重复提交同一输入会复用种子。LHE 支持 Part1–2147373；MiniAOD 当前支持 Part1–79，超过此范围需扩展 FullSim 种子映射。

LHE 每任务20M次生成尝试、20k×1000 warmup。MiniAOD 读完每份 LHE，执行 GEN/SIM→DIGI/L1/HLT→RECO→PAT，申请2 CPU、6000 MB内存、20 GiB磁盘和8小时。

## 输入、输出

EOS 父目录：

```text
/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV
```

```text
LHE_er3p0/Part<N>/
  chunk0.lhe ... chunk999.lhe
  chunk0_logs.tgz ... chunk999_logs.tgz

MiniAOD_<TUNE>_er3p0_GenJetFilter_v1/Part<N>/  # JET_FILTER=on
MiniAOD_<TUNE>_er3p0_GenJetMonitor_v1/Part<N>/ # JET_FILTER=off
  MiniAOD_<TUNE>_chunk0.root ... MiniAOD_<TUNE>_chunk999.root
  fullsim_<TUNE>_chunk0.tgz ... fullsim_<TUNE>_chunk999.tgz
```

- MiniAOD 固定读取 zhye 的 `LHE_er3p0`。其他账号的输出放在父目录的 `<账号>/` 下。
- 提交脚本把 `/tmp/x509up_u<UID>` 复制到 `$HOME/private/x509up_u<UID>`，由 Condor 传给 worker。
- 每个 Part 仅 Process=0 保留调度日志；`log/` 由提交脚本自动创建。每个任务的阶段配置、事件数、退出状态和运行日志随 tgz 回传。
- FullSim job ID 为 `(part-1)*1000+chunk+1`，同一输入的 CP2/CP5 使用相同种子映射。文件改名不改变事件身份或种子。
- LHE Part1–30 保留原种子；Part31起使用新的两阶段种子对分配。具体公式在 `ALPGEN/run_condor.sh` 中；不同版本若使用相同 Part/chunk，仍会复用种子。

## 生产配置

| 项目 | 设置 |
|---|---|
| ME | 13.6 TeV，六 parton，pT>25 GeV、\|η\|<3.0、ΔR>0.3 |
| PDF / 尺度 | NNPDF31_lo_as_0130；iqopt=1、qfac=1，μR=μF=√(ΣpT²/6) |
| Shower | CP2 / CP5，MPI off，保留 nominal 和 PS weights |
| Merging | CKKW-L，pp>jj，nJetMax=4，TMS=25 GeV，D=0.4，nQuarksMerge=5 |
| LHE 接口 | 直接读取原始 LHE |
| Gen-jet filter | 至少六个可见粒子 anti-kT R=0.4 gen jets，pT>20 GeV、\|η\|<3.0 |
| PU | 0≤μ<10 的 Run2026C profile；实际 Poisson 交互数不截断 |
| 软件 | ALPGEN 用EL7；GEN/SIM、MiniAOD用16_0_8，DIGI/HLT/RECO用16_0_6，FullSim用EL8 |
| Run3 配方 | Run3_2026，160X_mcRun3_2026_lowPU_v3，HLT:2026v11 |

`ALPGEN/alpgen/` 包含运行及重建程序需要的源码、可执行文件和库数据。FullSim 使用 CVMFS release 和容器，gen-jet filter 使用现有 CMSSW 插件，无需自行编译。

## 筛选与归一化

`JET_FILTER=on` 在 shower/hadronization 后、Geant4 前筛选；`off` 记录同样的筛选决定，但不拒绝事件。两种模式输出到不同目录。筛选不 veto 第七个 gen jet，不能代替最终重建级六喷注选择；物理接受度需要用 on/off 对照验证。

生产必须保留：

- `CentralGenJetFilter_cff.py`：筛选决定和 `GenFilterEfficiencyProducer`。
- `filter_summary.py`：导出每个 job 的筛选前后事件数、Σw、Σw²，包括零通过事件的任务。
- `central_jet_filter_summary.json`：随 `fullsim_*.tgz` 回传；先相加各 job 计数与权重，再计算总效率，不能平均各 job 的效率。
- MiniAOD 内的 `GenFilterInfo`、各 process 的 `TriggerResults` 和 nominal generator weight。

筛选统计的分母是 shower/merging 后的事件。若使用 filter 前截面，归一化为 `sigma_before * sumw_selected / sumw_before`；若截面已包含 filter，则不能再乘一次筛选效率。新的生成 cuts 需要独立确定截面，不能直接套用旧样本的94.7 nb。

## 运行与检查工具

```bash
# 本地 ALPGEN
./ALPGEN/run_alpgen.sh /tmp/alpgen_run
# FullSim：输入、工作目录、事件数（-1=全部）、阶段、job_id、tune、jet_filter
./ALPGEN/miniaod/fullsim.sh input.lhe /tmp/fullsim_run -1 all 1 CP2 on
# 统计某个 Part 的 LHE；也可传入 MiniAOD 目录统计 ROOT
./count_events_part.sh /eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV/LHE_er3p0/Part1 '*.lhe'
# 查询有自己任务的 schedd
./detect_active_schedds.sh
```
