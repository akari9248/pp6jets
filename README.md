# ALPGEN 六喷注 → Run2026C MiniAOD

## 提交

在 lxplus 上使用自己的工作区和 proxy：

```bash
git clone --branch run2026c https://github.com/akari9248/pp6jets.git
cd pp6jets
voms-proxy-init -voms cms -rfc -valid 192:00
```

修改对应 `submit.sh` 顶部的 `FIRST_PART` / `LAST_PART`，再运行。每个 part 固定1000任务，由用户手动提交。

| 阶段 | 提交入口 | 当前默认 |
|---|---|---|
| ALPGEN → LHE | `cd ALPGEN && ./submit.sh` | Part6–10 |
| LHE → MiniAOD | `cd ALPGEN/miniaod && ./submit.sh` | Part1–5，每个 part 同时跑 CP2/CP5 |

现有输入为新批次 `LHE/Part1`–`Part5`，共5000份文件；实际事件数用统计脚本检查。提交入口不预扫描输入。续产使用未占用的 part；重复提交会重复使用输入和种子。LHE 支持 Part1–30，MiniAOD 支持 Part1–79。

ALPGEN 每任务20M次生成尝试、20k×1000 warmup。MiniAOD 读完每份 LHE，执行 GEN/SIM→DIGI/L1/HLT→RECO→PAT，申请2 CPU、6000 MB内存、20 GiB磁盘和8小时。

## 输入、输出和账号

EOS 父目录：

```text
/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV
```

- LHE：`LHE/Part<N>/chunk<编号>.lhe` 及 `chunk<编号>_logs.tgz`。
- MiniAOD：`MiniAOD_CP2_AQCDUP_v1/Part<N>/` 和对应 CP5 目录；保留现有目录名。输出 ROOT 和 `fullsim_<TUNE>_chunk<编号>.tgz`。
- 全局编号为 `(part-1)*1000 + Process`；同一输入的 CP2/CP5 使用相同种子映射。
- MiniAOD 固定读取 zhye 的 LHE；其他账号的输出自动放在父目录的 `<账号>/` 下。各自使用自己的克隆、队列和 proxy。
- 提交脚本把 `/tmp/x509up_u<UID>` 复制到 `$HOME/private/x509up_u<UID>`，由 Condor 传给 worker。
- 每个 part 仅 Process=0 保存 Condor 日志；调度日志在本地 `log/`，out/err 在 EOS 对应 part 的 `log/`。每个任务的程序日志仍随输出归档。

不同账号需协调 part 分配；相同 tune/part 的重复输出不能作为独立样本合并。其他账号生成的 LHE 位于其自身输出目录，使用时需相应修改 MiniAOD 的输入来源。

## 当前配置

| 项目 | 设置 |
|---|---|
| ME | 13.6 TeV，六 parton，pT>25 GeV、\|η\|<5、ΔR>0.3 |
| PDF / 尺度 | NNPDF31_lo_as_0130；iqopt=1、qfac=1，μR=μF=√(ΣpT²/6) |
| Shower | CP2 / CP5，MPI off，保留 nominal 和 PS weights |
| Merging | CKKW-L，pp>jj，nJetMax=4，TMS=25 GeV，D=0.4，nQuarksMerge=5 |
| LHE 接口 | 直接读取原始 LHE |
| PU | 0≤μ<10 的 Run2026C profile；实际 Poisson 交互数不截断 |
| 软件 | ALPGEN 使用 EL7；GEN/SIM、MiniAOD 用16_0_8，DIGI/HLT/RECO用16_0_6 |
| Run3 配方 | Run3_2026，160X_mcRun3_2026_lowPU_v3，HLT:2026v11 |

worker 使用 CVMFS release 和容器；ALPGEN 可执行文件随仓库提供，无需新建 CMSSW。`ALPGEN/alpgen/` 保留配套源码、编译文件和库数据。

## 本地工具

```bash
# ALPGEN 小样本
./ALPGEN/run_alpgen.sh /tmp/alpgen_test
# FullSim：输入、工作目录、事件数、阶段、job_id、tune
./ALPGEN/miniaod/fullsim.sh input.lhe /tmp/fullsim_test 20 all 1 CP2
# 统计 LHE；或传入具体 MiniAOD Part 目录统计 ROOT
./count_events_part.sh all '*.lhe'
# 查询有自己任务的 schedd
./detect_active_schedds.sh
```

当前状态和待办见 [TODO_Run2026C.md](TODO_Run2026C.md)。工作区仅保留这两份说明文件。
