# ALPGEN 六喷注 → Run2026C MiniAOD

两个提交入口，每个 part 固定 1000 个任务。**只改脚本顶部的 FIRST_PART / LAST_PART，再运行脚本。** 助手不代交任务。

在 CERN lxplus 上克隆到自己的 AFS 工作目录并生成自己的 proxy：

```bash
git clone --branch run2026c https://github.com/akari9248/pp6jets.git
cd pp6jets
voms-proxy-init -voms cms -rfc -valid 192:00
```

不需要新建或编译 CMSSW；worker 使用 CVMFS 上的正式 release 和容器。ALPGEN 可执行文件随仓库提供。

当前 5000 份 LHE 已由 cluster 16823503 生产，对应 Part1–5，逐文件核验共 8,034,304 个事件。提交 CP2 和 CP5 MiniAOD：

```bash
cd ALPGEN/miniaod
./submit.sh
```

默认 Part1–5、CP2/CP5 各 5000 个任务。具体输入、输出见 [MiniAOD 说明](ALPGEN/miniaod/README.md)。不扫描输入数量，由你提交前检查。

后续生成新 LHE：

```bash
cd ALPGEN
./submit.sh
```

LHE 脚本默认 Part6–10，避免重复现有 Part1–5 的种子；输出 `LHE/Part<N>/chunk<全局编号>.lhe` 及统计日志。每任务固定 20M 次生成尝试、20k×1000 warmup；最终 unweighted 事件数需后续统计。现有 LHE 种子范围最多支持完整 Part30。脚本自动创建输出目录，继续生产使用未用过的 part；重复运行会重复提交。

物理设置保持：13.6 TeV，六 parton，pT>25、|eta|<5、ΔR>0.3，ME PDF NNPDF31_lo_as_0130；两版 shower 为 CP2/CP5、MPI=off、补真实 AQCDUP、0≤μ<10 PU。GEN/SIM 和 MiniAOD 用 CMSSW_16_0_8，DIGI/HLT/RECO 用 16_0_6。

本地 debug：`ALPGEN/run_alpgen.sh work/test`；FullSim 的逐步调试入口保留。当前待办见 [TODO_Run2026C.md](TODO_Run2026C.md)。历史资料已压缩归档到 EOS，AFS 的 archive 目录已清理；生产所需的 PU ROOT、fragment 和 MinBias 清单都随仓库提供。现有生产日志和 EOS 结果保留。

每个 part 只有第一个任务（Process=0）保存 Condor 的 log/out/err；其余999个任务均设为 `/dev/null`。第一个任务的调度 `.log` 保存在 AFS `log/`，`.out/.err` 随 `output_destination` 回传到该 part 的 EOS `log/`。CP2/CP5分别保留各自第一个任务。程序内部的物理/截面诊断归档仍随输出保存在 EOS。

提交脚本通过 `id` 和 `$HOME` 自动识别账号，将自己的 `/tmp/x509up_u<UID>` 复制到自己的 AFS `private/x509up_u<UID>`（权限600），JDL 用 `x509userproxy` 传给 worker。bigbird 无法直接读取 lxplus 的 `/tmp`。

## zhye / mitang 协作

两人分别在自己的 AFS 克隆仓库、生成自己的 proxy、运行 `ALPGEN/miniaod/submit.sh`。只改脚本顶部 FIRST_PART / LAST_PART，无需修改 JDL 或替换代码中的用户名。

- 输入固定复用 zhye 的 Run2026C LHE；无需复制约10 GB的输入文件。
- zhye 的输出路径保持原样。
- mitang 的输出为 `.../pp6j_25GeV/mitang/MiniAOD_CP2_AQCDUP_v1/Part<N>/` 和对应 CP5 目录；目录位于 zhye 的 EOS 空间，已单独授权 mitang。
- AFS 日志写入各自克隆目录；各自使用自己的 Condor 队列和 proxy。
- 每次提交都会跑指定 part 的 CP2 和 CP5。两人需分配不同 part 或明确接手旧任务；相同 tune/part 的重复输出使用相同 LHE 和种子，不能作为独立样本合并。现有队列不会因拉取新代码而取消。

默认 MiniAOD 的 Part6 以后仍读取 zhye 的 `LHE/Part<N>`。若使用 ALPGEN 入口由 mitang 生成新LHE，其输出在 `mitang/LHE/Part<N>`；需要使用这批新源时，应在 MiniAOD 提交脚本中明确修改 LHE 来源。LHE part 编号与种子需统一分配，不能按账号重新从头计数。
