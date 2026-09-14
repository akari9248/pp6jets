# LHE → CP2 / CP5 MiniAOD

只改 `submit.sh` 顶部：

```bash
FIRST_PART=1
LAST_PART=5
```

然后提交：

```bash
voms-proxy-init -voms cms -rfc -valid 192:00
./submit.sh
```

当前设置提交 Part1–5：每个 part、每种 tune 固定 1000 个任务，共 10000 个任务。脚本自动创建输出目录。提交前你自己检查输入完整性，脚本不扫描文件、不生成清单。重复运行会重复提交相同 part，续跑新批次时修改起止 part。

Part1–5 读取现有 `LHE_5M/chunk16823503_0.lhe` 至 `chunk16823503_4999.lhe`。后续 Part6 起读取 `LHE/Part<N>/chunk<全局编号>.lhe`，例如 Part6 是 chunk5000–5999。编号和种子自动按 part 分配，同一 part 的 CP2/CP5使用相同输入和种子映射。

EOS 父目录是 `/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV`。输出在 `MiniAOD_CP2_AQCDUP_v1/Part<N>/`、`MiniAOD_CP5_AQCDUP_v1/Part<N>/`；每任务只传回 `MiniAOD_<TUNE>_chunk<编号>.root` 和 `fullsim_<TUNE>_chunk<编号>.tgz` 日志包。

每任务固定读完全部 LHE，执行 LHE→GEN/SIM→DIGI/L1/HLT→RECO/RECOSIM→PAT/MiniAOD。MPI=off，补写原 ME 的 AQCDUP，PU 为 0≤μ<10；物理配置不变。固定请求 2 CPU、6000 MB 内存、20 GiB 磁盘、8 小时（workday）。本地小样本峰值约3 GB；完整任务的资源需求仍需实测。

proxy 使用当前账号的 `/tmp/x509up_u<UID>`；提交脚本自动复制到自己的 `$HOME/private/x509up_u<UID>`（权限600），JDL自动传到worker。worker使用Condor传入的 `X509_USER_PROXY`。MiniAOD 的现有种子布局支持完整 Part1–79，LHE 生成支持 Part1–30，其中 Part1–5 已保留给现有批次。

本地逐步 debug 仍可使用：`./fullsim.sh input.lhe workdir 20 all 1 CP2`。全链路物理配置和已完成验证摘要见仓库 TODO；详细历史日志在 zhye 本地 `archive/`，不随 GitHub 发布。

每个 part 只有第一个任务（Process=0）保存 Condor 的 log/out/err；其余999个任务均设为 `/dev/null`。第一个任务的调度 `.log` 保存在 AFS `log/`，`.out/.err` 随 `output_destination` 回传到该 part 的 EOS `log/`。CP2/CP5分别保留各自第一个任务。程序内部的物理/截面诊断归档仍随输出保存在 EOS。

协作提交：zhye 和 mitang 各自克隆 GitHub 的 `run2026c` 分支，在自己的克隆中执行本页命令即可。输入共用 zhye 的原始 LHE；mitang 的输出自动加 `/mitang/` 子目录，zhye 的输出不变。两人各自的日志和proxy互不共用。重复的 tune/part 不能当作独立样本合并。
