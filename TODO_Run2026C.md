# Run2026C ALPGEN 六喷注 SPS：当前 TODO

核对日期：2026-09-14。目标：恰好 6 个 AK4 PUPPI jets，pT>30 GeV、|eta|<5；SPS MC 对照 data-mixing TPS。用户已指定 MPI=off、PU 只使用 0≤mu<10，任务由用户手动提交。proxy统一使用提交主机 `/tmp/x509up_u170369`，用 `voms-proxy-init -voms cms -rfc -valid 192:00` 生成，不放工作区；经用户同意，提交脚本复制到 AFS `private/x509up_u170369`，JDL自动传给worker；worker使用Condor传入的 `X509_USER_PROXY`。

本文中的 `archive/` 路径指 zhye 本地历史证据，不随 GitHub 发布，生产无需这些文件。

这是当前清单，替代 `archive/TODO_Run2026C.md` 中迁移初期的状态和旧路径。实现配置、跑通小样本、完成物理验证分别记录。

**用户已授权两组各5000任务的CP2/CP5 MiniAOD生产，后续自行比较物理分布。两版AQCDUP补写和本地全链路已通过；按part展开的Condor dry-run通过（当前CP2/CP5各5×1000任务，及后续part接口），助手未提交。提交入口已简化为 submit.sh 顶部指定起止 part，每个 part 固定1000任务；无输入清单或预扫描，由用户检查源文件并运行脚本提交。**

## 已完成

- [x] 新工作区 CMSSW_16_0_8/src/pp6jets；移除 AMPLICOL，简化为 ALPGEN→MiniAOD。
- [x] ALPGEN 保留原 EL7/LHAPDF 环境；设置 ebeam=6800，重新生成 13.6 TeV LHE；没有仅修改旧 LHE 能量。
- [x] 六 parton 配置已明确：njets=6、ptjmin=25、etajmax=5、drjmin=0.3、ndns=315200、iqopt=1、qfac=1。其分析接受度仍待验证。
- [x] Run3 CP5 fragment 和 PS weights 已加载，comEnergy=13600；按用户决定保持 MPI=off；移除旧截面/效率常数。
- [x] nJetMax=4 已核对为 pp>jj 核心加四个额外 parton；nQuarksMerge=5 对应当前含 massless b 的样本。
- [x] 确认 CMSSW 将 merging weight 乘入 nominal generator weight；MiniAOD 中已实际读到非单位权重。分析端使用这些权重仍待验证。
- [x] 查到同一官方 low-PU 参考链并保存证据：GEN/SIM 16_0_8，DIGI/HLT/RECO 16_0_6，MiniAOD 16_0_8；EL8、Run3_2026、160X_mcRun3_2026_lowPU_v3、DB:Extended、DBrealistic、HLT:2026v11。
- [x] 替换 Run2 MinBias，保存 193 个官方 2026 GEN-SIM 文件；使用直接 mixing，25 ns、BX -5..3。
- [x] 查询 JetMET0–5 与 certification 的逐 LS 交集；去重并集 54,053 LS/48 runs。BRIL 校准排除 29 LS，计算使用 54,024 LS。
- [x] 用 BRIL 逐 bunch luminosity 和标准 makePileupJSON/pileupCalc 计算 true-mu profile；采用当前 hfet26v00、69.2 mb 参考假设。
- [x] 按用户决定截取 0≤mu<10，保留 0.1 分箱、重新归一化为 TH1F；平均 mu=4.38255、RMS=0.59640。实际 Poisson N 不作上限截断。
- [x] FullSim 已使用新 PU histogram，并开启 Poisson OOT；配置中的旧 averageNumber 已被整个替换。
- [x] RECO/PAT 使用 Run3 配方，移除 run2_miniAOD_UL；slimmedJetsPuppi 实际可读取。
- [x] 全链本地小样本：3 LHE→1 GEN-SIM→1 DIGI/HLT→1 RECO→1 MiniAOD；每阶段 FJR 无 FrameworkError。
- [x] 检查实际执行路径：GEN、SIM、DIGI、L1、DIGI2RAW、RAW2DIGI、L1Reco、RECO、RECOSIM 均执行成功；读取到 799 Geant4 SimTracks、810 条已执行的 HLT 路径及 MiniAOD 产品。
- [x] 新 PU 事例的 true mu=4.3174667，所有 BX 共享该 mu、实际 N 有 Poisson 波动。
- [x] Condor 容器、输入传输、EOS 目录、资源配置和 dry-run 已准备；新 PU ROOT/配置随任务传输。
- [x] 实现每 job/stage 的 seeds 和 job lumi、保留中间配置/FJR/日志、失败归档；较早的非零 job_id 小样本已验证 lumi 贯通。
- [x] 用户已提交 LHE cluster 16823503、5000 jobs。约 5M 是目标，尚非已核验的最终事件数。

证据：`archive/validation_Run2026C/`（旧本地测试产物已按用户要求清理，保留简短验证记录）；配方：`archive/migration/reference/production_recipe.json`；PU：`archive/pu2026/`。

## CP2/CP5 对照生产：当前执行状态

- [x] `ALPGEN/miniaod/submit.sh` 统一提交CP2/CP5，内部使用同一份 `condor.jdl`；每part固定1000任务，读完全部LHE、只传回MiniAOD和日志。输出目录、文件名、batch名分开。
- [x] `fill_aqcdup.py` 在worker派生副本中按原NNPDF31_lo_as_0130、member 0的alpha_s(SCALUP)补AQCDUP。核对输入参数和逐事件ME尺度；保留原LHE。已核对真实生产文件2380个事件只改AQCDUP，CP2/CP5得到相同修正文件哈希。
- [x] 同一生产LHE前20事件、同job ID=1：CP2各阶段20→2→2→2→2；CP5为20→7→7→7→7。两组MiniAOD的jets、nominal weights、PU、HLT及日志归档读取通过。
- [x] 两组使用相同seed映射；当前5000输入对应Part1–5，输出按tune和part分目录，提交脚本自动创建。后续LHE从Part6开始，种子按part自动分配。
- [ ] 用户自行检查源LHE数量和完整性，设置submit.sh的FIRST_PART/LAST_PART，运行脚本提交。已删除自动输入清单脚本和多余提交参数。
- [ ] 用户刷新有效proxy后运行miniaod/submit.sh，默认提交CP2/CP5各5个part。助手不提交。
- [ ] 用户后续比较两组的merging接受率、权重和六喷注/BDT形状。当前LO ME PDF + CP5的PDF/merging整体一致性仍需评估；AQCDUP补写不等于完整物理验证。

## 物理验证：随用户后续样本分析推进

- [ ] 核对exactly-six AK4 PUPPI pT>30、|eta|<5的生成边界迁移；检查pT6、pT7、第七jet veto、DeltaRmin、DJR、配对角度和平衡。
- [ ] 比较ME scale/PDF、PS/merging变化，并保证生成相空间覆盖；不要用本地20事件测试推断正式效率或哪版更好。
- [ ] 可用GEN级独立诊断定位差异；用户已决定先产两版MiniAOD，物理验证不再作为此轮配置工作的前置审批。

详细依据见 `archive/PHYSICS_REVIEW_Run2026C.md`。现有原始LHE无需因补AQCDUP重跑。

## 批任务运行与生产规模

- [ ] 用户提交两组各5000任务后检查早期任务、Hold原因和输出，记录足够的接受后事件数。
- [ ] 核验真实 Condor 排队/运行、proxy、远端 MinBias 读取、EOS stage-out、失败归档和输出完整性。dry-run 和本地 worker 运行不能替代实际批任务验证。
- [ ] 按实际完整LHE任务的耗时/内存/磁盘检查资源。20输入事件本地CP2 4:58、CP5 5:42，峰值RSS约2.91/2.97 GiB；启动开销大，不能线性外推全批。提交配置已调整为2 CPU、6000 MB内存、20GiB磁盘、8小时（用户指定workday）；原8000 MB被站点上调为3 CPU/9000 MB。完整任务内存仍待实测。
- [ ] 大样本再次检查 PU、顶点z/宽度、rho、jets、触发等分布；配置采用DBrealistic不等于这些分布已对数据验证。
- [x] 已固定读完全部事件、只保存MiniAOD和日志，并使用独立输出目录；用户提交后仍须避免重复提交同一组，不混用CP2/CP5事件。
- [ ] 冻结 LHE 派生版本、shower/PU 配置、软件版本和随机种子映射；新 source batch 使用不重叠的模拟 job ID 区间。

## P1：MiniAOD→ntuple→BDT 尚未贯通

- [ ] 建立并实际跑通2026 SPS MC ntuple配置。提供的2026RunC.json是data配置，isMC=false；不能原样用于MC。
- [ ] 保存/使用 nominal generator weight、LHE/PS weights、true PU；检查collection/process名和计数。不能因LHE unweighted而在BDT/模板中把全部MC权重设为1。
- [ ] 对齐 data/SPS/data-mixing TPS 的 AK4 PUPPI、pT>30、|eta|<5、恰好6个jets、overlap处理和第七jet veto。
- [ ] 确认2026适用的 AK4PFPuppi JEC/JER、jet ID、触发效率；只有使用b tagging时才需要对应校正。
- [ ] 按实际run核对目标HLT路径、L1 seeds、prescales、trigger objects及turn-on。菜单能够运行不等于数据触发效率已验证。
- [ ] 保留用户指定MPI=off，验证SPS控制区的jet/BDT形状；PU不补偿主碰撞soft MPI。检查BDT是否利用MC/data软活动或响应差异。
- [ ] 核对TPS data mixing的额外jets、配对规则及全局rho/NPV/MET定义；若物理模型只有SPS+TPS，核对DPS(2+4/3+3)的处理。此项不要求增加一个MC样本。

## P2：最终数据范围、PU与归一化

- [ ] 获取最终CRAB processedLumis，交叉检查成功处理的数据范围；目前提供的2026 CRAB results中仍未找到processedLumis.json。
- [ ] 生产PU按用户要求保持0≤mu<10；数据与TPS模板需要采用相应低PU定义，不能期待这份MC重加权覆盖mu>=10。
- [ ] 确认2026正式lumi/PU校准和推荐 minBiasXsec/系统变化。目前使用BRIL hfet26v00预备校准、69.2 mb参考假设，所选runs尚未进入当前normtag_PHYSICS。
- [ ] 最终数据范围确定后计算PU权重及系统变化，并用true PU/NPV/rho等检查。
- [ ] LHE任务结束后核验每个输出完整性和真实总事件数、失败/重试/重复种子情况；核实是否达到5M。
- [ ] 按用户之前安排，后续再汇总ALPGEN积分截面/误差、unweight效率、matching后的截面和sumw/sumw2；同一过程的分片截面不能直接相加。

## 生产快照（不同时间，不作同一快照相减）

- 2026-09-14较早查询：EOS扫描看到3056个该cluster的LHE及3056份日志归档；随后队列查询看到1932个Running、未看到Idle/Held。两个快照并非同一时刻，文件出现也不是完整性认证。
- 本轮已修改CP2/CP5的FullSim配置并跑完本地测试；助手没有提交或修改任何队列中的任务。

工作区整理：已删除旧本地测试目录、缓存和输入准备脚本；PU来源与物理审查移到archive。ALPGEN程序、输入卡、生产日志及EOS结果保留。后续只改提交脚本的起止part，无需生成清单。

2026-09-14批任务启动检查：MiniAOD clusters 16824020–16824029 使用提交机本地/tmp proxy，bigbird输入传输失败（Hold 13/2，NumJobStarts=0）。已改为private共享副本；已有任务需修改X509UserProxy后释放，尚待实际worker验证。

协作接口：提交脚本自动识别当前用户和UID；mitang使用自己的proxy和克隆，读取zhye LHE，输出到生产父目录下mitang子目录。用户要求发布GitHub run2026c分支，已准备协作权限与迁移说明；现有队列不自动变更。
