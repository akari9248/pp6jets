# Run2026C TODO

更新：2026-09-16。分析目标：恰好六个 AK4 PUPPI jets，pT>30 GeV，前六个 |η|<5。用户手动提交生产；物理参数和种子分配见 [README.md](README.md)。

## 当前状态

- [x] ALPGEN→CP2/CP5 MiniAOD 入口已就绪；此前小样本全链路及 Condor dry-run 通过。当前入口直接读取原始 LHE。
- [x] 新批次 LHE Part1–5 的5000个文件名已核对；默认续产 Part6–10。真实事件数和输出完整性以批次统计为准。
- [x] 已确认 nominal generator weight 包含 merging 权重，分析须保留它。
- [x] 沿用用户已确认的 Data 检查：666万选中事件，无重复，L1/HLT prescale均为1；Data 要求 HLT_QuadPFJet30_v，MC 不加 trigger；Data 用2 fb⁻¹/6归一化。
- [x] 已确认比较数值：CP2 共934831事件，原始Σw=1026279.880658；归一化总量94700 pb，六喷注6439.477 pb；Data 19991.043 pb，Data/MC=3.104。后续不重复排查这些已确认项目。
- [x] 取消旧补写步骤和生成端写出改动，删除废弃脚本；按用户要求结束该项工作，后续不再主动提起。

## 生产与分析待办

- [ ] 按当前批次检查队列、失败/重试、输出完整性、重复种子，以及各 tune 实际保留事件数；旧状态快照不代表当前队列。
- [ ] 记录实际资源消耗，并核对 proxy、远端 MinBias 读取和 EOS 回传。
- [ ] 对比 CP2/CP5 的接受率、nominal/PS weights、PU/NPV/rho、pT6/pT7、DJR、六喷注和 BDT 形状。
- [ ] 对齐 SPS/data/data-mixing 的 jet 选择、额外 jet veto、JEC/JER、jet ID及触发效率；检查 ntuple 中 nominal/LHE/PS weights 和 true PU 的保存。
- [ ] 评估一致的 ME/PS/merging 尺度变化与多 multiplicity 覆盖；分别报告绝对截面及固定总截面后的接受度。
- [ ] 保持用户指定 MPI off，检查 SPS 模板软活动和喷注响应；明确 data-mixing 中 DPS/TPS 的处理。
- [ ] 随分析需要更新最终 lumi/PU 校准与覆盖；现有 PU profile 为0≤μ<10，不能覆盖更高 μ。

## 已有物理检查结论

以下均为生成级诊断，不能直接代替重建级 PUPPI 闭合检验；尚未找到3.104的单一来源。

- CP2 的2380事件小样本：CKKW-L 非零权重保留率12.23%；关闭 merging 后，归一化六喷注接受度中心值增加约1.34倍，统计有限。固定总Σw后，不能再额外补除总体 merging 效率。
- 0130→0118 的未合并 LO 对照：仅 PDF 密度平均比1.029，连同 ME αs⁶ 为0.556；后者不是最终 CKKW-L 截面比。
- 同一50790个输入上的 ALPGEN-style MLM/CKKW-L 六 jet 绝对截面比0.738±0.047。MLM 接受度更高但总截面更低，不能沿用 CKKW-L 的94.7 nb宣称增强。
- 固定 TMS=25，只降 ME pTmin至20：新增区间的1061个带权相空间点全部被输入 merging cut 拒绝；没有新增六 jet 贡献。单独降低生成 pT 门槛不能检验这一段的 shower 迁移。

外部历史材料保留在 EOS，生产不依赖它们：

- 迁移与生产配方：生产父目录的 `config/pp6jets_history_20260914_113441.tgz`。
- MLM 诊断：`/eos/home-z/zhye/ALPGEN_diagnostics/mlm_20260915/`。
- pTmin 边界诊断：`/eos/home-z/zhye/ALPGEN_diagnostics/ckkwl_pt20_20260915/`。
