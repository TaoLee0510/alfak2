# alfak2 交接给 Claude

交接日期：2026-11-10（America/New_York）。实验结果完成于 2026-11-05；本次只做文件复核与交接，没有启动新的模型实验。

## 1. 接手目标与当前结论

用户的研究目标是提高**所有距离层级的 fitness 点估计及不确定性质量**，包括观测节点、1-hop、2-hop 和更远场，不能只改善某个距离或只缩窄区间。

上一轮具体目标已经完成：保持原统计目标，修复联合后验数值近似，完成解析/两节点验证、两个8节点独立后验核对、两个27节点资源预检，以及48例固定布局机制对照。**更大的“所有距离精度同时提高”目标没有完成。** 不要把48例运行成功解释为精度目标通过。

- 仓库：`/Users/4482173/Documents/GitHub/alfak2`
- 当前分支：`dev`
- 当前 HEAD：`8d82f03c50f5fbda2d3611bb078580e1b5afdbea`
- 上轮没有 commit、push、HPC 提交或生产默认切换。
- 当前公开入口 `fit_alfak2()` 仍默认 `pipeline="two_layer"`；`joint_count` 为显式选项。数值校正 v16 是验证目录中的独立原型，**没有成为通用包接口或默认引擎**。
- 本次进程查询未发现 `/alfak2/validation/observation_design/` 相关运行进程；这不代表检查或停止了其他项目/HPC任务。当前没有需要恢复的半途模型批次。
- 本轮接手应从新的模型/信息条件研究开始，不要重新从头调 scalar quadrature 或重跑已完成48例。

首先阅读[最终报告](/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/importance/REVIEW.md)，然后阅读[完整协议与版本变更](/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/importance/PROTOCOL.md)。协议前部是最初方案，后续修订累积到 v16；最终执行代码以冻结快照为准。

## 2. 2026-11-10 的当前文件复核

本次新审计保存在独立目录，未改写旧结果或旧审计：

- [审计摘要](/Users/4482173/Documents/GitHub/alfak2/validation/handoff_20261110/audit_summary.json)
- [逐文件审计](/Users/4482173/Documents/GitHub/alfak2/validation/handoff_20261110/integrity_audit.csv)
- [交接前 Git 状态](/Users/4482173/Documents/GitHub/alfak2/validation/handoff_20261110/git_status_before_handoff.txt)

| 复核范围 | 匹配 | 缺失 | 内容改变 |
| --- | ---: | ---: | ---: |
| 最终产物清单 | 413 | 0 | 0 |
| 唯一消耗源码/运行库路径 | 338 | 17 | 0 |
| 仓库内运行库备份 | 4 | 0 | 0 |
| 两个用户原始文件 | 2 | 0 | 0 |

17个缺失路径全部位于 `/private/tmp`，包括旧 native DLL、各版本 Rcpp 缓存和旧 observation TMB DLL。**最终结果没有丢失或改变；旧临时编译环境不能假定仍然存在。** 各历史缓存二进制并非全部另有备份，不要声称每个历史缓存都能位级恢复。

关键备份目录：`/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/runtime_archive`

其中 `alfak2.so` 是所需 native DLL 的完好副本；另存有最终 importance TMB、核心 Rcpp 和参考采样器二进制。参见该目录的 `README.md`、`backup_manifest.csv` 和 `session_info.txt`。这些是 Mac arm64 / R4.4 平台产物，不是跨平台安装包。

**不要直接运行旧 `check_integrity.R` 来“修复”审计。** 它会先覆盖 `importance_v1/consumed_source_integrity.csv`，再因缺失临时路径而失败；这会改写已经哈希归档的证据。应在新的审计目录记录当前状态，区分临时路径缺失、备份匹配和内容变化，保持旧清单原样。

## 3. 已完成的数值工作及其边界

最终计算快照：

`/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_source_v16`

不要编辑这个目录。根原型目录 `validation/observation_design/importance` 可用于理解，但确定旧结果的行为要读冻结版本。

v16 的主要流程：

1. 把 Gaussian landscape 系数解析边缘化，得到 `g ~ N(offset,K)`，白化为标准正态坐标；不重定义原先验。
2. 联合 mode/Hessian 只用于提议分布。初始16,384点 pilot 若表现不足，转为四条独立、只使用拟合数据的 ESS 训练链，每条2,000 warmup +60,000 retained。
3. 提议变换包含基线方向的自然样条条件均值、16个基线分组的完整条件协方差白化、边际 Gaussianization 与8层 PCA/边际变换；逆变换和 Jacobian 明确计入密度。
4. 复杂分支混合：0.30 Gaussian、0.20 沿 mean-log-birth 梯度上移1.5 SD的 Gaussian、0.30 独立坐标 t7、0.19 全局多变量 t7、0.01 原始标准正态先验。该偏移只用于覆盖绝对 birth 矩的尾部，不改变目标后验。简单 pilot 通过的分支仍使用原0.80/0.19/0.01混合。
5. 最终估计来自独立 iid 批次的**原始权重**；`loo::psis` 只提供诊断，没有对估计使用平滑权重。pilot/训练链不进入最终估计，也不使用参照 MCMC 或真值训练提议。
6. 八个批次逐档增加至每批最多131,072点，最多1,048,576点。16,384行分块计算矩和 MCSE，逐列计算加权分位数，控制内存。

门槛：每个普通完整案例外部300秒期限、峰值RSS≤6 GiB；原始权重ESS≥10,000；八个批次及合并样本的 k 全部≤0.5；均值MCSE≤0.005后验SD，SD的MCSE≤0.01后验SD。包括 log-birth、birth、相对/全局中心化/边对比、未知初始 logit 坐标的诊断。均值/SD MCSE同时用影响函数与独立批次交叉检查。

历史失败全部保留：简单提议尾部失败、Rcpp并发缓存初始化错误、Jacobian有限差分跨结点问题、绝对birth尾部MCSE失败、百万点诊断内存超限等。不能删除失败、放宽门槛或换掉难例来通过。

最终验证107/107通过。8节点参照各四条链，每条分别80,000/160,000 retained；关键对比均值最大误差为0.0136/0.0092参照SD，SD比约0.994–1.006。27节点两例预检全部通过。48/48最终案例均通过；最长200.224秒、峰值5.516 GiB。

**这证明的是已声明的数值检查通过，不是生成真值精度、频率学覆盖率或任意大图的保证。** `is_spec()` 明确只接受一个 culture、两列计数、时间0和T；不要直接拿它替换包内多时间点/多culture/可变beta/增长观测等通用路径。

## 4. 模型、数据与评价定义：不要悄悄改变

当前机制实验为3个染色体、CN1:3，共27节点；`beta=.005`，`T=1.5`，初始总群体1,200，typed cells总预算600，t0/T分别240/360。固定 opt_single 起始布局相对2.2.2的shell质量为 `[.45,0,.35,.20]`，组内均分；实际初始整数群体来自原模拟记录。

- 数据生成：有限种群、相关双子代的 exact Gillespie。
- 拟合似然：确定性均值传播后的 multinomial counts，尚未显式包含有限分枝过程协方差。
- `A = (2*P^T-I) diag(exp(g))`，`n(T)=exp(A*T)n(0)`；`P`为保留域内单子代核。域外损失不能重新分配、不能把 `P` 行归一化。
- 原初始组成先验：标准正态白化坐标经 `2*orthonormal Helmert` 得到 logit，再softmax；未知组成允许图上全部状态。
- 本实验固定先验尺度：baseline .5、main .07、curvature .04、interaction .04、residual .125、edge .25。读取冻结 context 中的原 offset/K，不要按新数据重算或重新居中 offset。
- fitness 指标使用 birth `b=exp(g)`，相对量为 `b - reference_weight^T b`。参照权重沿用原 nominal-A，不能改成新布局或新样本的bulk。
- 主距离分组沿用每个重复的原 control300 观测支持，不能随候选观测支持移动。当前域主要是距离0/1/2；没有完成任意更远场验证。
- Shape 对全27节点误差**统一去中心一次**；不能逐shell去中心。边按两个端点距离的最大值归组。
- 汇总：生物学重复→景观→家族→家族等权。6个独立景观、每景观2个重复；48例不是48个独立景观。
- 真值只进入生成/计分；`inference/` 中不含真实fitness。真实初始组成只允许出现在标注oracle的C/D组。

## 5. 48例机制实验与关键结论

| 组 | 生成过程 | 推断中的初始组成 |
| --- | --- | --- |
| A_finite_unknown | 有限分枝 | 未知 |
| B_mean_unknown | 确定性均值 | 未知 |
| C_finite_oracle | 有限分枝 | 真实组成，oracle |
| D_mean_oracle | 确定性均值 | 真实组成，oracle |

每个重复配对实际初始整数群体、t0读数、终点抽样uniforms、fitness参照和原控制距离。两个通过的A组预检以完全一致的源码/数据/seed复用；其余46例新拟合。**使用完整48行manifest，而非仅46行run_manifest进行总结。** A组两个结果目录实际在 `preflight_cases_v16`，应以manifest的 `output` 为准。

数值校正A组相对旧 opt_single：

| 原控制距离 | RMSE变化 | Shape变化 | IS95变化 | 新95%覆盖率 |
| --- | ---: | ---: | ---: | ---: |
| 0 | +0.180% | -0.277% | +3.584% | 87.15% |
| 1 | +0.110% | +0.148% | +1.466% | 82.68% |
| 2 | -0.090% | +0.386% | -0.124% | 89.31% |

原所有距离开发门槛全部未通过。数值校正没有实质解决生成真值误差，区间约窄2.1%，也没有修复欠覆盖。

机制线索：C−A 的RMSE在距离0/1/2下降7.22%/12.82%/8.49%；B−A下降3.85%/4.51%/2.58%。初始组成的信息条件值得优先研究，但这不是可部署改善。C在距离2的覆盖率86.67%，反而低于A的89.31%；D为91.67%。D的1-hop覆盖率也只有92.79%。

解释限制：B−A/D−C是有限过程替换为均值传播的综合影响，包含波动与均值传播近似之差；C−A/D−B是知道真实初始组成的综合影响，也包含真实结构零与连续先验支持的区别。不能分解成某种噪声的纯因果百分比。D的残余覆盖偏离不能唯一归因于似然错误：还涉及先验收缩、可识别性、少量景观及Bayesian区间的频率学性质。

## 6. 权威路径索引

以下均位于 `/Users/4482173/Documents/GitHub/alfak2`，完整入口如下：

| 用途 | 绝对路径 |
| --- | --- |
| 最终报告 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/importance/REVIEW.md` |
| 最终数值源码 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_source_v16` |
| 数据生成/计分/分析源码 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_analysis_source_v2` |
| 报告/图表源码 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_report_source_v3` |
| 48例manifest | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/manifest.csv` |
| 原context | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/context.rds` |
| 48例状态/资源 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/all_fit_status.csv` |
| 全部逐例指标 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/all_case_metrics.csv` |
| 家族/景观/总体指标 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism` 下的 `family_metrics.csv`、`landscape_metrics.csv`、`overall_metrics.csv` |
| 配对变化/交互作用 | 同上目录的 `paired_landscape_differences.csv`、`paired_interaction.csv`、`overall_risk_changes.csv` |
| 数值校正开发门槛 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/numerical_correction_development_gate.csv` |
| 小图参照与比较 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/references`、`/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/small_cases_v16_analysis` |
| 验证明细 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/validation_v16/validation_tests.csv` |
| 最终413产物哈希 | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/final_artifact_sha256.csv` |
| 图表PDF/PNG | `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/mechanism_summary.pdf`、同目录 `mechanism_summary.png` |

每例 `result.rds` 有矩、协方差、MCSE、加权分位数、提议和seed等信息。对于n>8，**`fit$samples`特意未重复存储**；完整原始样本及未归一化log权重在 `fit$sample_checkpoint` 指向的 `sample_checkpoint.rds`，按 `batches` 存储。不要把缺少 `fit$samples` 当成数据丢失。数据含计分真值的副本在manifest的 `scoring_file`，拟合输入在 `file`。

更早背景按需读取：

- `/Users/4482173/Documents/GitHub/alfak2/doc/accuracy_upgrade.md`：包入口及早期升级，部分算法描述早于最新验证原型。
- `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/constrained/REVIEW.md`：设计筛选与近场信息约束，无所有距离胜者。
- `/Users/4482173/Documents/GitHub/alfak2/validation/observation_design/mechanism/REVIEW.md`：小图参考、条件Laplace/基线证据权重定位；单纯加密一维积分节点不是主要修复方向。
- `/Users/4482173/Documents/GitHub/alfak2/validation/accuracy_round2/README.md`：新增开发/校准/正式测试边界。

## 7. 读取与复现入口

只读进度入口（现在应显示48通过、无active；不会重新拟合）：

```sh
python3 /Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_report_source_v3/progress.py /Users/4482173/Documents/GitHub/alfak2/validation/observation_design/results/importance_v1/mechanism/manifest.csv
```

环境基础：R4.4.2、Mac arm64；当时TMB1.9.19、Rcpp1.1.0、Matrix1.7.3、posterior1.6.1、loo2.8.0、mvtnorm1.3.3。当前重新运行前仍需核对本机版本。

临时native库已缺失，加载旧算法时先设置完好的备份路径；TMB/Rcpp可由冻结源码重新编译缓存：

```sh
export ALFAK2_REPO=/Users/4482173/Documents/GitHub/alfak2
export ALFAK2_IS_SOURCE="$ALFAK2_REPO/validation/observation_design/results/importance_source_v16"
export ALFAK2_REFERENCE_DLL="$ALFAK2_REPO/validation/observation_design/results/importance_v1/runtime_archive/alfak2.so"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```

`model.R` 的 `is_load(repo,src)` 会构建/加载新importance缓存与冻结参考采样器；内部有编译缓存目录锁。首次加载应串行预热；不要盲删锁，先确认owner进程是否仍存在。重新编译的二进制哈希未必与历史相同，应单独记录，不能篡改旧provenance。

需要另起**已经获准的新实验**时可复用 `case_controller.py --manifest NEW_MANIFEST --part 0/1`，manifest列至少为 `case_id,file,context,output,seed`。所有 `output` 必须是新的、不存在的目录，最多两个model进程。不要对旧manifest直接启动控制器，它会拒绝覆盖。

单例脚本 `run_case.R` 读取 `ALFAK2_IS_CASE`、`ALFAK2_IS_CONTEXT`、`ALFAK2_IS_CASE_OUTPUT`、`ALFAK2_IS_SEED`；它自身没有外部300秒强制终止器，应通过控制器运行普通完整案例。资源完整耗时以控制器的 `resource.csv` 为准；代码中300.3秒是轮询记账容差，不是扩大科学预算的许可。

`prepare_mechanism.R`、`analyze_mechanism.R`、`report.py`等旧脚本会写入目标目录。不要以旧 `ALFAK2_IS_OUTPUT` 盲目重跑；新分析也要另存，保留旧报告与哈希。阅读既有CSV/RDS不需要恢复全部临时库。

## 8. Claude 接下来应推进什么

先完成只读复核，接着把下一阶段范围和具体候选写成可确认的方案。上轮具体目标已经完成，**没有授予生物学模型变更、生产切换或HPC批次的无限授权**。

建议按以下顺序研究，而不是继续仅追求更紧积分容差：

1. **初始组成的现实信息约束。** 核对真实工作流中实际可得的播种构成、允许支持集、N0和初始抽样信息。区分“未检出”与“确实为零”，不要从 `y0=0` 自动断定真实组成零。不要把模拟真实pi0作为部署输入。先清楚说明新的先验/数据条件相对旧模型改变了什么。
2. **有限分枝过程的状态相关协方差。** 从相关双子代机制出发，推导/验证均值、协方差及抽样观测组合；保留域外损失定义。均值群体之比不等于随机频率的均值，不能把B−A全部称为纯噪声效应。先做解析/两节点/小图验证再扩到27节点。
3. **先验收缩、可识别性与分层校准。** D组仍有覆盖偏离，需要检查而非立即认定某个误设。总种群增长与typed sampling depth是不同数据；若加入增长信息，需要真实可获得的计数和新信息预算说明。
4. 用新增独立开发景观验证候选。保持共同fitness参照、原控制距离、完整案例可用性及家族安全检查。没有可信开发候选前，不启用预留正式测试。

原所有距离开发判据应继续显式报告：点RMSE不劣；各距离Shape至少改善5%；Edge不劣；IS50/80/95不劣；每个家族各距离的RMSE/Shape/Edge/IS95风险比≤1.05；全部预测/区间可用、数值与资源门槛通过。覆盖率和宽度单独报告，不能用“更窄”代替“更准”。不宜为了通过而更换指标、距离定义、参照或删除案例。

新的统计目标必须建立新的正确性参考；v16在旧目标上的证明不能自动转移到新似然或新先验。当前27节点诊断也不能替代更大图/更远场验证。

## 9. Git、结果与用户操作规则

工作树有大量未提交修改、未跟踪代码/验证目录，以及35个已删除的旧benchmark路径。它们在交接前已存在，不能一概归为上一轮新增，也不能擅自恢复、reset、clean或打包提交。本次交接只新增本文件和 `validation/handoff_20261110` 审计目录。

尤其保留用户原始修改：

- `/Users/4482173/Documents/GitHub/alfak2/R/benchmark-second-layer.R`
- `/Users/4482173/Documents/GitHub/alfak2/doc/alfak2_method.tex`

本次已核对二者仍与既存用户原始哈希一致。所有 `results/*_source_v*` 冻结快照及历史失败目录保持不动；修改原型后写新版本、新输出目录，不修改旧快照。

用户约定：

- 说“制定方案”：先复述需求，再提出方案，等待确认；确认后再执行该方案。
- 说“写commit”：只写当前未提交改动的commit文本，不提交。
- 说“commit”：commit并push；“push commit”是push到origin；“commit, hpc”还包括HPC pull/同步。
- 当前交接不授权上述Git/HPC动作。HPC、workstation与外部通信都不应自行启动；凭据不放入交接文件。
- “停止所有任务”指本会话提交的任务，不是HPC全体作业；停止前核对归属。
- 本机没有git-lfs；只有获准push时才处理相应hook，不要借机安装或改动全局配置。
- “写信”是给邮件内容，不操作邮件客户端发送。
- 回复使用中文；每次结束加America/New_York时间戳。提供可审阅的决策和验证摘要即可。

## 10. 可直接发送给 Claude 的启动指令

> 请接手 `/Users/4482173/Documents/GitHub/alfak2`，先阅读 `HANDOFF_CLAUDE_2026-11-10.md` 和其中的最终报告、2026-11-10审计。我的总目标仍是提高所有距离层级的fitness点估计和不确定性质量。上一轮107项验证和48例机制实验已完成，但所有距离改善门槛未通过。请保留旧数据、冻结源码、失败记录及工作树现有改动；特别注意17个临时库路径已缺失，但最终413产物与4个运行库备份完好。先核对状态，再复述需求并制定下一阶段可执行方案，优先研究实际可得的初始组成约束、有限分枝协方差及分层校准，等待我确认。不要重跑旧48例、覆盖旧审计、把oracle当部署方法、启用新预留测试、切换生产默认、提交HPC或commit/push。

时间戳：2026-11-10，America/New_York；精确复核时间见 `validation/handoff_20261110/audit_summary.json`。
