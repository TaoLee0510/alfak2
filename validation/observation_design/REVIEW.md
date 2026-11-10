# ALFAK2 固定预算观测设计：自适应积分完成检查点

原定60个开发案例现已全部数值可用，完整比较已完成；两种预先指定的600细胞候选均未达到“所有距离的fitness与不确定性同时改善”的开发条件。共同锚点布局改善了远场均值，但牺牲了观测层，不能据此升级默认方法。

## 本轮完成的数值改进

对原Gaussian baseline先验进行自适应Gauss–Kronrod积分，保留完整支持、原似然、原fitness和初始组成先验、固定参照、原始距离分组与统计阈值。直接prior-CDF坐标在解析压力测试中效率较差，因此实际使用 `a=sd*logit(u)` 并显式保留Gaussian密度和Jacobian；这是变量替换。每次通过前都强制细分两端面板，核对归一化、均值、总协方差和固定参照协方差。误差值是数值估计，不是严格的全局数学误差界。

新增65项解析断言全部通过；此前211项实现核对结果保留。8个4/8节点后验参照全部通过兼容性检查。相对旧收敛积分的最大均值变化为 9.75e-05 个MCMC参照SD，最大相对SD变化为 4.71e-05；条件Laplace近似对MCMC的最大均值误差仍约 0.122 SD，不能称为精确后验。

六个原不可用案例均重新独立拟合成功，最多同时运行两个实际拟合，均通过300秒／6 GiB门槛：

| 恢复案例 | 新运行墙钟秒 | 新峰值GiB | 计算节点 |
| --- | --- | --- | --- |
| design_interaction_L01_R1_composition_anchor600 | 257.84 | 0.545 | 195 |
| design_interaction_L01_R2_time600 | 143.79 | 0.535 | 195 |
| design_rough_L01_R1_composition_anchor600 | 265.06 | 0.508 | 195 |
| design_interaction_L02_R1_composition_anchor600 | 234.54 | 0.550 | 165 |
| design_rough_L01_R1_composition600 | 266.27 | 0.520 | 195 |
| design_interaction_L02_R1_composition600 | 238.63 | 0.554 | 165 |

另对smooth、rough、interaction各一个已完成27节点案例进行不使用生成真值的数值兼容性核对。三例全部通过，最大缩放均值变化 1.7e-06，最大协方差相对变化 1.06e-05，50/80/95%区间分布检查也通过。因此沿用54个已收敛旧后验并纳入6个新恢复结果；这不是全60例统一重算。三例抽查支持兼容性，不能证明其余所有案例完全一致。

60个后验均通过独立的混合均值／总协方差、固定参照、MSE分解、前向误差、收敛状态及已使用源码核对；新积分另通过完整区间分割和Gaussian先验／Jacobian权重重建。原始失败、旧原型、运行成本和源码快照全部保留。诊断布局composition600仍保留一例旧运行305秒的预算失败；旧批次部分案例无逐例峰值记录，不能声称所有60例都重新满足资源门槛。

## 完整60例结果

同一27节点、beta=.005、T=1.5、初始细胞总数1200。6个独立landscape（每个家族2个），每个2个生物重复，共12个配对单元、5个布局。除control300外均分型600个细胞。培养数和采样次数另计。各布局估计各自未知初始组成、共享fitness，没有输入真实组成、生长量或家族标签。

主距离固定为control300的实际观测支持所定义的0/1/2-hop，参照固定为名义A。先在landscape内平均生物重复，再在家族内平均，最后三个家族等权。Shape始终只做一次全图居中；Edge按两个端点中较远的原始距离分组。负百分比表示风险改善，IS为原始未校准的proper interval score。

| 设计（相对 depth600） | 原始距离 | RMSE变化% | Shape变化% | Edge变化% | IS50变化% | IS80变化% | IS95变化% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| composition600 | 0 | -4.59 | +4.15 | -3.10 | +1.21 | +8.50 | +3.90 |
| composition600 | 1 | -6.01 | -3.75 | -0.03 | -4.36 | -7.58 | -9.16 |
| composition600 | 2 | +0.27 | +4.55 | +0.33 | -3.57 | -0.47 | -11.13 |
| composition_anchor600 | 0 | +2.02 | +8.37 | +5.54 | +5.63 | +15.57 | +13.45 |
| composition_anchor600 | 1 | -7.29 | -4.76 | +3.39 | -6.02 | -6.26 | -3.90 |
| composition_anchor600 | 2 | -6.59 | -2.39 | -3.32 | -7.98 | -1.59 | -10.74 |
| control300 | 0 | +1.72 | +5.53 | +1.40 | +2.98 | +4.63 | +5.96 |
| control300 | 1 | -2.23 | -1.36 | +0.52 | -2.48 | -2.33 | +1.96 |
| control300 | 2 | -3.60 | -0.41 | -0.88 | -3.88 | -1.36 | -0.63 |
| time600 | 0 | +4.71 | +4.41 | +3.67 | +3.58 | +5.34 | +0.11 |
| time600 | 1 | -0.30 | -1.98 | +0.19 | -0.28 | +2.50 | +4.53 |
| time600 | 2 | -1.07 | -0.37 | -0.80 | -0.16 | +2.22 | +9.17 |

共同锚点相对depth600的1-hop／2-hop RMSE下降7.29%／6.59%，观测层RMSE上升2.02%、Shape上升8.37%、IS95上升13.45%。三时间点布局观测层RMSE上升4.71%，远场IS95也变差。无共享锚点的composition600只作预先披露的诊断对照，不能因某列较好临时替换主要候选；它同样没有全距离一致改善。

共同锚点2-hop的名义95%区间经验覆盖率为86.94%，depth600为89.03%。这些是同一小型开发集上的描述性覆盖率；远场IS改善不代表区间已经校准。

## 冻结的开发关卡

| 候选 | 各距离点风险不退步 | 各距离Shape改善≥5% | Edge不退步 | IS50/80/95不退步 | 家族安全界限 |
| --- | --- | --- | --- | --- | --- |
| time600 | FALSE | FALSE | FALSE | FALSE | FALSE |
| composition_anchor600 | FALSE | FALSE | FALSE | FALSE | FALSE |

time600有7个、composition_anchor600有10个家族／距离／风险单元超过1.05的描述性安全界限。共同锚点的smooth家族1-hop IS95风险比最高约1.333，不能用总体远场均值改善掩盖。两种主要候选的资源标记通过，但统计条件均未通过。这里只包含6个独立landscape；这些是开发集描述性比较，不是显著性、正式非劣效或频率覆盖率认证。

## 下一阶段建议

本轮已经从“数值不完整”推进到可判断的近场／远场权衡。建议下一轮固定600个分型细胞、相同图和观察时长，把保持原观测层信息设为硬约束，优化初始核型组成以及t0／终点或两培养组之间的分型分配。先在事先指定的参数点用包含未知初始组成的profiled Fisher规划可行区，再冻结候选，使用完整开发样本和相同数值门槛检验各距离RMSE、Shape、Edge及IS50/80/95。Fisher预测不能代替实际精度证据。该设计优化属于下一轮研究，尚未执行。

仍需检查有限分裂过程噪声对区间的影响：数据来自精确Gillespie模拟，拟合仍采用确定性均值轨迹加multinomial似然，未显式表示全部过程噪声和时间相关性。完成积分解决了计算瓶颈，没有证明当前模型足以实现所有距离的精度目标。

本轮未扩展81/64节点，未开启预留校准／正式测试，未替换生产默认，未提交HPC任务，未commit或push。原two_layer默认、joint_count opt-in、用户原有修改与旧脚本删除均保留。

## 可核查文件

- [完整风险变化表](results/adaptive_complete_v1/overall_risk_changes.csv)、[逐家族风险比](results/adaptive_complete_v1/family_risk_ratios.csv)、[全部landscape配对差异](results/adaptive_complete_v1/paired_landscape_differences.csv)
- [60例状态及原失败记录](results/adaptive_complete_v1/all_fit_status.csv)、[恢复耗时与内存](results/adaptive_complete_v1/recovery_summary.csv)、[开发关卡](results/adaptive_complete_v1/development_triage.csv)
- [后验独立核对](results/adaptive_complete_v1/posterior_integrity_audit.csv)、[实际使用源码核对](results/adaptive_complete_v1/consumed_source_integrity.csv)、[8例小图参照](results/adaptive_complete_v1/small_reference_comparisons.csv)、[3例新旧兼容性](results/adaptive_complete_v1/legacy_compatibility.csv)
- [自适应积分协议](adaptive/PROTOCOL.md)、[冻结的统计规则](ANALYSIS_RULES.md)、[54/60旧检查点](REVIEW_pre_adaptive.md)、[原始报告](REVIEW_v4.md)

生成时间：2026-11-04 18:28:29 EST
