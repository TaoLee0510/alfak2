# ALFAK2 固定预算观测设计：最终开发检查点

60 个原定案例均已尝试，54 个数值可用，6 个不可用。完整五设计精度排名被预设规则阻止，尚不能判断哪种布局能同时改善所有距离的 fitness 与不确定性。

这不是某种设计在统计精度上被证伪；当前的决定依据是数值／资源关卡未通过。不能删除失败案例后用剩余均值宣布胜出。

## 已完成的比较设置

同一 27 节点图、beta=.005、T=1.5、初始总细胞数 1,200；6 个独立 landscape，smooth/rough/interaction 各 2 个，每个有 2 个生物重复。数值重复计算不会增加独立样本量。

| 设计 | 分型细胞 | 培养组 | 采样次数 | 分配 |
| --- | --- | --- | --- | --- |
| control300 | 300 | 1 | 2 | t0/T 各150 |
| depth600 | 600 | 1 | 2 | t0/T 各300 |
| time600 | 600 | 1 | 3 | t0、T/2、T 各200 |
| composition_anchor600 | 600 | 2 | 4 | 每组t0/T各150，共享核型锚点 |
| composition600 | 600 | 2 | 4 | 每组t0/T各150，无共享初始锚点的诊断对照 |

各布局估计自己的未知初始组成，共享 fitness；所有状态／培养使用相同的、与观测支持无关的初始-logit SD=2 prior。未提供真实初始组成、家族标签或真实生长量。主比较固定使用 control300 支持定义的原始距离及同一名义 A 参照。培养数和操作数另计，600 个分型细胞不代表实验总成本相等。

## 完整数值可用性

| design | attempted | numerically_available | raw_numerical_failures | recovered | budget_failures |
| --- | --- | --- | --- | --- | --- |
| composition_anchor600 | 12 | 9 | 3 | 0 | 3 |
| composition600 | 12 | 10 | 2 | 0 | 3 |
| control300 | 12 | 12 | 0 | 0 | 0 |
| depth600 | 12 | 12 | 0 | 0 | 0 |
| time600 | 12 | 11 | 1 | 0 | 1 |

| 案例 | 原始状态 | 最终数值状态 | 原始耗时秒 |
| --- | --- | --- | --- |
| design_interaction_L01_R2_time600 | integration_error | integration_error | 237.57 |
| design_rough_L01_R1_composition_anchor600 | integration_error | integration_error | 319.58 |
| design_rough_L01_R1_composition600 | integration_error | integration_error | 345.51 |
| design_interaction_L01_R1_composition_anchor600 | external_time_limit | external_time_limit | 300.04 |
| design_interaction_L02_R1_composition_anchor600 | external_time_limit | external_time_limit | 300.11 |
| design_interaction_L02_R1_composition600 | external_time_limit | external_time_limit | 300.05 |

状态、预算与恢复尝试全部保留在 [all_fit_status.csv](results/analysis_input/all_fit_status.csv)；不可用案例在 [all_case_metrics.csv](results/analysis_input/all_case_metrics.csv) 中保留零预测／区间可用率与 NA 风险。原始报告见 [REVIEW_v4.md](REVIEW_v4.md)。

## 完整的预先指定对照：300→600 个分型细胞

control300 与 depth600 的24个拟合全部可用，因此可单独描述这个预先指定的预算翻倍对照：全部12对生物单元、6个独立landscape均纳入，没有删除失败样本。生物重复先在landscape内平均，再在家族内平均，最后三个家族等权。下表是平均风险变化，负值改善，正值退步。

| 原始距离 | RMSE变化% | Shape变化% | IS50变化% | IS80变化% | IS95变化% |
| --- | --- | --- | --- | --- | --- |
| 0 | -1.69 | -5.24 | -2.89 | -4.43 | -5.63 |
| 1 | +2.28 | +1.38 | +2.54 | +2.39 | -1.92 |
| 2 | +3.74 | +0.42 | +4.03 | +1.38 | +0.64 |

本开发集的观测层均值改善，1-hop与2-hop的点误差均值反而上升；加深分型未呈现所有距离同时改善。这里只是6个landscape的描述性结果，不是显著性或正式非劣效检验，也不能替代三个600细胞布局的完整比较。未据此推广默认方法或开启正式测试。

数据：[完整配对对照](results/analysis_input/complete_depth_control_diagnostic.csv)、[全部家族结果](results/analysis_input/complete_depth_control_family_metrics.csv)、[全部landscape结果](results/analysis_input/complete_depth_control_landscape_metrics.csv)。这一完整对照不对缺失的布局作排名。

## 已定位和测试的计算问题

原始 R 超时可能被 TMB 捕获并转成 NaN，已改为外层 R 子进程硬时限；中断前用时累计。固定尺度下没有外层自由参数，联合众数＋完整 Hessian 可计算同一个条件 Laplace 近似，避免不必要的嵌套求导。原始失败、较早原型的失败和所有源码快照均保留。

高增长核对先分段归一化并累计对数增长，避免直接构造极大总人口量；Taylor 前向运算也改为每步归一化，高尺度使用较高阶数。相邻积分节点复用已求得众数作为初值；全部实际积分节点及混合权重保留。似然、先验、数据、参照、距离分组和统计阈值均未修改。

当前实现通过 211 项核对：134 项模型／导数／数据检查与 77 项等价性、稳定性及积分初值复用检查，均无失败、错误或警告。这些检查不认证 27 节点的真实 fitness 准确率、完整后验精确性或频率覆盖率。

| 版本 | 秒 | 峰值GiB | 资源与数值关卡通过 |
| --- | --- | --- | --- |
| direct v3 | 79.22 | 0.654 | False |
| direct v4 | 76.63 | 0.482 | False |
| direct v5 | 220.47 | 0.537 | False |
| direct v6 | 251.29 | 0.764 | False |
| direct v7 | 300.05 | 1.148 | False |

最新最差案例尝试峰值为 1.15 GiB，但在 300.05 秒被硬时限终止；停止在 65 阶的极端 baseline 节点约 7.333。该案例原先曾达到约 9.06 GiB。内存下降不等于完整数值关卡通过，后续全批重跑没有启动。

详细记录见 [checkpoint_status.json](results/checkpoint_status.json)、[direct_smoke.log](results/direct_smoke.log)、[direct_validation_tests.csv](results/direct_validation_tests.csv) 与 [DEADLINE_AMENDMENT.md](DEADLINE_AMENDMENT.md)。

## 统计决定与下一步

固定的筛选条件仍为：所有原始距离点风险不退步、centered shape 至少改善5%、edge不退步、IS50/80/95不退步，以及各家族／距离风险比不超过1.05。由于存在不可用案例，本轮未计算完整风险排名、未认证这些条件，也未认证正式测试样本量。Fisher 信息只是规划诊断，不能替代实际精度证据。

下一步优先按证据与矩贡献自适应安排一维 baseline 积分，保留完整原先验与显式尾部误差核对。先恢复六个不可用案例，再统一比较全60个；具体建议见 [NEXT_NUMERICAL_PLAN.md](NEXT_NUMERICAL_PLAN.md)，尚未执行。

未扩展到81/64节点，未开启预留校准／正式测试，未替换生产默认，未提交HPC任务，未commit或push。原有 two_layer 默认及 joint_count 的 opt-in 状态保持；用户原有文件修改和旧脚本删除保持。

数据生成采用精确有限细胞分裂过程，拟合仍以确定性均值轨迹的 multinomial 似然近似，未显式建模全部过程噪声和跨时间相关性。数值正确、模型适配和新增观测信息仍须分别验证。

协议与时间线：[PROTOCOL.md](PROTOCOL.md)、[ANALYSIS_RULES.md](ANALYSIS_RULES.md)、[POSTPROCESSING_AMENDMENT.md](POSTPROCESSING_AMENDMENT.md)、[DIRECT_JOINT_AMENDMENT.md](DIRECT_JOINT_AMENDMENT.md)。内部分析规则在拟合启动后、精度审阅前固定，非外部预注册。

生成时间：2026-11-04 17:47:52 EST
