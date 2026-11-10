# 固定细胞预算的观测设计：数值与资源关卡未通过

60 个预定案例均已尝试；54 个数值可用，6 个仍不可用。完整五设计精度排名不成立，不能声称所有距离的 fitness 和不确定性均获改善。

没有丢弃失败 seed、家族或距离，没有用有限值均值替代缺失结果。未扩展到 81/64 节点，未开启预留校准/测试。

| design | attempted | numerically_available | raw_numerical_failures | recovered | budget_failures |
| --- | --- | --- | --- | --- | --- |
| composition_anchor600 | 12 | 9 | 3 | 0 | 3 |
| composition600 | 12 | 10 | 2 | 0 | 3 |
| control300 | 12 | 12 | 0 | 0 | 0 |
| depth600 | 12 | 12 | 0 | 0 | 0 |
| time600 | 12 | 11 | 1 | 0 | 1 |

## 不可用案例

| case_id | status | elapsed_seconds | recovery |
| --- | --- | --- | --- |
| design_interaction_L01_R2_time600 | integration_error | 237.57 | none |
| design_rough_L01_R1_composition_anchor600 | integration_error | 319.58 | none |
| design_rough_L01_R1_composition600 | integration_error | 345.51 | none |
| design_interaction_L01_R1_composition_anchor600 | external_time_limit | 300.04 | none |
| design_interaction_L02_R1_composition_anchor600 | external_time_limit | 300.11 | none |
| design_interaction_L02_R1_composition600 | external_time_limit | 300.05 | none |

原始状态、计算预算和恢复标签见 [all_fit_status.csv](results/analysis_input/all_fit_status.csv)；不可用案例的逐距离 prediction/interval availability 为零、风险为 NA，保留在 [all_case_metrics.csv](results/analysis_input/all_case_metrics.csv)。

R 内部超时可能被 TMB 捕获，原批次出现超过 300 秒的案例。外层独立 R 子进程的硬时限已修复控制问题；中断前耗时累计，最多两例并行。详见 [DEADLINE_AMENDMENT.md](DEADLINE_AMENDMENT.md)。已有完整缓存可离线重构，但不能补造缺失拟合或使超预算变成通过。

已完成 134 项模型/数据核对。共同先验、参照、目标节点和原始距离分组保持固定；具体设计与修订见 [PROTOCOL.md](PROTOCOL.md)、[DESIGN_AMENDMENT.md](DESIGN_AMENDMENT.md)。

证据值的独立联合众数/Hessian 核对、原始失败与不进行新拟合的恢复规则见 [EVIDENCE_RECONSTRUCTION.md](EVIDENCE_RECONSTRUCTION.md) 和 [POSTPROCESSING_AMENDMENT.md](POSTPROCESSING_AMENDMENT.md)。分析规则在拟合启动后、精度审阅前固定，属于内部开发规则，非外部预注册。

下一阶段应先解决原型在计算预算内的积分效率和证据值可靠性，再完成同一冻结比较。当前数据不足以推荐一个替换默认方法的设计。生产默认、用户已有修改保持不变；未 commit、push 或提交 HPC 任务。

生成时间： 2026-11-04 17:06:44 EST
