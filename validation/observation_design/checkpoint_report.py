#!/usr/bin/env python3
"""Write the completed resource-checkpoint report without ranking missing cases."""
from pathlib import Path
import csv
import json
from datetime import datetime
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(data, fields):
    return ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"] + [
        "| " + " | ".join(str(row.get(field, "")) for field in fields) + " |" for row in data]


def main():
    availability = rows(OUT / "analysis_input/availability_by_design.csv")
    status = rows(OUT / "analysis_input/all_fit_status.csv")
    assert len(status) == 60 and len({row["case_id"] for row in status}) == 60
    missing = [row for row in status if row["status"] != "passed_quadrature_refinement"]
    resources = [row for row in status if row["resource_pass"].lower() != "true"]
    assert missing, "This checkpoint generator must not replace a complete statistical report"
    baseline_tests = rows(OUT / "validation/validation_tests.csv")
    direct_tests = rows(OUT / "direct_validation_tests.csv")
    for row in baseline_tests + direct_tests:
        assert float(row["failed"]) == 0 and row["error"].lower() == "false" and float(row["warning"]) == 0
    assertions = int(sum(float(row["passed"]) for row in baseline_tests + direct_tests))
    latest = json.loads((OUT / "direct_smoke_resources.json").read_text())
    assert not latest["passed_6GiB_300sec"]
    experiments = []
    for version in (3, 4, 5, 6, 7):
        path = OUT / (f"v{version}_direct_smoke_resources.json" if version < 7 else "direct_smoke_resources.json")
        if path.exists():
            record = json.loads(path.read_text())
            experiments.append({"版本": f"direct v{version}", "秒": round(record["child_seconds"], 2),
                "峰值GiB": round(record["maximum_rss_bytes"] / 1024**3, 3),
                "资源与数值关卡通过": record["passed_6GiB_300sec"]})
    with (OUT / "checkpoint_resource_summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(experiments[0]))
        writer.writeheader()
        writer.writerows(experiments)
    metadata = {"timestamp_utc": datetime.now(ZoneInfo("UTC")).isoformat(),
        "original_attempts": 60, "numerically_available": 60-len(missing), "unavailable": len(missing),
        "resource_failures": len(resources), "validation_assertions_passed": assertions,
        "latest_smoke": latest, "statistical_ranking": "WITHHELD_NUMERICAL_GATE_FAILED",
        "larger_domains_run": False, "reserved_calibration_or_test_opened": False,
        "production_default_changed": False, "commit_or_push_performed": False}
    (OUT / "checkpoint_status.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    unavailable = [{"案例": row["case_id"], "原始状态": row["raw_status"],
        "最终数值状态": row["status"], "原始耗时秒": round(float(row["elapsed_seconds"]), 2)} for row in missing]
    diagnostic = rows(OUT / "analysis_input/complete_depth_control_diagnostic.csv")
    diagnostic_table = [{"原始距离": int(float(row["distance"])),
        "RMSE变化%": f"{float(row['rmse_percent_change']):+.2f}",
        "Shape变化%": f"{float(row['shape_rmse_percent_change']):+.2f}",
        "IS50变化%": f"{float(row['score50_percent_change']):+.2f}",
        "IS80变化%": f"{float(row['score80_percent_change']):+.2f}",
        "IS95变化%": f"{float(row['score95_percent_change']):+.2f}"}
        for row in diagnostic if float(row["distance"]) >= 0]
    columns = ["design", "attempted", "numerically_available", "raw_numerical_failures", "recovered", "budget_failures"]
    text = ["# ALFAK2 固定预算观测设计：最终开发检查点", "",
        f"60 个原定案例均已尝试，{60-len(missing)} 个数值可用，{len(missing)} 个不可用。完整五设计精度排名被预设规则阻止，尚不能判断哪种布局能同时改善所有距离的 fitness 与不确定性。", "",
        "这不是某种设计在统计精度上被证伪；当前的决定依据是数值／资源关卡未通过。不能删除失败案例后用剩余均值宣布胜出。", "",
        "## 已完成的比较设置", "",
        "同一 27 节点图、beta=.005、T=1.5、初始总细胞数 1,200；6 个独立 landscape，smooth/rough/interaction 各 2 个，每个有 2 个生物重复。数值重复计算不会增加独立样本量。", "",
        "| 设计 | 分型细胞 | 培养组 | 采样次数 | 分配 |", "| --- | --- | --- | --- | --- |",
        "| control300 | 300 | 1 | 2 | t0/T 各150 |", "| depth600 | 600 | 1 | 2 | t0/T 各300 |",
        "| time600 | 600 | 1 | 3 | t0、T/2、T 各200 |",
        "| composition_anchor600 | 600 | 2 | 4 | 每组t0/T各150，共享核型锚点 |",
        "| composition600 | 600 | 2 | 4 | 每组t0/T各150，无共享初始锚点的诊断对照 |", "",
        "各布局估计自己的未知初始组成，共享 fitness；所有状态／培养使用相同的、与观测支持无关的初始-logit SD=2 prior。未提供真实初始组成、家族标签或真实生长量。主比较固定使用 control300 支持定义的原始距离及同一名义 A 参照。培养数和操作数另计，600 个分型细胞不代表实验总成本相等。", "",
        "## 完整数值可用性", "", *table(availability, columns), "",
        *table(unavailable, ["案例", "原始状态", "最终数值状态", "原始耗时秒"]), "",
        "状态、预算与恢复尝试全部保留在 [all_fit_status.csv](results/analysis_input/all_fit_status.csv)；不可用案例在 [all_case_metrics.csv](results/analysis_input/all_case_metrics.csv) 中保留零预测／区间可用率与 NA 风险。原始报告见 [REVIEW_v4.md](REVIEW_v4.md)。", "",
        "## 完整的预先指定对照：300→600 个分型细胞", "",
        "control300 与 depth600 的24个拟合全部可用，因此可单独描述这个预先指定的预算翻倍对照：全部12对生物单元、6个独立landscape均纳入，没有删除失败样本。生物重复先在landscape内平均，再在家族内平均，最后三个家族等权。下表是平均风险变化，负值改善，正值退步。", "",
        *table(diagnostic_table, ["原始距离", "RMSE变化%", "Shape变化%", "IS50变化%", "IS80变化%", "IS95变化%"]), "",
        "本开发集的观测层均值改善，1-hop与2-hop的点误差均值反而上升；加深分型未呈现所有距离同时改善。这里只是6个landscape的描述性结果，不是显著性或正式非劣效检验，也不能替代三个600细胞布局的完整比较。未据此推广默认方法或开启正式测试。", "",
        "数据：[完整配对对照](results/analysis_input/complete_depth_control_diagnostic.csv)、[全部家族结果](results/analysis_input/complete_depth_control_family_metrics.csv)、[全部landscape结果](results/analysis_input/complete_depth_control_landscape_metrics.csv)。这一完整对照不对缺失的布局作排名。", "",
        "## 已定位和测试的计算问题", "",
        "原始 R 超时可能被 TMB 捕获并转成 NaN，已改为外层 R 子进程硬时限；中断前用时累计。固定尺度下没有外层自由参数，联合众数＋完整 Hessian 可计算同一个条件 Laplace 近似，避免不必要的嵌套求导。原始失败、较早原型的失败和所有源码快照均保留。", "",
        "高增长核对先分段归一化并累计对数增长，避免直接构造极大总人口量；Taylor 前向运算也改为每步归一化，高尺度使用较高阶数。相邻积分节点复用已求得众数作为初值；全部实际积分节点及混合权重保留。似然、先验、数据、参照、距离分组和统计阈值均未修改。", "",
        f"当前实现通过 {assertions} 项核对：134 项模型／导数／数据检查与 {assertions-134} 项等价性、稳定性及积分初值复用检查，均无失败、错误或警告。这些检查不认证 27 节点的真实 fitness 准确率、完整后验精确性或频率覆盖率。", "",
        *table(experiments, ["版本", "秒", "峰值GiB", "资源与数值关卡通过"]), "",
        f"最新最差案例尝试峰值为 {latest['maximum_rss_bytes']/1024**3:.2f} GiB，但在 {latest['child_seconds']:.2f} 秒被硬时限终止；停止在 65 阶的极端 baseline 节点约 7.333。该案例原先曾达到约 9.06 GiB。内存下降不等于完整数值关卡通过，后续全批重跑没有启动。", "",
        "详细记录见 [checkpoint_status.json](results/checkpoint_status.json)、[direct_smoke.log](results/direct_smoke.log)、[direct_validation_tests.csv](results/direct_validation_tests.csv) 与 [DEADLINE_AMENDMENT.md](DEADLINE_AMENDMENT.md)。", "",
        "## 统计决定与下一步", "",
        "固定的筛选条件仍为：所有原始距离点风险不退步、centered shape 至少改善5%、edge不退步、IS50/80/95不退步，以及各家族／距离风险比不超过1.05。由于存在不可用案例，本轮未计算完整风险排名、未认证这些条件，也未认证正式测试样本量。Fisher 信息只是规划诊断，不能替代实际精度证据。", "",
        "下一步优先按证据与矩贡献自适应安排一维 baseline 积分，保留完整原先验与显式尾部误差核对。先恢复六个不可用案例，再统一比较全60个；具体建议见 [NEXT_NUMERICAL_PLAN.md](NEXT_NUMERICAL_PLAN.md)，尚未执行。", "",
        "未扩展到81/64节点，未开启预留校准／正式测试，未替换生产默认，未提交HPC任务，未commit或push。原有 two_layer 默认及 joint_count 的 opt-in 状态保持；用户原有文件修改和旧脚本删除保持。", "",
        "数据生成采用精确有限细胞分裂过程，拟合仍以确定性均值轨迹的 multinomial 似然近似，未显式建模全部过程噪声和跨时间相关性。数值正确、模型适配和新增观测信息仍须分别验证。", "",
        "协议与时间线：[PROTOCOL.md](PROTOCOL.md)、[ANALYSIS_RULES.md](ANALYSIS_RULES.md)、[POSTPROCESSING_AMENDMENT.md](POSTPROCESSING_AMENDMENT.md)、[DIRECT_JOINT_AMENDMENT.md](DIRECT_JOINT_AMENDMENT.md)。内部分析规则在拟合启动后、精度审阅前固定，非外部预注册。", "",
        "生成时间：" + datetime.now(ZoneInfo("America/New_York")).strftime("%Y-%m-%d %H:%M:%S %Z")]
    (ROOT / "REVIEW.md").write_text("\n".join(text) + "\n")
    print(json.dumps({key: metadata[key] for key in ["original_attempts", "numerically_available", "unavailable", "resource_failures", "validation_assertions_passed"]}))


if __name__ == "__main__":
    main()
