#!/usr/bin/env python3
"""Summarize the frozen sixty-case comparison without changing its decisions."""
import csv
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    out = results / "adaptive_complete_v1"
    audit = read(out / "posterior_integrity_audit.csv")
    assert len(audit) == 60
    assert len({r["case_id"] for r in audit}) == 60
    compat = read(results / "adaptive_compatibility_v1/compatibility.csv")
    assert len(compat) == 3 and all(r["passed"] == "TRUE" for r in compat)
    references = [read(p)[0] for p in sorted(
        (results / "adaptive_stage_v1/reference").glob("*/reference_comparison.csv"))]
    assert len(references) == 8 and all(r["passed"] == "TRUE" for r in references)
    write(out / "small_reference_comparisons.csv", references)
    write(out / "legacy_compatibility.csv", compat)
    old = {r["case_id"]: r for r in read(results / "analysis_input/all_fit_status.csv")}
    recovery = []
    for folder in ["adaptive_primary_smoke_v1", "adaptive_repairs_v1"]:
        for part in range(2):
            for r in read(results / folder / f"primary_part_{part}/status.csv"):
                assert r["validation_pass"] == "TRUE" and r["resource_pass"] == "True"
                recovery.append(dict(case_id=r["case_id"], original_status=old[r["case_id"]]["status"],
                    original_seconds=old[r["case_id"]]["elapsed_seconds"],
                    fresh_fit_seconds=float(r["elapsed_seconds"]),
                    fresh_wall_seconds=float(r["child_wall_seconds"]),
                    fresh_peak_gib=float(r["maximum_rss_bytes"]) / 1024**3,
                    evaluated_nodes=int(r["nodes"])))
    assert len(recovery) == 6
    write(out / "recovery_summary.csv", recovery)
    overall = read(out / "baseline_overall_metrics.csv")
    keys = ["rmse", "shape_rmse", "edge_rmse", "score50", "score80", "score95"]
    changes = []
    for ref in ["control300", "depth600"]:
        for design in sorted({r["design"] for r in overall} - {ref}):
            for d in range(3):
                x = next(r for r in overall if r["design"] == design and int(r["distance"]) == d)
                b = next(r for r in overall if r["design"] == ref and int(r["distance"]) == d)
                changes.append(dict(design=design, reference=ref, distance=d,
                    **{k + "_change_percent": 100 * (float(x[k]) / float(b[k]) - 1) for k in keys}))
    write(out / "overall_risk_changes.csv", changes)
    triage = read(out / "development_triage.csv")
    stage = (out / "stage_status.txt").read_text().strip()
    status = dict(stage=stage, required_cases=60, available_cases=60,
        independent_landscapes=6, biological_units=12, repaired_cases=6,
        legacy_posteriors_retained=54, representative_compatibility_cases=3,
        analytic_assertions_passed=65, previous_implementation_assertions_passed=211,
        small_graph_references_passed=8,
        max_recovery_wall_seconds=max(r["fresh_wall_seconds"] for r in recovery),
        max_recovery_peak_gib=max(r["fresh_peak_gib"] for r in recovery),
        max_compatibility_mean_change=max(float(r["mean_change"]) for r in compat),
        max_compatibility_covariance_change=max(float(r["covariance_change"]) for r in compat),
        max_small_reference_mean_error_sd=max(float(r["new_mean_error_reference_sd"]) for r in references),
        legacy_resource_failures_retained=len(read(out / "resource_failures.csv")),
        new_calibration_or_test_opened=False, production_default_changed=False,
        hpc_jobs_submitted=False, commit_or_push=False)
    status["generated_at"] = datetime.now(ZoneInfo("America/New_York")).isoformat()
    (out / "adaptive_checkpoint_status.json").write_text(json.dumps(status, indent=2) + "\n")
    headers = ["设计（相对 depth600）", "原始距离", "RMSE变化%", "Shape变化%", "Edge变化%", "IS50变化%", "IS80变化%", "IS95变化%"]
    rows = [[r["design"], r["distance"]] + [f'{r[k + "_change_percent"]:+.2f}' for k in keys]
            for r in changes if r["reference"] == "depth600"]
    rec_rows = [[r["case_id"], f'{r["fresh_wall_seconds"]:.2f}', f'{r["fresh_peak_gib"]:.3f}', r["evaluated_nodes"]]
                for r in recovery]
    flags = [[r["design"]] + [r[k] for k in ["point_nonworse", "shape_5percent_all_distances", "edge_nonworse", "interval_scores_nonworse", "family_safety"]]
             for r in triage]
    text = f'''# ALFAK2 固定预算观测设计：自适应积分完成检查点

原定60个开发案例现已全部数值可用，完整比较已完成；两种预先指定的600细胞候选均未达到“所有距离的fitness与不确定性同时改善”的开发条件。共同锚点布局改善了远场均值，但牺牲了观测层，不能据此升级默认方法。

## 本轮完成的数值改进

对原Gaussian baseline先验进行自适应Gauss–Kronrod积分，保留完整支持、原似然、原fitness和初始组成先验、固定参照、原始距离分组与统计阈值。直接prior-CDF坐标在解析压力测试中效率较差，因此实际使用 `a=sd*logit(u)` 并显式保留Gaussian密度和Jacobian；这是变量替换。每次通过前都强制细分两端面板，核对归一化、均值、总协方差和固定参照协方差。误差值是数值估计，不是严格的全局数学误差界。

新增65项解析断言全部通过；此前211项实现核对结果保留。8个4/8节点后验参照全部通过兼容性检查。相对旧收敛积分的最大均值变化为 {max(float(r['mean_change_reference_sd']) for r in references):.3g} 个MCMC参照SD，最大相对SD变化为 {max(float(r['relative_sd_change']) for r in references):.3g}；条件Laplace近似对MCMC的最大均值误差仍约 {status['max_small_reference_mean_error_sd']:.3f} SD，不能称为精确后验。

六个原不可用案例均重新独立拟合成功，最多同时运行两个实际拟合，均通过300秒／6 GiB门槛：

{table(['恢复案例', '新运行墙钟秒', '新峰值GiB', '计算节点'], rec_rows)}

另对smooth、rough、interaction各一个已完成27节点案例进行不使用生成真值的数值兼容性核对。三例全部通过，最大缩放均值变化 {status['max_compatibility_mean_change']:.3g}，最大协方差相对变化 {status['max_compatibility_covariance_change']:.3g}，50/80/95%区间分布检查也通过。因此沿用54个已收敛旧后验并纳入6个新恢复结果；这不是全60例统一重算。三例抽查支持兼容性，不能证明其余所有案例完全一致。

60个后验均通过独立的混合均值／总协方差、固定参照、MSE分解、前向误差、收敛状态及已使用源码核对；新积分另通过完整区间分割和Gaussian先验／Jacobian权重重建。原始失败、旧原型、运行成本和源码快照全部保留。诊断布局composition600仍保留一例旧运行305秒的预算失败；旧批次部分案例无逐例峰值记录，不能声称所有60例都重新满足资源门槛。

## 完整60例结果

同一27节点、beta=.005、T=1.5、初始细胞总数1200。6个独立landscape（每个家族2个），每个2个生物重复，共12个配对单元、5个布局。除control300外均分型600个细胞。培养数和采样次数另计。各布局估计各自未知初始组成、共享fitness，没有输入真实组成、生长量或家族标签。

主距离固定为control300的实际观测支持所定义的0/1/2-hop，参照固定为名义A。先在landscape内平均生物重复，再在家族内平均，最后三个家族等权。Shape始终只做一次全图居中；Edge按两个端点中较远的原始距离分组。负百分比表示风险改善，IS为原始未校准的proper interval score。

{table(headers, rows)}

共同锚点相对depth600的1-hop／2-hop RMSE下降7.29%／6.59%，观测层RMSE上升2.02%、Shape上升8.37%、IS95上升13.45%。三时间点布局观测层RMSE上升4.71%，远场IS95也变差。无共享锚点的composition600只作预先披露的诊断对照，不能因某列较好临时替换主要候选；它同样没有全距离一致改善。

共同锚点2-hop的名义95%区间经验覆盖率为86.94%，depth600为89.03%。这些是同一小型开发集上的描述性覆盖率；远场IS改善不代表区间已经校准。

## 冻结的开发关卡

{table(['候选', '各距离点风险不退步', '各距离Shape改善≥5%', 'Edge不退步', 'IS50/80/95不退步', '家族安全界限'], flags)}

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

生成时间：{datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')}
'''
    archive = root / "REVIEW_pre_adaptive.md"
    if not archive.exists():
        archive.write_bytes((root / "REVIEW.md").read_bytes())
    (root / "REVIEW.md").write_text(text)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
