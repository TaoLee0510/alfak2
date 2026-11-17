#!/usr/bin/env python3
"""Present the frozen constrained-design experiment with its full denominator."""
import csv
from datetime import datetime
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo


def read(p):
    with p.open() as stream:
        return list(csv.DictReader(stream))


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    root = Path(os.environ["ALFAK2_REPO"]) / "validation/observation_design"
    out = Path(os.environ["ALFAK2_CONSTRAINED_OUTPUT"])
    status = read(out / "all_fit_status.csv")
    selected = read(out / "selected_designs.csv")
    noise = read(out / "noise/variance_decomposition.csv")
    stage = (out / "stage_status.txt").read_text().strip()
    complete = all(r["status"] == "passed_quadrature_refinement" and r["resource_pass"].lower() == "true" for r in status)
    assert len(status) == 12 * (1 + len(selected))
    costs = [float(r["child_wall_seconds"]) for r in status if r.get("child_wall_seconds") not in ["", "NA"]]
    peaks = [float(r["maximum_rss_bytes"]) / 1024**3 for r in status if r.get("maximum_rss_bytes") not in ["", "NA"]]
    summary = dict(stage=stage, complete_comparison=complete, planned_cases=len(status),
        numeric_available=sum(r["status"] == "passed_quadrature_refinement" for r in status),
        resource_pass=sum(r["resource_pass"].lower() == "true" for r in status),
        independent_landscapes=6, biological_units=12, candidates=len(selected),
        noise_trajectories=300, max_wall_seconds=max(costs, default=None), max_peak_gib=max(peaks, default=None),
        local_only=True, reserved_calibration_test_opened=False, production_default_changed=False)
    now = datetime.now(ZoneInfo("America/New_York"))
    summary["generated_at"] = now.isoformat()
    (out / "checkpoint_status.json").write_text(json.dumps(summary, indent=2) + "\n")
    if complete:
        changes = read(out / "overall_risk_changes.csv")
        triage = read(out / "development_triage.csv")
        keys = ["rmse", "shape_rmse", "edge_rmse", "score50", "score80", "score95"]
        risk_table = table(["候选", "原始距离", "RMSE变化%", "Shape变化%", "Edge变化%", "IS50变化%", "IS80变化%", "IS95变化%"],
            [[r["design"], r["distance"]] + [f'{float(r[k]):+.2f}' for k in keys] for r in changes])
        flags = table(["候选", "点风险不退步", "Shape各层≥5%", "Edge不退步", "IS不退步", "家族安全", "全部通过"],
            [[r["design"]] + [r[k] for k in ["point_nonworse", "shape_5percent_all_distances", "edge_nonworse", "interval_scores_nonworse", "family_safety", "all_distance_promising"]] for r in triage])
        promising = [r["design"] for r in triage if r["all_distance_promising"] == "TRUE"]
        conclusion = ("通过全部开发关卡的候选：" + "、".join(promising) + "。仍需后续独立及更大图验证。") if promising else "两个候选均未通过全部开发关卡，本轮没有满足所有距离目标的布局。"
    else:
        risk_table = "存在数值、资源失败或未启动案例，完整精度比较被关卡阻止；全部案例保留。"
        flags = table(["案例", "状态", "资源通过"], [[r["case_id"], r["status"], r["resource_pass"]] for r in status if r["status"] != "passed_quadrature_refinement" or r["resource_pass"].lower() != "true"])
        conclusion = "数值／资源关卡未完成，不能宣布布局胜出。"
    noise_table = table(["规划点", "原计划距离", "总读数方差／确定性multinomial", "轨迹bootstrap区间"],
        [[r["profile"], r["planned_distance"], f'{float(r["variance_ratio"]):.3f}',
          f'{float(r["bootstrap_lo"]):.3f}–{float(r["bootstrap_hi"]):.3f}'] for r in noise])
    text = f'''# 保持近场信息的观测设计：开发结果

{conclusion}

主比较原定{len(status)}例，现有{summary['numeric_available']}例数值可用，{summary['resource_pass']}例通过资源门槛。使用6个独立landscape、每个2个生物重复；这是12个开发单元，不是12个独立landscape。所有主拟合采用同一个冻结的自适应积分器。最大实际墙钟为{max(costs, default=0):.2f}秒，最大峰值{max(peaks, default=0):.3f} GiB。预检计入主比较并复用，最多同时运行两个实际拟合。

## 筛选与冻结的候选

对5304个布局在9个固定参数点、两种规划平滑尺度1e-4/1e-6下计算局部信息。纳入原fitness先验和未知初始组成先验，通过联合精度的Schur补边缘化初始组成。每个规划点及七个原计划观测核型的相对fitness方差都不能高于depth600（只允许1e-6数值容差）。3026/5250个单培养布局、36/54个双培养布局满足该规划约束；每类只选择预先定义的最差远场方差比最小者。

单培养候选：中心45%、第一壳层0%、第二壳层35%、第三壳层20%，每组内均分；初始1200细胞，t0/终点分型240/360。它确实没有在第一壳层播种；该参数来自冻结的Fisher规则，不能因后续精度表现手动更改。最差观测节点规划方差比0.99850，最差远场组比0.84816。

双培养候选：各培养初始600细胞。A保持原组成、分型400（t0/终点160/240）；B分型200（100/100），初始45%在共同中心锚点、两个极端各20%、六个极端邻居合计15%。最差观测节点规划方差比0.98492，最差远场组比0.94194。

候选选择仅用原先固定的prior/smooth/interaction参数点及0.8/1/1.2速率尺度，未用六个生成真值或其拟合精度挑选。局部方差规划不包含全部非线性偏差和先验失配，不能保证真实误差不退步。已另外保存不含初始组成先验的profile-information敏感性结果。规划代码通过14项导数、Schur补和预算检查。

## 完整精度比较

共同条件：27节点，beta=.005，T=1.5，初始总细胞1200，分型600。保留原12个depth600数据文件并重新统一拟合；两个候选各生成12个开发案例，种子提前冻结。配对依据相同landscape和重复编号，不意味着不同布局使用了同一个培养过程。实际模拟不添加规划用的小概率平滑或伪细胞。

主评价固定使用原control300的实际支持距离及原名义A参照，防止扩大覆盖后重新划分距离掩盖误差。原计划距离和新实际距离仅作辅助结果。生物重复先在landscape内平均，再在家族内平均，最后三个家族等权。Shape误差始终只做一次全图居中；负百分比表示改善。IS50/80/95为原始未校准proper interval score。

{risk_table}

冻结的关卡：每个原始距离RMSE和Edge不退步、Shape至少下降5%、IS50/80/95不退步；每个家族／距离的RMSE、Shape、Edge、IS95比值<=1.05，并且全部案例及区间完整且通过资源门槛。

{flags}

这些是6个开发landscape上的描述性判断，没有显著性、正式非劣效、频率覆盖率或最终测试样本量认证。全部家族结果和六个landscape配对差异均保留，不从成功案例子集计算胜出结果。

## 独立的过程噪声诊断

3个固定规划点各100条精确Gillespie轨迹，统一从原A的固定整数初始人口出发，排除初始人口随机分配这一来源。对终点300细胞频率读数，用全协方差公式分解期望multinomial采样协方差与过程频率协方差。下表比较总频率读数方差与确定性均值轨迹的multinomial预测；它不是fitness后验方差比或区间校准系数。

{noise_table}

远场稀有状态的100条轨迹估计较不稳定，bootstrap区间用于展示Monte Carlo敏感性，不能当成严格尾部误差界或fitness区间覆盖率。诊断没有新增300次后验拟合，没有修改当前似然或区间。候选的实际区间覆盖和宽度另见完整指标表。

## 可复核文件

- [冻结协议](PROTOCOL.md)、[所有Fisher候选](../results/{out.name}/fisher_all_candidates.csv)、[选中参数](../results/{out.name}/selected_designs.csv)、[规划敏感性](../results/{out.name}/selected_fisher_diagnostics.csv)
- [数据与SHA256](../results/{out.name}/dataset_manifest.csv)、[开发种子](../results/{out.name}/development_seed_ledger.csv)、[全部运行状态](../results/{out.name}/all_fit_status.csv)、[全部风险指标](../results/{out.name}/all_case_metrics.csv)
- [逐家族风险比](../results/{out.name}/family_risk_ratios.csv)、[逐landscape差异](../results/{out.name}/paired_landscape_differences.csv)、[开发关卡](../results/{out.name}/development_triage.csv)
- [后验汇总核对](../results/{out.name}/posterior_integrity_audit.csv)、[实际使用源码核对](../results/{out.name}/consumed_source_integrity.csv)、[过程方差分解](../results/{out.name}/noise/variance_decomposition.csv)

历史记录：最初规划执行在完成网格计算后遇到R的换行else语法错误；原constrained_v1输出和源码完整保留，constrained_source_v2只修复该语法并重算相同网格与阈值。它没有改变统计规则或候选范围。所有主比较均来自constrained_v2。

本轮仅本地开发验证。生产默认、预留校准／正式测试及既有用户文件保持；未提交HPC任务，未commit或push。

生成时间：{now.strftime('%Y-%m-%d %H:%M:%S %Z')}
'''
    (root / "constrained/REVIEW.md").write_text(text)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
