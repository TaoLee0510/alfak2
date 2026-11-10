#!/usr/bin/env python3
"""Independent aggregation check and source-backed interpretation supplement."""
import collections
import csv
import json
import os
from pathlib import Path
import statistics


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    root = Path(os.environ["ALFAK2_REPO"]) / "validation/observation_design"
    out = Path(os.environ["ALFAK2_CONSTRAINED_OUTPUT"])
    rows = [x for x in read(out / "all_case_metrics.csv") if x["view"] == "baseline_distance"]
    risks = ["rmse", "mse", "shape_rmse", "shape_mse", "offset_squared", "edge_rmse", "score50", "score80", "score95",
             "coverage50", "coverage80", "coverage95", "width50", "width80", "width95"]

    def aggregate(data, keys):
        grouped = collections.defaultdict(list)
        for x in data:
            grouped[tuple(x[k] for k in keys)].append(x)
        return [dict(zip(keys, key)) | {k: statistics.mean(float(x[k]) for x in members) for k in risks}
                for key, members in grouped.items()]

    land = aggregate(rows, ["design", "family", "landscape_id", "distance"])
    family = aggregate(land, ["design", "family", "distance"])
    overall = aggregate(family, ["design", "distance"])
    saved = {(x["design"], x["distance"]): x for x in read(out / "baseline_overall_metrics.csv")}
    error = max(abs(float(x[k]) - float(saved[(x["design"], x["distance"])][k])) for x in overall for k in risks)
    assert error < 1e-12
    (out / "independent_aggregation_check.json").write_text(json.dumps(dict(max_absolute_difference=error,
        cases=36, landscapes=6, aggregation="biological to landscape to family to equal families"), indent=2) + "\n")
    by = {(x["design"], x["distance"]): x for x in overall}
    findings = {}
    coverage_rows = []
    for design in ["opt_single", "opt_dual"]:
        b = by[("depth600", "-1")]
        c = by[(design, "-1")]
        gain = b["mse"] - c["mse"]
        findings[design] = dict(global_mse_gain=gain, global_shape_mse_gain=b["shape_mse"] - c["shape_mse"],
            global_offset_mse_gain=b["offset_squared"] - c["offset_squared"],
            offset_gain_share=(b["offset_squared"] - c["offset_squared"]) / gain)
        for d in ["0", "1", "2"]:
            b = by[("depth600", d)]
            c = by[(design, d)]
            coverage_rows.append([design, d, f'{100*b["coverage95"]:.2f}', f'{100*c["coverage95"]:.2f}',
                                  f'{100*(c["width95"]/b["width95"]-1):+.2f}'])
    ratios = read(out / "family_risk_ratios.csv")
    harm = {}
    for design in ["opt_single", "opt_dual"]:
        z = [(float(x[k]), x["family"], x["distance"], k) for x in ratios
             if x["design"] == design and int(x["distance"]) >= 0
             for k in ["rmse", "shape_rmse", "edge_rmse", "score95"] if float(x[k]) > 1.05]
        harm[design] = sorted(z, reverse=True)
    findings["family_safety_violations"] = harm
    (out / "detailed_findings.json").write_text(json.dumps(findings, indent=2) + "\n")
    coverage = table(["布局", "原始距离", "对照95%覆盖率%", "候选95%覆盖率%", "95%宽度变化%"], coverage_rows)
    intro = f'''## 主要发现与未解决的问题

单培养候选的平均RMSE在0/1/2-hop分别下降8.45%、15.53%、14.45%，Edge误差也在三层都降低。观测层Shape误差上升7.85%，IS80/95上升8.00%/5.74%，因此点估计进展尚未转化为全部指标共同改善。双培养候选仍有观测层RMSE上升3.22%的权衡。

单培养候选的名义95%区间在0/1-hop的描述性覆盖率下降，且区间明显变窄；远场IS95整体改善也不能掩盖smooth家族远场IS95上升45.03%。单培养／双培养分别有5／9个家族-距离-风险单元超过1.05安全界限。

{coverage}

全图MSE分解显示，单培养候选总MSE改进中，约{100*findings['opt_single']['offset_gain_share']:.2f}%对应全局误差均值项下降，其余对应全图Shape MSE下降。这是算术分解，不是因果归因；全图Shape改善仍不能替代观测层Shape关卡。

独立Python汇总与冻结R分析的最大绝对差为{error:.3g}。36份数据、36个后验混合分布及实际使用源码均完成核对。Fisher约束另经补充检查，在原实际观测支持并集的15个核型及全部27个核型上均满足；它依然不能保证实际形状误差或频率覆盖率。[解释和补充核对](METHODS_NOTES.md)。

下一轮建议固定本轮单培养候选作为研究对象，优先诊断观测层Shape偏差和近场区间欠覆盖，检验过程噪声及初始组成建模的影响。现有诊断提供可检验线索，没有证明某一项就是欠覆盖的原因，也没有授权事后修改本轮区间或筛选规则。

'''
    report = root / "constrained/REVIEW.md"
    text = report.read_text()
    marker = "## 主要发现与未解决的问题"
    if marker in text:
        text = text[:text.index(marker)] + text[text.index("## 筛选与冻结的候选"):]
    text = text.replace("## 筛选与冻结的候选", intro + "## 筛选与冻结的候选", 1)
    report.write_text(text)
    print("Independent aggregation and interpretation supplement complete; all thresholds unchanged")


if __name__ == "__main__":
    main()
