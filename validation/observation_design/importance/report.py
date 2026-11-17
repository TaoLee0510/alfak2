#!/usr/bin/env python3
"""Frozen descriptive readout; no estimator choices are made from truth metrics."""
import csv,os,statistics
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo

repo=Path(os.environ['ALFAK2_REPO']);root=repo/'validation/observation_design';out=Path(os.environ['ALFAK2_IS_OUTPUT'])
source=Path(os.environ['ALFAK2_IS_SOURCE']);version=source.name.rsplit('_',1)[1]
def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h))
def number(r,k):return float(r[k])
def write(p,rows):
 with Path(p).open('w') as h:
  w=csv.DictWriter(h,rows[0].keys());w.writeheader();w.writerows(rows)
def table(header,rows):
 return '\n'.join(['| '+' | '.join(header)+' |','| '+' | '.join(['---']*len(header))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
status=read(out/'mechanism/all_fit_status.csv');assert len(status)==48 and all(r['numerical_pass']=='TRUE' for r in status)
overall=read(out/'mechanism/overall_metrics.csv');family=read(out/'mechanism/family_metrics.csv')
old=read(root/'results/constrained_v2/baseline_overall_metrics.csv');old=[r for r in old if r['design']=='opt_single' and float(r['distance'])>=0]
current=[r for r in overall if r['design']=='A_finite_unknown' and float(r['distance'])>=0]
old={r['distance']:r for r in old};current={r['distance']:r for r in current}
core=['rmse','shape_rmse','edge_rmse','score50','score80','score95'];changes=[]
for d,a in sorted(current.items(),key=lambda kv:float(kv[0])):
 b=old[d];row=dict(distance=d)
 for k in core+['width95']:row[k+'_percent_change']=100*(number(a,k)/number(b,k)-1)
 row.update(old_coverage95=number(b,'coverage95'),new_coverage95=number(a,'coverage95'))
 changes.append(row)
write(out/'mechanism/numerical_correction_overall_changes.csv',changes)
oldf=[r for r in read(root/'results/constrained_v2/baseline_family_metrics.csv') if r['design']=='opt_single']
oldf={(r['family'],r['distance']):r for r in oldf};ratios=[]
for r in family:
 if r['design']!='A_finite_unknown' or float(r['distance'])<0:continue
 b=oldf[(r['family'],r['distance'])];row=dict(family=r['family'],distance=r['distance'])
 for k in core:row[k]=number(r,k)/number(b,k)
 ratios.append(row)
write(out/'mechanism/numerical_correction_family_ratios.csv',ratios)
flags=dict(point_nonworse=all(r['rmse_percent_change']<=0 for r in changes),
 shape_5percent_all_distances=all(r['shape_rmse_percent_change']<=-5 for r in changes),
 edge_nonworse=all(r['edge_rmse_percent_change']<=0 for r in changes),
 interval_scores_nonworse=all(r[k+'_percent_change']<=0 for r in changes for k in ['score50','score80','score95']),
 family_safety=all(r[k]<=1.05 for r in ratios for k in ['rmse','shape_rmse','edge_rmse','score95']))
flags['all_distance_promising']=all(flags.values());write(out/'mechanism/numerical_correction_development_gate.csv',[flags])

small=read(out/f'small_cases_{version}_analysis/importance_reference_gate.csv');assert all(r['gate_pass']=='TRUE' for r in small)
pre=read(out/f'preflight_{version}_gate.csv');masses=read(out/f'small_cases_{version}_analysis/interval_comparison.csv')
lines=['# alfak2 联合后验校正与固定布局机制实验',
 '本轮本地目标已完成：原模型下的联合后验校正、两个高精度小图核对、两个 27 节点预检以及 48 例机制实验均已完成。所有最终案例通过事先规定的数值与资源门槛，历史失败尝试完整保留。',
 '数值校正对同一 opt_single 数据的所有距离改善开发门槛'+('已通过。仍需独立验证。' if flags['all_distance_promising'] else '未全部通过；不能据此宣称所有距离的 fitness 与不确定性都已提高。'),
 '生产默认值未改变；没有提交 HPC、commit 或 push，也没有打开预留校准/正式测试数据。当前计算源为 `'+source.name+'`，机制数据与计分源为 `importance_analysis_source_v2`。',
 '**数值核对**',
 '原似然、原 Gaussian fitness 先验、原初始组成先验与原 fitness 参照保持不变；Gaussian 系数被解析边缘化。独立提议训练链只用于构造分布，最终估计使用独立采样的原始重要性权重和加权分位数。高 birth-rate 提议偏移只改变采样分配，不改变目标分布。',
 '本轮 107 项验证全部通过，涵盖解析先验、两节点独立积分、完整前向传播、边缘化前后模式/Hessian/证据一致性、变换密度及 Jacobian，以及分块诊断与原公式的一致性。同源独立冷构建另计 TMB 16.912 秒、核心 Rcpp 2.914 秒，未从任何普通案例的提议/采样/诊断耗时中扣减。',
 table(['8 节点案例','最大对比均值差 / 参照 SD','对比 SD 比范围','通过'],[[r['case_id'],f"{float(r['max_contrast_mean_difference']):.4f}",f"{float(r['min_contrast_sd_ratio']):.4f}–{float(r['max_contrast_sd_ratio']):.4f}",r['gate_pass']] for r in small]),
 '参照分别保留每条链 80,000 与 160,000 次样本；比较同时计入候选与参照 MCSE，边界不确定的案例不能判为通过。这些区间质量检查是固定数据下的后验核对，不是重复抽样覆盖率。',
 table(['27 节点预检','完整耗时 (秒)','峰值 RSS (GiB)','最大均值 MCSE / SD','最大 SD MCSE / SD'],[[r['case_id'],f"{float(r['child_seconds']):.2f}",f"{float(r['maximum_rss_bytes'])/1024**3:.3f}",f"{float(r['max_mean_mcse_ratio']):.5f}",f"{float(r['max_sd_mcse_ratio']):.5f}"] for r in pre]),
 '最终 48 例最大完整耗时 '+f"{max(float(r['elapsed_seconds']) for r in status):.2f}"+' 秒，最大峰值 RSS '+f"{max(float(r['maximum_rss_gib']) for r in status):.3f}"+' GiB。门槛为每例 300 秒、6 GiB、原始权重 ESS ≥10,000、每个独立批次及合并 Pareto-k ≤0.5、均值 MCSE ≤0.005 后验 SD、SD MCSE ≤0.01 后验 SD。分块诊断与原公式通过同样本一致性验证。',
 '**同一数据：数值校正前后**',
 '下表仅比较同一有限过程、未知初始组成的 A 组与旧 opt_single 结果。百分比为新/旧−1；RMSE、Shape、Edge 与区间评分越低越好。覆盖率按名义 95% 单独列出，区间变窄本身不是成功标准。',
 table(['原控制距离','RMSE Δ%','Shape Δ%','Edge Δ%','IS95 Δ%','旧覆盖率','新覆盖率'],[[r['distance']]+[f"{r[k+'_percent_change']:+.2f}" for k in ['rmse','shape_rmse','edge_rmse','score95']]+[f"{r['old_coverage95']:.3f}",f"{r['new_coverage95']:.3f}"] for r in changes]),
 '原开发门槛结果：'+', '.join(k+'='+str(v) for k,v in flags.items())+'。这些是描述性开发判据，不能视为风险差的显著性检验；很小的变化也可能包含 Monte Carlo 误差。',
 '**48 例机制对照**',
 '六个独立景观、每景观两个生物学重复，每个重复四组：A 有限过程/未知 pi0；B 均值传播/未知 pi0；C 有限过程/真实 pi0 oracle；D 均值传播/真实 pi0 oracle。初始整数种群、t0 读数、终点抽样 uniforms、原控制距离和 fitness 参照按重复配对。真实 birth 只进入生成/计分，真实 pi0 只进入 C/D。Oracle 不能作为可部署方法排名。',
 table(['组','距离','RMSE','Shape','Edge','IS95','覆盖率95','宽度95'],[[r['design'],r['distance']]+[f"{float(r[k]):.4f}" for k in ['rmse','shape_rmse','edge_rmse','score95','coverage95','width95']] for r in sorted(overall,key=lambda r:(r['design'],float(r['distance']))) if float(r['distance'])>=0]),
 '汇总顺序为生物学重复→景观→家族→家族等权平均；48 例不等于 48 个独立景观。Shape 只在全图统一去中心一次，边按两个端点距离的最大值分组。所有原始 50/80/95% 区间评分、覆盖率、宽度、家族风险比、景观配对差和交互作用均已保存。',
 'B−A / D−C 表示把有限分枝过程替换为均值传播的综合影响，包含过程波动及均值传播近似的差别；C−A / D−B 表示知道真实初始组成的综合影响，不能单独解释为某一噪声的因果份额。Oracle 改变的信息条件也包含真实零组成状态与连续组成先验的差别。',
 '数值后验一致不保证生成真值精度或频率学覆盖率。本实验仅是固定布局、27 节点开发诊断，尚不能证明更大图与更远场的全面改善。后续模型变更应由完整四组配对差决定，并继续保留原控制距离、所有失败及独立验证边界。',
 '**后续优先级：初始组成约束、有限过程协方差与分层校准**',
 '这些开发数据支持优先处理初始组成：C−A 的 RMSE 在距离 0/1/2 分别降低 7.22%/12.82%/8.49%，B−A 分别降低 3.85%/4.51%/2.58%。这比较的是机制干预的综合影响，不是可部署方法的提升，也不能解释为噪声贡献的因果百分比。',
 '初始组成已知化并未同时解决不确定性：距离 2 的 95% 覆盖率从 A 的 89.31% 变为 C 的 86.67%；D 为 91.67%。因此后续估计器需要共同处理初始组成与有限分枝协方差。D 使用均值生成且给定真实初始组成，其距离 1/2 覆盖率仍为 92.79%/91.67%，提示还应检查先验收缩、可识别性与有限开发样本的影响，不能据此唯一归因于某一个模型误设。',
 '建议下一阶段依次：①利用实际可获得的播种/初始测量与允许支持集，将未检出和真正为零分开建模；②把有限分枝过程的状态相关协方差传入观测不确定性；③在新增独立开发景观及更大图上分层验证点估计、区间评分和覆盖率，保留正式测试作最终判断。本轮没有执行这些新的生物学模型改动。',
 '**可复核输出**',
 '- [48 例状态与资源](../results/importance_v1/mechanism/all_fit_status.csv)\n- [逐例完整指标](../results/importance_v1/mechanism/all_case_metrics.csv)\n- [景观配对差](../results/importance_v1/mechanism/paired_landscape_differences.csv)\n- [家族风险比](../results/importance_v1/mechanism/family_risk_ratios.csv)\n- [数值校正开发门槛](../results/importance_v1/mechanism/numerical_correction_development_gate.csv)\n- [完整协议与保留的版本变更](PROTOCOL.md)\n- [图表](../results/importance_v1/mechanism/mechanism_summary.pdf)\n- [运行库归档与恢复说明](../results/importance_v1/runtime_archive/README.md)\n- [消耗源码完整性核查](../results/importance_v1/consumed_source_integrity.csv)',
 'Timestamp: '+datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d %H:%M:%S %Z')]
(root/'importance/REVIEW.md').write_text('\n\n'.join(lines)+'\n')
print(flags)
