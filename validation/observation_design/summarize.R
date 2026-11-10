#!/usr/bin/env Rscript
# Human-readable development report; never used to choose or tune the model.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
root<-file.path(repo,"validation/observation_design")
out<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(root,"results"))
status<-read.csv(file.path(out,"all_fit_status.csv"))
metrics<-read.csv(file.path(out,"all_case_metrics.csv"))
availability<-read.csv(file.path(out,"availability_by_design.csv"))
format_table<-function(x,digits=4) {
  x<-as.data.frame(x,stringsAsFactors=FALSE)
  for(j in seq_along(x)) {
    if(is.logical(x[[j]]))x[[j]]<-ifelse(x[[j]],"是","否")
    else if(is.numeric(x[[j]]))x[[j]]<-ifelse(is.finite(x[[j]]),formatC(x[[j]],digits=digits,format="f"),"NA")
    else x[[j]]<-as.character(x[[j]])
  }
  c(paste0("| ",paste(names(x),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),
    apply(x,1,function(z)paste0("| ",paste(z,collapse=" | ")," |")))
}
relative_out<-substring(normalizePath(out),nchar(normalizePath(root))+2L)
if(any(status$status!="passed_quadrature_refinement")) {
  unavailable<-status[status$status!="passed_quadrature_refinement",c("case_id","status","elapsed_seconds","recovery")]
  lines<-c("# 固定细胞预算的观测设计：数值与资源关卡未通过","",
    sprintf("60 个预定案例均已尝试；%d 个数值可用，%d 个仍不可用。完整五设计精度排名不成立，不能声称所有距离的 fitness 和不确定性均获改善。",sum(status$status=="passed_quadrature_refinement"),nrow(unavailable)),"",
    "没有丢弃失败 seed、家族或距离，没有用有限值均值替代缺失结果。未扩展到 81/64 节点，未开启预留校准/测试。","",
    format_table(availability,0),"","## 不可用案例","",format_table(unavailable,2),"",
    "原始状态、计算预算和恢复标签见 [all_fit_status.csv](results/all_fit_status.csv)；不可用案例的逐距离 prediction/interval availability 为零、风险为 NA，保留在 [all_case_metrics.csv](results/all_case_metrics.csv)。","",
    "R 内部超时可能被 TMB 捕获，原批次出现超过 300 秒的案例。外层独立 R 子进程的硬时限已修复控制问题；中断前耗时累计，最多两例并行。详见 [DEADLINE_AMENDMENT.md](DEADLINE_AMENDMENT.md)。已有完整缓存可离线重构，但不能补造缺失拟合或使超预算变成通过。","",
    "已完成 134 项模型/数据核对。共同先验、参照、目标节点和原始距离分组保持固定；具体设计与修订见 [PROTOCOL.md](PROTOCOL.md)、[DESIGN_AMENDMENT.md](DESIGN_AMENDMENT.md)。","",
    "证据值的独立联合众数/Hessian 核对、原始失败与不进行新拟合的恢复规则见 [EVIDENCE_RECONSTRUCTION.md](EVIDENCE_RECONSTRUCTION.md) 和 [POSTPROCESSING_AMENDMENT.md](POSTPROCESSING_AMENDMENT.md)。分析规则在拟合启动后、精度审阅前固定，属于内部开发规则，非外部预注册。","",
    "下一阶段应先解决原型在计算预算内的积分效率和证据值可靠性，再完成同一冻结比较。当前数据不足以推荐一个替换默认方法的设计。生产默认、用户已有修改保持不变；未 commit、push 或提交 HPC 任务。","",
    paste("生成时间：",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z",tz="America/New_York")))
  lines<-gsub("](results/",paste0("](",relative_out,"/"),lines,fixed=TRUE)
  writeLines(lines,file.path(root,"REVIEW.md"));cat("Wrote incomplete-case review; full ranking withheld\n")
  quit(save="no",status=0)
}
overall<-read.csv(file.path(out,"baseline_overall_metrics.csv"))
family<-read.csv(file.path(out,"baseline_family_metrics.csv"))
ratios<-read.csv(file.path(out,"family_risk_ratios.csv"))
triage<-read.csv(file.path(out,"development_triage.csv"))
audit<-read.csv(file.path(out,"posterior_integrity_audit.csv"))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
stopifnot(nrow(status)==60L,all(status$status=="passed_quadrature_refinement"),nrow(audit)==60L)
designs<-c("control300","depth600","time600","composition_anchor600","composition600")
risks<-c("rmse","mse","shape_rmse","shape_mse","offset_squared","edge_rmse",
  "score50","score80","score95","coverage50","coverage80","coverage95","width50","width80","width95")
# These diagnostic views cannot replace the frozen primary groups.
for(view in c("planned_distance","actual_distance")) {
  x<-metrics[metrics$view==view,]
  l<-aggregate(x[,risks],x[,c("design","family","landscape_id","distance")],mean,na.action=na.pass)
  f<-aggregate(l[,risks],l[,c("design","family","distance")],mean,na.action=na.pass)
  o<-aggregate(f[,risks],f[,c("design","distance")],mean,na.action=na.pass)
  write.csv(l,file.path(out,paste0(view,"_landscape_metrics.csv")),row.names=FALSE)
  write.csv(f,file.path(out,paste0(view,"_family_metrics.csv")),row.names=FALSE)
  write.csv(o,file.path(out,paste0(view,"_overall_metrics.csv")),row.names=FALSE)
}
denominators<-do.call(rbind,lapply(split(metrics,paste(metrics$design,metrics$view,metrics$distance)),function(x)
  data.frame(design=x$design[1],view=x$view[1],distance=x$distance[1],biological_units=nrow(x),
    landscapes=length(unique(x$landscape_id)),families=length(unique(x$family)),
    minimum_nodes=min(x$n_nodes),maximum_nodes=max(x$n_nodes),edge_available_units=sum(is.finite(x$edge_rmse)))))
write.csv(denominators,file.path(out,"distance_risk_denominators.csv"),row.names=FALSE)
support<-do.call(rbind,lapply(manifest$file,function(p) {
  c<-readRDS(p);data.frame(case_id=c$id,design=c$design,family=c$family,
    baseline_observed_nodes=sum(c$baseline_distance==0),actual_observed_nodes=sum(c$actual_distance==0),
    baseline_max_distance=max(c$baseline_distance),actual_max_distance=max(c$actual_distance))
}))
write.csv(support,file.path(out,"support_by_case.csv"),row.names=FALSE)
fisher<-read.csv(file.path(out,"fisher_planning.csv"))
fisher_average<-aggregate(fisher[,c("effective_directions","relative_mean_variance","shape_mean_variance")],
  fisher[,"design",drop=FALSE],mean)
write.csv(fisher_average,file.path(out,"fisher_average_planning.csv"),row.names=FALSE)
percent<-list()
for(reference in c("control300","depth600"))for(design in setdiff(designs,reference)) {
  a<-overall[overall$design==reference,];b<-overall[overall$design==design,]
  j<-match(a$distance,b$distance);z<-data.frame(design=design,reference=reference,distance=a$distance)
  for(metric in c("rmse","shape_rmse","edge_rmse","score50","score80","score95"))
    z[[paste0(metric,"_percent_change")]]<-100*(b[[metric]][j]/a[[metric]]-1)
  percent[[paste(reference,design)]]<-z
}
percent<-do.call(rbind,percent)
write.csv(percent,file.path(out,"overall_risk_percent_changes.csv"),row.names=FALSE)
harm<-ratios[ratios$reference=="depth600" & ratios$design%in%c("time600","composition_anchor600") & ratios$distance>=0,]
harm<-harm[apply(harm[,c("rmse","shape_rmse","edge_rmse","score95")],1,max)>1.05,]
write.csv(harm,file.path(out,"primary_family_safety_flags.csv"),row.names=FALSE)
order_rows<-function(x) x[order(x$distance,match(x$design,designs)),]
primary<-order_rows(overall[overall$distance>=0,c("design","distance","rmse","shape_rmse","edge_rmse","score95","coverage95")])
cost<-unique(manifest[,c("design","typed_cells","initial_cells","cultures","observations")])
cost<-cost[match(designs,cost$design),]
cost$median_seconds<-vapply(cost$design,function(d)median(status$elapsed_seconds[status$design==d]),numeric(1))
cost$maximum_seconds<-vapply(cost$design,function(d)max(status$elapsed_seconds[status$design==d]),numeric(1))
benefits<-order_rows(percent[percent$reference=="depth600"&percent$design%in%c("time600","composition_anchor600")&percent$distance>=0,
  c("design","distance","rmse_percent_change","shape_rmse_percent_change","edge_rmse_percent_change","score95_percent_change")])
gate<-if(any(triage$all_distance_promising))
  "开发阶段至少一个设计通过预设的所有距离筛选；这仅允许考虑更大域的独立资源与验证关卡，尚不支持替换默认方法。" else
  "本轮没有设计通过预设的所有距离筛选。因此不扩展到 81/64 节点，不开启预留校准或正式测试，也不替换默认方法。"
lines<-c("# 固定细胞预算的观测设计：开发阶段审阅", "",gate,"",
  "本轮 60 个预定案例均已尝试，且都有数值可用的后验近似；部分原始失败通过既有完整缓存的离线证据重构恢复，但超预算事实保留。6 个独立 landscape（smooth/rough/interaction 各 2 个），每个有 2 个生物重复，比较 5 种布局。所有风险先在 landscape 内平均生物重复，再在家族内平均，最后给三个家族相同权重。12 个生物单元不视为 12 个独立 landscape。", "",
  format_table(availability,0),"",
  "## 设计、成本与推断", "",format_table(cost,2),"",
  "主要等预算比较为 time600、composition_anchor600 相对 depth600；control300 用于量化增加分型细胞数的效果。composition600 是未共享初始锚点的诊断对照，完整保留其结果。600 个分型细胞并不等于相同的培养或操作成本。所有布局初始总细胞数均为 1,200，终点 T=1.5，无生长量观测。", "",
  "共同 fitness 使用原有 proper Gaussian prior；各培养分别估计未知初始组成，初始 logits 独立 N(0,2²)，没有把真实或名义初始组成提供给推断。对原有高斯 baseline 系数进行积分，保留其先验而不新增信息。仍有条件 Laplace 近似；积分收敛不等于完整后验精确。", "",
  "生成模型使用精确 Gillespie 有限细胞分裂过程；读数为不消耗细胞的 multinomial 采样。拟合采用确定性均值轨迹的 multinomial 似然，未显式建模分支过程噪声与跨时间的额外相关性。区间结果因此是对这一固定估计方法的开发检验。", "",
  "## 原始距离层级的主要结果", "",
  "每个配对单元沿用 control300 实际观测支持计算的距离；同一单元的五种设计评估同一组节点。所有结果使用固定的名义 A 参照差值。Shape 只去除一次全域误差均值，绝不逐距离重新居中。Edge 评估相邻节点的 fitness 差值误差。RMSE、shape、edge 与 IS 越小越好；覆盖率须结合宽度和 proper interval score 阅读。", "",
  format_table(primary),"",
  "### 相对相同 600 细胞深度布局的变化（%）", "",
  "负值为风险下降；这不是置信界限或正式非劣效检验。", "",format_table(benefits,2),"",
  "### 预先固定的开发筛选", "",format_table(triage),"",
  "筛选要求各原始距离的点风险不劣、centered shape 至少改善 5%、edge 不劣、IS50/80/95 不劣，并且每个家族/距离的 RMSE、shape、edge、IS95 比值不超过 1.05。所有条件须同时满足，候选及深度参照还须通过计算预算。6 个 landscape 仅提供描述性证据，不能认证显著性、正式覆盖率或最终测试样本量。规则在拟合启动后、精度结果审阅前固定，属于内部开发规则，非外部预注册；时间与数值修订见 [POSTPROCESSING_AMENDMENT.md](POSTPROCESSING_AMENDMENT.md)。", "",
  "### 家族安全界限触发项", "",
  if(nrow(harm))format_table(harm)else"未触发 1.05 家族风险比值界限。", "",
  "全部六个 landscape 的配对差值见 [paired_landscape_differences.csv](results/paired_landscape_differences.csv)；完整家族风险比值见 [family_risk_ratios.csv](results/family_risk_ratios.csv)。", "",
  "## 信息规划与距离重新分类", "",format_table(fisher_average),"",
  "Fisher 值是在 9 个预先指定参数点的局部规划平均，profile 掉初始组成的 nuisance 参数，未加入 eta prior。它衡量局部信息，不是观测到的准确率；不能用它覆盖上面的实际拟合风险。倍增深度导致相同参数点的数据 Fisher 恰好倍增的核对已通过。", "",
  "原始两个初始组成的 15 节点支持并集在 beta=0 时有两个分量（birth-contrast rank 13）；加入共同 2.2.2 锚点后为一个分量（rank 14）。这是在主要精度比较前完成的拓扑修订，详见 [DESIGN_AMENDMENT.md](DESIGN_AMENDMENT.md)。beta=.005 下迁移可弱连接原始支持，不能据此断言不可辨识。", "",
  "每种设计自己的观测距离结果见 [actual_distance_overall_metrics.csv](results/actual_distance_overall_metrics.csv)；固定名义 A 支持结果见 [planned_distance_overall_metrics.csv](results/planned_distance_overall_metrics.csv)。实际支持不同会重新分类节点，不能替代上面的固定距离比较；某距离在个别样本不存在时，其均值也对应不同样本集合。逐样本支持数量见 [support_by_case.csv](results/support_by_case.csv)，每种视图的样本、landscape、家族和节点分母见 [distance_risk_denominators.csv](results/distance_risk_denominators.csv)。", "",
  "## 数值与完整性", "",
  sprintf("134 项单元核对通过；60 个保存的混合后验近似均通过全部正权重节点、前向误差与 9/17/33/65 阶积分细化检查。最终均值最大变化 %.6g，协方差最大相对变化 %.6g，参照差值协方差最大变化 %.6g。",max(audit$mean_change),max(audit$covariance_change),max(audit$contrast_covariance_change)),"",
  sprintf("独立重构保存的 lognormal 分量矩与混合总协方差：最大均值差 %.3g，最大协方差差 %.3g；全域 MSE=shape MSE+offset² 的最大残差 %.3g。所有输入与实际消费的源码 SHA256 未变。",max(audit$mixture_mean_error),max(audit$mixture_covariance_error),max(audit$global_mse_decomposition_error)),"",
  "最初单例资源检查用 6.798 GiB，超过 6 GiB 门槛，已完整保留。独立原型改用正 Taylor 前向运算并主动释放 AD tape 后，同一数据用 2.670 GiB / 109.16 秒通过；均值及协方差变化小于数值阈值。此为计算修复，不改变生物模型或生产代码，详见 [RESOURCE_AMENDMENT.md](RESOURCE_AMENDMENT.md)。原批次的 R 超时随后被证实可能被 TMB 内部捕获，出现超预算案例；后续改为外层独立子进程硬时限，累计中断前耗时，最多两例并行。原批次只报告单例 smoke 的峰值与抽查 RSS，后续硬控制案例保存逐例峰值。完整耗时与资源失败见 [all_fit_status.csv](results/all_fit_status.csv)、[resource_failures.csv](results/resource_failures.csv) 和 [DEADLINE_AMENDMENT.md](DEADLINE_AMENDMENT.md)。", "",
  "NaN 证据值通过保存的联合众数与 Hessian 独立核对并离线重构，原始失败未覆盖，详见 [EVIDENCE_RECONSTRUCTION.md](EVIDENCE_RECONSTRUCTION.md)。重构不执行优化，也不使计算预算失败变成通过。", "",
  "本轮 fit/source 冻结副本在 results/computation_source_v4；数据、数值状态、后验、全部距离视图和失败的早期资源检查均保留。没有按结果丢弃 seed、节点、家族或设计。生产函数、默认路径及用户已有修改保持本轮启动时状态；未 commit、push 或提交 HPC 任务。", "",
  "## 证据边界与后续", "",gate,"",
  "无论某个布局的平均风险是否下降，本开发集只覆盖 27 节点、固定 T、beta、初始群体规模和一个分型预算。不能外推到所有图域、真实数据或任意 fitness 景观。后续更改应明确区分新增观测信息、结构先验失配、有限群体过程噪声和残余条件 Laplace 误差；不得把一个来源的改善解释为其他来源已解决。", "",
  paste("生成时间：",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z",tz="America/New_York")))
lines<-gsub("](results/",paste0("](",relative_out,"/"),lines,fixed=TRUE)
writeLines(lines,file.path(root,"REVIEW.md"))
cat("Wrote",file.path(root,"REVIEW.md"),"\n")
