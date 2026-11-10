#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
gate<-read.csv(file.path(out,"reference_gate.csv"));mom<-read.csv(file.path(out,"moment_comparison.csv"))
stopifnot(nrow(gate)==2L,all(gate$reference_certified),!all(gate$gate_pass))
primary_source<-file.path(root,"results/mechanism_reference_source_v2")
localization<-file.path(out,"localization")
for(i in 1:2) {
  resource<-read.csv(file.path(localization,paste0("controller_",i,"_resource.csv")))
  stopifnot(resource$exit_code==0,tolower(resource$memory_pass)=="true")
}
ids<-gate$case_id;quads<-do.call(rbind,lapply(ids,function(id)read.csv(file.path(localization,id,"quadrature_tolerance_check.csv"))))
conditional<-do.call(rbind,lapply(ids,function(id)read.csv(file.path(localization,id,"conditional_moment_comparison.csv"))))
stopifnot(all(conditional$reference_certified),!any(conditional$material_mean|conditional$material_sd))
baseline<-read.csv(file.path(out,"baseline_identity_moments.csv"));density<-read.csv(file.path(out,"baseline_identity_density.csv"))
target_equal<-identical(readRDS(file.path(root,"results/mechanism_v1/frozen_cases.rds")),readRDS(file.path(out,"frozen_cases.rds")))
stopifnot(target_equal)
validation<-read.csv(file.path(out,"validation_tests.csv"));stopifnot(sum(validation$passed)==133,all(validation$failed==0),!any(validation$error))
case_summary<-do.call(rbind,lapply(ids,function(id) {
  z<-mom[mom$case_id==id&mom$type=="relative",]
  data.frame(case_id=id,max_mean_error_reference_sd=max(abs(z$standardized_mean_difference)),
    min_sd_ratio=min(z$sd_ratio),max_sd_ratio=max(z$sd_ratio),
    material_relative_means=sum(z$material_mean),material_relative_sds=sum(z$material_sd))
}))
write.csv(case_summary,file.path(out,"case_comparison_summary.csv"),row.names=FALSE)
# Verify every consumed input/source, proposal-only pilot and frozen analysis.
manifests<-c(file.path(out,"consumed_source_integrity.csv"),file.path(out,"validation_source_sha256.csv"),
  file.path(out,"proposal_freeze_sha256.csv"),file.path(out,"baseline_identity_source_sha256.csv"),
  file.path(root,"results/mechanism_v1/validation_source_sha256.csv"),
  file.path(root,"results/mechanism_v1",ids,"source_sha256.csv"),
  file.path(localization,ids,"source_sha256.csv"),
  file.path(root,"results",c("mechanism_reference_source_v1_sha256.csv","mechanism_reference_source_v2_sha256.csv",
    "mechanism_analysis_source_v1_sha256.csv","mechanism_analysis_source_v2_sha256.csv",
    "mechanism_localization_source_v1_sha256.csv","mechanism_localization_source_v2_sha256.csv","mechanism_identity_source_v1_sha256.csv")))
ledger<-do.call(rbind,lapply(manifests,function(p){z<-read.csv(p);data.frame(file=z$file,sha256=z$sha256,manifest=p)}))
hash<-round2_sha256(unique(ledger$file));ledger$current_sha256<-unname(hash[ledger$file]);stopifnot(all(ledger$sha256==ledger$current_sha256))
write.csv(ledger,file.path(out,"final_consumed_source_integrity.csv"),row.names=FALSE)
mdtable<-function(headers,rows)paste(c(paste0("| ",paste(headers,collapse=" | ")," |"),
  paste0("| ",paste(rep("---",length(headers)),collapse=" | ")," |"),
  apply(rows,1,function(x)paste0("| ",paste(x,collapse=" | ")," |"))),collapse="\n")
summary_table<-mdtable(c("案例","最大相对均值误差／参考SD","相对SD比范围","门槛"),
  cbind(case_summary$case_id,sprintf("%.5f",case_summary$max_mean_error_reference_sd),
    sprintf("%.4f–%.4f",case_summary$min_sd_ratio,case_summary$max_sd_ratio),ifelse(gate$gate_pass,"通过","触发均值材料阈值")))
precision_table<-mdtable(c("案例","每链保留","最大Rhat","最小bulk/tail ESS","均值／SD最大MCSE比"),
  cbind(gate$case_id,gate$retained_per_chain,sprintf("%.6f",gate$max_rhat),sprintf("%.0f / %.0f",gate$min_bulk_ess,gate$min_tail_ess),
    sprintf("%.4f / %.4f",gate$max_mean_mcse_sd,gate$max_sd_mcse_sd)))
quad_table<-mdtable(c("案例","原／加密节点","均值变化／参考SD","协方差相对变化"),
  cbind(quads$case_id,paste(quads$default_nodes,quads$tight_nodes,sep=" / "),format(quads$mean_change_reference_sd,digits=3),format(quads$covariance_relative_change,digits=3)))
conditional_table<-mdtable(c("案例","固定a0=0最大相对均值误差／参考SD","材料标记"),
  do.call(rbind,lapply(ids,function(id){z<-conditional[conditional$case_id==id&conditional$type=="relative",];c(id,sprintf("%.5f",max(abs(z$standardized_mean_difference))),"无")})))
baseline_table<-mdtable(c("案例","MCMC baseline均值","近似baseline均值","差异／参考SD","差异／MCSE"),
  cbind(baseline$case_id,sprintf("%.6f",baseline$reference_mean),sprintf("%.6f",baseline$quadrature_mean),
    sprintf("%.5f",baseline$mean_error_reference_sd),sprintf("%.2f",baseline$mean_error_mcse)))
old_status<-do.call(rbind,lapply(ids,function(id)read.csv(file.path(root,"results/mechanism_v1",id,"status.csv"))))
status<-list(stage="REFERENCE_GATE_FAILED_ATTRIBUTION_BATCH_NOT_STARTED",reference_cases=2L,certified_references=2L,
  passing_cases=sum(gate$gate_pass),failing_cases=sum(!gate$gate_pass),reference_validation_assertions=133L,
  biological_target_unchanged=target_equal,tight_quadrature_max_mean_change_reference_sd=max(quads$mean_change_reference_sd),
  conditional_references_certified=2L,conditional_material_errors=sum(conditional$material_mean|conditional$material_sd),
  attribution_new_fits_started=0L,production_default_changed=FALSE,reserved_calibration_test_opened=FALSE,
  hpc_jobs_submitted=FALSE,commit_or_push=FALSE)
jsonlite::write_json(status,file.path(out,"checkpoint_status.json"),auto_unbox=TRUE,pretty=TRUE)
time<-format(Sys.time(),tz="America/New_York",format="%Y-%m-%d %H:%M:%S %Z")
text<-paste0('# 固定布局机制研究：参照门槛及数值定位\n\n',
  '两个8节点MCMC参照均已认证；一个案例通过近似核对，另一个相对均值单元触发预设材料阈值。按用户确认的“参照检查通过后，再运行配对机制对照”顺序，本轮没有启动36例机制归因拟合。\n\n',
  '触发值为0.05408个参考SD，略高于0.05的点估计阈值，差异约7.72个均值MCSE。它证明了小幅非零差异，但不是对“真实误差严格超过0.05 SD”的正式检验。保留原规则，不用事后修改阈值获得通过。\n\n',
  summary_table,'\n\n',
  '## 固定案例与参照认证\n\n',
  '两个案例均为CN1–2的三染色体8节点图、beta=.005、T=1.5、初始1200细胞、分型240/360。初始壳层总质量为.45/0/.35/.20，使用事先固定的prior-offset及interaction参数点，不按拟合结果挑选数据。每例都有三个初始零计数、终点非零的状态。参照始终是原名义A的8节点对应权重。\n\n',
  'MCMC使用相同Gaussian log-birth和SD=2初始-logit先验，不使用真实初始组成、真实birth或外部生长观测。新独立前向实现对应完整均值生成元，已核对当前TMB似然和独立R均值传播。两个采样版本各通过133项检查，检查存在重叠，不能计为266项独立断言。\n\n',
  precision_table,'\n\n',
  '原blocked ESS在每链16万次上限时仍未达到全部MCSE要求：两例均值MCSE约为后验SD的2.82%/2.93%，门槛为2.5%；该次参照没有被认证。原始链、收敛记录、耗时和源码均保留。\n\n',
  '随后用旧链的固定pilot构造Gaussian提议Q，采样时显式加入log phi(z)-log Q(z)修正，目标分布未变。Gaussian提议在正式保留采样阶段保持固定，未使用Laplace比较器的均值／协方差或生成真值来调提议。四条链各5000次warmup、正式每链2万次即通过全部Rhat、ESS和MCSE门槛；所有初始-logit二阶乘积也纳入检查。\n\n',
  '新参照各耗时约31/32秒（含各自adaptive比较器），峰值约1.72/1.59 GiB。此耗时不包含旧失败尝试和生成pilot的成本，不能当作从零开始的总运行成本或性能提升倍数。所有过程均在本地，最多两个模型进程，链在案例内串行。\n\n',
  '## 标量积分与条件后验定位\n\n',
  '把归一化、均值、协方差和参照协方差的积分容差全部收紧10倍，并保留原支持与尾部复核，得到：\n\n',
  quad_table,'\n\n',
  '当前差异不是通过继续增加这一标量积分的分辨率就能消除的。误差数值是本组案例的验证结果，不是全局数学误差界。\n\n',
  '再固定原独立intercept a0=0，保持原offset不变，严格从原K中减去sd(a0)^2*11\'，得到原模型的条件先验。未知初始组成及其原先验保持。两个独立条件MCMC参照均认证，条件Laplace比较未触发材料标记：\n\n',
  conditional_table,'\n\n',
  '这是单个a0位置的诊断，不能证明所有条件节点或其证据权重都准确。条件参照及pilot的约33秒/例成本也单独保存。\n\n',
  '## 对baseline混合权重的进一步核对\n\n',
  '利用原Gaussian层级模型的恒等式a|g,Y,eta = a|g，可从认证MCMC的g样本恢复baseline后验，而不改变先验或加入信息。具体条件均值为s²1\'K^-1(g-mu)，条件方差为s²-s⁴1\'K^-1 1。\n\n',
  baseline_table,'\n\n',
  '在prior案例，a=-0.5处的近似baseline密度比Rao–Blackwell参照低约7.55%，差异约4.70个MCSE；所有五个预先声明的位置和两个案例的密度值均保存。该核对是事后数值定位，不是新的晋级条件。\n\n',
  '加密积分不改变结果，baseline混合后验仍有偏移，数值证据因此支持优先检查剩余的条件Laplace证据／混合权重近似。单点条件检查不足以把误差全部归于权重，也不能把8节点结果直接当成27节点Shape偏差或区间欠覆盖的原因。\n\n',
  '## 本轮决定与后续方向\n\n',
  '本轮完成了参照认证和未通过门槛后的数值定位。尚未执行的工作明确为：均值轨迹／有限过程与未知／真实初始组成的2×2机制对照，共36个新增常规拟合；oracle初始组成分支未启动，未据此宣布过程噪声或初始组成处理是近场问题的根因。\n\n',
  '下一项应先在相同目标下验证对条件Laplace证据及混合权重的校正，复核这两个固定案例；达到原门槛后，再恢复已拟定的机制对照。这个建议不等于已实现校正，也不自动替换生产方法。正式测试、生产默认、HPC和Git发布状态保持。\n\n',
  '## 可复核文件\n\n',
  '- [确认协议](PROTOCOL.md)、[数值提议修正](PRECONDITIONING.md)、[条件定位说明](LOCALIZATION.md)\n',
  '- [参照门槛](../results/mechanism_v2/reference_gate.csv)、[全部均值与SD比较](../results/mechanism_v2/moment_comparison.csv)、[协方差比较](../results/mechanism_v2/covariance_comparison.csv)、[区间与后验质量](../results/mechanism_v2/interval_comparison.csv)\n',
  '- [baseline均值恒等式核对](../results/mechanism_v2/baseline_identity_moments.csv)、[全部密度位置](../results/mechanism_v2/baseline_identity_density.csv)、[输入案例](../results/mechanism_v2/case_manifest.csv)\n',
  '- [全部实际使用源码与输入完整性](../results/mechanism_v2/final_consumed_source_integrity.csv)、[主验证检查](../results/mechanism_v2/validation_tests.csv)、[阶段状态](../results/mechanism_v2/checkpoint_status.json)\n\n',
  '实现记录：首次分析脚本在比较前因Python布尔字符串被R读取为文本而停止，后续仅规范该资源标记的解析；第一次条件定位控制器在启动R前有导入调用错误，修复后运行。原失败源码和记录均保留，统计规则与输入未因这些问题改变。\n\n',
  '生成时间：',time,'\n')
writeLines(text,file.path(root,"mechanism/REVIEW.md"))
print(status)
