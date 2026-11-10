#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT"));src<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_SOURCE"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
source(file.path(root,"results/computation_source_v4/score.R"))
source(file.path(root,"results/adaptive_source_v1/quadrature.R"))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"));selection<-readRDS(file.path(out,"selected_designs.rds"))
expected_n<-12L*(length(selection)+1L);stopifnot(nrow(manifest)==expected_n,!anyDuplicated(manifest$case_id))
ss<-list.files(file.path(out,"workers"),"^status[.]csv$",recursive=TRUE,full.names=TRUE)
statuses<-do.call(rbind,lapply(ss,read.csv));stopifnot(!anyDuplicated(statuses$case_id))
status<-merge(manifest,statuses,by="case_id",all.x=TRUE,sort=FALSE)
status$status[is.na(status$status)]<-"not_started_after_failed_preflight"
status$resource_pass<-!is.na(status$resource_pass)&tolower(as.character(status$resource_pass))=="true"
status$result_path<-file.path(out,"fits",status$case_id,paste0(status$case_id,"_result.rds"))
write.csv(status,file.path(out,"all_fit_status.csv"),row.names=FALSE)
all_metrics<-audits<-list();source_files<-character()
for(i in seq_len(nrow(manifest))) {
  row<-manifest[i,];s<-status[match(row$case_id,status$case_id),];case<-readRDS(row$file)
  stopifnot(unname(round2_sha256(row$file))==row$sha256,sum(case$y)==600,case$total_initial_cells==1200)
  if(file.exists(s$result_path)&&s$status=="passed_quadrature_refinement") {
    result<-readRDS(s$result_path);fit<-result$fit
    stopifnot(result$key$dataset==row$sha256,fit$status==s$status)
    n<-length(fit$mean);w<-fit$weights
    stopifnot(all(w>=0),abs(sum(w)-1)<1e-12,all(is.finite(fit$mean)),all(is.finite(fit$covariance)))
    means<-lapply(fit$components,function(c)exp(c$latent_mean+diag(c$latent_covariance)/2))
    Cs<-Map(function(m,c)tcrossprod(m)*expm1(c$latent_covariance),means,fit$components)
    mu<-Reduce(`+`,Map(`*`,means,w));C<-Reduce(`+`,Map(function(m,V,p)p*(V+tcrossprod(m-mu)),means,Cs,w))
    A<-diag(n)-matrix(case$reference_weight,n,n,byrow=TRUE)
    mean_error<-max(abs(mu-fit$mean));cov_error<-max(abs(C-fit$covariance))
    ref_error<-max(abs(result$scored$nodes$pred_relative-as.numeric(A%*%mu)))
    stopifnot(mean_error<1e-9,cov_error<1e-8,ref_error<1e-10,
      all(fit$error_estimates<=c(normalizer=1e-4,mean=2e-4,covariance=.001,contrast=.001)),
      tail(fit$history$passed,1),tail(fit$history$tail_confirmed,1),tail(fit$tail_history$passed,1),
      all(fit$conditional_status$status=="trusted_direct_joint_hessian"),
      all(fit$conditional_status$gradient<=1e-6),all(fit$conditional_status$forward_error<=1e-5),
      all(fit$conditional_status$backend_error<=1e-8)|is.null(fit$conditional_status$backend_error))
    partition<-fit$partition[order(fit$partition$lower),]
    stopifnot(partition$lower[1]==0,tail(partition$upper,1)==1,
      max(abs(head(partition$upper,-1)-tail(partition$lower,-1)))<1e-14)
    rule<-adaptive_gk15_rule();u<-q<-numeric()
    for(j in seq_len(nrow(partition))) {
      half<-(partition$upper[j]-partition$lower[j])/2
      u<-c(u,(partition$upper[j]+partition$lower[j])/2+half*rule$x);q<-c(q,half*rule$k)
    }
    coord<-vapply(fit$components,`[[`,numeric(1),"integration_coordinate")
    idx<-match(sprintf("%.17g",coord),sprintf("%.17g",u));stopifnot(!anyNA(idx))
    a<-vapply(fit$components,`[[`,numeric(1),"intercept");ell<-vapply(fit$components,`[[`,numeric(1),"log_evidence")
    jac<-dnorm(a,0,fit$prior_sd,log=TRUE)+log(fit$prior_sd)-log(coord)-log1p(-coord)
    lw<-log(q[idx])+ell+jac;ww<-exp(lw-max(lw));ww<-ww/sum(ww)
    weight_error<-max(abs(ww-w));stopifnot(weight_error<1e-12)
    full<-result$scored$metrics[result$scored$metrics$distance==-1,]
    decomposition_error<-max(abs(full$mse-full$shape_mse-full$offset_squared));stopifnot(decomposition_error<1e-12)
    audits[[row$case_id]]<-data.frame(case_id=row$case_id,mean_error=mean_error,covariance_error=cov_error,
      reference_error=ref_error,jacobian_weight_error=weight_error,mse_decomposition_error=decomposition_error,
      max_gradient=max(fit$conditional_status$gradient),max_forward_error=max(fit$conditional_status$forward_error))
    metrics<-result$scored$metrics
    source_files<-c(source_files,file.path(out,"fits",row$case_id,"source_sha256.csv"))
  } else metrics<-design_score(case,list(status=s$status))$metrics
  expected<-unlist(lapply(c("baseline_distance","planned_distance","actual_distance"),function(view)paste(view,c(-1L,sort(unique(case[[view]]))))))
  stopifnot(setequal(paste(metrics$view,metrics$distance),expected))
  all_metrics[[row$case_id]]<-metrics
}
metrics<-do.call(rbind,all_metrics);write.csv(metrics,file.path(out,"all_case_metrics.csv"),row.names=FALSE)
if(length(audits))write.csv(do.call(rbind,audits),file.path(out,"posterior_integrity_audit.csv"),row.names=FALSE)
sources<-c(source_files,file.path(out,"generation_source_sha256.csv"),file.path(out,"planning_consumed_source_sha256.csv"),
  file.path(out,"selection_freeze_sha256.csv"),file.path(out,"noise/source_sha256.csv"),
  file.path(root,"results/constrained_analysis_source_v1_sha256.csv"))
rows<-do.call(rbind,lapply(sources,function(p){z<-read.csv(p);z$manifest_path<-p;z}))
hash<-round2_sha256(unique(rows$file));rows$current_sha256<-unname(hash[rows$file]);stopifnot(all(rows$sha256==rows$current_sha256))
write.csv(rows,file.path(out,"consumed_source_integrity.csv"),row.names=FALSE)
if(any(status$status!="passed_quadrature_refinement")||any(!status$resource_pass)) {
  writeLines("INCOMPLETE_NUMERICAL_OR_RESOURCE_GATE: full accuracy comparison withheld",file.path(out,"stage_status.txt"))
  stop("Incomplete numerical or resource gate; all cases and unavailable risks retained")
}
stopifnot(nrow(status)==expected_n,all(metrics$prediction_fraction==1),all(metrics$interval_fraction50==1),
  all(metrics$interval_fraction80==1),all(metrics$interval_fraction95==1))
primary<-metrics[metrics$view=="baseline_distance",]
risks<-c("rmse","mse","shape_rmse","shape_mse","offset_squared","edge_rmse","score50","score80","score95",
  "coverage50","coverage80","coverage95","width50","width80","width95")
stopifnot(all(is.finite(as.matrix(primary[,risks]))),all(primary$n_edges>0))
land<-aggregate(primary[,risks],primary[,c("design","family","landscape_id","distance")],mean)
family<-aggregate(land[,risks],land[,c("design","family","distance")],mean)
overall<-aggregate(family[,risks],family[,c("design","distance")],mean)
write.csv(land,file.path(out,"baseline_landscape_metrics.csv"),row.names=FALSE)
write.csv(family,file.path(out,"baseline_family_metrics.csv"),row.names=FALSE)
write.csv(overall,file.path(out,"baseline_overall_metrics.csv"),row.names=FALSE)
ratios<-changes<-differences<-flags<-list()
core<-c("rmse","shape_rmse","edge_rmse","score50","score80","score95")
ref<-overall[overall$design=="depth600"&overall$distance>=0,]
for(design in setdiff(unique(overall$design),"depth600")) {
  z<-overall[overall$design==design&overall$distance>=0,];z<-z[match(ref$distance,z$distance),]
  change<-data.frame(design=design,distance=ref$distance)
  for(k in core)change[[k]]<-100*(z[[k]]/ref[[k]]-1)
  changes[[design]]<-change
  b<-family[family$design=="depth600",];f<-family[family$design==design,];ii<-match(paste(f$family,f$distance),paste(b$family,b$distance))
  rr<-f[,c("family","distance")];rr$design<-design
  for(k in core)rr[[k]]<-f[[k]]/b[[k]][ii]
  ratios[[design]]<-rr
  base<-land[land$design=="depth600",];x<-land[land$design==design,];ii<-match(paste(x$landscape_id,x$distance),paste(base$landscape_id,base$distance))
  dd<-x[,c("family","landscape_id","distance")];dd$design<-design
  for(k in core)dd[[k]]<-x[[k]]-base[[k]][ii]
  differences[[design]]<-dd
  f<-c(point_nonworse=all(z$rmse<=ref$rmse),shape_5percent_all_distances=all(z$shape_rmse<=.95*ref$shape_rmse),
    edge_nonworse=all(z$edge_rmse<=ref$edge_rmse),interval_scores_nonworse=all(as.matrix(z[,c("score50","score80","score95")])<=as.matrix(ref[,c("score50","score80","score95")])),
    family_safety=all(as.matrix(rr[rr$distance>=0,c("rmse","shape_rmse","edge_rmse","score95")])<=1.05))
  flags[[design]]<-cbind(data.frame(design=design,all_distance_promising=all(f),resource_budget_pass=TRUE),as.data.frame(as.list(f)))
}
write.csv(do.call(rbind,ratios),file.path(out,"family_risk_ratios.csv"),row.names=FALSE)
write.csv(do.call(rbind,changes),file.path(out,"overall_risk_changes.csv"),row.names=FALSE)
write.csv(do.call(rbind,differences),file.path(out,"paired_landscape_differences.csv"),row.names=FALSE)
triage<-do.call(rbind,flags);write.csv(triage,file.path(out,"development_triage.csv"),row.names=FALSE)
writeLines(if(any(triage$all_distance_promising))"DEVELOPMENT_PROMISING_REQUIRES_INDEPENDENT_LARGER_DOMAIN_VALIDATION"else
  "NO_ALL_DISTANCE_WINNER_UNDER_FROZEN_NEAR_INFORMATION_CONSTRAINT",file.path(out,"stage_status.txt"))
print(do.call(rbind,changes),row.names=FALSE);print(triage,row.names=FALSE)
