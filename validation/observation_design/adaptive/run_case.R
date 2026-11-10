#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_ADAPTIVE_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_ADAPTIVE_OUTPUT"))
cpp<-file.path(root,"results/direct_forward_source_v2")
Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp)
source(file.path(cpp,"model.R"));source(file.path(root,"results/direct_source_v7/direct_joint.R"))
source(file.path(src,"quadrature.R"));source(file.path(src,"model.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"));source(file.path(root,"results/computation_source_v4/score.R"))
runtime<-design_load(repo);invisible(TMB::openmp(1,DLL="observation_model"))
mode<-Sys.getenv("ALFAK2_ADAPTIVE_MODE");id<-Sys.getenv("ALFAK2_ADAPTIVE_CASE")
reference_dir<-file.path(repo,"validation/posterior_reference_4to8/results")
extra<-character()
if(mode=="primary") {
  manifest<-read.csv(file.path(root,"results/dataset_manifest.csv"))
  row<-manifest[manifest$case_id==id,];stopifnot(nrow(row)==1L,unname(round2_sha256(row$file))==row$sha256)
  case<-readRDS(row$file);context<-readRDS(file.path(root,"results/context.rds"));data_files<-c(row$file,file.path(root,"results/context.rds"))
} else {
  cases<-readRDS(file.path(reference_dir,"frozen_cases.rds"))
  old<-cases[[which(vapply(cases,function(c)c$id==id,logical(1)))]];stopifnot(!old$growth)
  context<-design_context(old$graph$karyotypes,old$beta,runtime)
  stopifnot(max(abs(context$covariance-old$covariance))<1e-10)
  case<-list(id=id,y=old$counts,times=c(0,old$specification$dt),trajectory=c(0L,0L),reference_weight=old$weight)
  data_files<-c(file.path(reference_dir,"frozen_cases.rds"),file.path(reference_dir,paste0(id,"_summary.rds")),
    file.path(repo,"validation/baseline_quadrature/results",paste0(id,"_integration.rds")))
}
destination<-file.path(out,mode,id);dir.create(destination,recursive=TRUE,showWarnings=FALSE)
files<-c(file.path(src,c("PROTOCOL.md","quadrature.R","model.R","run_case.R","controller.py")),
  file.path(cpp,c("model.R","observation_model.cpp","positive_taylor_forward.h")),
  file.path(root,"results/direct_source_v7/direct_joint.R"),file.path(root,"results/computation_source_v4/score.R"),
  file.path(repo,"validation/baseline_quadrature/quadrature.R"),file.path(repo,"validation/accuracy_round2/common.R"),
  runtime$design_binary,runtime$native,data_files,list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE))
sha<-round2_sha256(files);key<-round2_hash_object(sha)
write.csv(data.frame(file=names(sha),sha256=unname(sha)),file.path(destination,"source_sha256.csv"),row.names=FALSE)
started<-proc.time()[["elapsed"]]
fit<-tryCatch(adaptive_design_integrate(case,context,runtime,trace=TRUE),
  error=function(e)list(status="integration_error",message=conditionMessage(e)))
fit$elapsed_seconds<-proc.time()[["elapsed"]]-started
saveRDS(list(case_id=id,key=key,fit=fit),file.path(destination,"posterior.rds"))
if(!is.null(fit$history))write.csv(fit$history,file.path(destination,"refinement.csv"),row.names=FALSE)
if(!is.null(fit$tail_history))write.csv(fit$tail_history,file.path(destination,"tail_checks.csv"),row.names=FALSE)
if(!is.null(fit$conditional_status))write.csv(fit$conditional_status,file.path(destination,"nodes.csv"),row.names=FALSE)
validation_pass<-identical(fit$status,"passed_quadrature_refinement")
if(mode=="reference"&&validation_pass) {
  ref<-readRDS(data_files[2]);old_fit<-readRDS(data_files[3]);stopifnot(ref$reference$convergence$passed,old_fit$status=="passed_quadrature_refinement")
  n<-context$n;A<-diag(n)-matrix(case$reference_weight,n,n,byrow=TRUE)
  Cref<-ref$reference$covariance$relative
  mean_change<-max(abs(A%*%(fit$mean-old_fit$mean))/sqrt(diag(Cref)))
  oldC<-A%*%old_fit$covariance%*%t(A);newC<-A%*%fit$covariance%*%t(A)
  sd_change<-max(abs(sqrt(diag(newC))/sqrt(diag(oldC))-1))
  cov_change<-norm(newC-oldC,"F")/norm(oldC,"F")
  new_draw<-A%*%design_draws(fit,160000L,825001L+context$n)
  old_draw<-A%*%design_draws(old_fit,160000L,825001L+context$n)
  probs<-c(.025,.1,.25,.5,.75,.9,.975);cdf_rows<-list()
  for(p in probs) {
    q<-apply(old_draw,1,quantile,p)
    actual<-rowMeans(new_draw<=q)
    tolerance<-.01+3*sqrt(2*p*(1-p)/ncol(new_draw))
    cdf_rows[[as.character(p)]]<-data.frame(node=seq_len(n),probability=p,new_cdf=actual,
      difference=actual-p,tolerance=tolerance,passed=abs(actual-p)<=tolerance)
  }
  cdf<-do.call(rbind,cdf_rows);write.csv(cdf,file.path(destination,"interval_cdf_compatibility.csv"),row.names=FALSE)
  rsummary<-ref$reference$summary;rr<-rsummary[match(paste0("relative_",seq_len(n)),rsummary$variable),]
  comparison<-data.frame(case_id=id,mean_change_reference_sd=mean_change,relative_sd_change=sd_change,
    covariance_change=cov_change,max_cdf_difference=max(abs(cdf$difference)),
    new_mean_error_reference_sd=max(abs(as.numeric(A%*%fit$mean)-rr$mean)/rr$sd),
    old_mean_error_reference_sd=max(abs(as.numeric(A%*%old_fit$mean)-rr$mean)/rr$sd),
    min_sd_ratio_reference=min(sqrt(diag(newC))/rr$sd),max_sd_ratio_reference=max(sqrt(diag(newC))/rr$sd),
    passed=mean_change<=.01&&sd_change<=.01&&cov_change<=.01&&all(cdf$passed))
  write.csv(comparison,file.path(destination,"reference_comparison.csv"),row.names=FALSE)
  validation_pass<-comparison$passed;print(comparison,row.names=FALSE)
}
if(mode=="primary") {
  scored<-design_score(case,fit)
  saveRDS(list(case_id=id,key=key,fit=fit,scored=scored),file.path(destination,paste0(id,"_result.rds")))
  write.csv(scored$metrics,file.path(destination,"case_metrics.csv"),row.names=FALSE)
}
status<-data.frame(case_id=id,mode=mode,status=fit$status,validation_pass=validation_pass,
  elapsed_seconds=fit$elapsed_seconds,nodes=if(!is.null(fit$evaluated))nrow(fit$evaluated)else NA_integer_,
  message=if(is.null(fit$message))""else fit$message)
write.csv(status,file.path(destination,"status.csv"),row.names=FALSE)
cat("Adaptive result",id,fit$status,"validated",validation_pass,"seconds",fit$elapsed_seconds,"\n")
