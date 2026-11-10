#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT"))
cpp<-file.path(root,"results/direct_forward_source_v2");adaptive<-file.path(root,"results/adaptive_source_v1")
Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp)
source(file.path(cpp,"model.R"));source(file.path(root,"results/direct_source_v7/direct_joint.R"))
source(file.path(adaptive,"quadrature.R"));source(file.path(adaptive,"model.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"));source(file.path(root,"results/computation_source_v4/score.R"))
runtime<-design_load(repo);invisible(TMB::openmp(1,DLL="observation_model"))
id<-Sys.getenv("ALFAK2_CONSTRAINED_CASE");manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
row<-manifest[manifest$case_id==id,];stopifnot(nrow(row)==1L,unname(round2_sha256(row$file))==row$sha256)
case<-readRDS(row$file);context<-readRDS(file.path(root,"results/context.rds"))
destination<-file.path(out,"fits",id);dir.create(destination,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(destination,"status.csv")))
files<-c(file.path(src,c("run_case.R","controller.py","launch.py")),
  file.path(cpp,c("model.R","observation_model.cpp","positive_taylor_forward.h")),
  file.path(adaptive,c("quadrature.R","model.R")),file.path(root,"results/direct_source_v7/direct_joint.R"),
  file.path(root,"results/computation_source_v4/score.R"),file.path(repo,"validation/baseline_quadrature/quadrature.R"),
  file.path(repo,"validation/accuracy_round2/common.R"),runtime$design_binary,runtime$native,
  file.path(root,"results/control_source_v1/bounded_run.py"),file.path(out,c("dataset_manifest.csv","queue.json","selection_freeze_sha256.csv")),
  row$file,file.path(root,"results/context.rds"),list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE))
sha<-round2_sha256(files);key<-round2_hash_object(sha)
write.csv(data.frame(file=names(sha),sha256=unname(sha)),file.path(destination,"source_sha256.csv"),row.names=FALSE)
infer_case<-case[c("id","y","times","trajectory","reference_weight")]
started<-proc.time()[["elapsed"]]
fit<-tryCatch(adaptive_design_integrate(infer_case,context,runtime,trace=TRUE),
  error=function(e)list(status="integration_error",message=conditionMessage(e)))
fit$elapsed_seconds<-proc.time()[["elapsed"]]-started
saveRDS(list(case_id=id,key=key,fit=fit),file.path(destination,"posterior.rds"))
filespec<-c(history="refinement",tail_history="tail_checks",conditional_status="nodes")
for(name in names(filespec)) {
  if(!is.null(fit[[name]]))write.csv(fit[[name]],file.path(destination,paste0(filespec[[name]],".csv")),row.names=FALSE)
}
scored<-design_score(case,fit)
saveRDS(list(case_id=id,key=list(source=key,dataset=row$sha256),fit=fit,scored=scored),file.path(destination,paste0(id,"_result.rds")))
write.csv(scored$metrics,file.path(destination,"case_metrics.csv"),row.names=FALSE)
status<-data.frame(case_id=id,status=fit$status,validation_pass=identical(fit$status,"passed_quadrature_refinement"),
  elapsed_seconds=fit$elapsed_seconds,nodes=if(!is.null(fit$evaluated))nrow(fit$evaluated)else NA_integer_,
  message=if(is.null(fit$message))""else fit$message)
write.csv(status,file.path(destination,"status.csv"),row.names=FALSE)
cat("Constrained design case",id,fit$status,"seconds",fit$elapsed_seconds,"\n")
