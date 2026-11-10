#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
source_root<-normalizePath(Sys.getenv("ALFAK2_DIRECT_BATCH_SOURCE"))
out<-normalizePath(Sys.getenv("ALFAK2_DIRECT_BATCH_OUTPUT"))
cpp_source<-Sys.getenv("ALFAK2_DIRECT_CPP_SOURCE",file.path(root,"results/direct_forward_source_v1"))
Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp_source)
source(file.path(cpp_source,"model.R"))
source(file.path(source_root,"direct_joint.R"));source(file.path(root,"results/computation_source_v4/score.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
runtime<-design_load(repo);invisible(TMB::openmp(1,DLL="observation_model"))
context<-readRDS(file.path(root,"results/context.rds"))
manifest<-read.csv(file.path(root,"results/dataset_manifest.csv"))
row<-manifest[manifest$case_id==Sys.getenv("ALFAK2_DIRECT_CASE"),];stopifnot(nrow(row)==1L)
case<-readRDS(row$file);stopifnot(unname(round2_sha256(row$file))==row$sha256)
destination<-file.path(out,"cases",case$id);dir.create(destination,recursive=TRUE,showWarnings=FALSE)
files<-c(file.path(source_root,c("direct_joint.R","run_direct.R","bounded_direct.py")),
  file.path(cpp_source,c("model.R","observation_model.cpp","positive_taylor_forward.h")),
  file.path(root,"results/computation_source_v4/score.R"),
  file.path(root,"results",c("context.rds","dataset_manifest.csv")),runtime$design_binary,runtime$native,
  file.path(repo,"validation/baseline_quadrature/quadrature.R"),
  file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"),
  file.path(repo,"validation/accuracy_round2/common.R"),
  list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE))
sha<-round2_sha256(files);source_key<-round2_hash_object(sha);key<-list(source=source_key,dataset=row$sha256)
write.csv(data.frame(file=names(sha),sha256=unname(sha)),file.path(destination,"source_sha256.csv"),row.names=FALSE)
cache<-file.path("/private/tmp/alfak2_observation_design_direct",substr(source_key,1,12));dir.create(cache,recursive=TRUE,showWarnings=FALSE)
cat("Direct fitting PID",Sys.getpid(),case$id,"\n")
started<-proc.time()[["elapsed"]]
fit<-tryCatch(direct_joint_integrate(case,context,runtime,cache,key),
  error=function(e)list(status="integration_error",message=conditionMessage(e)))
fit$elapsed_seconds<-proc.time()[["elapsed"]]-started
scored<-design_score(case,fit)
path<-file.path(destination,paste0(case$id,"_result.rds"))
saveRDS(list(case_id=case$id,key=key,fit=fit,scored=scored),path)
status<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,family=case$family,design=case$design,
  status=fit$status,raw_status=fit$status,elapsed_seconds=fit$elapsed_seconds,
  final_order=if(!is.null(fit$history))tail(fit$history$order,1)else NA_integer_,
  conditional_nodes=if(!is.null(fit$conditional_status))nrow(fit$conditional_status)else 0L,
  message=if(is.null(fit$message))""else fit$message,result_path=path,recovery="none")
write.csv(status,file.path(destination,"status.csv"),row.names=FALSE)
write.csv(scored$metrics,file.path(destination,"case_metrics.csv"),row.names=FALSE)
cat("Direct result",fit$status,"seconds",fit$elapsed_seconds,"\n")
