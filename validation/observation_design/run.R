#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
src<-Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
source(file.path(src,"model.R"));source(file.path(src,"score.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
output<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(repo,"validation/observation_design/results"))
part<-Sys.getenv("ALFAK2_DESIGN_PART","all")
destination<-file.path(output,paste0("fits_",part));dir.create(destination,recursive=TRUE,showWarnings=FALSE)
validation<-read.csv(file.path(output,"validation/validation_tests.csv"))
if(any(validation$failed>0|validation$error|validation$warning>0))stop("Joint model validation not passed")
runtime<-design_load(repo)
invisible(TMB::openmp(1,DLL="observation_model"))
context<-readRDS(file.path(output,"context.rds"))
manifest<-read.csv(file.path(output,"dataset_manifest.csv"))
files<-c(file.path(src,c("PROTOCOL.md","DESIGN_AMENDMENT.md","add_anchor.R","model.R","data.R","run.R","score.R","observation_model.cpp","positive_taylor_forward.h")),
  file.path(output,"context.rds"),file.path(output,"dataset_manifest.csv"),runtime$design_binary,runtime$native,
  file.path(repo,"validation/baseline_quadrature/quadrature.R"),
  file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"),
  list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE))
sha<-round2_sha256(files);source_key<-round2_hash_object(sha)
write.csv(data.frame(file=names(sha),sha256=unname(sha)),file.path(destination,"source_sha256.csv"),row.names=FALSE)
cache<-file.path("/private/tmp/alfak2_observation_design_fits",substr(source_key,1,12))
dir.create(cache,recursive=TRUE,showWarnings=FALSE)
only<-Sys.getenv("ALFAK2_DESIGN_ONLY","")
if(nzchar(only)) {
  manifest<-manifest[manifest$case_id==only,]
} else if(part %in% c("0","1")) {
  units<-unique(manifest$unit_id);selected<-units[seq_along(units)%%2==as.integer(part)]
  manifest<-manifest[manifest$unit_id%in%selected,]
}
if(!nrow(manifest))stop("No design cases selected")
manifest<-manifest[order(match(manifest$unit_id,unique(manifest$unit_id)),
  match(manifest$design,c("control300","depth600","time600","composition_anchor600","composition600"))),]
cat("Owned process PID",Sys.getpid(),"cases",nrow(manifest),"\n")
statuses<-metrics<-list()
for(i in seq_len(nrow(manifest))) {
  row<-manifest[i,]
  if(unname(round2_sha256(row$file))!=row$sha256)stop("Frozen dataset SHA mismatch")
  case<-readRDS(row$file)
  stopifnot(identical(as.character(context$graph$labels),as.character(case$target_labels)))
  key<-list(source=source_key,dataset=row$sha256)
  filename<-file.path(cache,paste0(case$id,"_complete.rds"))
  if(file.exists(filename)) {
    saved<-readRDS(filename)
    if(!identical(saved$key,key))stop("Complete design cache mismatch")
    fit<-saved$fit
  } else {
    started<-proc.time()[["elapsed"]]
    fit<-tryCatch({
      setTimeLimit(cpu=Inf,elapsed=300,transient=TRUE)
      answer<-design_integrate(case,context,runtime,cache=cache,key=key)
      setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE)
      answer
    },error=function(e){setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE)
      list(status="integration_error",message=conditionMessage(e))})
    fit$elapsed_seconds<-proc.time()[["elapsed"]]-started
    saveRDS(list(key=key,fit=fit),filename)
  }
  scored<-design_score(case,fit)
  saveRDS(list(case_id=case$id,key=key,fit=fit,scored=scored),
    file.path(destination,paste0(case$id,"_result.rds")))
  statuses[[case$id]]<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,
    family=case$family,design=case$design,status=fit$status,elapsed_seconds=fit$elapsed_seconds,
    final_order=if(!is.null(fit$history))tail(fit$history$order,1)else NA_integer_,
    conditional_nodes=if(!is.null(fit$conditional_status))nrow(fit$conditional_status)else 0L,
    message=if(is.null(fit$message))""else fit$message)
  metrics[[case$id]]<-scored$metrics
  write.csv(do.call(rbind,statuses),file.path(destination,"status.csv"),row.names=FALSE)
  write.csv(do.call(rbind,metrics),file.path(destination,"case_metrics.csv"),row.names=FALSE)
  cat(sprintf("%d/%d %s %s order=%s %.2fs\n",i,nrow(manifest),case$id,fit$status,
    statuses[[case$id]]$final_order,fit$elapsed_seconds))
}
