#!/usr/bin/env Rscript
# Reconcile all attempted cases without refitting or deleting failed outcomes.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
root<-file.path(repo,"validation/observation_design");raw<-file.path(root,"results")
destination<-file.path(raw,"analysis_input");dir.create(destination,showWarnings=FALSE)
source(file.path(raw,"computation_source_v4/score.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
manifest<-read.csv(file.path(raw,"dataset_manifest.csv"))
statuses<-lapply(0:1,function(p)read.csv(file.path(raw,paste0("fits_bounded_",p),"status.csv")))
status<-do.call(rbind,statuses)
stopifnot(nrow(status)==60L,!anyDuplicated(status$case_id),setequal(status$case_id,manifest$case_id),
  all(file.exists(file.path(raw,paste0("fits_bounded_",0:1),"COMPLETE"))))
reference_path<-status$result_path[which(!is.na(status$result_path)&nzchar(status$result_path))[1]]
source_key<-readRDS(reference_path)$key$source
cache<-file.path("/private/tmp/alfak2_observation_design_fits",substr(source_key,1,12))
file.copy(file.path(raw,"dataset_manifest.csv"),file.path(destination,"dataset_manifest.csv"),overwrite=FALSE)
recovery_rows<-list()
for(part in 0:1) {
  target<-file.path(destination,paste0("fits_",part));dir.create(target,showWarnings=FALSE)
  file.copy(file.path(raw,paste0("fits_",part),"source_sha256.csv"),file.path(target,"source_sha256.csv"),overwrite=FALSE)
  ss<-statuses[[part+1L]];mm<-list()
  for(i in seq_len(nrow(ss))) {
    s<-ss[i,];mr<-manifest[match(s$case_id,manifest$case_id),]
    stopifnot(mr$sha256==unname(round2_sha256(mr$file)))
    case<-readRDS(mr$file);key<-list(source=source_key,dataset=mr$sha256)
    result<-NULL;recovery<-"none"
    if(!is.na(s$result_path)&&nzchar(s$result_path)&&file.exists(s$result_path))result<-readRDS(s$result_path)
    if(is.null(result)||result$fit$status!="passed_quadrature_refinement") {
      complete<-file.path(cache,paste0(s$case_id,"_complete.rds"))
      cached<-if(file.exists(complete))tryCatch(readRDS(complete),error=function(e)NULL)else NULL
      if(!is.null(cached)&&identical(cached$key,key)&&cached$fit$status=="passed_quadrature_refinement") {
        result<-list(case_id=s$case_id,key=key,fit=cached$fit,scored=NULL)
        recovery<-"completed_cache_after_external_cut"
      } else {
        descriptor<-list(case_id=s$case_id,key=key,fit=list(status=as.character(s$status),
          elapsed_seconds=s$elapsed_seconds,message=as.character(s$message)),synthetic_failure_descriptor=TRUE)
        input<-file.path(target,paste0(s$case_id,"_failure_descriptor.rds"));saveRDS(descriptor,input)
        old_case<-Sys.getenv("ALFAK2_DESIGN_RECONSTRUCT_CASE",unset=NA_character_)
        old_input<-Sys.getenv("ALFAK2_DESIGN_RECONSTRUCT_INPUT",unset=NA_character_)
        Sys.setenv(ALFAK2_DESIGN_RECONSTRUCT_CASE=s$case_id,ALFAK2_DESIGN_RECONSTRUCT_INPUT=input)
        log<-file.path(target,paste0(s$case_id,"_reconstruction.log"))
        rc<-suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),
          shQuote(file.path(root,"evidence_reconstruct.R")),stdout=log,stderr=log))
        if(is.na(old_case))Sys.unsetenv("ALFAK2_DESIGN_RECONSTRUCT_CASE")else Sys.setenv(ALFAK2_DESIGN_RECONSTRUCT_CASE=old_case)
        if(is.na(old_input))Sys.unsetenv("ALFAK2_DESIGN_RECONSTRUCT_INPUT")else Sys.setenv(ALFAK2_DESIGN_RECONSTRUCT_INPUT=old_input)
        repaired<-file.path(raw,"evidence_reconstruction",paste0(s$case_id,"_result.rds"))
        if(rc==0L&&file.exists(repaired)) {
          result<-readRDS(repaired);recovery<-"joint_mode_hessian_evidence_reconstruction"
        } else result<-descriptor
      }
    }
    stopifnot(identical(result$key,key))
    result$fit$elapsed_seconds<-s$elapsed_seconds
    result$fit$resource_pass<-tolower(as.character(s$resource_pass))=="true"
    result$original_execution<-as.list(s)
    result$recovery<-recovery
    result$scored<-design_score(case,result$fit)
    path<-file.path(target,paste0(s$case_id,"_result.rds"));saveRDS(result,path)
    ss$raw_status[i]<-as.character(s$status);ss$status[i]<-result$fit$status
    ss$raw_final_order[i]<-s$final_order
    ss$final_order[i]<-if(!is.null(result$fit$history))tail(result$fit$history$order,1)else NA_integer_
    ss$conditional_nodes[i]<-if(!is.null(result$fit$conditional_status))nrow(result$fit$conditional_status)else 0L
    ss$raw_result_path[i]<-as.character(s$result_path);ss$result_path[i]<-normalizePath(path)
    ss$recovery[i]<-recovery
    mm[[s$case_id]]<-result$scored$metrics
    recovery_rows[[s$case_id]]<-data.frame(case_id=s$case_id,design=s$design,raw_status=s$status,
      numerical_status=result$fit$status,recovery=recovery,elapsed_seconds=s$elapsed_seconds,
      resource_pass=result$fit$resource_pass)
  }
  write.csv(ss,file.path(target,"status.csv"),row.names=FALSE)
  write.csv(do.call(rbind,mm),file.path(target,"case_metrics.csv"),row.names=FALSE)
}
recovery<-do.call(rbind,recovery_rows)
write.csv(recovery,file.path(destination,"recovery_status.csv"),row.names=FALSE)
cat("All 60 attempts preserved; numerically available",sum(recovery$numerical_status=="passed_quadrature_refinement"),
  "; resource failures",sum(!recovery$resource_pass),"; no additional fitting\n")
