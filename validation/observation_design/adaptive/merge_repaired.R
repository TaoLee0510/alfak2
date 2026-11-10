#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design");r<-file.path(root,"results")
out<-file.path(r,"adaptive_complete_v1");dir.create(out,showWarnings=FALSE)
source(file.path(repo,"validation/accuracy_round2/common.R"))
compat<-read.csv(file.path(r,"adaptive_compatibility_v1/compatibility.csv"))
stopifnot(nrow(compat)==3L,all(compat$passed))
manifest<-read.csv(file.path(r,"dataset_manifest.csv"));old<-read.csv(file.path(r,"analysis_input/all_fit_status.csv"))
failed<-old$case_id[old$status!="passed_quadrature_refinement"]
new_status<-do.call(rbind,lapply(c("adaptive_primary_smoke_v1","adaptive_repairs_v1"),function(folder)
  do.call(rbind,lapply(0:1,function(part)read.csv(file.path(r,folder,paste0("primary_part_",part),"status.csv"))))))
stopifnot(nrow(new_status)==6L,!anyDuplicated(new_status$case_id),setequal(failed,new_status$case_id),
  all(new_status$status=="passed_quadrature_refinement"),all(tolower(new_status$resource_pass)=="true"))
file.copy(file.path(r,"dataset_manifest.csv"),file.path(out,"dataset_manifest.csv"),overwrite=FALSE)
units<-unique(manifest$unit_id)
index<-list()
for(part in 0:1) {
  dest<-file.path(out,paste0("fits_",part));dir.create(dest,showWarnings=FALSE)
  selected<-units[seq_along(units)%%2==part]
  ss<-old[old$unit_id%in%selected,];metrics<-list()
  for(i in seq_len(nrow(ss))) {
    id<-ss$case_id[i];mr<-manifest[match(id,manifest$case_id),]
    stopifnot(unname(round2_sha256(mr$file))==mr$sha256)
    raw<-as.list(ss[i,]);fresh<-id%in%failed
    if(fresh) {
      folder<-if(file.exists(file.path(r,"adaptive_primary_smoke_v1/primary",id,paste0(id,"_result.rds"))))
        "adaptive_primary_smoke_v1" else "adaptive_repairs_v1"
      path<-file.path(r,folder,"primary",id,paste0(id,"_result.rds"));result<-readRDS(path)
      source_file<-file.path(dirname(path),"source_sha256.csv");source<-read.csv(source_file)
      stopifnot(source$sha256[match(mr$file,source$file)]==mr$sha256)
      result$original_source_key<-result$key
      result$key<-list(source=result$key,dataset=mr$sha256)
      current<-new_status[match(id,new_status$case_id),]
      ss$status[i]<-result$fit$status;ss$elapsed_seconds[i]<-current$elapsed_seconds
      ss$resource_pass[i]<-current$resource_pass;ss$maximum_rss_bytes[i]<-current$maximum_rss_bytes
      ss$recovery[i]<-"adaptive_same_target_fresh_fit"
      ss$conditional_nodes[i]<-nrow(result$fit$conditional_status);ss$final_order[i]<-NA_integer_
      ss$message[i]<-"Original failure preserved; fresh adaptive result passed"
      method<-result$fit$computational_method
    } else {
      path<-ss$result_path[i];result<-readRDS(path)
      source_file<-file.path(r,paste0("fits_",part),"source_sha256.csv")
      method<-"legacy_converged_gauss_hermite"
    }
    stopifnot(result$fit$status=="passed_quadrature_refinement",result$key$dataset==mr$sha256)
    result$original_execution<-raw;result$recovery<-ss$recovery[i]
    dest_path<-file.path(dest,paste0(id,"_result.rds"));saveRDS(result,dest_path)
    ss$result_path[i]<-normalizePath(dest_path);ss$numerical_method[i]<-method
    ss$raw_status[i]<-raw$raw_status;ss$original_elapsed_seconds[i]<-raw$elapsed_seconds
    ss$memory_evidence[i]<-if(is.na(ss$maximum_rss_bytes[i])||!nzchar(as.character(ss$maximum_rss_bytes[i])))
      "legacy_smoke_only" else "per_case_wait4"
    metrics[[id]]<-result$scored$metrics
    index[[id]]<-data.frame(case_id=id,result_path=normalizePath(dest_path),source_manifest=normalizePath(source_file),
      numerical_method=method,original_status=raw$status,fresh_recovery=fresh)
  }
  write.csv(ss,file.path(dest,"status.csv"),row.names=FALSE)
  write.csv(do.call(rbind,metrics),file.path(dest,"case_metrics.csv"),row.names=FALSE)
}
write.csv(do.call(rbind,index),file.path(out,"result_index.csv"),row.names=FALSE)
cat("All sixty estimates restored under the unchanged statistical model; raw failures retained\n")
