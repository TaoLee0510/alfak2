#!/usr/bin/env Rscript
# Input/source integrity is checked even if statistical availability fails.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
root<-file.path(repo,"validation/observation_design");out<-file.path(root,"results")
source(file.path(repo,"validation/accuracy_round2/common.R"))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
manifest$current_sha256<-unname(round2_sha256(manifest$file))
manifest$unchanged<-manifest$sha256==manifest$current_sha256
write.csv(manifest,file.path(out,"final_dataset_integrity.csv"),row.names=FALSE)
stopifnot(nrow(manifest)==60L,all(manifest$unchanged))
consumed<-c(file.path(out,paste0("fits_",0:1),"source_sha256.csv"),
  file.path(list.files(out,"^fits_bounded_child_",full.names=TRUE),"source_sha256.csv"))
consumed<-consumed[file.exists(consumed)]
rows<-do.call(rbind,lapply(consumed,function(p){z<-read.csv(p);z$manifest_path<-p;z}))
hash<-round2_sha256(unique(rows$file));rows$current_sha256<-unname(hash[rows$file])
rows$unchanged<-rows$sha256==rows$current_sha256
write.csv(rows,file.path(out,"final_consumed_source_integrity.csv"),row.names=FALSE)
stopifnot(all(rows$unchanged))
frozen<-list(list(csv="frozen_source_v4_sha256.csv",directory="computation_source_v4"),
  list(csv="analysis_source_v1_sha256.csv",directory=NULL),
  list(csv="control_source_v1_sha256.csv",directory=NULL),
  list(csv="analysis_source_v2_sha256.csv",directory=NULL))
extra<-list.files(out,"^(direct_source_v[0-9]+|direct_forward_source_v[0-9]+|direct_batch_source_v[0-9]+|checkpoint_source_v[0-9]+|depth_diagnostic_source_v[0-9]+|adaptive_source_v[0-9]+|adaptive_compatibility_source_v[0-9]+|constrained(_[a-z]+)?_source_v[0-9]+|mechanism(_[a-z]+)?_source_v[0-9]+|importance(_[a-z]+)?_source_v[0-9]+)_sha256[.]csv$")
extra<-c(extra,list.files(out,"source_sha256[.]csv$",recursive=TRUE)[grepl("adaptive_postprocess_source_v",list.files(out,"source_sha256[.]csv$",recursive=TRUE))])
frozen<-c(frozen,lapply(extra,function(p)list(csv=p,directory=NULL)))
checks<-list()
for(spec in frozen) {
  z<-read.csv(file.path(out,spec$csv))
  paths<-if(is.null(spec$directory))z$file else file.path(out,spec$directory,z$file)
  z$current_sha256<-unname(round2_sha256(paths));z$unchanged<-z$sha256==z$current_sha256
  z$manifest<-spec$csv;checks[[spec$csv]]<-z[,c("file","sha256","current_sha256","unchanged","manifest")]
}
checks<-do.call(rbind,checks)
write.csv(checks,file.path(out,"final_frozen_source_integrity.csv"),row.names=FALSE)
stopifnot(all(checks$unchanged))
original<-jsonlite::read_json(file.path(repo,"validation/accuracy_upgrade/results/initial_source_manifest.json"),simplifyVector=TRUE)
names<-c("R/benchmark-second-layer.R","doc/alfak2_method.tex")
user<-data.frame(file=names,original_sha256=as.character(original[names]),
  current_sha256=unname(round2_sha256(file.path(repo,names))))
user$unchanged<-user$original_sha256==user$current_sha256
write.csv(user,file.path(out,"user_original_file_integrity.csv"),row.names=FALSE)
stopifnot(all(user$unchanged))
cat("60 data hashes, all consumed primary sources, frozen stages and the two user-original files unchanged\n")
