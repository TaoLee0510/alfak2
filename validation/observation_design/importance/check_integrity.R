#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));files<-list.files(out,"source_sha256[.]csv$",recursive=TRUE,full.names=TRUE)
rows<-do.call(rbind,lapply(files,function(f){z<-read.csv(f);z$manifest<-f;z}))
paths<-unique(rows$file);hash<-vapply(paths,function(p)if(file.exists(p))digest::digest(file=p,algo="sha256")else NA_character_,character(1))
rows$current_sha256<-unname(hash[rows$file]);rows$unchanged<-!is.na(rows$current_sha256)&rows$sha256==rows$current_sha256
write.csv(rows,file.path(out,"consumed_source_integrity.csv"),row.names=FALSE)
stopifnot(all(rows$unchanged))
manifestfile<-file.path(out,"mechanism/manifest.csv")
if(file.exists(manifestfile)) {
 m<-read.csv(manifestfile);m$current_sha256<-vapply(m$file,digest::digest,character(1),algo="sha256",file=TRUE)
 m$current_scoring_sha256<-vapply(m$scoring_file,digest::digest,character(1),algo="sha256",file=TRUE)
 m$unchanged<-m$sha256==m$current_sha256&m$scoring_sha256==m$current_scoring_sha256
 write.csv(m,file.path(out,"mechanism/final_dataset_integrity.csv"),row.names=FALSE);stopifnot(all(m$unchanged))
}
cat(length(unique(rows$file)),"consumed files unchanged, including all failed attempts\n")
