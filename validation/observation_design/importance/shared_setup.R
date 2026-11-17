#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));src<-normalizePath(Sys.getenv("ALFAK2_IS_SOURCE"))
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));dest<-tempfile("alfak2_importance_cold_compile_",tmpdir="/private/tmp")
dir.create(dest);files<-file.path(src,c("importance_model.cpp","positive_taylor_forward.h","importance.cpp"))
stopifnot(all(file.copy(files,dest)));start<-proc.time()[["elapsed"]]
TMB::compile(file.path(dest,"importance_model.cpp"),flags="-O2 -g0")
tmb<-proc.time()[["elapsed"]]-start;start<-proc.time()[["elapsed"]]
env<-new.env();Rcpp::sourceCpp(file.path(dest,"importance.cpp"),env=env,cacheDir=file.path(dest,"rcpp"),showOutput=FALSE)
core<-proc.time()[["elapsed"]]-start
write.csv(data.frame(setup="independent_cold_build_of_identical_source",tmb_compile_seconds=tmb,core_compile_seconds=core,
 source=src,temporary_directory=dest),file.path(out,"shared_cold_compile.csv"),row.names=FALSE)
files<-c(files,TMB::dynlib(file.path(dest,"importance_model")))
write.csv(data.frame(file=files,sha256=vapply(files,digest::digest,character(1),file=TRUE,algo="sha256")),file.path(out,"shared_cold_compile_sha256.csv"),row.names=FALSE)
