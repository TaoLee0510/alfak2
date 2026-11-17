#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));i<-as.integer(Sys.getenv("ALFAK2_IS_INDEX"))
src<-file.path(root,"results/mechanism_reference_source_v2")
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src)
cases<-readRDS(file.path(root,"results/mechanism_v2/frozen_cases.rds"));case<-cases[[i]]
case$proposal<-readRDS(file.path(root,"results/mechanism_v2",paste0("proposal_",i,".rds")))
old<-readRDS(file.path(root,"results/mechanism_v2",case$id,"reference.rds"))
chains<-old$chains;retained<-nrow(chains[[1]]$draws)
dest<-file.path(out,"references",case$id);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(dest,"COMPLETE")))
history<-list();started<-proc.time()[["elapsed"]]
for(length in c(80000L,160000L)) {
  chains<-reference_chains(case,runtime,length-retained,warmup=0L,seed=10100000L+i*1000L+length,previous=chains,cores=1L)
  retained<-length;d<-reference_diagnostics(chains,case)
  precision<-d$convergence$passed&&d$convergence$max_mean_mcse_sd<=.005&&d$convergence$max_sd_mcse_sd<=.01
  history[[as.character(length)]]<-cbind(d$convergence,precision_pass=precision)
  saveRDS(list(case_id=case$id,chains=chains,reference=d),file.path(dest,paste0("checkpoint_",length,".rds")))
  write.csv(do.call(rbind,history),file.path(dest,"convergence_history.csv"),row.names=FALSE)
  cat(case$id,"retained",length,"mean/SD MCSE",d$convergence$max_mean_mcse_sd,d$convergence$max_sd_mcse_sd,"precision",precision,"\n");flush.console()
  if(precision)break
}
saveRDS(list(case_id=case$id,chains=chains,reference=d,precision_pass=precision,
  forward_pass=old$forward_pass,elapsed_seconds=proc.time()[["elapsed"]]-started),file.path(dest,"reference.rds"))
write.csv(d$summary,file.path(dest,"reference_diagnostics.csv"),row.names=FALSE)
files<-c(file.path(root,"results/mechanism_v2",case$id,"reference.rds"),
  file.path(root,"results/mechanism_v2",paste0("proposal_",i,".rds")),
  file.path(root,"results/mechanism_v2/frozen_cases.rds"),
  normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
write.csv(mechanism_source_manifest(repo,src,out,files,runtime),file.path(dest,"source_sha256.csv"),row.names=FALSE)
writeLines(if(precision)"High-precision reference passed"else"High-precision reference uncertified at cap",file.path(dest,"COMPLETE"))
