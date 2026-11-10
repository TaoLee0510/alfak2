#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src,FALSE)
cases<-readRDS(file.path(out,"frozen_cases.rds"));rows<-covariances<-intervals<-statuses<-list()
source_files<-character()
for(i in seq_along(cases)) {
  case<-cases[[i]];dest<-file.path(out,case$id);n<-case$n
  resource<-read.csv(file.path(out,paste0("controller_",i,"_resource.csv")))
  memory_pass<-tolower(as.character(resource$memory_pass))=="true"
  stopifnot(resource$exit_code==0,memory_pass)
  fit<-readRDS(file.path(dest,"adaptive.rds"));ref<-readRDS(file.path(dest,"reference.rds"));d<-ref$reference
  stopifnot(fit$status=="passed_quadrature_refinement")
  E<-matrix(0,nrow(case$edges),n);E[cbind(seq_len(nrow(E)),case$edges[,1])]<- -1
  E[cbind(seq_len(nrow(E)),case$edges[,2])]<-1
  transform<-list(log_birth=diag(n),birth=diag(n),relative=diag(n)-matrix(case$weight,n,n,byrow=TRUE),
    centered=diag(n)-matrix(1/n,n,n),edge=E)
  obs<-lapply(ref$chains,function(ch)reference_observables(ch$draws,case))
  draws<-design_draws(fit,160000L,9600000L+i)
  material<-FALSE
  for(type in names(transform)) {
    A<-transform[[type]];m<-if(type=="log_birth")fit$latent_mean else as.numeric(A%*%fit$mean)
    C<-if(type=="log_birth")fit$latent_covariance else A%*%fit$covariance%*%t(A)
    names<-paste0(type,"_",seq_len(nrow(A)));r<-d$summary[match(names,d$summary$variable),]
    sd<-sqrt(pmax(0,diag(C)));dm<-m-r$mean;ds<-sd-r$sd
    z<-data.frame(case_id=case$id,type=type,index=seq_along(m),reference_mean=r$mean,approximate_mean=m,
      reference_sd=r$sd,approximate_sd=sd,mcse_mean=r$mcse_mean,mcse_sd=r$mcse_sd,
      standardized_mean_difference=dm/r$sd,sd_ratio=sd/r$sd,
      mean_difference_in_mcse=dm/r$mcse_mean,sd_difference_in_mcse=ds/r$mcse_sd,
      material_mean=abs(dm)>.05*r$sd & abs(dm)>4*r$mcse_mean,
      material_sd=abs(ds)>.1*r$sd & abs(ds)>4*r$mcse_sd,
      reference_certified=d$convergence$passed & ref$forward_pass)
    rows[[paste(case$id,type)]]<-z
    if(type%in%c("relative","centered","edge"))material<-material||any(z$material_mean|z$material_sd)
    RC<-d$covariance[[type]]
    covariances[[paste(case$id,type)]]<-data.frame(case_id=case$id,type=type,
      relative_frobenius_error=norm(C-RC,"F")/norm(RC,"F"),
      reference_covariance_relative_mcse=d$covariance_relative_mcse[[type]])
    if(type=="log_birth")next
    approx<-A%*%draws
    for(level in c(.5,.8,.95)) {
      q<-t(apply(approx,1,quantile,c((1-level)/2,(1+level)/2)))
      for(j in seq_len(nrow(A))) {
        v<-do.call(cbind,lapply(obs,function(o)o[,names[j]]));event<-v>=q[j,1]&v<=q[j,2]
        qq<-quantile(as.numeric(v),c((1-level)/2,(1+level)/2));qse<-posterior::mcse_quantile(v,probs=c((1-level)/2,(1+level)/2))
        intervals[[length(intervals)+1L]]<-data.frame(case_id=case$id,type=type,index=j,level=level,
          reference_low=qq[1],reference_high=qq[2],reference_low_mcse=qse[1],reference_high_mcse=qse[2],
          approximate_low=q[j,1],approximate_high=q[j,2],reference_mass=mean(event),mass_mcse=posterior::mcse_mean(event*1),
          interpretation="fixed-data posterior mass; not repeated-sampling coverage")
      }
    }
  }
  certified<-isTRUE(d$convergence$passed)&&isTRUE(ref$forward_pass)&&memory_pass
  statuses[[case$id]]<-data.frame(case_id=case$id,reference_certified=certified,material_contrast_error=material,
    gate_pass=certified&&!material,retained_per_chain=d$convergence$retained_per_chain,
    max_rhat=d$convergence$max_rhat,min_bulk_ess=d$convergence$min_bulk_ess,min_tail_ess=d$convergence$min_tail_ess,
    max_mean_mcse_sd=d$convergence$max_mean_mcse_sd,max_sd_mcse_sd=d$convergence$max_sd_mcse_sd,
    adaptive_seconds=fit$elapsed_seconds,reference_seconds=ref$elapsed_seconds,
    maximum_rss_gib=resource$maximum_rss_bytes/1024^3)
  source_files<-c(source_files,file.path(dest,"source_sha256.csv"))
}
status<-do.call(rbind,statuses);moments<-do.call(rbind,rows)
write.csv(status,file.path(out,"reference_gate.csv"),row.names=FALSE)
write.csv(moments,file.path(out,"moment_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,covariances),file.path(out,"covariance_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,intervals),file.path(out,"interval_comparison.csv"),row.names=FALSE)
ledger<-do.call(rbind,lapply(c(source_files,file.path(out,"input_freeze_sha256.csv")),function(p){z<-read.csv(p);z$manifest_path<-p;z}))
hash<-round2_sha256(unique(ledger$file));ledger$current_sha256<-unname(hash[ledger$file])
stopifnot(all(ledger$sha256==ledger$current_sha256))
write.csv(ledger,file.path(out,"consumed_source_integrity.csv"),row.names=FALSE)
writeLines(if(all(status$gate_pass))"REFERENCE_GATE_PASSED_REVIEW_DISTRIBUTION_BEFORE_ATTRIBUTION"else
  "REFERENCE_GATE_FAILED_ATTRIBUTION_BATCH_NOT_STARTED",file.path(out,"stage_status.txt"))
print(status,row.names=FALSE)
print(aggregate(cbind(material_mean,material_sd)~case_id+type,moments,sum),row.names=FALSE)
