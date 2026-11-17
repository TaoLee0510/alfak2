#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));cases<-readRDS(file.path(root,"results/mechanism_v2/frozen_cases.rds"))
folder<-Sys.getenv("ALFAK2_IS_SMALL_FOLDER","small_cases")
source(file.path(root,"results/mechanism_reference_source_v2/diagnostics.R"))
reportdir<-file.path(out,paste0(folder,"_analysis"));dir.create(reportdir,showWarnings=FALSE)
rows<-gates<-covs<-intervals<-identities<-list()
for(case in cases) {
 dest<-file.path(out,folder,case$id);fit<-readRDS(file.path(dest,"result.rds"))
 ref<-readRDS(file.path(out,"references",case$id,"reference.rds"));d<-ref$reference
 resource<-read.csv(file.path(dest,"resource.csv"));m<-fit$moments
 r<-d$summary[match(paste0(m$type,"_",m$index),d$summary$variable),];stopifnot(!anyNA(r$variable))
 se<-sqrt(m$mcse_mean^2+r$mcse_mean^2);sse<-sqrt(m$mcse_sd^2+r$mcse_sd^2)
 dm<-m$mean-r$mean;ds<-m$sd-r$sd
 z<-data.frame(case_id=case$id,type=m$type,index=m$index,reference_mean=r$mean,importance_mean=m$mean,
 reference_sd=r$sd,importance_sd=m$sd,joint_mean_mcse=se,joint_sd_mcse=sse,
 standardized_mean_difference=dm/r$sd,sd_ratio=m$sd/r$sd,
 material_mean=abs(dm)>.05*r$sd&abs(dm)>4*se,material_sd=abs(ds)>.1*r$sd&abs(ds)>4*sse,
 mean_within_tolerance=abs(dm)+4*se<=.05*r$sd,sd_within_tolerance=abs(ds)+4*sse<=.1*r$sd)
 rows[[case$id]]<-z
 obs<-lapply(ref$chains,function(ch)reference_observables(ch$draws,case))
 for(type in c("log_birth","birth","relative","centered","edge")) {
  C<-fit$covariances[[type]];RC<-d$covariance[[type]]
  covs[[paste(case$id,type)]]<-data.frame(case_id=case$id,type=type,
   relative_frobenius_error=norm(C-RC,"F")/norm(RC,"F"),
   reference_covariance_relative_mcse=d$covariance_relative_mcse[[type]])
  if(type=="log_birth")next
  for(level in c(.5,.8,.95)) {
   cols<-switch(as.character(level),"0.5"=c(3L,5L),"0.8"=c(2L,6L),"0.95"=c(1L,7L))
   q<-fit$quantiles[[type]][,cols,drop=FALSE]
   for(j in seq_len(nrow(q))) {
    v<-sapply(obs,function(o)o[,paste0(type,"_",j)]);e<-1*(v>=q[j,1]&v<=q[j,2])
    qq<-quantile(as.numeric(v),c((1-level)/2,(1+level)/2))
    intervals[[length(intervals)+1L]]<-data.frame(case_id=case$id,type=type,index=j,level=level,
     reference_low=qq[1],reference_high=qq[2],importance_low=q[j,1],importance_high=q[j,2],
     reference_mass=mean(e),reference_mass_mcse=posterior::mcse_mean(e),
     interpretation="fixed-data posterior mass; not repeated-sampling coverage")
   }
  }
 }
 # Integrate a | g analytically under the unchanged original Gaussian prior.
 n<-case$n;s2<-.25;v<-solve(case$covariance,rep(1,n));conditional_variance<-s2-s2*s2*sum(v)
 g<-sweep(fit$samples[,1:n,drop=FALSE]%*%t(case$specification$lower),2,case$specification$offset,"+")
 mm<-as.numeric((sweep(g,2,case$specification$offset,"-"))%*%v)*s2
 refm<-sapply(ref$chains,function(ch){gg<-sweep(ch$draws[,1:n,drop=FALSE]%*%t(case$specification$lower),2,case$specification$offset,"+");as.numeric(sweep(gg,2,case$specification$offset,"-")%*%v)*s2})
 mu<-sum(fit$weights*mm);se<-sqrt(sum(fit$weights^2*(mm-mu)^2))
 identities[[length(identities)+1L]]<-data.frame(case_id=case$id,quantity="baseline_mean",point=NA,
  importance=mu,reference=mean(refm),importance_mcse=se,reference_mcse=posterior::mcse_mean(refm))
 for(a in c(-1,-.5,0,.5,1)) {
  value<-dnorm(a,mm,sqrt(conditional_variance));rm<-matrix(dnorm(a,as.numeric(refm),sqrt(conditional_variance)),nrow(refm))
  mu<-sum(fit$weights*value);se<-sqrt(sum(fit$weights^2*(value-mu)^2))
  identities[[length(identities)+1L]]<-data.frame(case_id=case$id,quantity="baseline_density",point=a,
   importance=mu,reference=mean(rm),importance_mcse=se,reference_mcse=posterior::mcse_mean(rm))
 }
 sel<-z$type%in%c("relative","centered","edge")
 resource_pass<-tolower(as.character(resource$resource_pass))=="true"
 gate<-ref$precision_pass&&ref$forward_pass&&all(fit$flags)&&resource_pass&&
  all(z$mean_within_tolerance[sel])&&all(z$sd_within_tolerance[sel])
 gates[[case$id]]<-data.frame(case_id=case$id,reference_precision_pass=ref$precision_pass,
 importance_certified=all(fit$flags),resource_pass=resource_pass,
 max_contrast_mean_difference=max(abs(z$standardized_mean_difference[sel])),
 min_contrast_sd_ratio=min(z$sd_ratio[sel]),max_contrast_sd_ratio=max(z$sd_ratio[sel]),
 material_errors=sum(z$material_mean[sel]|z$material_sd[sel]),
 inconclusive=sum(!z$mean_within_tolerance[sel]|!z$sd_within_tolerance[sel]),gate_pass=gate)
}
write.csv(do.call(rbind,rows),file.path(reportdir,"importance_reference_moments.csv"),row.names=FALSE)
write.csv(do.call(rbind,covs),file.path(reportdir,"covariance_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,intervals),file.path(reportdir,"interval_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,identities),file.path(reportdir,"baseline_identity.csv"),row.names=FALSE)
g<-do.call(rbind,gates);write.csv(g,file.path(reportdir,"importance_reference_gate.csv"),row.names=FALSE);print(g,row.names=FALSE)
writeLines(if(all(g$gate_pass))"SMALL_REFERENCE_GATE_PASSED"else"SMALL_REFERENCE_GATE_FAILED_OR_INCONCLUSIVE",file.path(reportdir,"small_gate.txt"))
