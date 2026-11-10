#!/usr/bin/env Rscript
# Post-hoc numerical localization, not a replacement promotion/accuracy gate.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src)
cases<-readRDS(file.path(out,"frozen_cases.rds"));context<-readRDS(file.path(out,"context.rds"))
gate<-read.csv(file.path(out,"reference_gate.csv"));stopifnot(all(gate$reference_certified))
summary<-density<-list();files<-character()
for(i in seq_along(cases)) {
  case<-cases[[i]];n<-case$n;dest<-file.path(out,case$id)
  fit<-readRDS(file.path(dest,"adaptive.rds"));ref<-readRDS(file.path(dest,"reference.rds"))
  col<-which(colnames(context$prior$design)=="intercept");s2<-context$prior$coef_sd[col]^2
  c<-rep(s2,n);h<-as.numeric(solve(case$covariance,c));v<-s2-sum(c*h);stopifnot(v>0)
  # a | g,Y,eta = a | g ~ N(c'K^-1(g-mu), s2-c'K^-1*c).
  direction<-as.numeric(t(case$specification$lower)%*%h)
  M<-do.call(cbind,lapply(ref$chains,function(ch)as.numeric(ch$draws[,seq_len(n),drop=FALSE]%*%direction)))
  refmean<-mean(M);refsd<-sqrt(var(as.numeric(M))+v);se<-posterior::mcse_mean(M)
  a<-vapply(fit$components,`[[`,numeric(1),"intercept");w<-fit$weights
  amean<-sum(w*a);asd<-sqrt(sum(w*(a-amean)^2))
  summary[[case$id]]<-data.frame(case_id=case$id,conditional_prior_variance=v,
    reference_mean=refmean,reference_sd=refsd,reference_mean_mcse=se,
    quadrature_mean=amean,quadrature_sd=asd,mean_error_reference_sd=(amean-refmean)/refsd,
    mean_error_mcse=(amean-refmean)/se,
    exact_identity_mean_from_approximate_g=sum(h*(fit$latent_mean-case$prior$offset)),
    exact_identity_sd_from_approximate_g=sqrt(v+as.numeric(crossprod(h,fit$latent_covariance%*%h))),
    reference_exp_a_mean=mean(exp(M+v/2)),reference_exp_a_mcse=posterior::mcse_mean(exp(M+v/2)),
    quadrature_exp_a_mean=sum(w*exp(a)))
  infer<-case[c("id","y","times","trajectory","reference_weight")]
  for(value in c(-1,-.5,0,.5,1)) {
    p<-dnorm(value,M,sqrt(v));dim(p)<-dim(M)
    component<-direct_joint_fit(infer,context,runtime,shift=value)
    approximation<-exp(-component$diagnostics$objective+dnorm(value,0,sqrt(s2),log=TRUE)-fit$log_evidence)
    truth<-mean(p);mcse<-posterior::mcse_mean(p)
    density[[length(density)+1L]]<-data.frame(case_id=case$id,intercept=value,reference_density=truth,
      reference_density_mcse=mcse,laplace_mixture_density=approximation,
      relative_density_difference=(approximation-truth)/truth,density_difference_mcse=(approximation-truth)/mcse,
      interpretation="post-hoc Rao-Blackwell conditional-prior identity; not a new gate")
  }
  files<-c(files,file.path(dest,c("adaptive.rds","reference.rds")))
}
write.csv(do.call(rbind,summary),file.path(out,"baseline_identity_moments.csv"),row.names=FALSE)
write.csv(do.call(rbind,density),file.path(out,"baseline_identity_density.csv"),row.names=FALSE)
own<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
ledger<-mechanism_source_manifest(repo,src,out,c(files,own,file.path(out,c("frozen_cases.rds","context.rds","reference_gate.csv"))),runtime)
write.csv(ledger,file.path(out,"baseline_identity_source_sha256.csv"),row.names=FALSE)
print(do.call(rbind,summary),row.names=FALSE)
