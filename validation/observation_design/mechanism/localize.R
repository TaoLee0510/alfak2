#!/usr/bin/env Rscript
# Run only after a certified primary reference finds material approximation error.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src)
i<-as.integer(Sys.getenv("ALFAK2_MECHANISM_INDEX"));cases<-readRDS(file.path(out,"frozen_cases.rds"));case<-cases[[i]]
gate<-read.csv(file.path(out,"reference_gate.csv"));stopifnot(all(gate$reference_certified),any(gate$material_contrast_error))
context<-readRDS(file.path(out,"context.rds"));infer<-case[c("id","y","times","trajectory","reference_weight")]
dest<-file.path(out,"localization",case$id);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(dest,"COMPLETE")))
original<-readRDS(file.path(out,case$id,"adaptive.rds"));primary<-readRDS(file.path(out,case$id,"reference.rds"))
orig_integrator<-adaptive_prior_integrate
adaptive_prior_integrate<-local({fun<-orig_integrator;function(evaluate,prior_sd,weight,trace=FALSE,...)
  fun(evaluate,prior_sd,weight,tolerances=c(normalizer=1e-5,mean=2e-5,covariance=1e-4,contrast=1e-4),
    max_panels=64L,max_evaluations=2000L,trace=trace)})
fine<-adaptive_design_integrate(infer,context,runtime,trace=TRUE)
adaptive_prior_integrate<-orig_integrator
saveRDS(fine,file.path(dest,"tight_quadrature.rds"));stopifnot(fine$status=="passed_quadrature_refinement")
n<-case$n;A<-diag(n)-matrix(case$weight,n,n,byrow=TRUE)
refrows<-primary$reference$summary[match(paste0("relative_",1:n),primary$reference$summary$variable),]
quad<-data.frame(case_id=case$id,mean_change_reference_sd=max(abs(as.numeric(A%*%(fine$mean-original$mean)))/refrows$sd),
  covariance_relative_change=norm(A%*%(fine$covariance-original$covariance)%*%t(A),"F")/norm(A%*%fine$covariance%*%t(A),"F"),
  default_nodes=nrow(original$evaluated),tight_nodes=nrow(fine$evaluated),
  default_mean_error_reference_sd=max(abs(as.numeric(A%*%original$mean)-refrows$mean)/refrows$sd),
  tight_mean_error_reference_sd=max(abs(as.numeric(A%*%fine$mean)-refrows$mean)/refrows$sd))
write.csv(quad,file.path(dest,"quadrature_tolerance_check.csv"),row.names=FALSE)
# Condition on the original independent intercept exactly at a0=0. The mean-one
# offset remains unchanged; only its independent variance contribution is removed.
column<-which(colnames(context$prior$design)=="intercept");stopifnot(length(column)==1L)
s2<-context$prior$coef_sd[column]^2
conditional<-case;conditional$id<-paste0(case$id,"_conditional_a0")
conditional$covariance<-case$covariance-matrix(s2,n,n)
conditional$specification$lower<-t(chol(conditional$covariance))
mode<-direct_joint_fit(infer,context,runtime,shift=0)
stopifnot(mode$diagnostics$covariance_status=="trusted_direct_joint_hessian")
saveRDS(list(case=conditional,intercept=0,original_intercept_variance=s2,laplace=mode),file.path(dest,"conditional_target.rds"))
old_helper<-new.env(parent=globalenv())
old_helper_file<-file.path(root,"results/mechanism_reference_source_v1/diagnostics.R")
sys.source(old_helper_file,old_helper)
pilot<-old_helper$reference_chains(conditional,runtime,10000L,warmup=1000L,seed=9710000L+i*1000L,cores=1L)
z<-do.call(rbind,lapply(pilot,function(ch)tail(ch$draws,5000L)));S<-2*cov(z)+diag(1e-4,ncol(z))
conditional$proposal<-list(mean=colMeans(z),covariance=S,lower=t(chol(S)),precision=solve(S))
saveRDS(conditional$proposal,file.path(dest,"conditional_proposal.rds"))
chains<-NULL;retained<-0L;history<-list()
for(length in c(20000L,40000L,80000L,160000L)) {
  chains<-reference_chains(conditional,runtime,length-retained,warmup=5000L,
    seed=9700000L+i*1000L+length,previous=chains,cores=1L)
  retained<-length;d<-reference_diagnostics(chains,conditional);history[[as.character(length)]]<-d$convergence
  write.csv(do.call(rbind,history),file.path(dest,"conditional_convergence.csv"),row.names=FALSE)
  cat(conditional$id,"retained",length,"passed",d$convergence$passed,"Rhat",d$convergence$max_rhat,
    "ESS",d$convergence$min_bulk_ess,d$convergence$min_tail_ess,"\n");flush.console()
  if(d$convergence$passed)break
}
saveRDS(list(case_id=conditional$id,reference=d,chains=chains),file.path(dest,"conditional_reference.rds"))
E<-matrix(0,nrow(case$edges),n);E[cbind(seq_len(nrow(E)),case$edges[,1])]<- -1;E[cbind(seq_len(nrow(E)),case$edges[,2])]<-1
transforms<-list(relative=A,centered=diag(n)-matrix(1/n,n,n),edge=E);rows<-list()
for(type in names(transforms)) {
  T<-transforms[[type]];m<-as.numeric(T%*%mode$summary$fitness_mean)
  sd<-sqrt(diag(T%*%mode$fitness_covariance%*%t(T)))
  r<-d$summary[match(paste0(type,"_",seq_along(m)),d$summary$variable),]
  rows[[type]]<-data.frame(case_id=case$id,type=type,index=seq_along(m),reference_mean=r$mean,laplace_mean=m,
    reference_sd=r$sd,laplace_sd=sd,mcse_mean=r$mcse_mean,mcse_sd=r$mcse_sd,
    standardized_mean_difference=(m-r$mean)/r$sd,sd_ratio=sd/r$sd,
    material_mean=abs(m-r$mean)>.05*r$sd & abs(m-r$mean)>4*r$mcse_mean,
    material_sd=abs(sd-r$sd)>.1*r$sd & abs(sd-r$sd)>4*r$mcse_sd,
    reference_certified=d$convergence$passed)
}
write.csv(do.call(rbind,rows),file.path(dest,"conditional_moment_comparison.csv"),row.names=FALSE)
own<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
files<-c(file.path(out,c("reference_gate.csv","frozen_cases.rds","context.rds")),
  file.path(out,case$id,c("adaptive.rds","reference.rds")),own,old_helper_file)
ledger<-mechanism_source_manifest(repo,src,out,files,runtime)
write.csv(ledger,file.path(dest,"source_sha256.csv"),row.names=FALSE)
writeLines("Single fixed-intercept localization complete; does not certify all conditional nodes.",file.path(dest,"COMPLETE"))
print(quad,row.names=FALSE)
print(aggregate(cbind(material_mean,material_sd)~type,do.call(rbind,rows),sum),row.names=FALSE)
