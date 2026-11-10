#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src)
test<-read.csv(file.path(out,"validation_tests.csv"));stopifnot(all(test$failed==0),!any(test$error),all(test$warning==0))
cases<-readRDS(file.path(out,"frozen_cases.rds"));i<-as.integer(Sys.getenv("ALFAK2_MECHANISM_INDEX"));case<-cases[[i]];case$proposal<-readRDS(file.path(out,paste0("proposal_",i,".rds")))
context<-readRDS(file.path(out,"context.rds"));dest<-file.path(out,case$id)
dir.create(dest,showWarnings=FALSE);stopifnot(!file.exists(file.path(dest,"COMPLETE")))
ledger<-mechanism_source_manifest(repo,src,out,file.path(out,c("frozen_cases.rds","context.rds","input_freeze_sha256.csv","validation_tests.csv",paste0("proposal_",i,".rds"))),runtime)
write.csv(ledger,file.path(dest,"source_sha256.csv"),row.names=FALSE)
infer<-case[c("id","y","times","trajectory","reference_weight")]
t0<-proc.time()[["elapsed"]];fit<-adaptive_design_integrate(infer,context,runtime,trace=TRUE)
fit$elapsed_seconds<-proc.time()[["elapsed"]]-t0;stopifnot(fit$status=="passed_quadrature_refinement",fit$elapsed_seconds<=300)
saveRDS(fit,file.path(dest,"adaptive.rds"));write.csv(fit$history,file.path(dest,"adaptive_refinement.csv"),row.names=FALSE)
chains<-NULL;retained<-0L;history<-list();started<-proc.time()[["elapsed"]]
for(length in c(20000L,40000L,80000L,160000L)) {
  chains<-reference_chains(case,runtime,length-retained,warmup=5000L,
    seed=9800000L+1000L*i+length,previous=chains,cores=1L)
  retained<-length;d<-reference_diagnostics(chains,case);history[[as.character(length)]]<-d$convergence
  saveRDS(list(case_id=case$id,chains=chains,reference=d,history=do.call(rbind,history)),file.path(dest,paste0("checkpoint_",length,".rds")))
  write.csv(do.call(rbind,history),file.path(dest,"convergence_history.csv"),row.names=FALSE)
  cat(case$id,"retained",length,"Rhat",d$convergence$max_rhat,"ESS bulk/tail",d$convergence$min_bulk_ess,d$convergence$min_tail_ess,
    "passed",d$convergence$passed,"\n");flush.console()
  if(d$convergence$passed)break
}
forward<-list();sp<-case$specification;n<-case$n
for(ch in seq_along(chains)) {
  z<-chains[[ch]]$draws;pick<-unique(round(seq(1,nrow(z),length.out=24L)))
  for(j in pick) {
    state<-z[j,];g<-as.numeric(sp$offset+sp$lower%*%state[1:n]);eta<-as.numeric(sp$eta_basis%*%state[n+seq_len(n-1L)])
    p<-exp(eta-max(eta));p<-p/sum(p)
    A<-sweep(2*t(sp$probability)-diag(n),2,exp(g),"*")
    mass<-runtime$independent$process_deterministic_mean(A,p,1.5);q<-mass/sum(mass)
    exact_ll<-sum(sp$y0*log(p+1e-14))+sum(sp$y1*log(q+1e-14))
    ref<-runtime$cpp$ref_likelihood_cpp(state,sp)
    prod<-direct_forward_r(p,exp(g),context$operator,1.5,15L,8L)
    prod_ll<-sum(sp$y0*log(p+1e-14))+sum(sp$y1*log(prod$pi1+1e-14))
    forward[[length(forward)+1L]]<-data.frame(chain=ch,draw=j,pi_error=max(abs(q-ref$pi1)),
      reference_ll_error=abs(exact_ll-ref$log_likelihood),production_ll_error=abs(prod_ll-ref$log_likelihood))
  }
}
forward<-do.call(rbind,forward);write.csv(forward,file.path(dest,"independent_forward_checks.csv"),row.names=FALSE)
forward_pass<-all(forward$pi_error<=1e-10&forward$reference_ll_error<=1e-7&forward$production_ll_error<=1e-5)
saveRDS(list(case_id=case$id,reference=d,chains=chains,elapsed_seconds=proc.time()[["elapsed"]]-started,
  forward_pass=forward_pass),file.path(dest,"reference.rds"))
write.csv(d$summary,file.path(dest,"reference_diagnostics.csv"),row.names=FALSE)
write.csv(d$product_summary,file.path(dest,"reference_product_diagnostics.csv"),row.names=FALSE)
status<-cbind(d$convergence,data.frame(forward_pass=forward_pass,adaptive_seconds=fit$elapsed_seconds,
  mcmc_seconds=proc.time()[["elapsed"]]-started))
write.csv(status,file.path(dest,"status.csv"),row.names=FALSE)
writeLines("All declared sampling stages attempted until convergence or cap; no cases dropped.",file.path(dest,"COMPLETE"))
print(status,row.names=FALSE)
