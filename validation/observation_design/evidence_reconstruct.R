#!/usr/bin/env Rscript
# Deterministic Laplace-evidence verification at saved modes; no optimization.
# Original worker failures and caches are never overwritten.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
root<-file.path(repo,"validation/observation_design")
out<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(root,"results"))
Sys.setenv(ALFAK2_DESIGN_SOURCE=file.path(out,"computation_source_v4"))
source(file.path(root,"model.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
case_id<-Sys.getenv("ALFAK2_DESIGN_RECONSTRUCT_CASE","design_rough_L01_R1_composition_anchor600")
input<-Sys.getenv("ALFAK2_DESIGN_RECONSTRUCT_INPUT","")
paths<-if(nzchar(input))input else file.path(out,paste0("fits_",0:1),paste0(case_id,"_result.rds"))
paths<-paths[file.exists(paths)];stopifnot(length(paths)==1L)
original<-readRDS(paths);stopifnot(original$fit$status!="passed_quadrature_refinement")
manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
row<-manifest[manifest$case_id==case_id,];stopifnot(nrow(row)==1L)
stopifnot(row$sha256==unname(round2_sha256(row$file)),identical(original$key$dataset,row$sha256))
case<-readRDS(row$file);context<-readRDS(file.path(out,"context.rds"))
cache<-file.path("/private/tmp/alfak2_observation_design_fits",substr(original$key$source,1,12))
runtime<-design_load(repo);invisible(TMB::openmp(1,DLL="observation_model"))
destination<-file.path(out,"evidence_reconstruction");dir.create(destination,showWarnings=FALSE)
started<-proc.time()[["elapsed"]]
# Importance proposal can be recovered exactly from its stored affine nodes;
# the proposal is not a model parameter and does not need a new fit.
rule<-oracle_normal_quadrature(9L)
nodes<-lapply(seq_along(rule$nodes),function(j)readRDS(file.path(cache,sprintf("%s_order9_node%d.rds",case_id,j))))
values<-vapply(nodes,function(z)z$key$shift,numeric(1))
proposal_mean<-mean(values);proposal_sd<-sum((values-proposal_mean)*rule$nodes)/sum(rule$nodes^2)
stopifnot(is.finite(proposal_mean),is.finite(proposal_sd),proposal_sd>0,
  max(abs(values-proposal_mean-proposal_sd*rule$nodes))<1e-10)
proposal<-list(mean=proposal_mean,sd=proposal_sd,
  prior_sd=context$prior$coef_sd[which(colnames(context$prior$design)=="intercept")])
previous<-NULL;history<-audits<-conditional<-list();passed<-FALSE
A<-diag(context$n)-matrix(case$reference_weight,context$n,context$n,byrow=TRUE)
for(order in c(9L,17L,33L,65L)) {
  rule<-oracle_normal_quadrature(order);components<-list();lw<-numeric()
  paths<-file.path(cache,sprintf("%s_order%d_node%d.rds",case_id,order,which(rule$weights>0)))
  if(!all(file.exists(paths)))stop("Incomplete cached order: reconstruction cannot create missing fits")
  for(j in which(rule$weights>0)) {
    saved<-readRDS(file.path(cache,sprintf("%s_order%d_node%d.rds",case_id,order,j)))
    stopifnot(identical(saved$key$key,original$key),
      abs(saved$key$shift-proposal$mean-proposal$sd*rule$nodes[j])<1e-10)
    fit<-saved$fit
    stopifnot(fit$diagnostics$covariance_status=="TMB_sdreport",
      fit$diagnostics$forward_error<=1e-5,fit$diagnostics$backend_error<=1e-8)
    data<-design_tmb_data(case,context,saved$key$shift,max_step=fit$diagnostics$max_step)
    obj<-TMB::MakeADFun(data,fit$parameter_mode,random=NULL,DLL="observation_model",silent=TRUE)
    value<-obj$fn(obj$par);gradient<-obj$gr(obj$par);H<-obj$he(obj$par)
    stopifnot(is.finite(value),all(is.finite(gradient)),max(abs(gradient))<=.001,all(is.finite(H)))
    R<-chol((H+t(H))/2)
    objective<-value+sum(log(diag(R)))-.5*length(obj$par)*log(2*pi)
    J<-cbind(matrix(0,context$n,length(fit$parameter_mode$zeta)),diag(context$n),data$design)
    Cg<-J%*%chol2inv(R)%*%t(J)
    covariance_error<-max(abs(Cg-fit$latent_fitness_covariance))
    stored<-fit$diagnostics$objective
    difference<-if(is.finite(stored))abs(stored-objective)else NA_real_
    stopifnot(is.finite(objective),covariance_error<1e-8,
      !is.finite(stored)||difference<1e-6)
    audits[[length(audits)+1L]]<-data.frame(case_id=case_id,order=order,node=j,shift=saved$key$shift,
      stored_objective=stored,reconstructed_objective=objective,finite_objective_difference=difference,
      max_joint_gradient=max(abs(gradient)),latent_covariance_error=covariance_error)
    TMB::FreeADFun(obj);rm(obj);invisible(gc(full=TRUE))
    conditional[[length(conditional)+1L]]<-data.frame(order=order,node=j,shift=saved$key$shift,
      status=fit$diagnostics$covariance_status,forward_error=fit$diagnostics$forward_error,
      refinements=fit$diagnostics$refinements)
    components[[length(components)+1L]]<-list(mean=fit$summary$fitness_mean,covariance=fit$fitness_covariance,
      latent_mean=fit$latent_fitness_mean,latent_covariance=fit$latent_fitness_covariance,intercept=saved$key$shift)
    lw<-c(lw,log(rule$weights[j])-objective+dnorm(saved$key$shift,0,proposal$prior_sd,log=TRUE)-
      dnorm(saved$key$shift,proposal$mean,proposal$sd,log=TRUE))
  }
  weight<-exp(lw-max(lw));weight<-weight/sum(weight)
  current<-baseline_combine(components,weight)
  stopifnot(all(is.finite(current$mean)),all(is.finite(current$covariance)))
  dm<-dc<-dr<-Inf
  if(!is.null(previous)) {
    dm<-max(abs(current$mean-previous$mean)/pmax(1,abs(current$mean)))
    dc<-norm(current$covariance-previous$covariance,"F")/norm(current$covariance,"F")
    dr<-norm(A%*%(current$covariance-previous$covariance)%*%t(A),"F")/norm(A%*%current$covariance%*%t(A),"F")
  }
  passed<-isTRUE(dm<.001&&dc<.005&&dr<.005)
  history[[length(history)+1L]]<-data.frame(order=order,mean_change=dm,covariance_change=dc,
    contrast_covariance_change=dr,passed=passed)
  if(passed)break
  previous<-current
}
audit<-do.call(rbind,audits)
fit<-c(current,list(status=if(passed)"passed_quadrature_refinement"else"quadrature_not_converged",
  full=NULL,history=do.call(rbind,history),conditional_status=do.call(rbind,conditional),proposal=proposal,
  elapsed_seconds=original$fit$elapsed_seconds,original_status=original$fit$status,
  original_message=original$fit$message,evidence_reconstructed=TRUE,
  reconstruction_seconds=proc.time()[["elapsed"]]-started,
  resource_pass=original$fit$elapsed_seconds<=300))
write.csv(audit,file.path(destination,paste0(case_id,"_node_audit.csv")),row.names=FALSE)
saveRDS(list(case_id=case_id,key=original$key,fit=fit,scored=NULL,audit=audit),
  file.path(destination,paste0(case_id,"_result.rds")))
print(fit$history,row.names=FALSE)
cat("Recovered numerical status",fit$status,"; original resource gate",fit$resource_pass,
  "; reconstruction seconds",fit$reconstruction_seconds,"; no optimization performed\n")
