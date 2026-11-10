adaptive_design_integrate<-function(case,context,runtime,cache=NULL,key=NULL,trace=FALSE) {
  full<-direct_joint_fit(case,context,runtime)
  j<-which(colnames(context$prior$design)=="intercept")
  stopifnot(length(j)==1L,all(abs(context$prior$design[,j]-1)<1e-12))
  pool<-list(list(shift=full$parameter_mode$coefficients[j],fit=full,full=TRUE))
  status<-list()
  evaluate<-function(value) {
    nearest<-pool[[which.min(vapply(pool,function(p)abs(p$shift-value),numeric(1)))]]
    initial<-nearest$fit$parameter_mode
    if(!isTRUE(nearest$full))initial$u<-initial$u+nearest$shift-value
    node_key<-list(source_case=key,intercept=value)
    filename<-if(is.null(cache))NULL else file.path(cache,paste0(round2_hash_object(node_key),".rds"))
    if(!is.null(filename)&&file.exists(filename)) {
      saved<-readRDS(filename);stopifnot(identical(saved$key,node_key));fit<-saved$fit
    } else {
      fit<-direct_joint_fit(case,context,runtime,shift=value,initial=initial)
      if(!is.null(filename))saveRDS(list(key=node_key,fit=fit),filename)
    }
    stopifnot(fit$diagnostics$covariance_status=="trusted_direct_joint_hessian",
      is.finite(fit$diagnostics$objective),fit$diagnostics$inner_gradient<=.001)
    pool[[length(pool)+1L]]<<-list(shift=value,fit=fit,full=FALSE)
    status[[length(status)+1L]]<<-data.frame(intercept=value,status=fit$diagnostics$covariance_status,
      objective=fit$diagnostics$objective,forward_error=fit$diagnostics$forward_error,
      backend_error=fit$diagnostics$backend_error,gradient=fit$diagnostics$inner_gradient)
    list(mean=fit$summary$fitness_mean,covariance=fit$fitness_covariance,
      latent_mean=fit$latent_fitness_mean,latent_covariance=fit$latent_fitness_covariance,
      log_evidence=-as.numeric(fit$diagnostics$objective))
  }
  result<-adaptive_prior_integrate(evaluate,context$prior$coef_sd[j],case$reference_weight,trace=trace)
  result$full<-full;result$conditional_status<-do.call(rbind,status)
  result
}
