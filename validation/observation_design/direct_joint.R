# Equivalent conditional Laplace computation when there are no free scales.
# Optimize the joint density once, then use its joint Hessian directly. This
# avoids differentiating a Laplace objective that has zero outer parameters.
# All statistical data/prior/forward helpers remain the frozen v4 functions.
direct_joint_fit <- function(case,context,runtime,shift=NULL,initial=NULL,max_step=.1,
                             refinements=0L,check_forward=TRUE) {
  data<-design_tmb_data(case,context,shift,max_step)
  parameters<-design_parameters(case,context,shift,initial)
  obj<-TMB::MakeADFun(data,parameters,random=NULL,DLL="observation_model",silent=TRUE)
  on.exit(TMB::FreeADFun(obj),add=TRUE)
  sparse_hessian<-function(par) {
    HH<-obj$he(par);d<-nrow(HH)
    Matrix::forceSymmetric(Matrix::sparseMatrix(i=rep(seq_len(d),d),
      j=rep(seq_len(d),each=d),x=as.numeric(HH),dims=c(d,d)),uplo="L")
  }
  newton<-get("newton",asNamespace("TMB"))
  attempt<-try(newton(obj$par,obj$fn,function(p)as.numeric(obj$gr(p)),sparse_hessian,
    trace=0,maxit=2000L,tol=1e-10,grad.tol=1e-10,step.tol=1e-12,tol10=0,
    env=new.env(parent=globalenv())),silent=TRUE)
  if(!inherits(attempt,"try-error")&&is.finite(attempt$value)&&
     all(is.finite(attempt$gradient))&&max(abs(attempt$gradient))<=.001) {
    opt<-list(par=attempt$par,convergence=0L,evaluations=attempt$iterations)
    optimization_method<-"joint_newton"
  } else {
    opt<-nlminb(obj$par,obj$fn,obj$gr,
      control=list(eval.max=5000L,iter.max=3000L,rel.tol=1e-12,x.tol=1e-10))
    optimization_method<-"nlminb_fallback"
  }
  value<-obj$fn(opt$par);gradient<-as.numeric(obj$gr(opt$par))
  if(!is.finite(value)||any(!is.finite(gradient)))stop("Nonfinite direct joint mode/evidence")
  if(max(abs(gradient))>.001||opt$convergence!=0L) {
    opt<-nlminb(opt$par,obj$fn,obj$gr,
      control=list(eval.max=8000L,iter.max=5000L,rel.tol=1e-13,x.tol=1e-11))
    value<-obj$fn(opt$par);gradient<-as.numeric(obj$gr(opt$par))
  }
  # The trust criterion is stationarity plus positive curvature, as for the
  # zero-outer-parameter nested fit; preserve the optimizer flag separately.
  if(!is.finite(value)||any(!is.finite(gradient))||max(abs(gradient))>.001)
    stop(sprintf("Untrusted direct joint optimization at shift %s: gradient %.6g, convergence %s",
      if(is.null(shift))"full"else format(shift,digits=16),max(abs(gradient)),opt$convergence))
  par<-opt$par;polishing_steps<-0L
  for(polish in seq_len(8L)) {
    if(max(abs(gradient))<=1e-8)break
    HH<-obj$he(par);RR<-chol((HH+t(HH))/2)
    direction<-backsolve(RR,forwardsolve(t(RR),gradient))
    accepted<-FALSE
    for(alpha in 2^-(0:20)) {
      candidate<-par-alpha*direction
      candidate_value<-obj$fn(candidate);candidate_gradient<-as.numeric(obj$gr(candidate))
      if(is.finite(candidate_value)&&all(is.finite(candidate_gradient))&&candidate_value<=value+1e-9&&
         max(abs(candidate_gradient))<max(abs(gradient))) {
        par<-candidate;value<-candidate_value;gradient<-candidate_gradient
        accepted<-TRUE;polishing_steps<-polishing_steps+1L;break
      }
    }
    if(!accepted)break
  }
  opt$par<-par
  if(max(abs(gradient))>1e-6)stop("Direct joint Newton polishing did not reach stationarity")
  H<-obj$he(opt$par);H<-(H+t(H))/2
  if(any(!is.finite(H)))stop("Nonfinite direct joint Hessian")
  R<-chol(H)
  objective<-value+sum(log(diag(R)))-.5*length(opt$par)*log(2*pi)
  mode<-obj$env$parList(opt$par);report<-obj$report(opt$par)
  J<-cbind(matrix(0,context$n,length(mode$zeta)),diag(context$n),data$design)
  Cg<-J%*%chol2inv(R)%*%t(J);Cg<-(Cg+t(Cg))/2
  mu<-as.numeric(report$f_latent)
  b<-exp(mu+.5*diag(Cg));Cb<-tcrossprod(b)*expm1(Cg)
  if(!is.finite(objective)||any(!is.finite(b))||any(!is.finite(Cb)))stop("Nonfinite direct posterior moments")
  forward_error<-backend_error<-0
  if(check_forward)for(o in which(case$times>0)) {
    pi<-report$pi0[,case$trajectory[o]+1L];birth<-as.numeric(report$f_report)
    coarse<-direct_forward_r(pi,birth,context$operator,case$times[o],data$steps[o],data$taylor_order)
    fine<-direct_forward_r(pi,birth,context$operator,case$times[o],2L*data$steps[o],data$taylor_order)
    exact<-direct_exact_forward(coarse$generator,pi,case$times[o],runtime)
    exact_pi<-exact$pi1;exact_growth<-exact$log_total_growth
    forward_error<-max(forward_error,max(abs(coarse$pi1-fine$pi1)),abs(coarse$log_total_growth-fine$log_total_growth),
      max(abs(coarse$pi1-exact_pi)),abs(coarse$log_total_growth-exact_growth))
    backend_error<-max(backend_error,max(abs(coarse$pi1-report$predicted[,o])))
  }
  if(!is.finite(forward_error)||!is.finite(backend_error))stop("Nonfinite direct forward check")
  if(backend_error>1e-8)stop("Direct forward backend disagreement")
  if(forward_error>1e-5) {
    if(refinements>=3L)stop("Direct forward accuracy untrusted")
    TMB::FreeADFun(obj);on.exit(NULL);rm(obj);invisible(gc(full=TRUE))
    return(direct_joint_fit(case,context,runtime,shift,mode,max_step/2,refinements+1L,check_forward))
  }
  list(summary=data.frame(karyotype=context$graph$labels,fitness_mean=b),
    fitness_covariance=Cb,latent_fitness_mean=mu,latent_fitness_covariance=Cg,
    parameter_mode=mode,posterior_predictive=report,
    diagnostics=list(objective=objective,joint_objective=value,
      covariance_status="trusted_direct_joint_hessian",inner_gradient=max(abs(gradient)),
      forward_error=forward_error,backend_error=backend_error,max_step=max_step,refinements=refinements,
      optimization_convergence=opt$convergence,optimization_evaluations=opt$evaluations,
      newton_polishing_steps=polishing_steps,taylor_order=data$taylor_order,
      optimization_method=optimization_method))
}

direct_joint_integrate <- function(case,context,runtime,cache=NULL,key=NULL) {
  full<-direct_joint_fit(case,context,runtime)
  proposal<-baseline_proposal(c(context,list(weight=case$reference_weight)),full)
  previous<-NULL;history<-conditional<-list()
  baseline_index<-which(colnames(context$prior$design)=="intercept")
  pool<-list(list(shift=full$parameter_mode$coefficients[baseline_index],fit=full,full=TRUE))
  A<-diag(context$n)-matrix(case$reference_weight,context$n,context$n,byrow=TRUE)
  for(order in c(9L,17L,33L,65L)) {
    rule<-oracle_normal_quadrature(order);components<-vector("list",length(rule$nodes));lw<-rep(NA_real_,length(rule$nodes))
    live<-which(rule$weights>0)
    visit<-live[order(abs(rule$nodes[live]))]
    for(j in visit) {
      value<-proposal$mean+proposal$sd*rule$nodes[j]
      nearest<-pool[[which.min(vapply(pool,function(p)abs(p$shift-value),numeric(1)))]]
      initial<-nearest$fit$parameter_mode
      if(!isTRUE(nearest$full))initial$u<-initial$u+nearest$shift-value
      if(Sys.getenv("ALFAK2_DIRECT_TRACE","")=="1")cat("Conditional order",order,"node",j,"shift",value,"\n")
      filename<-if(is.null(cache))NULL else file.path(cache,sprintf("%s_order%d_node%d.rds",case$id,order,j))
      node_key<-list(key=key,shift=value)
      if(!is.null(filename)&&file.exists(filename)) {
        saved<-readRDS(filename);stopifnot(identical(saved$key,node_key));fit<-saved$fit
      } else {
        fit<-tryCatch(direct_joint_fit(case,context,runtime,shift=value,initial=initial),
          error=function(e)stop(sprintf("%s order %d node %d shift %.16g: %s",case$id,order,j,value,conditionMessage(e)),call.=FALSE))
        if(!is.null(filename))saveRDS(list(key=node_key,fit=fit),filename)
      }
      pool[[length(pool)+1L]]<-list(shift=value,fit=fit,full=FALSE)
      stopifnot(fit$diagnostics$covariance_status=="trusted_direct_joint_hessian",
        is.finite(fit$diagnostics$objective),all(is.finite(fit$summary$fitness_mean)),
        all(is.finite(fit$fitness_covariance)))
      conditional[[length(conditional)+1L]]<-data.frame(order=order,node=j,shift=value,
        status=fit$diagnostics$covariance_status,forward_error=fit$diagnostics$forward_error,
        refinements=fit$diagnostics$refinements)
      components[[j]]<-list(mean=fit$summary$fitness_mean,covariance=fit$fitness_covariance,
        latent_mean=fit$latent_fitness_mean,latent_covariance=fit$latent_fitness_covariance,intercept=value)
      lw[j]<-log(rule$weights[j])-fit$diagnostics$objective+dnorm(value,0,proposal$prior_sd,log=TRUE)-
        dnorm(value,proposal$mean,proposal$sd,log=TRUE)
    }
    components<-components[live];lw<-lw[live]
    weight<-exp(lw-max(lw));weight<-weight/sum(weight)
    current<-baseline_combine(components,weight)
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
  c(current,list(status=if(passed)"passed_quadrature_refinement"else"quadrature_not_converged",
    full=full,history=do.call(rbind,history),conditional_status=do.call(rbind,conditional),proposal=proposal,
    computational_method="direct_joint_mode_hessian"))
}
