design_load <- function(repo, compile=TRUE) {
  helper <- new.env(parent=globalenv())
  sys.source(file.path(repo,"validation/posterior_reference_4to8/reference.R"),helper)
  runtime <- helper$reference_load(repo,compile=FALSE)
  module <- Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
  files <- file.path(module,c("observation_model.cpp","positive_taylor_forward.h"))
  sha <- paste(tools::md5sum(files),collapse="")
  cache <- file.path("/private/tmp/alfak2_observation_design_model",sha)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  for(f in files)file.copy(f,file.path(cache,basename(f)),overwrite=FALSE)
  binary <- TMB::dynlib(file.path(cache,"observation_model"))
  if(!file.exists(binary)&&compile)TMB::compile(file.path(cache,"observation_model.cpp"),flags="-O2 -g0")
  if(!file.exists(binary))stop("Joint design model not compiled")
  dyn.load(binary)
  runtime$design_binary <- binary
  runtime
}

design_context <- function(karyotypes,beta,runtime) {
  labels <- apply(karyotypes,1,paste,collapse=".")
  graph <- runtime$pkg$build_karyotype_graph(labels,beta=beta,shell_depth=0,
    min_cn=min(karyotypes),max_cn=max(karyotypes),max_nodes=nrow(karyotypes))
  prior <- runtime$pkg$build_joint_landscape_prior(graph,baseline_sd=.5,main_sd=.07,
    curvature_sd=.04,interaction_sd=.04,residual_sd=.125,edge_sd=.25,interaction_support="isotropic")
  K <- solve(as.matrix(prior$precision))+tcrossprod(sweep(prior$design,2,prior$coef_sd,"*"))
  basis <- contr.helmert(nrow(karyotypes));basis <- sweep(basis,2,sqrt(colSums(basis^2)),"/")
  list(n=nrow(karyotypes),graph=graph,prior=prior,covariance=K,
    eta_basis=2*basis,operator=runtime$pkg$division_transition_data(graph))
}

design_tmb_data <- function(case,context,shift=NULL,max_step=.1) {
  n <- context$n
  stopifnot(nrow(case$y)==n,length(case$trajectory)==ncol(case$y),
    length(case$times)==ncol(case$y),all(case$y>=0),all(colSums(case$y)>0),
    all(case$times>=0),all(case$trajectory>=0),!anyNA(case$y))
  prior <- context$prior
  if(!is.null(shift)) {
    j<-which(colnames(prior$design)=="intercept")
    prior$design<-prior$design[,-j,drop=FALSE];prior$coef_sd<-prior$coef_sd[-j]
    prior$offset<-prior$offset+shift
  }
  op <- context$operator
  list(y=case$y,trajectory=as.integer(case$trajectory),times=as.numeric(case$times),
    steps=as.integer(pmax(4,ceiling(case$times/max_step))),eta_basis=context$eta_basis,
    design=prior$design,offset=prior$offset,precision=prior$precision,logdet=prior$logdet,
    coef_sd=prior$coef_sd,trans_from=as.integer(op$from0),trans_to=as.integer(op$to0),
    trans_weight=as.numeric(op$weight))
}

design_parameters <- function(case,context,shift=NULL,initial=NULL) {
  n<-context$n; nr<-max(case$trajectory)+1L
  nc<-ncol(context$prior$design)-as.integer(!is.null(shift))
  if(!is.null(initial)) {
    par<-initial
    if(!is.null(shift)&&length(par$coefficients)==ncol(context$prior$design)) {
      j<-which(colnames(context$prior$design)=="intercept")
      par$u<-par$u+par$coefficients[j]-shift
      par$coefficients<-par$coefficients[-j]
    }
    return(par)
  }
  zeta<-matrix(0,n-1L,nr)
  for(r in seq_len(nr)) {
    j<-which(case$trajectory==r-1L & case$times==0)
    if(length(j)!=1L)stop("Every trajectory must have one initial count observation")
    p<-(case$y[,j]+.5)/sum(case$y[,j]+.5)
    zeta[,r]<-as.numeric(crossprod(context$eta_basis,log(p)-mean(log(p))))/4
  }
  list(zeta=zeta,u=rep(0,n),coefficients=rep(0,nc))
}

design_forward_r <- function(pi0,birth,operator,dt,steps) {
  n<-length(pi0);P<-matrix(0,n,n)
  P[cbind(operator$from0+1L,operator$to0+1L)]<-operator$weight
  A<-sweep(2*t(P)-diag(n),2,birth,"*")
  shift<-min(diag(A));B<-A-diag(shift,n)
  state<-as.numeric(pi0);growth<-0;h<-dt/steps
  for(s in seq_len(steps)) {
    term<-total<-state
    for(order in 1:8){term<-as.numeric(B%*%term)*h/order;total<-total+term}
    state<-total;growth<-growth+h*shift
    if(s%%16L==0L||s==steps){z<-sum(state);growth<-growth+log(z);state<-state/z}
  }
  list(pi1=state,log_total_growth=growth,generator=A)
}

design_fit <- function(case,context,runtime,shift=NULL,initial=NULL,max_step=.1,
                        refinements=0L,keep_object=FALSE,check_forward=TRUE) {
  data<-design_tmb_data(case,context,shift,max_step)
  parameters<-design_parameters(case,context,shift,initial)
  obj<-TMB::MakeADFun(data,parameters,random=c("zeta","u","coefficients"),
    DLL="observation_model",silent=TRUE,inner.control=list(maxit=2000,tol=1e-10))
  attempt<-runtime$pkg$run_local_tmb_attempt(obj,obj$par,list(eval.max=1000,iter.max=1000),context$n,require_latent=TRUE)
  status<-runtime$pkg$assess_local_covariance(attempt,.001)
  if(status!="TMB_sdreport") {
    attempt<-runtime$pkg$run_local_tmb_attempt(obj,attempt$opt$par,list(eval.max=2000,iter.max=2000),context$n,require_latent=TRUE)
    status<-runtime$pkg$assess_local_covariance(attempt,.001)
  }
  mode<-runtime$pkg$local_tmb_parameter_list(obj,attempt$opt$par)
  report<-obj$report()
  forward_error<-0;backend_error<-0
  if(check_forward)for(o in which(case$times>0)) {
    pi<-report$pi0[,case$trajectory[o]+1L];birth<-as.numeric(report$f_report)
    coarse<-design_forward_r(pi,birth,context$operator,case$times[o],data$steps[o])
    fine<-design_forward_r(pi,birth,context$operator,case$times[o],2L*data$steps[o])
    exact_mass<-runtime$independent$process_deterministic_mean(coarse$generator,pi,case$times[o])
    exact_pi<-exact_mass/sum(exact_mass);exact_growth<-log(sum(exact_mass))
    forward_error<-max(forward_error,max(abs(coarse$pi1-fine$pi1)),abs(coarse$log_total_growth-fine$log_total_growth))
    forward_error<-max(forward_error,max(abs(coarse$pi1-exact_pi)),abs(coarse$log_total_growth-exact_growth))
    backend_error<-max(backend_error,max(abs(coarse$pi1-report$predicted[,o])))
  }
  if(backend_error>1e-8)status<-"untrusted_backend_disagreement"
  if(forward_error>1e-5) {
    if(refinements<3L) {
      TMB::FreeADFun(obj)
      rm(obj,attempt)
      invisible(gc(full=TRUE))
      return(design_fit(case,context,runtime,shift,mode,max_step/2,
        refinements+1L,keep_object,check_forward))
    }
    status<-"untrusted_forward_accuracy"
  }
  Cg<-runtime$pkg$extract_local_fitness_covariance(attempt$sdrep,attempt$latent_frep,status)
  mu<-if(!is.null(attempt$latent_frep))attempt$latent_frep$mean else as.numeric(report$f_latent)
  b<-if(!is.null(Cg))exp(mu+.5*diag(Cg))else rep(NA_real_,context$n)
  Cb<-if(!is.null(Cg))tcrossprod(b)*expm1(Cg)else NULL
  answer<-list(summary=data.frame(karyotype=context$graph$labels,fitness_mean=b),
    fitness_covariance=Cb,latent_fitness_mean=mu,latent_fitness_covariance=Cg,
    parameter_mode=mode,posterior_predictive=report,
    diagnostics=list(objective=attempt$opt$objective,covariance_status=status,
      inner_gradient=attempt$conditional_gradient_norm,forward_error=forward_error,
      backend_error=backend_error,max_step=max_step,refinements=refinements))
  if(keep_object)answer$obj<-obj
  else {
    TMB::FreeADFun(obj)
    rm(obj,attempt)
    invisible(gc(full=TRUE))
  }
  answer
}

design_integrate <- function(case,context,runtime,cache=NULL,key=NULL) {
  full<-design_fit(case,context,runtime)
  if(full$diagnostics$covariance_status!="TMB_sdreport")return(list(status="untrusted_full_fit",full=full))
  scalar_case<-c(context,list(weight=case$reference_weight))
  proposal<-baseline_proposal(scalar_case,full)
  previous<-NULL;history<-list();conditional<-list()
  A<-diag(context$n)-matrix(case$reference_weight,context$n,context$n,byrow=TRUE)
  for(order in c(9L,17L,33L,65L)) {
    rule<-oracle_normal_quadrature(order);components<-list();lw<-numeric()
    for(j in which(rule$weights>0)) {
      value<-proposal$mean+proposal$sd*rule$nodes[j]
      filename<-if(is.null(cache))NULL else file.path(cache,sprintf("%s_order%d_node%d.rds",case$id,order,j))
      node_key<-list(key=key,shift=value)
      if(!is.null(filename)&&file.exists(filename)) {
        saved<-readRDS(filename)
        if(!identical(saved$key,node_key))stop("Joint conditional cache mismatch")
        fit<-saved$fit
      } else {
        fit<-design_fit(case,context,runtime,shift=value,initial=full$parameter_mode)
        if(!is.null(filename))saveRDS(list(key=node_key,fit=fit),filename)
      }
      conditional[[length(conditional)+1L]]<-data.frame(order=order,node=j,shift=value,
        status=fit$diagnostics$covariance_status,forward_error=fit$diagnostics$forward_error,
        refinements=fit$diagnostics$refinements)
      if(fit$diagnostics$covariance_status!="TMB_sdreport")return(list(status="untrusted_conditional_node",
        full=full,conditional_status=do.call(rbind,conditional)))
      components[[length(components)+1L]]<-list(mean=fit$summary$fitness_mean,covariance=fit$fitness_covariance,
        latent_mean=fit$latent_fitness_mean,latent_covariance=fit$latent_fitness_covariance,intercept=value)
      lw<-c(lw,log(rule$weights[j])-fit$diagnostics$objective+
        dnorm(value,0,proposal$prior_sd,log=TRUE)-dnorm(value,proposal$mean,proposal$sd,log=TRUE))
    }
    weight<-exp(lw-max(lw));weight<-weight/sum(weight)
    current<-baseline_combine(components,weight)
    dm<-dc<-dr<-Inf
    if(!is.null(previous)) {
      dm<-max(abs(current$mean-previous$mean)/pmax(1,abs(current$mean)))
      dc<-norm(current$covariance-previous$covariance,"F")/norm(current$covariance,"F")
      dr<-norm(A%*%(current$covariance-previous$covariance)%*%t(A),"F")/norm(A%*%current$covariance%*%t(A),"F")
    }
    passed<-dm<.001&&dc<.005&&dr<.005
    history[[length(history)+1L]]<-data.frame(order=order,mean_change=dm,covariance_change=dc,
      contrast_covariance_change=dr,passed=passed)
    if(passed)break
    previous<-current
  }
  c(current,list(status=if(passed)"passed_quadrature_refinement"else"quadrature_not_converged",
    full=full,history=do.call(rbind,history),conditional_status=do.call(rbind,conditional),proposal=proposal))
}
