is_load<-function(repo,src,compile=TRUE) {
  lock<-"/private/tmp/alfak2_importance_setup.lock";started<-proc.time()[["elapsed"]]
  while(!dir.create(lock,showWarnings=FALSE)) {
    if(proc.time()[["elapsed"]]-started>120)stop("Shared compiled-cache lock timeout")
    Sys.sleep(.05)
  }
  writeLines(as.character(Sys.getpid()),file.path(lock,"owner"));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  helper<-new.env(parent=globalenv());sys.source(file.path(repo,"validation/posterior_reference_4to8/reference.R"),helper)
  r<-helper$reference_load(repo,compile=FALSE)
  files<-file.path(src,c("importance_model.cpp","positive_taylor_forward.h"));key<-paste(tools::md5sum(files),collapse="")
  cache<-file.path("/private/tmp/alfak2_importance_model",key);dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  for(f in files)file.copy(f,file.path(cache,basename(f)),overwrite=FALSE)
  binary<-TMB::dynlib(file.path(cache,"importance_model"))
  if(!file.exists(binary)&&compile)TMB::compile(file.path(cache,"importance_model.cpp"),flags="-O2 -g0")
  stopifnot(file.exists(binary));dyn.load(binary);invisible(TMB::openmp(1,DLL="importance_model"))
  r$importance_binary<-binary;r$core<-new.env(parent=globalenv())
  Rcpp::sourceCpp(file.path(src,"importance.cpp"),env=r$core,
    cacheDir=file.path("/private/tmp/alfak2_importance_core",unname(tools::md5sum(file.path(src,"importance.cpp")))),showOutput=FALSE)
  r$pilot_core<-new.env(parent=globalenv())
  Rcpp::sourceCpp(file.path(repo,"validation/observation_design/results/mechanism_reference_source_v2/reference.cpp"),
    env=r$pilot_core,cacheDir="/private/tmp/alfak2_mechanism_reference_cpp_v2",showOutput=FALSE)
  r
}

is_spec<-function(case,context) {
  y<-if(is.null(case$y))case$counts else case$y;n<-nrow(y)
  stopifnot(ncol(y)==2L,length(case$times)==2L,case$times[1]==0,case$times[2]>0,all(case$trajectory==0))
  known<-if(is.null(case$known_pi0))numeric()else as.numeric(case$known_pi0)
  stopifnot(length(known)%in%c(0L,n));if(length(known))stopifnot(all(known>=0),abs(sum(known)-1)<1e-12)
  op<-context$operator;P<-matrix(0,n,n);P[cbind(op$from0+1L,op$to0+1L)]<-op$weight
  list(n=n,q=as.integer(n+if(length(known))0L else n-1L),dt=case$times[2],offset=context$prior$offset,
    lower=t(chol(context$covariance)),eta_basis=context$eta_basis,known_pi0=known,probability=P,
    beta_zero=all(abs(P-diag(n))<1e-14),y0=as.numeric(y[,1]),y1=as.numeric(y[,2]))
}

is_mode_proposal<-function(spec,runtime,max_step=.1,initial=NULL,refinement=0L) {
  n<-spec$n;P<-spec$probability;edge<-which(P>0,arr.ind=TRUE)
  data<-list(y0=spec$y0,y1=spec$y1,dt=spec$dt,steps=as.integer(max(4,ceiling(spec$dt/max_step))),
    beta_zero=as.integer(spec$beta_zero),offset=spec$offset,lower=spec$lower,eta_basis=spec$eta_basis,
    known_pi0=spec$known_pi0,trans_from=as.integer(edge[,1]-1L),trans_to=as.integer(edge[,2]-1L),trans_weight=as.numeric(P[edge]))
  if(is.null(initial)) {
    p<-(spec$y0+.5)/sum(spec$y0+.5)
    ze<-if(spec$q>n)as.numeric(crossprod(spec$eta_basis,log(p)-mean(log(p))))/4 else numeric()
    initial<-list(z_g=numeric(n),z_eta=ze)
  }
  obj<-TMB::MakeADFun(data,initial,random=NULL,DLL="importance_model",silent=TRUE)
  on.exit(TMB::FreeADFun(obj),add=TRUE)
  hessian<-function(par)Matrix::forceSymmetric(Matrix::Matrix(obj$he(par),sparse=TRUE),uplo="L")
  attempt<-try(get("newton",asNamespace("TMB"))(obj$par,obj$fn,function(p)as.numeric(obj$gr(p)),hessian,
    trace=0,maxit=2000L,tol=1e-10,grad.tol=1e-10,step.tol=1e-12,tol10=0,env=new.env(parent=globalenv())),silent=TRUE)
  if(!inherits(attempt,"try-error")&&is.finite(attempt$value)&&max(abs(attempt$gradient))<=.001)par<-attempt$par
  else par<-nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000L,iter.max=3000L,rel.tol=1e-12,x.tol=1e-10))$par
  value<-obj$fn(par);gradient<-as.numeric(obj$gr(par))
  for(j in 1:8) {
    if(max(abs(gradient))<=1e-8)break
    H<-obj$he(par);R<-chol((H+t(H))/2);step<-backsolve(R,forwardsolve(t(R),gradient))
    accepted<-FALSE
    for(alpha in 2^-(0:20)) {
      nextpar<-par-alpha*step;v<-obj$fn(nextpar);g<-as.numeric(obj$gr(nextpar))
      if(is.finite(v)&&all(is.finite(g))&&v<=value+1e-9&&max(abs(g))<max(abs(gradient))) {
        par<-nextpar;value<-v;gradient<-g;accepted<-TRUE;break
      }
    }
    if(!accepted)break
  }
  stopifnot(all(is.finite(gradient)),max(abs(gradient))<=1e-6,is.finite(value))
  H<-obj$he(par);R<-chol((H+t(H))/2);C<-chol2inv(R);report<-obj$report(par)
  target<-runtime$core$is_forward_cpp(as.numeric(report$p0),as.numeric(report$b),P,spec$dt,spec$beta_zero)
  err<-max(abs(target$pi1-as.numeric(report$p1)))
  exactll<-sum(spec$y0*log(report$p0+1e-14))+sum(spec$y1*log(target$pi1+1e-14))
  approxll<-sum(spec$y0*log(report$p0+1e-14))+sum(spec$y1*log(report$p1+1e-14))
  if(err>1e-8||abs(exactll-approxll)>1e-5) {
    if(refinement>=3L)stop("Proposal forward accuracy not certified")
    parameters<-obj$env$parList(par);TMB::FreeADFun(obj);on.exit(NULL)
    return(is_mode_proposal(spec,runtime,max_step/2,parameters,refinement+1L))
  }
  list(mean=as.numeric(par),covariance=C,lower=t(chol(C)),mode_objective=value,
    gradient=max(abs(gradient)),forward_error=err,forward_loglik_error=abs(exactll-approxll),
    steps=data$steps,parameter_mode=obj$env$parList(par),report=report,refinement=refinement)
}

is_draw_proposal<-function(proposal,N,seed) {
  if(!is.null(proposal$transport))return(is_transport_draw(proposal,N,seed))
  set.seed(seed);d<-length(proposal$mean);component<-findInterval(runif(N),c(.80,.99))+1L
  z<-matrix(rnorm(N*d),N,d)
  ii<-component!=3L
  if(any(ii))z[ii,]<-z[ii,,drop=FALSE]%*%t(proposal$lower)
  jj<-component==2L
  if(any(jj))z[jj,]<-z[jj,,drop=FALSE]*sqrt(5/rchisq(sum(jj),7))
  if(any(ii))z[ii,]<-sweep(z[ii,,drop=FALSE],2,proposal$mean,"+")
  z
}

is_logq<-function(z,proposal) {
  if(!is.null(proposal$transport))return(is_transport_logq(z,proposal))
  d<-ncol(z);D<-forwardsolve(proposal$lower,t(sweep(z,2,proposal$mean,"-")));q<-colSums(D*D)
  ld<-sum(log(diag(proposal$lower)))
  normal<--.5*(d*log(2*pi)+q)-ld
  student<-lgamma((7+d)/2)-lgamma(7/2)-d/2*log(7*pi)-ld-d/2*log(5/7)-(7+d)/2*log1p(q/5)
  prior<--.5*(d*log(2*pi)+rowSums(z*z))
  a<-normal+log(.80);b<-student+log(.19);c<-prior+log(.01);m<-pmax(a,b,c)
  m+log(exp(a-m)+exp(b-m)+exp(c-m))
}

is_weights<-function(logw) {w<-exp(logw-max(logw));w/sum(w)}
is_khat<-function(logw)as.numeric(loo::psis(as.numeric(logw),r_eff=1,cores=1)$diagnostics$pareto_k)

is_pilot_proposal<-function(spec,proposal,runtime,seed) {
  z<-is_draw_proposal(proposal,16384L,seed)
  ll<-runtime$core$is_loglikelihood_cpp(z,spec)
  lw<-ll-.5*(spec$q*log(2*pi)+rowSums(z*z))-is_logq(z,proposal)
  w<-is_weights(lw);ess<-1/sum(w*w);k<-is_khat(lw);p<-proposal
  if(ess<8192||!is.finite(k)||k>.4) {
    p<-is_transport_fit(spec,proposal,runtime,seed+1000L)
    p$pilot$initial_importance<-list(ess=ess,pareto_k=k,samples=nrow(z),seed=seed)
    return(p)
  }
  if(ess>=max(1000,20*spec$q)) {
    mu<-colSums(z*w);X<-sweep(z,2,mu,"-");C<-crossprod(X*w,X)
    p$mean<-.9*mu+.1*proposal$mean;p$covariance<-.9*C+.1*proposal$covariance+diag(1e-6,spec$q)
    p$lower<-t(chol(p$covariance));p$adapted<-TRUE
  } else p$adapted<-FALSE
  p$pilot<-list(ess=ess,pareto_k=k,samples=nrow(z),seed=seed)
  p
}

is_transport_design<-function(t,p) {
 x<-(t-p$center)/p$scale
 if(!is.null(p$knots))cbind(1,splines::ns(x,knots=p$knots,Boundary.knots=p$boundary))
 else cbind(1,x,x*x,x*x*x)
}
is_transport_to_y<-function(z,p) {
 y<-z%*%p$rotation;y[,-1]<-y[,-1,drop=FALSE]-is_transport_design(y[,1],p)%*%p$regression
 jac<-rep(0,nrow(y))
 if(!is.null(p$conditional_covariance)) {
  bin<-findInterval(y[,1],p$conditional_covariance$cuts,all.inside=TRUE)
  for(k in seq_along(p$conditional_covariance$lower)) {
   ii<-which(bin==k);if(!length(ii))next
   L<-p$conditional_covariance$lower[[k]]
   y[ii,-1]<-t(forwardsolve(L,t(y[ii,-1,drop=FALSE])));jac[ii]<-jac[ii]-sum(log(diag(L)))
  }
 }
 if(!is.null(p$marginal))for(j in seq_len(ncol(y))) {
  m<-is_piecewise(y[,j],p$marginal[[j]]$x,p$marginal[[j]]$u);y[,j]<-m$value;jac<-jac+m$log_jacobian
 }
 if(!is.null(p$layers))for(layer in p$layers) {
  y<-y%*%layer$rotation
  for(j in seq_len(ncol(y))) {
   m<-is_piecewise(y[,j],layer$marginal[[j]]$x,layer$marginal[[j]]$u);y[,j]<-m$value;jac<-jac+m$log_jacobian
  }
 }
 attr(y,"log_jacobian")<-jac
 y
}
is_transport_to_z<-function(y,p) {
 if(!is.null(p$layers))for(layer in rev(p$layers)) {
  for(j in seq_len(ncol(y)))y[,j]<-is_piecewise(y[,j],layer$marginal[[j]]$u,layer$marginal[[j]]$x)$value
  y<-y%*%t(layer$rotation)
 }
 if(!is.null(p$marginal))for(j in seq_len(ncol(y)))y[,j]<-is_piecewise(y[,j],p$marginal[[j]]$u,p$marginal[[j]]$x)$value
 if(!is.null(p$conditional_covariance)) {
  bin<-findInterval(y[,1],p$conditional_covariance$cuts,all.inside=TRUE)
  for(k in seq_along(p$conditional_covariance$lower)) {
   ii<-which(bin==k);if(!length(ii))next
   y[ii,-1]<-y[ii,-1,drop=FALSE]%*%t(p$conditional_covariance$lower[[k]])
  }
 }
 y[,-1]<-y[,-1,drop=FALSE]+is_transport_design(y[,1],p)%*%p$regression;y%*%t(p$rotation)
}
is_piecewise<-function(x,from,to) {
 k<-pmax(1L,pmin(length(from)-1L,findInterval(x,from)))
 slope<-(to[k+1L]-to[k])/(from[k+1L]-from[k])
 list(value=to[k]+slope*(x-from[k]),log_jacobian=log(slope))
}
is_transport_logjac<-function(z,p) {
 attr(is_transport_to_y(z,p),"log_jacobian")
}
is_gaussian_logq<-function(y,p,student=FALSE,product=FALSE) {
 d<-ncol(y);D<-forwardsolve(p$lower,t(sweep(y,2,p$mean,"-")));q<-colSums(D*D);ld<-sum(log(diag(p$lower)))
 if(product)colSums(matrix(dt(as.numeric(D)/sqrt(5/7),df=7,log=TRUE),d))-d/2*log(5/7)-ld
 else if(student)lgamma((7+d)/2)-lgamma(7/2)-d/2*log(7*pi)-ld-d/2*log(5/7)-(7+d)/2*log1p(q/5)
 else -.5*(d*log(2*pi)+q)-ld
}
is_transport_logq<-function(z,p) {
 y<-is_transport_to_y(z,p$transport)
 jac<-attr(y,"log_jacobian")
 K<-length(p$components)
 V<-cbind(vapply(p$components,function(c)is_gaussian_logq(y,c)+log(.30/K),numeric(nrow(z))),
  vapply(p$components,function(c)is_gaussian_logq(y,c,product=TRUE)+log(.30/K),numeric(nrow(z))))
 V<-cbind(V,is_gaussian_logq(y,p$global,TRUE)+log(.19),is_gaussian_logq(y,p$tilted)+log(.20))+jac
 V<-cbind(V,-.5*(ncol(z)*log(2*pi)+rowSums(z*z))+log(.01))
 mx<-apply(V,1,max);mx+log(rowSums(exp(V-mx)))
}
is_transport_draw<-function(p,N,seed) {
 set.seed(seed);d<-length(p$mean);K<-length(p$components);u<-runif(N)
 c<-ifelse(u<.3,pmin(K,floor(u/.3*K)+1L),ifelse(u<.5,2L*K+2L,ifelse(u<.8,K+pmin(K,floor((u-.5)/.3*K)+1L),ifelse(u<.99,2L*K+1L,2L*K+3L))))
 y<-matrix(rnorm(N*d),N,d);z<-y
 for(j in 1:(2L*K+2L)) {
  ii<-which(c==j);if(!length(ii))next
  g<-if(j<=2L*K)p$components[[(j-1L)%%K+1L]]else if(j==2L*K+1L)p$global else p$tilted
  v<-y[ii,,drop=FALSE]
  if(j>K&&j<=2L*K)v<-v*matrix(sqrt(5/rchisq(length(ii)*d,7)),length(ii),d)
  v<-v%*%t(g$lower)
  if(j==2L*K+1L)v<-v*sqrt(5/rchisq(length(ii),7))
  z[ii,]<-is_transport_to_z(sweep(v,2,g$mean,"+"),p$transport)
 }
 z
}
is_transport_fit<-function(spec,mode,runtime,seed) {
 # Independent data-only chains train a proposal; their draws never enter estimates.
 sp<-c(spec,list(steps=15L,growth=0,growth_sd=0));d<-spec$q
 Q<-2*mode$covariance;lower<-t(chol(Q));precision<-solve(Q);chains<-list();evals<-0
 for(j in 1:4) {
  set.seed(seed+j);initial<-as.numeric(mode$mean+lower%*%rnorm(d,sd=1.5))
  set.seed(seed+100L+j)
  ch<-runtime$pilot_core$ref_preconditioned_ess_chain_cpp(sp,initial,2000L,60000L,seed+100L+j,mode$mean,lower,precision)
  chains[[j]]<-ch$draws;evals<-evals+ch$likelihood_evaluations
 }
 z<-do.call(rbind,chains)
 v<-c(colMeans(spec$lower),rep(0,d-spec$n));v<-v/sqrt(sum(v*v))
 rotation<-qr.Q(qr(cbind(v,diag(d))))[,seq_len(d),drop=FALSE]
 y<-z%*%rotation;p<-list(rotation=rotation,center=mean(y[,1]),scale=sd(y[,1]))
 p$knots<-as.numeric(quantile((y[,1]-p$center)/p$scale,c(.1,.3,.5,.7,.9)))
 p$boundary<-as.numeric(quantile((y[,1]-p$center)/p$scale,c(.001,.999)))
 D<-is_transport_design(y[,1],p);p$regression<-qr.coef(qr(D),y[,-1,drop=FALSE])
 y[,-1]<-y[,-1,drop=FALSE]-D%*%p$regression
 cuts<-quantile(y[,1],seq(0,1,length.out=17));bin<-findInterval(y[,1],cuts,all.inside=TRUE)
 Cglobal<-cov(y[,-1,drop=FALSE])
 lowers<-lapply(1:16,function(k)t(chol(.98*cov(y[bin==k,-1,drop=FALSE])+.02*Cglobal+diag(1e-6,d-1L))))
 p$conditional_covariance<-list(cuts=cuts,lower=lowers)
 for(k in 1:16) {
  ii<-which(bin==k);y[ii,-1]<-t(forwardsolve(lowers[[k]],t(y[ii,-1,drop=FALSE])))
 }
 probs<-sort(unique(c(.0001,pnorm(seq(-3.5,3.5,length.out=43)),.9999)))
 p$marginal<-lapply(seq_len(d),function(j)list(x=as.numeric(quantile(y[,j],probs)),u=qnorm(probs)))
 for(j in seq_len(d))y[,j]<-is_piecewise(y[,j],p$marginal[[j]]$x,p$marginal[[j]]$u)$value
 p$layers<-list()
 for(layer in 1:8) {
  R<-eigen(cov(y),symmetric=TRUE)$vectors;y<-y%*%R
  marginal<-lapply(seq_len(d),function(j)list(x=as.numeric(quantile(y[,j],probs)),u=qnorm(probs)))
  for(j in seq_len(d))y[,j]<-is_piecewise(y[,j],marginal[[j]]$x,marginal[[j]]$u)$value
  p$layers[[layer]]<-list(rotation=R,marginal=marginal)
 }
 C<-cov(y);mu<-colMeans(y);cuts<-quantile(y[,1],seq(0,1,length.out=9));group<-findInterval(y[,1],cuts,all.inside=TRUE)
 S<-1.1*C+diag(1e-7,d)
 components<-list(list(mean=mu,lower=t(chol(S)),covariance=S))
 p0<-mode;p0$transport<-p;p0$components<-components;p0$global<-list(mean=mu,lower=t(chol(C)),covariance=C)
 # A fixed proposal-only tilt covers the birth and birth-squared integrands.
 h<-.01;plus<-minus<-matrix(mu,d,d,byrow=TRUE);diag(plus)<-diag(plus)+h;diag(minus)<-diag(minus)-h
 raw_direction<-c(colMeans(spec$lower),rep(0,d-spec$n))
 gradient<-as.numeric((is_transport_to_z(plus,p)-is_transport_to_z(minus,p))%*%raw_direction)/(2*h)
 direction<-as.numeric(S%*%gradient)/sqrt(as.numeric(crossprod(gradient,S%*%gradient)))
 p0$tilted<-list(mean=mu+1.5*direction,lower=t(chol(S)),covariance=S,shift_sd=1.5)
 p0$adapted<-TRUE;p0$pilot<-list(method="independent_data_only_ESS_transport",samples=nrow(z),
  ess=NA_real_,pareto_k=NA_real_,seed=seed,likelihood_evaluations=evals,
  baseline_rhat=posterior::rhat(sapply(chains,function(a)as.numeric(a%*%v))))
 p0
}
