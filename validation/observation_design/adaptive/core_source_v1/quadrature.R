# Positive vector-valued Gauss-Kronrod integration in the original prior CDF.
# QK15 constants: QUADPACK, as published in SciPy's __quadpack.c.
adaptive_gk15_rule<-function() {
  x<-c(.991455371120812639206854697526329,.949107912342758524526189684047851,
    .864864423359769072789712788640926,.741531185599394439863864773280788,
    .586087235467691130294144838258730,.405845151377397166906606412076961,
    .207784955007898467600689403773245)
  k<-c(.022935322010529224963732008058970,.063092092629978553290700663189204,
    .104790010322250183839876322541518,.140653259715525918745189590510238,
    .169004726639267902826583426598550,.190350578064785409913256402421014,
    .204432940075298892414161999234649)
  g<-numeric(15L)
  g[c(2,4,6,8,10,12,14)]<-c(.129484966168869693270611432679082,
    .279705391489276667901467771423780,.381830050505118944950369775488975,
    .417959183673469387755102040816327,.381830050505118944950369775488975,
    .279705391489276667901467771423780,.129484966168869693270611432679082)
  list(x=c(-x,0,rev(x)),k=c(k,.209482141084727828012999174891714,rev(k)),g=g)
}

adaptive_error<-function(values,k,g) {
  values<-as.matrix(values)
  ik<-colSums(values*k);ig<-colSums(values*g)
  absolute<-colSums(abs(values)*k)
  asc<-colSums(abs(sweep(values,2,ik/sum(k),"-"))*k)
  error<-abs(ik-ig)
  ii<-asc>0 & error>0
  error[ii]<-asc[ii]*pmin(1,(200*error[ii]/asc[ii])^1.5)
  pmax(error,50*.Machine$double.eps*absolute)
}

adaptive_prior_integrate<-function(evaluate,prior_sd,weight,
    tolerances=c(normalizer=1e-4,mean=2e-4,covariance=.001,contrast=.001),
    max_panels=64L,max_evaluations=1000L,trace=FALSE) {
  stopifnot(prior_sd>0,abs(sum(weight)-1)<1e-10)
  n<-length(weight);A<-diag(n)-matrix(weight,n,n,byrow=TRUE)
  rule<-adaptive_gk15_rule();cache<-new.env(parent=emptyenv())
  evaluated<-list();history<-list();tail_history<-list();evaluations<-0L
  panel<-function(lower,upper) {
    h<-(upper-lower)/2;u<-(lower+upper)/2+h*rule$x
    if(any(u<=0|u>=1))stop("Prior-CDF nodes reached floating-point boundary; no clipping allowed")
    keys<-sprintf("%.17g",u)
    for(j in order(abs(qnorm(u))))if(!exists(keys[j],cache,inherits=FALSE)) {
      if(evaluations>=max_evaluations)stop("Adaptive node limit reached without a certificate")
      value<-prior_sd*qnorm(u[j]);piece<-evaluate(value)
      stopifnot(length(piece$mean)==n,all(is.finite(piece$mean)),
        all(is.finite(piece$covariance)),is.finite(piece$log_evidence))
      piece$intercept<-value;piece$prior_cdf<-u[j]
      assign(keys[j],piece,cache);evaluations<<-evaluations+1L
      evaluated[[evaluations]]<<-data.frame(node=evaluations,u=u[j],intercept=value,
        log_evidence=piece$log_evidence)
      if(trace)cat("Adaptive node",evaluations,"a",format(value,digits=8),"\n")
    }
    list(lower=lower,upper=upper,keys=keys,k=h*rule$k,g=h*rule$g)
  }
  summarize<-function(panels) {
    all_keys<-unlist(lapply(panels,`[[`,"keys"))
    pieces<-lapply(all_keys,get,envir=cache,inherits=FALSE)
    ll<-vapply(pieces,`[[`,numeric(1),"log_evidence");scale<-max(ll);L<-exp(ll-scale)
    k<-unlist(lapply(panels,`[[`,"k"));wk<-k*L;Z<-sum(wk)
    current<-baseline_combine(pieces,wk/Z)
    centered<-lapply(pieces,function(p)p$mean-current$mean)
    X<-Map(function(p,d)p$covariance+tcrossprod(d),pieces,centered)
    RX<-lapply(X,function(x)A%*%x%*%t(A))
    dmat<-do.call(rbind,centered);Xmat<-do.call(rbind,lapply(X,as.numeric));Rmat<-do.call(rbind,lapply(RX,as.numeric))
    errors<-list();offset<-0L
    for(i in seq_along(panels)) {
      p<-panels[[i]];ii<-offset+seq_along(p$keys);offset<-max(ii)
      eZ<-adaptive_error(matrix(L[ii],ncol=1),p$k,p$g)
      em<-adaptive_error(dmat[ii,,drop=FALSE]*L[ii],p$k,p$g)
      ec<-matrix(adaptive_error(Xmat[ii,,drop=FALSE]*L[ii],p$k,p$g),n,n)
      er<-matrix(adaptive_error(Rmat[ii,,drop=FALSE]*L[ii],p$k,p$g),n,n)
      errors[[i]]<-list(Z=eZ,mean=em,covariance=ec,contrast=er)
    }
    combine_error<-function(indices) {
      eZ<-sum(vapply(errors[indices],`[[`,numeric(1),"Z"))
      em<-Reduce(`+`,lapply(errors[indices],`[[`,"mean"))
      ec<-Reduce(`+`,lapply(errors[indices],`[[`,"covariance"))
      er<-Reduce(`+`,lapply(errors[indices],`[[`,"contrast"))
      denominator<-Z-eZ
      if(denominator<=0)return(c(normalizer=Inf,mean=Inf,covariance=Inf,contrast=Inf))
      em<-em/denominator
      ec<-(ec+abs(current$covariance)*eZ)/denominator+tcrossprod(em)
      RC<-A%*%current$covariance%*%t(A);re<-as.numeric(abs(A)%*%em)
      er<-(er+abs(RC)*eZ)/denominator+tcrossprod(re)
      c(normalizer=eZ/Z,mean=max(em/pmax(1,abs(current$mean))),
        covariance=norm(ec,"F")/norm(current$covariance,"F"),contrast=norm(er,"F")/norm(RC,"F"))
    }
    err<-combine_error(seq_along(panels))
    endpoint<-which(vapply(panels,function(p)p$lower==0||p$upper==1,logical(1)))
    tail_error<-combine_error(endpoint)
    priority<-vapply(seq_along(panels),function(i)max(combine_error(i)/tolerances),numeric(1))
    current$log_evidence<-scale+log(Z)
    list(current=current,error=err,tail_error=tail_error,priority=priority,
      endpoint=endpoint,tail_prior_mass=sum(vapply(panels[endpoint],function(p)p$upper-p$lower,numeric(1))),
      tail_posterior_mass=sum(vapply(endpoint,function(i) {
        ii<-((i-1L)*15L)+seq_len(15L);sum(wk[ii])
      },numeric(1)))/Z)
  }
  panels<-list(panel(.05,.95),panel(0,.05),panel(.95,1))
  tail_previous<-NULL;tail_confirmed<-FALSE;passed<-FALSE
  for(iteration in seq_len(max_panels)) {
    stopifnot(abs(sum(vapply(panels,function(p)p$upper-p$lower,numeric(1)))-1)<1e-14)
    s<-summarize(panels)
    tm<-tc<-tr<-Inf
    if(!is.null(tail_previous)) {
      tm<-max(abs(s$current$mean-tail_previous$mean)/pmax(1,abs(s$current$mean)))
      tc<-norm(s$current$covariance-tail_previous$covariance,"F")/norm(s$current$covariance,"F")
      tr<-norm(A%*%(s$current$covariance-tail_previous$covariance)%*%t(A),"F")/
        norm(A%*%s$current$covariance%*%t(A),"F")
      tail_confirmed<-isTRUE(tm<.001&&tc<.005&&tr<.005)
      tail_history[[length(tail_history)+1L]]<-data.frame(iteration=iteration,mean_change=tm,
        covariance_change=tc,contrast_covariance_change=tr,passed=tail_confirmed)
      tail_previous<-NULL
    }
    estimated_pass<-all(s$error<=tolerances)
    passed<-estimated_pass&&tail_confirmed
    history[[iteration]]<-data.frame(iteration=iteration,evaluations=evaluations,active_nodes=15L*length(panels),
      panels=length(panels),normalizer_error=s$error["normalizer"],mean_error=s$error["mean"],
      covariance_error=s$error["covariance"],contrast_covariance_error=s$error["contrast"],
      tail_normalizer_error=s$tail_error["normalizer"],tail_mean_error=s$tail_error["mean"],
      tail_covariance_error=s$tail_error["covariance"],tail_contrast_covariance_error=s$tail_error["contrast"],
      tail_prior_mass=s$tail_prior_mass,tail_posterior_mass=s$tail_posterior_mass,
      tail_confirmed=tail_confirmed,passed=passed)
    if(trace)cat("Adaptive iteration",iteration,"nodes",evaluations,"errors",format(s$error,digits=3),"tail",tail_confirmed,"\n")
    if(passed)break
    split<-if(estimated_pass&&!tail_confirmed)s$endpoint else which.max(s$priority)
    if(length(panels)+length(split)>max_panels)break
    if(estimated_pass&&!tail_confirmed)tail_previous<-s$current
    replacement<-list()
    for(i in seq_along(panels)) {
      p<-panels[[i]]
      if(i%in%split) {
        mid<-(p$lower+p$upper)/2
        replacement<-c(replacement,list(panel(p$lower,mid),panel(mid,p$upper)))
      } else replacement<-c(replacement,list(p))
    }
    panels<-replacement
  }
  current<-s$current
  order<-order(vapply(current$components,`[[`,numeric(1),"intercept"),decreasing=TRUE)
  current$components<-current$components[order];current$weights<-current$weights[order]
  c(current,list(status=if(passed)"passed_quadrature_refinement"else"adaptive_not_converged",
    history=do.call(rbind,history),tail_history=if(length(tail_history))do.call(rbind,tail_history)else data.frame(),
    evaluated=do.call(rbind,evaluated),partition=do.call(rbind,lapply(panels,function(p)data.frame(lower=p$lower,upper=p$upper))),
    error_estimates=s$error,tail_error_estimates=s$tail_error,
    computational_method="adaptive_prior_cdf_gk15",prior_sd=prior_sd))
}
