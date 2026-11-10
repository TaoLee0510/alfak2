design_draws <- function(fit,n_draws,seed) {
  set.seed(seed);u<-runif(n_draws)
  index<-pmin(length(fit$weights),findInterval(u,cumsum(fit$weights))+1L)
  n<-length(fit$mean);Z<-matrix(rnorm(n*n_draws),nrow=n)
  answer<-matrix(NA_real_,n,n_draws)
  for(k in unique(index)) {
    j<-which(index==k);c<-fit$components[[k]]
    E<-eigen((c$latent_covariance+t(c$latent_covariance))/2,symmetric=TRUE)
    if(min(E$values)< -1e-7)stop("Indefinite design component")
    L<-sweep(E$vectors,2,sqrt(pmax(0,E$values)),"*")
    answer[,j]<-exp(sweep(L%*%Z[,j,drop=FALSE],1,c$latent_mean,"+"))
  }
  answer
}

design_score <- function(case,fit,n_draws=40000L) {
  n<-length(case$truth);w<-case$reference_weight
  valid<-identical(fit$status,"passed_quadrature_refinement")
  pred<-se<-rep(NA_real_,n)
  intervals<-list()
  if(valid) {
    A<-diag(n)-matrix(w,n,n,byrow=TRUE)
    pred<-as.numeric(A%*%fit$mean)
    se<-sqrt(pmax(0,diag(A%*%fit$covariance%*%t(A))))
    draws<-design_draws(fit,n_draws,case$sampling_seeds[1]+61000000L)
    draws<-sweep(draws,2,as.numeric(crossprod(w,draws)),"-")
    for(level in c(.5,.8,.95))intervals[[as.character(level)]]<-t(apply(draws,1,quantile,c((1-level)/2,(1+level)/2)))
  }
  nodes<-data.frame(case_id=case$id,karyotype=case$target_labels,truth_relative=case$truth_relative,
    pred_relative=pred,sd_relative=se,baseline_distance=case$baseline_distance,
    planned_distance=case$planned_distance,actual_distance=case$actual_distance)
  error<-pred-case$truth_relative;centered<-error-mean(error)
  edges<-case$target_edges;rows<-list()
  cor_safe<-function(a,b)if(length(a)>=5L&&all(is.finite(a))&&sd(a)>1e-12&&sd(b)>1e-12)cor(a,b,method="spearman")else NA_real_
  for(view in c("baseline_distance","planned_distance","actual_distance")) {
    d<-case[[view]];edge_d<-pmax(d[edges[,1]],d[edges[,2]])
    for(distance in c(-1L,sort(unique(d)))) {
      ii<-if(distance<0)seq_len(n)else which(d==distance)
      ee<-if(distance<0)seq_len(nrow(edges))else which(edge_d==distance)
      gradient<-error[edges[ee,2]]-error[edges[ee,1]]
      x<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,
        family=case$family,replicate=case$biological_replicate,design=case$design,
        view=view,distance=distance,n_nodes=length(ii),n_edges=length(ee),
        prediction_fraction=mean(is.finite(pred[ii])),
        rmse=sqrt(mean(error[ii]^2)),mse=mean(error[ii]^2),
        shape_rmse=sqrt(mean(centered[ii]^2)),shape_mse=mean(centered[ii]^2),
        offset_squared=mean(error)^2,edge_rmse=if(length(ee))sqrt(mean(gradient^2))else NA_real_,
        spearman=cor_safe(pred[ii],case$truth_relative[ii]),
        amplitude_ratio=if(length(ii)>1&&sd(case$truth_relative[ii])>0)sd(pred[ii])/sd(case$truth_relative[ii])else NA_real_)
      for(level in c(.5,.8,.95)) {
        suffix<-sprintf("%02d",round(100*level));q<-intervals[[as.character(level)]]
        lo<-hi<-rep(NA_real_,length(ii))
        if(!is.null(q)){lo<-q[ii,1];hi<-q[ii,2]}
        alpha<-1-level;t<-case$truth_relative[ii]
        x[[paste0("interval_fraction",suffix)]]<-mean(is.finite(lo)&is.finite(hi))
        x[[paste0("coverage",suffix)]]<-mean(t>=lo&t<=hi)
        x[[paste0("score",suffix)]]<-mean(hi-lo+2/alpha*pmax(lo-t,0)+2/alpha*pmax(t-hi,0))
        x[[paste0("width",suffix)]]<-mean(hi-lo)
      }
      rows[[length(rows)+1L]]<-x
    }
  }
  list(nodes=nodes,metrics=do.call(rbind,rows))
}
