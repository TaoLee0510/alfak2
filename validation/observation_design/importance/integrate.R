is_quantile<-function(x,w,probs=c(.025,.1,.25,.5,.75,.9,.975)) {
  o<-order(x);cw<-pmin(1,cumsum(w[o]));cw[length(cw)]<-1
  x[o[pmax(1L,pmin(length(x),findInterval(probs,c(0,cw),left.open=TRUE)))]]
}

is_observables<-function(z,spec,weight,edges) {
  n<-spec$n;g<-sweep(z[,seq_len(n),drop=FALSE]%*%t(spec$lower),2,spec$offset,"+");b<-exp(g)
  out<-list(log_birth=g,birth=b,relative=b-as.numeric(b%*%weight),centered=b-rowMeans(b),
    edge=b[,edges[,2],drop=FALSE]-b[,edges[,1],drop=FALSE])
  if(spec$q>n)out$eta_contrast<-2*z[,n+seq_len(spec$q-n),drop=FALSE]
  out
}

is_summarize_legacy<-function(batches,spec,weight,edges,include_quantiles=FALSE) {
  z<-do.call(rbind,lapply(batches,`[[`,"z"));lw<-unlist(lapply(batches,`[[`,"logw"));w<-is_weights(lw)
  group<-rep(seq_along(batches),vapply(batches,function(b)nrow(b$z),integer(1)))
  n<-spec$n;g<-sweep(z[,seq_len(n),drop=FALSE]%*%t(spec$lower),2,spec$offset,"+");birth<-exp(g)
  rows<-list();covariances<-list();cov_se<-list();quantiles<-list()
  types<-c("log_birth","birth","relative","centered","edge",if(spec$q>n)"eta_contrast"else character())
  for(type in types) {
    x<-switch(type,log_birth=g,birth=birth,relative=birth-as.numeric(birth%*%weight),centered=birth-rowMeans(birth),
      edge=birth[,edges[,2],drop=FALSE]-birth[,edges[,1],drop=FALSE],eta_contrast=2*z[,n+seq_len(spec$q-n),drop=FALSE])
    mu<-colSums(x*w);X<-sweep(x,2,mu,"-")
    C<-crossprod(X*w,X);sd<-sqrt(pmax(0,diag(C)))
    mse<-sqrt(colSums(X*X*(w*w)))
    psi<-sweep(X*X,2,sd*sd,"-");ii<-sd>1e-12
    psi[,ii]<-sweep(psi[,ii,drop=FALSE],2,2*sd[ii],"/");psi[,!ii]<-0
    sse<-sqrt(colSums(psi*psi*(w*w)))
    bmean<-bsd<-matrix(NA_real_,length(batches),ncol(x))
    for(j in seq_along(batches)) {
      sel<-group==j;ww<-w[sel]/sum(w[sel]);bmean[j,]<-colSums(x[sel,,drop=FALSE]*ww)
      bsd[j,]<-sqrt(colSums(sweep(x[sel,,drop=FALSE],2,bmean[j,],"-")^2*ww))
    }
    mse<-pmax(mse,apply(bmean,2,stats::sd)/sqrt(length(batches)))
    sse<-pmax(sse,apply(bsd,2,stats::sd)/sqrt(length(batches)))
    rows[[type]]<-data.frame(type=type,index=seq_along(mu),mean=mu,sd=sd,mcse_mean=mse,mcse_sd=sse,
      mean_mcse_ratio=ifelse(sd>1e-12,mse/sd,0),sd_mcse_ratio=ifelse(sd>1e-12,sse/sd,0))
    covariances[[type]]<-C
    if(type%in%c("birth","relative","log_birth")) {
      V<-crossprod((X*X)*(w*w),X*X)-2*C*crossprod(X*(w*w),X)+C*C*sum(w*w)
      V[V<0]<-0;cov_se[[type]]<-sqrt(V)
    }
    if(include_quantiles&&type%in%c("birth","relative","centered","edge"))
      quantiles[[type]]<-t(vapply(seq_len(ncol(x)),function(j)is_quantile(x[,j],w),numeric(7)))
    rm(list=intersect(c("x","X","psi","bmean","bsd","V"),ls()));invisible(gc())
  }
  moments<-do.call(rbind,rows);allk<-c(vapply(batches,function(b)is_khat(b$logw),numeric(1)),is_khat(lw))
  ess<-1/sum(w*w);finite<-all(is.finite(moments$mean))&&all(is.finite(moments$sd))
  flags<-c(weight_ess=ess>=10000,pareto_tail=all(is.finite(allk))&&max(allk)<=.5,
    mean_mcse=max(moments$mean_mcse_ratio)<=.005,sd_mcse=max(moments$sd_mcse_ratio)<=.01,finite=finite)
  b<-moments[moments$type=="birth",];g<-moments[moments$type=="log_birth",]
  list(status=if(all(flags))"passed_joint_importance"else"importance_not_certified",flags=flags,
    mean=b$mean,covariance=covariances$birth,latent_mean=g$mean,latent_covariance=covariances$log_birth,
    weights=w,log_weights=lw,samples=z,moments=moments,covariances=covariances,covariance_mcse=cov_se,
    quantiles=quantiles,ess=ess,pareto_k=allk,max_mean_mcse_ratio=max(moments$mean_mcse_ratio),
    max_sd_mcse_ratio=max(moments$sd_mcse_ratio),n_samples=nrow(z),
    log_evidence=max(lw)+log(mean(exp(lw-max(lw)))))
}

is_summarize<-function(batches,spec,weight,edges,include_quantiles=FALSE) {
  z<-do.call(rbind,lapply(batches,`[[`,"z"));lw<-unlist(lapply(batches,`[[`,"logw"));w<-is_weights(lw)
  n<-spec$n;N<-nrow(z);g<-sweep(z[,seq_len(n),drop=FALSE]%*%t(spec$lower),2,spec$offset,"+");birth<-exp(g)
  reference<-as.numeric(birth%*%weight);center<-rowMeans(birth)
  getx<-function(ii,type)switch(type,log_birth=g[ii,,drop=FALSE],birth=birth[ii,,drop=FALSE],
    relative=birth[ii,,drop=FALSE]-reference[ii],centered=birth[ii,,drop=FALSE]-center[ii],
    edge=birth[ii,edges[,2],drop=FALSE]-birth[ii,edges[,1],drop=FALSE],
    eta_contrast=2*z[ii,n+seq_len(spec$q-n),drop=FALSE])
  getcolumn<-function(j,type)switch(type,birth=birth[,j],relative=birth[,j]-reference,
    centered=birth[,j]-center,edge=birth[,edges[j,2]]-birth[,edges[j,1]])
  sizes<-vapply(batches,function(b)nrow(b$z),integer(1));starts<-c(0L,cumsum(sizes));J<-length(batches)
  blocks<-lapply(seq_len(J),function(j)split(starts[j]+seq_len(sizes[j]),ceiling(seq_len(sizes[j])/16384L)))
  rows<-covariances<-cov_se<-quantiles<-list();sw2<-sum(w*w)
  types<-c("log_birth","birth","relative","centered","edge",if(spec$q>n)"eta_contrast"else character())
  for(type in types) {
    d<-ncol(getx(1L,type));bmean<-matrix(0,J,d);wb<-numeric(J)
    for(j in seq_len(J))for(ii in blocks[[j]]) {
      x<-getx(ii,type);wb[j]<-wb[j]+sum(w[ii]);bmean[j,]<-bmean[j,]+colSums(x*w[ii])
    }
    mu<-colSums(bmean);bmean<-bmean/wb
    C<-V1<-V2<-matrix(0,d,d);E2<-E4<-numeric(d);bsd<-matrix(0,J,d)
    for(j in seq_len(J))for(ii in blocks[[j]]) {
      x<-getx(ii,type);X<-sweep(x,2,mu,"-");sq<-X*X;ww<-w[ii];ww2<-ww*ww
      C<-C+crossprod(X*ww,X);E2<-E2+colSums(sq*ww2);E4<-E4+colSums(sq*sq*ww2)
      bsd[j,]<-bsd[j,]+colSums(sweep(x,2,bmean[j,],"-")^2*ww)
      if(type%in%c("birth","relative","log_birth")) {
        V1<-V1+crossprod(sq*ww2,sq);V2<-V2+crossprod(X*ww2,X)
      }
    }
    sd<-sqrt(pmax(0,diag(C)));bsd<-sqrt(bsd/wb)
    mse<-pmax(sqrt(E2),apply(bmean,2,stats::sd)/sqrt(J))
    sse<-numeric(d);live<-sd>1e-12
    sse[live]<-sqrt(pmax(0,E4[live]-2*sd[live]^2*E2[live]+sd[live]^4*sw2))/(2*sd[live])
    sse<-pmax(sse,apply(bsd,2,stats::sd)/sqrt(J))
    rows[[type]]<-data.frame(type=type,index=seq_len(d),mean=mu,sd=sd,mcse_mean=mse,mcse_sd=sse,
      mean_mcse_ratio=ifelse(live,mse/sd,0),sd_mcse_ratio=ifelse(live,sse/sd,0))
    covariances[[type]]<-C
    if(type%in%c("birth","relative","log_birth")) {
      V<-V1-2*C*V2+C*C*sw2;V[V<0]<-0;cov_se[[type]]<-sqrt(V)
    }
    if(include_quantiles&&type%in%c("birth","relative","centered","edge")) {
      q<-matrix(NA_real_,d,7)
      for(k in seq_len(d)){q[k,]<-is_quantile(getcolumn(k,type),w);if(k%%8==0)invisible(gc())}
      quantiles[[type]]<-q
    }
    rm(x,X,sq);invisible(gc())
  }
  moments<-do.call(rbind,rows);allk<-c(vapply(batches,function(b)is_khat(b$logw),numeric(1)),is_khat(lw))
  ess<-1/sum(w*w);finite<-all(is.finite(moments$mean))&&all(is.finite(moments$sd))
  flags<-c(weight_ess=ess>=10000,pareto_tail=all(is.finite(allk))&&max(allk)<=.5,
    mean_mcse=max(moments$mean_mcse_ratio)<=.005,sd_mcse=max(moments$sd_mcse_ratio)<=.01,finite=finite)
  b<-moments[moments$type=="birth",];gg<-moments[moments$type=="log_birth",]
  list(status=if(all(flags))"passed_joint_importance"else"importance_not_certified",flags=flags,
    mean=b$mean,covariance=covariances$birth,latent_mean=gg$mean,latent_covariance=covariances$log_birth,
    weights=w,log_weights=lw,samples=z,moments=moments,covariances=covariances,covariance_mcse=cov_se,
    quantiles=quantiles,ess=ess,pareto_k=allk,max_mean_mcse_ratio=max(moments$mean_mcse_ratio),
    max_sd_mcse_ratio=max(moments$sd_mcse_ratio),n_samples=N,log_evidence=max(lw)+log(mean(exp(lw-max(lw)))))
}

is_integrate<-function(case,context,runtime,seed,trace=TRUE,save_directory=NULL) {
  started<-proc.time()[["elapsed"]];spec<-is_spec(case,context)
  weight<-case$reference_weight
  edges<-if(!is.null(case$target_edges))case$target_edges else case$edges
  if(is.null(edges)){D<-as.matrix(dist(context$graph$karyotypes,method="manhattan"));edges<-which(D==1&upper.tri(D),arr.ind=TRUE)}
  proposal<-is_mode_proposal(spec,runtime)
  proposal<-is_pilot_proposal(spec,proposal,runtime,seed+1L)
  if(!is.null(save_directory))saveRDS(proposal,file.path(save_directory,"proposal.rds"),compress=FALSE)
  if(trace)cat("Pilot ESS",proposal$pilot$ess,"k",proposal$pilot$pareto_k,"adapted",proposal$adapted,"\n")
  batches<-replicate(8,list(z=matrix(numeric(),0,spec$q),logw=numeric()),simplify=FALSE)
  history<-list();stage<-0L;oldN<-0L
  for(N in c(4096L,8192L,16384L,32768L,65536L,131072L)) {
    stage<-stage+1L
    for(j in seq_along(batches)) {
      s<-seed+j*100L+stage+10L;z<-is_draw_proposal(proposal,N-oldN,s)
      ll<-runtime$core$is_loglikelihood_cpp(z,spec)
      lw<-ll-.5*(spec$q*log(2*pi)+rowSums(z*z))-is_logq(z,proposal)
      stopifnot(all(is.finite(lw)))
      batches[[j]]$z<-rbind(batches[[j]]$z,z);batches[[j]]$logw<-c(batches[[j]]$logw,lw)
    }
    oldN<-N;if(exists("current",inherits=FALSE))rm(current);invisible(gc())
    current<-is_summarize(batches,spec,weight,edges,FALSE)
    history[[stage]]<-data.frame(stage=stage,samples_per_batch=N,total_samples=current$n_samples,ess=current$ess,
      max_pareto_k=max(current$pareto_k),max_mean_mcse_ratio=current$max_mean_mcse_ratio,
      max_sd_mcse_ratio=current$max_sd_mcse_ratio,passed=all(current$flags),elapsed=proc.time()[["elapsed"]]-started)
    if(trace){print(history[[stage]],row.names=FALSE);flush.console()}
    if(!is.null(save_directory)) {
      write.csv(do.call(rbind,history),file.path(save_directory,"refinement.csv"),row.names=FALSE)
      saveRDS(list(batches=batches,proposal=proposal,specification=spec,history=history),file.path(save_directory,"sample_checkpoint.rds"),compress=FALSE)
    }
    if(all(current$flags))break
  }
  rm(current);invisible(gc());current<-is_summarize(batches,spec,weight,edges,TRUE)
  current$proposal<-proposal;current$history<-do.call(rbind,history);current$specification<-spec
  current$elapsed_seconds<-proc.time()[["elapsed"]]-started
  current$numerical_seed<-seed;current$computational_method<-"raw_joint_importance_original_target"
  current
}
