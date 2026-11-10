#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design");r<-file.path(root,"results")
out<-file.path(r,"adaptive_complete_v1")
source(file.path(repo,"validation/accuracy_round2/common.R"));source(file.path(r,"adaptive_source_v1/quadrature.R"))
index<-read.csv(file.path(out,"result_index.csv"));manifest<-read.csv(file.path(out,"dataset_manifest.csv"));rows<-list()
for(i in seq_len(nrow(index))) {
  item<-index[i,];saved<-readRDS(item$result_path);fit<-saved$fit;n<-length(fit$mean)
  case_row<-manifest[match(item$case_id,manifest$case_id),];case<-readRDS(case_row$file)
  stopifnot(fit$status=="passed_quadrature_refinement",saved$key$dataset==case_row$sha256,
    unname(round2_sha256(case_row$file))==case_row$sha256,all(is.finite(fit$mean)),all(is.finite(fit$covariance)))
  means<-lapply(fit$components,function(c)exp(c$latent_mean+diag(c$latent_covariance)/2))
  covs<-Map(function(m,c)tcrossprod(m)*expm1(c$latent_covariance),means,fit$components)
  w<-fit$weights;stopifnot(all(w>=0),abs(sum(w)-1)<1e-12)
  mu<-Reduce(`+`,Map(`*`,means,w))
  covariance<-Reduce(`+`,Map(function(m,C,p)p*(C+tcrossprod(m-mu)),means,covs,w))
  me<-max(abs(mu-fit$mean));ce<-max(abs(covariance-fit$covariance))
  stopifnot(me<1e-9,ce<1e-8)
  A<-diag(n)-matrix(case$reference_weight,n,n,byrow=TRUE)
  stopifnot(max(abs(saved$scored$nodes$pred_relative-as.numeric(A%*%mu)))<1e-10)
  full<-saved$scored$metrics[saved$scored$metrics$distance==-1,]
  de<-max(abs(full$mse-full$shape_mse-full$offset_squared));stopifnot(de<1e-12)
  h<-tail(fit$history,1);stopifnot(h$passed)
  adaptive<-identical(item$fresh_recovery,TRUE)
  prior_weight_error<-NA_real_
  if(adaptive) {
    stopifnot(all(fit$error_estimates<=c(normalizer=1e-4,mean=2e-4,covariance=.001,contrast=.001)),
      h$tail_confirmed,tail(fit$tail_history$passed,1),all(fit$conditional_status$status=="trusted_direct_joint_hessian"),
      all(fit$conditional_status$gradient<=.001))
    p<-fit$partition[order(fit$partition$lower),]
    stopifnot(p$lower[1]==0,tail(p$upper,1)==1,max(abs(head(p$upper,-1)-tail(p$lower,-1)))<1e-14)
    rule<-adaptive_gk15_rule();u<-q<-numeric()
    for(j in seq_len(nrow(p))) {
      half<-(p$upper[j]-p$lower[j])/2
      u<-c(u,(p$upper[j]+p$lower[j])/2+half*rule$x);q<-c(q,half*rule$k)
    }
    coords<-vapply(fit$components,`[[`,numeric(1),"integration_coordinate")
    order<-match(sprintf("%.17g",coords),sprintf("%.17g",u));stopifnot(!anyNA(order))
    ell<-vapply(fit$components,`[[`,numeric(1),"log_evidence")
    a<-vapply(fit$components,`[[`,numeric(1),"intercept")
    jac<-dnorm(a,0,fit$prior_sd,log=TRUE)+log(fit$prior_sd)-log(coords)-log1p(-coords)
    lw<-log(q[order])+ell+jac;ww<-exp(lw-max(lw));ww<-ww/sum(ww)
    prior_weight_error<-max(abs(ww-w));stopifnot(prior_weight_error<1e-12)
  } else stopifnot(h$mean_change<.001,h$covariance_change<.005,h$contrast_covariance_change<.005,
    all(fit$conditional_status$status=="TMB_sdreport"))
  stopifnot(all(fit$conditional_status$forward_error<=1e-5))
  rows[[i]]<-data.frame(case_id=item$case_id,method=item$numerical_method,mean_error=me,covariance_error=ce,
    reference_error=max(abs(saved$scored$nodes$pred_relative-as.numeric(A%*%mu))),
    global_mse_decomposition_error=de,prior_jacobian_weight_error=prior_weight_error)
}
audit<-do.call(rbind,rows);write.csv(audit,file.path(out,"posterior_integrity_audit.csv"),row.names=FALSE)
source_files<-unique(index$source_manifest)
source_rows<-do.call(rbind,lapply(source_files,function(p){z<-read.csv(p);z$manifest_path<-p;z}))
hash<-round2_sha256(unique(source_rows$file));source_rows$current_sha256<-unname(hash[source_rows$file])
stopifnot(all(source_rows$sha256==source_rows$current_sha256))
write.csv(source_rows,file.path(out,"consumed_source_integrity.csv"),row.names=FALSE)
cat("60 complete posterior distributions passed independent moment, contrast, prior-Jacobian and source checks\n")
