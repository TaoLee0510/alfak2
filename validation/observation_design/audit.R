#!/usr/bin/env Rscript
# Post-batch checks of saved numerical distributions and immutable inputs.
# This performs no fits and does not change the prespecified analysis rules.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
out<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(repo,"validation/observation_design/results"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
status<-read.csv(file.path(out,"all_fit_status.csv"))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
stopifnot(nrow(status)==60L,all(status$status=="passed_quadrature_refinement"))
rows<-list()
for(i in seq_len(nrow(status))) {
  s<-status[i,]
  if("result_path"%in%names(status))path<-s$result_path else {
    paths<-file.path(out,paste0("fits_",0:1),paste0(s$case_id,"_result.rds"))
    paths<-paths[file.exists(paths)];stopifnot(length(paths)==1L);path<-paths
  }
  saved<-readRDS(path)
  fit<-saved$fit;w<-fit$weights;components<-fit$components
  stopifnot(all(is.finite(w)),all(w>=0),any(w>0),abs(sum(w)-1)<1e-12,length(w)==length(components))
  means<-lapply(components,function(c)exp(c$latent_mean+diag(c$latent_covariance)/2))
  covs<-Map(function(m,c)tcrossprod(m)*expm1(c$latent_covariance),means,components)
  mu<-Reduce(`+`,Map(`*`,means,w))
  second<-Reduce(`+`,Map(function(m,C,p)p*(C+tcrossprod(m)),means,covs,w))
  covariance<-second-tcrossprod(mu)
  mean_error<-max(abs(mu-fit$mean));cov_error<-max(abs(covariance-fit$covariance))
  stopifnot(mean_error<1e-9,cov_error<1e-8)
  component_error<-max(vapply(seq_along(w),function(j)
    max(abs(means[[j]]-components[[j]]$mean),abs(covs[[j]]-components[[j]]$covariance)),numeric(1)))
  stopifnot(component_error<1e-9)
  case_row<-manifest[match(s$case_id,manifest$case_id),]
  stopifnot(unname(round2_sha256(case_row$file))==case_row$sha256,
    identical(saved$key$dataset,case_row$sha256))
  case<-readRDS(case_row$file)
  A<-diag(length(mu))-matrix(case$reference_weight,length(mu),length(mu),byrow=TRUE)
  expected<-as.numeric(A%*%mu)
  stopifnot(max(abs(expected-saved$scored$nodes$pred_relative))<1e-10)
  all_rows<-saved$scored$metrics[saved$scored$metrics$distance==-1L,]
  decomposition_error<-max(abs(all_rows$mse-all_rows$shape_mse-all_rows$offset_squared))
  stopifnot(decomposition_error<1e-12)
  h<-tail(fit$history,1)
  trusted_status<-if(identical(fit$computational_method,"direct_joint_mode_hessian"))
    "trusted_direct_joint_hessian" else "TMB_sdreport"
  stopifnot(h$passed,h$mean_change<.001,h$covariance_change<.005,
    h$contrast_covariance_change<.005,
    all(fit$conditional_status$status==trusted_status),
    all(fit$conditional_status$forward_error<=1e-5))
  if(isTRUE(fit$evidence_reconstructed)) {
    stopifnot(nrow(saved$audit)==nrow(fit$conditional_status),
      all(saved$audit$max_joint_gradient<=.001),all(saved$audit$latent_covariance_error<1e-8),
      all(is.finite(saved$audit$reconstructed_objective)))
  } else stopifnot(!is.null(fit$full),fit$full$diagnostics$backend_error<=1e-8)
  rows[[s$case_id]]<-data.frame(case_id=s$case_id,design=s$design,
    mixture_mean_error=mean_error,mixture_covariance_error=cov_error,
    component_moment_error=component_error,global_mse_decomposition_error=decomposition_error,
    final_order=h$order,mean_change=h$mean_change,covariance_change=h$covariance_change,
    contrast_covariance_change=h$contrast_covariance_change,
    maximum_forward_error=max(fit$conditional_status$forward_error),
    minimum_component_weight=min(w),maximum_component_weight=max(w),
    evidence_reconstructed=isTRUE(fit$evidence_reconstructed),
    elapsed_seconds=s$elapsed_seconds,within_300_seconds=s$elapsed_seconds<=300)
}
audit<-do.call(rbind,rows)
write.csv(audit,file.path(out,"posterior_integrity_audit.csv"),row.names=FALSE)
source_files<-c(file.path(out,paste0("fits_",0:1),"source_sha256.csv"),
  file.path(list.files(file.path(out,"cases"),full.names=TRUE),"source_sha256.csv"),
  file.path(list.files(file.path(repo,"validation/observation_design/results"),"^fits_bounded_child_",full.names=TRUE),"source_sha256.csv"))
source_files<-unique(normalizePath(source_files[file.exists(source_files)]))
source_rows<-do.call(rbind,lapply(source_files,function(p) {
  z<-read.csv(p);z$manifest_path<-p;z
}))
current_hash<-round2_sha256(unique(source_rows$file))
source_rows$current_sha256<-unname(current_hash[source_rows$file])
source_rows$unchanged<-source_rows$sha256==source_rows$current_sha256
write.csv(source_rows,file.path(out,"consumed_source_integrity.csv"),row.names=FALSE)
stopifnot(all(source_rows$unchanged))
cat(nrow(audit),"saved posteriors passed mixture, reference, refinement and input-integrity checks\n")
