#!/usr/bin/env Rscript
# Numerical compatibility only: no generating truths are used here.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-Sys.getenv("ALFAK2_ADAPTIVE_OUTPUT")
source(file.path(root,"results/computation_source_v4/score.R"))
old_status<-read.csv(file.path(root,"results/analysis_input/all_fit_status.csv"))
manifest<-read.csv(file.path(root,"results/dataset_manifest.csv"))
paths<-list.files(file.path(out,"primary"),"_result[.]rds$",recursive=TRUE,full.names=TRUE)
rows<-list()
for(path in paths) {
  new<-readRDS(path);i<-match(new$case_id,old_status$case_id)
  stopifnot(old_status$status[i]=="passed_quadrature_refinement",new$fit$status=="passed_quadrature_refinement")
  old<-readRDS(old_status$result_path[i])$fit;fit<-new$fit
  case<-readRDS(manifest$file[match(new$case_id,manifest$case_id)])
  n<-length(fit$mean);A<-diag(n)-matrix(case$reference_weight,n,n,byrow=TRUE)
  dm<-max(abs(fit$mean-old$mean)/pmax(1,abs(fit$mean)))
  dc<-norm(fit$covariance-old$covariance,"F")/norm(fit$covariance,"F")
  dr<-norm(A%*%(fit$covariance-old$covariance)%*%t(A),"F")/norm(A%*%fit$covariance%*%t(A),"F")
  nd<-A%*%design_draws(fit,160000L,835001L);od<-A%*%design_draws(old,160000L,835001L)
  cdf_ok<-TRUE;maxcdf<-0
  for(p in c(.025,.1,.25,.5,.75,.9,.975)) {
    quantile<-apply(od,1,stats::quantile,p)
    difference<-abs(rowMeans(nd<=quantile)-p)
    maxcdf<-max(maxcdf,difference)
    cdf_ok<-cdf_ok&&all(difference<=.01+3*sqrt(2*p*(1-p)/ncol(nd)))
  }
  rows[[new$case_id]]<-data.frame(case_id=new$case_id,mean_change=dm,covariance_change=dc,
    contrast_covariance_change=dr,max_cdf_difference=maxcdf,passed=dm<.001&&dc<.005&&dr<.005&&cdf_ok)
}
result<-do.call(rbind,rows);stopifnot(nrow(result)==3L)
write.csv(result,file.path(out,"compatibility.csv"),row.names=FALSE);print(result,row.names=FALSE)
if(!all(result$passed))stop("Numerical compatibility not established; require uniform recomputation")
