#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_IS_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"))
source(file.path(src,"model.R"));source(file.path(src,"integrate.R"));runtime<-is_load(repo,src)
oldsrc<-file.path(root,"results/mechanism_reference_source_v2")
source(file.path(oldsrc,"common.R"));oldruntime<-mechanism_load(repo,oldsrc)
cases<-readRDS(file.path(root,"results/mechanism_v2/frozen_cases.rds"))
context<-readRDS(file.path(root,"results/mechanism_v2/context.rds"))
library(testthat);reporter<-ListReporter$new();details<-list()
with_reporter(reporter,{
 test_that("proposal densities and weighted inverse CDF are normalized as specified",{
  set.seed(10500001);d<-5L;L<-t(chol(crossprod(matrix(rnorm(d*d),d,d))+diag(d)))
  p<-list(mean=rnorm(d),lower=L,covariance=tcrossprod(L));z<-matrix(rnorm(100*d),100,d)
  a<-.80*mvtnorm::dmvnorm(z,p$mean,p$covariance)+.19*mvtnorm::dmvt(sweep(z,2,p$mean,"-"),sigma=p$covariance*5/7,df=7,log=FALSE)+.01*mvtnorm::dmvnorm(z)
  expect_equal(is_logq(z,p),log(a),tolerance=1e-12)
  expect_equal(is_quantile(c(30,10,20),c(.2,.3,.5),c(.1,.3,.31,.8,.81)),c(10,10,20,20,30))
 })
 test_that("independent positive mean action matches exact and extreme two-node actions",{
  for(case in cases) {
   sp<-is_spec(case,context)
   for(j in 1:12) {
    set.seed(10501000+j);z<-rnorm(sp$q,sd=.75);g<-sp$offset+sp$lower%*%z[1:sp$n]
    eta<-sp$eta_basis%*%z[sp$n+seq_len(sp$n-1L)];p<-exp(eta-max(eta));p<-as.numeric(p/sum(p))
    A<-sweep(2*t(sp$probability)-diag(sp$n),2,exp(g),"*")
    exact<-oldruntime$independent$process_deterministic_mean(A,p,sp$dt)
    f<-runtime$core$is_forward_cpp(p,as.numeric(exp(g)),sp$probability,sp$dt,FALSE)
    expect_equal(f$pi1,as.numeric(exact/sum(exact)),tolerance=1e-11)
    expect_equal(f$log_growth,log(sum(exact)),tolerance=1e-10)
    expect_equal(runtime$core$is_loglikelihood_cpp(matrix(z,1),sp),oldruntime$cpp$ref_likelihood_cpp(z,case$specification)$log_likelihood,tolerance=1e-9)
   }
  }
  for(rates in list(c(.1,1),c(10,2),c(exp(12),exp(-12)),c(exp(-12),exp(12)))) {
   P<-matrix(c(.99,.01,.015,.985),2,2,byrow=TRUE);A<-sweep(2*t(P)-diag(2),2,rates,"*")
   tr<-sum(diag(A))/2;D<-A-diag(tr,2);delta<-sqrt(D[1,1]^2+D[1,2]*D[2,1]);t<-1.5
   # exp(A*t), divided by its largest eigenvalue exponential.
   M<-(diag(2)+D/delta)/2+exp(-2*delta*t)*(diag(2)-D/delta)/2
   v<-as.numeric(M%*%c(.3,.7));expected<-v/sum(v)
   f<-runtime$core$is_forward_cpp(c(.3,.7),rates,P,t,FALSE)
   expect_equal(f$pi1,expected,tolerance=1e-10)
   expect_equal(f$log_growth,(tr+delta)*t+log(sum(v)),tolerance=1e-9)
  }
  f<-runtime$core$is_forward_cpp(c(1,0),c(.4,1e6),diag(2),1.5,TRUE)
  expect_equal(f$pi1,c(1,0));expect_equal(f$log_growth,.6)
 })
 test_that("collapsed Gaussian mode, Hessian and evidence equal the original expanded model",{
  for(case in cases) {
   p<-is_mode_proposal(is_spec(case,context),runtime)
   f<-direct_joint_fit(case,context,oldruntime)
   C<-p$covariance[1:case$n,1:case$n];L<-t(chol(context$covariance))
   expect_equal(as.numeric(p$report$g),f$latent_fitness_mean,tolerance=1e-6)
   expect_equal(L%*%C%*%t(L),f$latent_fitness_covariance,tolerance=1e-6,ignore_attr=TRUE)
   evidence<-p$mode_objective-sum(log(diag(p$lower)))-length(p$mean)/2*log(2*pi)
   expect_equal(evidence,f$diagnostics$objective,tolerance=1e-7)
  }
 })
 test_that("raw joint importance recovers analytic Gaussian and lognormal prior",{
  c<-cases[[1]];c$y[]<-0;c$counts[]<-0
  fit<-is_integrate(c,context,runtime,10510000L)
  m<-fit$moments[fit$moments$type=="birth",];expected<-exp(context$prior$offset+diag(context$covariance)/2)
  C<-tcrossprod(expected)*expm1(context$covariance)
  expect_equal(fit$status,"passed_joint_importance")
  expect_true(all(abs(m$mean-expected)<=6*m$mcse_mean))
  expect_true(all(abs(fit$covariance-C)<=6*fit$covariance_mcse$birth))
  idx<-split(seq_len(fit$n_samples),rep(1:8,each=fit$n_samples/8))
  batches<-lapply(idx,function(ii)list(z=fit$samples[ii,,drop=FALSE],logw=fit$log_weights[ii]))
  legacy<-is_summarize_legacy(batches,fit$specification,c$reference_weight,c$edges,TRUE)
  expect_equal(fit$moments,legacy$moments,tolerance=1e-11)
  expect_equal(fit$covariances,legacy$covariances,tolerance=1e-11)
  expect_equal(fit$covariance_mcse,legacy$covariance_mcse,tolerance=1e-10)
  expect_equal(fit$quantiles,legacy$quantiles,tolerance=1e-12)
  details$prior<<-fit
 })
 test_that("transport density uses an exact inverse and Jacobian",{
  sp<-is_spec(cases[[1]],context);p<-is_transport_fit(sp,is_mode_proposal(sp,runtime),runtime,10515000L)
  z<-matrix(seq(-.5,.5,length.out=10*sp$q),10,sp$q);y<-is_transport_to_y(z,p$transport)
  expect_equal(unname(is_transport_to_z(y,p$transport)),unname(z),tolerance=1e-12,ignore_attr=TRUE)
  expect_equal(crossprod(p$transport$rotation),diag(sp$q),tolerance=1e-12)
  f<-function(a)as.numeric(is_transport_to_y(matrix(a,1),p$transport));h<-1e-7
  point<-as.numeric(is_transport_to_z(matrix(p$global$mean,1),p$transport))
  J<-sapply(1:sp$q,function(j){a<-point;b<-a;a[j]<-a[j]+h;b[j]<-b[j]-h;(f(a)-f(b))/(2*h)})
  expect_equal(log(abs(det(J))),is_transport_logjac(matrix(point,1),p$transport),tolerance=1e-8)
  terms<-sapply(p$components,function(c){D<-forwardsolve(c$lower,t(sweep(y,2,c$mean,"-")));(.3*mvtnorm::dmvnorm(y,c$mean,c$covariance)+.3*exp(colSums(matrix(dt(as.numeric(D)/sqrt(5/7),7,log=TRUE),sp$q))-sp$q/2*log(5/7)-sum(log(diag(c$lower)))))/length(p$components)})
  dens<-(rowSums(terms)+.19*mvtnorm::dmvt(sweep(y,2,p$global$mean,"-"),sigma=p$global$covariance*5/7,df=7,log=FALSE)+.2*mvtnorm::dmvnorm(y,p$tilted$mean,p$tilted$covariance))*exp(is_transport_logjac(z,p$transport))+.01*mvtnorm::dmvnorm(z)
  expect_equal(is_logq(z,p),log(dens),tolerance=1e-12)
 })
 test_that("raw joint importance matches independent two-node quadrature",{
  old<-readRDS(file.path(repo,"validation/posterior_reference_4to8/inputs/two_node_reference.rds"))
  ctx<-list(n=2L,prior=old$prior,covariance=old$aggregate_prior_covariance,
   eta_basis=matrix(c(1,-1)*old$eta_prior_sd/sqrt(2),2,1),
   operator=list(from0=0:1,to0=0:1,weight=c(1,1)))
  c<-list(id="two_node_count_validation",y=old$counts,times=c(0,old$dt),trajectory=c(0L,0L),
   reference_weight=c(.5,.5),edges=matrix(c(1L,2L),1,2))
  fit<-is_integrate(c,ctx,runtime,10520000L)
  exact<-tail(old$results$count_only$reference,1)[[1]]
  m<-fit$moments[c(which(fit$moments$type=="birth"),which(fit$moments$type=="relative")[1]),]
  expect_true(all(abs(m$mean-c(exact$mean,exact$contrast_mean))<=6*m$mcse_mean+3e-5))
  expect_true(all(abs(m$sd-c(sqrt(diag(exact$covariance)),exact$contrast_sd))<=6*m$mcse_sd+3e-5))
  expect_equal(fit$status,"passed_joint_importance");details$two_node<<-fit
  c$known_pi0<-c(1,0);sp<-is_spec(c,ctx);z<-matrix(c(.1,-.2),1)
  expect_equal(sp$q,2L)
  ll<-runtime$core$is_loglikelihood_cpp(z,sp)
  expect_equal(ll,sum(rowSums(c$y)*log(c(1,0)+1e-14)),tolerance=1e-12)
  expect_true(is.finite(is_mode_proposal(sp,runtime)$mode_objective))
 })
})
saveRDS(details,file.path(out,"validation_details.rds"));res<-as.data.frame(reporter$get_results())
saveRDS(res,file.path(out,"validation_test_details.rds"));res$result<-NULL
write.csv(res,file.path(out,"validation_tests.csv"),row.names=FALSE);print(colSums(res[,c("nb","passed","failed","error","warning")]))
if(any(res$failed>0|res$error|res$warning>0))stop("Joint importance validation failed")
writeLines("All analytic, forward, collapsed-target and two-node validations passed",file.path(out,"VALIDATED"))
