#!/usr/bin/env Rscript
source("validation/baseline_quadrature/quadrature.R")
source("validation/observation_design/adaptive/quadrature.R")
library(testthat);reporter<-ListReporter$new()
with_reporter(reporter,{
  test_that("embedded positive rules have the stated polynomial precision",{
    rule<-adaptive_gk15_rule()
    expect_true(all(rule$k>0));expect_true(all(rule$g>=0))
    for(p in 0:23)expect_equal(sum(rule$k*rule$x^p),if(p%%2)0 else 2/(p+1),tolerance=1e-13)
    for(p in 0:13)expect_equal(sum(rule$g*rule$x^p),if(p%%2)0 else 2/(p+1),tolerance=1e-13)
  })
  test_that("normalizer and complete Gaussian/lognormal moments retain the prior and tails",{
    for(tilt in c(0,4)) {
      n<-4L;sd0<-.5;K0<-diag(.02,n)+matrix(.01,n,n);offset<-rep(-.5*(sd0^2+.03),n)
      evaluate<-function(a) {
        mu<-offset+a;birth<-exp(mu+diag(K0)/2)
        list(mean=birth,covariance=tcrossprod(birth)*expm1(K0),latent_mean=mu,latent_covariance=K0,
          log_evidence=tilt*a)
      }
      fit<-adaptive_prior_integrate(evaluate,sd0,rep(1/n,n))
      truth_mean<-exp(offset+tilt*sd0^2+.5*(sd0^2+diag(K0)))
      truth_cov<-tcrossprod(truth_mean)*expm1(K0+sd0^2)
      expect_identical(fit$status,"passed_quadrature_refinement")
      expect_lte(abs(fit$log_evidence-tilt^2*sd0^2/2),1e-3)
      expect_lte(max(abs(fit$mean-truth_mean)/pmax(1,truth_mean)),1e-3)
      expect_lte(norm(fit$covariance-truth_cov,"F")/norm(truth_cov,"F"),.005)
      expect_true(all(fit$weights>=0));expect_equal(sum(fit$weights),1,tolerance=1e-12)
      expect_true(tail(fit$tail_history$passed,1))
      expect_equal(sum(fit$partition$upper-fit$partition$lower),1,tolerance=1e-14)
      expect_true(any(fit$partition$lower==0));expect_true(any(fit$partition$upper==1))
      expect_equal(tail(fit$history$evaluations,1),nrow(fit$evaluated))
      expect_false(anyDuplicated(fit$evaluated$u)>0)
      cat("Analytic tilt",tilt,"nodes",nrow(fit$evaluated),"panels",nrow(fit$partition),"\n")
    }
  })
  test_that("a failed positive node cannot be discarded",{
    expect_error(adaptive_prior_integrate(function(a)stop("deliberate failure"),.5,c(.5,.5)),"deliberate failure")
  })
})
results<-as.data.frame(reporter$get_results());saveRDS(results,"validation/observation_design/adaptive/core_test_details.rds")
results$result<-NULL;write.csv(results,"validation/observation_design/adaptive/core_tests.csv",row.names=FALSE)
print(colSums(results[,c("nb","failed","error","warning","passed")]))
if(any(results$failed>0|results$error|results$warning>0))stop("Adaptive core validation failed")
