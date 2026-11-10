#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_MECHANISM_OUTPUT"))
source(file.path(src,"common.R"));runtime<-mechanism_load(repo,src)
cases<-readRDS(file.path(out,"frozen_cases.rds"));context<-readRDS(file.path(out,"context.rds"))
library(testthat);reporter<-ListReporter$new();validation_results<-list()
with_reporter(reporter,{
  test_that("new reference targets the current full-generator likelihood and priors",{
    for(case in cases) {
      n<-case$n;sp<-case$specification;infer<-case[c("id","y","times","trajectory","reference_weight")]
      data<-design_tmb_data(infer,context);param<-design_parameters(infer,context)
      obj<-TMB::MakeADFun(data,param,random=NULL,DLL="observation_model",silent=TRUE)
      expect_equal(tcrossprod(sp$lower),context$covariance,tolerance=1e-12)
      expect_equal(crossprod(sp$eta_basis)/4,diag(n-1L),tolerance=1e-12)
      expect_equal(unname(tcrossprod(sp$eta_basis)/4),diag(n)-matrix(1/n,n,n),tolerance=1e-12)
      for(j in 1:12) {
        set.seed(9400000L+j);z<-rnorm(sp$q,sd=.75)
        g<-as.numeric(sp$offset+sp$lower%*%z[1:n]);eta<-as.numeric(sp$eta_basis%*%z[n+seq_len(n-1L)])
        p<-exp(eta-max(eta));p<-p/sum(p)
        A<-sweep(2*t(sp$probability)-diag(n),2,exp(g),"*")
        mass<-runtime$independent$process_deterministic_mean(A,p,sp$dt)
        q<-mass/sum(mass);lg<-log(sum(mass));ll<-sum(sp$y0*log(p+1e-14))+sum(sp$y1*log(q+1e-14))
        ref<-runtime$cpp$ref_likelihood_cpp(z,sp)
        par<-obj$par;par[names(par)=="zeta"]<-z[n+seq_len(n-1L)]
        par[names(par)=="u"]<-g-sp$offset;par[names(par)=="coefficients"]<-0
        report<-obj$report(par)
        tmbll<-sum(sp$y0*log(report$pi0[,1]+1e-14))+sum(sp$y1*log(report$predicted[,2]+1e-14))
        expect_lt(max(abs(q-ref$pi1)),1e-11)
        expect_lt(abs(lg-ref$log_growth),1e-10)
        expect_lt(abs(ll-ref$log_likelihood),1e-8)
        expect_lt(max(abs(q-report$predicted[,2])),1e-8)
        expect_lt(abs(tmbll-ref$log_likelihood),1e-6)
      }
      TMB::FreeADFun(obj)
    }
  })
  test_that("density-corrected Gaussian-reference ESS recovers complete Gaussian/lognormal prior moments",{
    case<-cases[[1]];case$specification$y0[]<-case$specification$y1[]<-0
    chains<-reference_chains(case,runtime,10000L,warmup=1000L,seed=9401000L,cores=1L)
    d<-reference_diagnostics(chains,case)
    m<-d$summary[match(paste0("birth_",1:8),d$summary$variable),]
    expected<-exp(case$prior$offset+diag(case$covariance)/2)
    C<-tcrossprod(expected)*expm1(case$covariance)
    expect_true(all(abs(m$mean-expected)<=6*m$mcse_mean))
    expect_true(all(abs(d$covariance$birth-C)<=6*d$covariance_mcse$birth))
    expect_true(d$convergence$passed)
    expect_true("eta_contrast"%in%names(d$covariance))
    validation_results$prior<<-d
  })
  test_that("density-corrected Gaussian-reference ESS matches the independent two-node quadrature posterior",{
    old<-readRDS(file.path(repo,"validation/posterior_reference_4to8/inputs/two_node_reference.rds"))
    sp<-list(n=2L,q=3L,steps=15L,dt=old$dt,offset=old$prior$offset,
      lower=t(chol(old$aggregate_prior_covariance)),eta_basis=matrix(c(1,-1)*old$eta_prior_sd/sqrt(2),2,1),
      known_pi0=numeric(),probability=diag(2),beta_zero=TRUE,y0=old$counts[,1],y1=old$counts[,2],growth=0,growth_sd=0)
    case<-list(id="two_node_count_validation",n=2L,specification=sp,weight=c(.5,.5),edges=matrix(c(1L,2L),1,2))
    old_helper<-new.env(parent=globalenv())
    sys.source(file.path(repo,"validation/posterior_reference_4to8/reference.R"),old_helper)
    pilot<-old_helper$reference_chains(case,runtime,10000L,warmup=1000L,seed=9410000L,cores=1L)
    z<-do.call(rbind,lapply(pilot,function(ch)tail(ch$draws,5000L)));S<-2*cov(z)+diag(1e-4,3)
    case$proposal<-list(mean=colMeans(z),lower=t(chol(S)),precision=solve(S))
    chains<-reference_chains(case,runtime,30000L,warmup=5000L,seed=9402000L,cores=1L)
    d<-reference_diagnostics(chains,case);exact<-tail(old$results$count_only$reference,1)[[1]]
    m<-d$summary[match(c("birth_1","birth_2","relative_1"),d$summary$variable),]
    expect_true(all(abs(m$mean-c(exact$mean,exact$contrast_mean))<=6*m$mcse_mean+3e-5))
    expect_true(all(abs(m$sd-c(sqrt(diag(exact$covariance)),exact$contrast_sd))<=6*m$mcse_sd+3e-5))
    expect_true(d$convergence$passed)
    validation_results$two_node<<-d
  })
})
saveRDS(validation_results,file.path(out,"validation_reference_details.rds"))
res<-as.data.frame(reporter$get_results());saveRDS(res,file.path(out,"validation_test_details.rds"));res$result<-NULL
write.csv(res,file.path(out,"validation_tests.csv"),row.names=FALSE);print(colSums(res[,c("nb","passed","failed","error","warning")]))
if(any(res$failed>0|res$error|res$warning>0))stop("Reference validation failed")
ledger<-mechanism_source_manifest(repo,src,out,file.path(out,c("frozen_cases.rds","context.rds")),runtime)
write.csv(ledger,file.path(out,"validation_source_sha256.csv"),row.names=FALSE)
cat("New target, unchanged sampler, prior/gauge and two-node references passed\n")
