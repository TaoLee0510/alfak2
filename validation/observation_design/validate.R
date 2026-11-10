#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
src<-Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
source(file.path(src,"model.R"));source(file.path(src,"data.R"))
source(file.path(repo,"validation/posterior_reference_4to8/reference.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
runtime<-design_load(repo)
output<-file.path(repo,"validation/observation_design/results/validation")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
library(testthat)
reporter<-ListReporter$new()
with_reporter(reporter,{
  test_that("joint model agrees with the original single-trajectory likelihood and prior",{
    for(n in c(4L,8L))for(beta in c(0,.005)) {
      old_case<-reference_case(n,"balanced",beta,FALSE,runtime)
      context<-design_context(old_case$graph$karyotypes,beta,runtime)
      case<-list(id=paste0("equivalence",n,beta),y=old_case$counts,
        times=c(0,1.5),trajectory=c(0L,0L),reference_weight=old_case$weight)
      original<-runtime$pkg$fit_local_posterior(
        runtime$pkg$prepare_alfak2_data(old_case$counts,dt=1.5,beta=beta),
        graph=old_case$graph,structured_prior=old_case$prior,forward_model="division_process",
        forward_step=.025,eta_prior_sd=2,eta_borrowed_prior_mean=0,eta_borrowed_prior_sd=2,
        eta_distance_penalty=0,sdreport_mode="fitness_only",conditional_warmstart="none")
      joint<-design_fit(case,context,runtime,keep_object=TRUE,check_forward=FALSE)
      expect_identical(joint$diagnostics$covariance_status,"TMB_sdreport")
      expect_equal(as.numeric(joint$latent_fitness_mean),as.numeric(original$latent_fitness_mean),tolerance=1e-6)
      expect_equal(unname(joint$latent_fitness_covariance),unname(original$latent_fitness_covariance),tolerance=1e-6)
      expect_equal(joint$diagnostics$objective,original$diagnostics$objective,tolerance=1e-6)
      expect_equal(joint$posterior_predictive$predicted[,2],as.numeric(original$posterior_predictive$pi1),tolerance=1e-6)
      released<-design_fit(case,context,runtime,check_forward=FALSE)
      expect_equal(released$latent_fitness_mean,joint$latent_fitness_mean,tolerance=1e-10)
    }
  })
  test_that("multi-observation objective, derivatives, and trajectory permutation are correct",{
    k<-as.matrix(expand.grid(1:2,1:2));context<-design_context(k,.005,runtime)
    case<-list(y=matrix(c(20,10,8,2,14,12,10,4,3,7,11,19,4,9,12,15),4),
      times=c(0,1.5,0,1.5),trajectory=c(0L,0L,1L,1L))
    d<-design_tmb_data(case,context);p<-design_parameters(case,context)
    set.seed(812001)
    p$u<-rnorm(4,sd=.04);p$coefficients<-rnorm(length(p$coefficients))*d$coef_sd
    p$zeta[]<-rnorm(length(p$zeta),sd=.5)
    obj<-TMB::MakeADFun(d,p,DLL="observation_model",silent=TRUE)
    g<-d$offset+p$u+as.numeric(d$design%*%p$coefficients);b<-exp(g)
    eta<-d$eta_basis%*%p$zeta
    pi0<-apply(eta,2,function(x){x<-exp(x-max(x));x/sum(x)})
    ll<-0
    for(o in 1:4) {
      mass<-if(case$times[o]==0)pi0[,case$trajectory[o]+1L]else
        runtime$independent$process_deterministic_mean(
          sweep(2*t(runtime$independent$process_independent_operator(context$graph$karyotypes,b,.005)$probability)-diag(4),2,b,"*"),
          pi0[,case$trajectory[o]+1L],case$times[o])
      pr<-mass/sum(mass)
      ll<-ll+sum(case$y[,o]*log(pr+1e-14))
    }
    prior_nll<-.5*as.numeric(crossprod(p$u,d$precision%*%p$u))-.5*d$logdet+2*log(2*pi)-
      sum(dnorm(p$coefficients,0,d$coef_sd,log=TRUE))-sum(dnorm(p$zeta,0,1,log=TRUE))
    expect_equal(obj$fn(obj$par),prior_nll-ll,tolerance=1e-10)
    eps<-1e-6
    fd<-vapply(seq_along(obj$par),function(i){delta<-numeric(length(obj$par));delta[i]<-eps
      (obj$fn(obj$par+delta)-obj$fn(obj$par-delta))/(2*eps)},numeric(1))
    expect_lt(max(abs(obj$gr(obj$par)-fd)),1e-4)
    swap<-d;swap$trajectory<-1L-d$trajectory
    pp<-p;pp$zeta<-p$zeta[,2:1,drop=FALSE]
    obj2<-TMB::MakeADFun(swap,pp,DLL="observation_model",silent=TRUE)
    expect_equal(obj$fn(obj$par),obj2$fn(obj2$par),tolerance=1e-10)
    fake<-case;fake$truth<-rep(999,4);fake$family<-"fake";fake$latent_pi0<-matrix(1,4,2)
    expect_identical(design_tmb_data(case,context),design_tmb_data(fake,context))
    # Re-evaluate a tape after the minimum diagonal changes, not just at its
    # construction point. This detects a frozen AD comparison branch.
    p2<-p;p2$u<-c(1.5,-1.5,.7,-.7);v2<-unlist(p2,use.names=FALSE)
    g2<-d$offset+p2$u+as.numeric(d$design%*%p2$coefficients)
    b2<-exp(g2);ll2<-0
    for(o in 1:4) {
      pr<-if(case$times[o]==0)pi0[,case$trajectory[o]+1L]else
        design_forward_r(pi0[,case$trajectory[o]+1L],b2,context$operator,case$times[o],d$steps[o])$pi1
      ll2<-ll2+sum(case$y[,o]*log(pr+1e-14))
    }
    prior2<-.5*as.numeric(crossprod(p2$u,d$precision%*%p2$u))-.5*d$logdet+2*log(2*pi)-
      sum(dnorm(p2$coefficients,0,d$coef_sd,log=TRUE))-sum(dnorm(p2$zeta,0,1,log=TRUE))
    expect_equal(obj$fn(v2),prior2-ll2,tolerance=1e-10)
  })
  test_that("snapshot simulator has correct pure-birth moments and nested counts",{
    k<-matrix(2L,1,1);N<-100L;b<-.8;t<-.5
    x<-vapply(1:400,function(i)design_simulate_path(k,b,0,1,N,815000L+i,times=c(0,.25,t))$populations[1,3],numeric(1))
    expected<-N*exp(b*t);variance<-N*exp(b*t)*(exp(b*t)-1)
    expect_lt(abs(mean(x)-expected)/sqrt(variance/400),4)
    expect_lt(abs(var(x)-variance)/variance,.25)
    u<-seq(.001,.999,length.out=300)
    a<-design_count_uniform(c(40,30,20,10),u,150)
    z<-design_count_uniform(c(40,30,20,10),u,300)
    expect_equal(sum(a),150);expect_equal(sum(z),300);expect_true(all(z>=a))
    expect_equal(sum(design_count_uniform(c(1,1,1,0,0),u,300)),300)
    expect_equal(design_count_uniform(c(0,0,3,0),u,300),c(0,0,300,0))
    short<-design_simulate_path(k,b,0,1,N,818000L,times=c(0,t))
    middle<-design_simulate_path(k,b,0,1,N,818000L,times=c(0,.25,t))
    expect_identical(short$populations[,2],middle$populations[,3])
  })
  if(file.exists(file.path(repo,"validation/observation_design/results/dataset_manifest.csv")))
    test_that("design data have equal budgets and paired truth, references and strata",{
      manifest<-read.csv(file.path(repo,"validation/observation_design/results/dataset_manifest.csv"))
      expect_equal(nrow(manifest),60L)
      for(unit in unique(manifest$unit_id)) {
        m<-manifest[manifest$unit_id==unit,]
        cases<-lapply(m$file,readRDS);names(cases)<-m$design
        expect_setequal(names(cases),c("control300","depth600","time600","composition600","composition_anchor600"))
        expect_equal(vapply(cases,function(x)sum(x$y),numeric(1)),
          c(control300=300,depth600=600,time600=600,composition600=600,composition_anchor600=600)[names(cases)])
        expect_true(all(vapply(cases,function(x)identical(x$truth,cases[[1]]$truth),logical(1))))
        expect_true(all(vapply(cases,function(x)identical(x$reference_weight,cases[[1]]$reference_weight),logical(1))))
        expect_true(all(vapply(cases,function(x)identical(x$baseline_distance,cases[[1]]$baseline_distance),logical(1))))
        expect_true(all(cases$depth600$y>=cases$control300$y))
        expect_true(all(cases$depth600$y>=cases$time600$y[,c(1,3)]))
        expect_identical(cases$composition600$y[,1:2],cases$composition_anchor600$y[,1:2])
      }
    })
})
results<-as.data.frame(reporter$get_results());results$result<-NULL
write.csv(results,file.path(output,"validation_tests.csv"),row.names=FALSE)
print(colSums(results[,c("nb","failed","error","warning","passed")]))
if(any(results$failed>0|results$error|results$warning>0))stop("Joint design validation failed")
cat("Joint-model, derivative, gauge, no-oracle and simulator checks passed\n")
