#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
cpp_source<-Sys.getenv("ALFAK2_DIRECT_CPP_SOURCE",file.path(root,"results/direct_forward_source_v1"))
Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp_source)
source(file.path(cpp_source,"model.R"))
source(Sys.getenv("ALFAK2_DIRECT_SOURCE",file.path(root,"direct_joint.R")))
source(file.path(repo,"validation/posterior_reference_4to8/reference.R"))
runtime<-design_load(repo)
library(testthat);reporter<-ListReporter$new()
with_reporter(reporter,{
  test_that("direct joint Hessian matches nested Laplace on fixed single/multi-trajectory models",{
    for(n in c(4L,8L))for(beta in c(0,.005)) {
      old_case<-reference_case(n,"balanced",beta,FALSE,runtime)
      context<-design_context(old_case$graph$karyotypes,beta,runtime)
      for(two in c(FALSE,TRUE)) {
        y<-old_case$counts
        case<-list(id=paste0("direct_equivalence_",n,"_",beta,"_",two),y=if(two)cbind(y,y[n:1,,drop=FALSE])else y,
          times=if(two)c(0,1.5,0,1.5)else c(0,1.5),trajectory=if(two)c(0L,0L,1L,1L)else c(0L,0L),
          reference_weight=old_case$weight)
        a<-design_fit(case,context,runtime)
        b<-direct_joint_fit(case,context,runtime)
        expect_identical(b$diagnostics$covariance_status,"trusted_direct_joint_hessian")
        expect_lte(max(abs(b$latent_fitness_mean-a$latent_fitness_mean)),1e-5)
        expect_lte(max(abs(b$latent_fitness_covariance-a$latent_fitness_covariance)),1e-5)
        expect_lte(abs(as.numeric(b$diagnostics$objective)-as.numeric(a$diagnostics$objective)),1e-6)
        expect_lte(max(abs(b$posterior_predictive$predicted-a$posterior_predictive$predicted)),1e-5)
        expect_lte(b$diagnostics$inner_gradient,.001)
      }
    }
  })
  test_that("normalized independent mean action is stable at high common growth",{
    k<-as.matrix(expand.grid(1:3,1:3,1:3));context<-design_context(k,.005,runtime)
    initial<-rep(1/27,27)
    for(rate in c(1,1000,5000)) {
      trial<-direct_forward_r(initial,rep(rate,27),context$operator,1.5,120L,32L)
      exact<-direct_exact_forward(trial$generator,initial,1.5,runtime)
      expect_true(all(is.finite(exact$pi1)));expect_true(is.finite(exact$log_total_growth))
      expect_equal(sum(exact$pi1),1,tolerance=1e-12)
      if(rate==1) {
        old<-runtime$independent$process_deterministic_mean(trial$generator,initial,1.5)
        expect_lte(max(abs(exact$pi1-old/sum(old))),1e-12)
        expect_lte(abs(exact$log_total_growth-log(sum(old))),1e-12)
      }
    }
  })
  test_that("high-order TMB and R actions agree without changing normalized target",{
    context<-design_context(as.matrix(expand.grid(1:2,1:2)),.005,runtime)
    case<-list(y=matrix(c(20,10,8,2,14,12,10,4),4),times=c(0,1.5),trajectory=c(0L,0L))
    for(shift in c(4,8)) {
      d<-design_tmb_data(case,context,shift=shift);p<-design_parameters(case,context,shift=shift)
      if(shift==8)p$u<-c(0,.6,-.7,.3)
      obj<-TMB::MakeADFun(d,p,random=NULL,DLL="observation_model",silent=TRUE)
      report<-obj$report(obj$par)
      answer<-direct_forward_r(report$pi0[,1],as.numeric(report$f_report),context$operator,1.5,d$steps[2],d$taylor_order)
      expect_equal(d$taylor_order,if(shift==4)32L else 64L)
      expect_true(all(is.finite(report$predicted)))
      expect_true(all(is.finite(report$log_growth)))
      expect_lte(max(abs(answer$pi1-report$predicted[,2])),1e-10)
      expect_lte(abs(answer$log_total_growth-report$log_growth[2]),1e-10)
      TMB::FreeADFun(obj)
    }
  })
  test_that("nearest-node continuation preserves the integrated posterior",{
    cold<-new.env(parent=globalenv())
    sys.source(file.path(root,"results/direct_source_v5/direct_joint.R"),cold)
    source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
    source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
    for(n in c(4L,8L)) {
      old<-reference_case(n,"balanced",.005,FALSE,runtime)
      context<-design_context(old$graph$karyotypes,.005,runtime)
      case<-list(id=paste0("warm_equivalence",n),y=cbind(old$counts,old$counts[n:1,,drop=FALSE]),
        times=c(0,1.5,0,1.5),trajectory=c(0L,0L,1L,1L),reference_weight=old$weight)
      a<-cold$direct_joint_integrate(case,context,runtime)
      b<-direct_joint_integrate(case,context,runtime)
      A<-diag(n)-matrix(case$reference_weight,n,n,byrow=TRUE)
      expect_identical(b$status,"passed_quadrature_refinement")
      expect_lte(max(abs(a$mean-b$mean)),1e-5)
      expect_lte(norm(a$covariance-b$covariance,"F")/norm(a$covariance,"F"),1e-4)
      expect_lte(norm(A%*%(a$covariance-b$covariance)%*%t(A),"F")/norm(A%*%a$covariance%*%t(A),"F"),1e-4)
    }
  })
})
saveRDS(reporter$get_results(),file.path(root,"results/direct_validation_details.rds"))
answer<-as.data.frame(reporter$get_results());answer$result<-NULL
write.csv(answer,file.path(root,"results/direct_validation_tests.csv"),row.names=FALSE)
print(colSums(answer[,c("nb","failed","error","warning","passed")]))
if(any(answer$failed>0|answer$error|answer$warning>0))stop("Direct joint-Hessian equivalence not passed")
