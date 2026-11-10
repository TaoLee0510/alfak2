#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
cpp_source<-Sys.getenv("ALFAK2_DIRECT_CPP_SOURCE",file.path(root,"results/direct_forward_source_v1"))
Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp_source)
source(file.path(cpp_source,"model.R"))
source(Sys.getenv("ALFAK2_DIRECT_SOURCE",file.path(root,"direct_joint.R")))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/baseline_quadrature/quadrature.R"))
runtime<-design_load(repo)
case<-readRDS(file.path(root,"results/data/design_interaction_L01_R1_composition_anchor600.rds"))
context<-readRDS(file.path(root,"results/context.rds"))
started<-proc.time()[["elapsed"]]
fit<-direct_joint_integrate(case,context,runtime)
fit$elapsed_seconds<-proc.time()[["elapsed"]]-started
saveRDS(fit,file.path(root,"results/direct_smoke_fit.rds"))
print(fit$history,row.names=FALSE)
cat("Direct smoke",fit$status,"seconds",fit$elapsed_seconds,"\n")
if(fit$status!="passed_quadrature_refinement")stop("Integrated smoke did not pass numerical refinement")
