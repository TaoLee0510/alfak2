mechanism_load<-function(repo,src,compile_reference=TRUE) {
  root<-file.path(repo,"validation/observation_design")
  cpp<-file.path(root,"results/direct_forward_source_v2")
  Sys.setenv(ALFAK2_DESIGN_SOURCE=cpp)
  sys.source(file.path(cpp,"model.R"),.GlobalEnv)
  sys.source(file.path(root,"results/direct_source_v7/direct_joint.R"),.GlobalEnv)
  sys.source(file.path(root,"results/adaptive_source_v1/quadrature.R"),.GlobalEnv)
  sys.source(file.path(root,"results/adaptive_source_v1/model.R"),.GlobalEnv)
  sys.source(file.path(repo,"validation/baseline_quadrature/quadrature.R"),.GlobalEnv)
  sys.source(file.path(repo,"validation/accuracy_round2/common.R"),.GlobalEnv)
  sys.source(file.path(root,"results/computation_source_v4/score.R"),.GlobalEnv)
  sys.source(file.path(src,"diagnostics.R"),.GlobalEnv)
  runtime<-design_load(repo);invisible(TMB::openmp(1,DLL="observation_model"))
  runtime$cpp<-new.env(parent=globalenv())
  if(compile_reference)Rcpp::sourceCpp(file.path(src,"reference.cpp"),env=runtime$cpp,
    cacheDir="/private/tmp/alfak2_mechanism_reference_cpp_v2",showOutput=FALSE)
  runtime
}

mechanism_source_manifest<-function(repo,src,out,case_files,runtime) {
  root<-file.path(repo,"validation/observation_design")
  dlls<-vapply(getLoadedDLLs(),function(x)x[["path"]],character(1))
  sampler_dlls<-dlls[grepl("alfak2_mechanism_reference_cpp",dlls,fixed=TRUE)]
  files<-c(list.files(src,full.names=TRUE),
    file.path(root,"results/direct_forward_source_v2",c("model.R","observation_model.cpp","positive_taylor_forward.h")),
    file.path(root,"results/direct_source_v7/direct_joint.R"),
    file.path(root,"results/adaptive_source_v1",c("quadrature.R","model.R")),
    file.path(repo,"validation/baseline_quadrature/quadrature.R"),
    file.path(root,"results/computation_source_v4",c("score.R","data.R")),
    file.path(repo,"validation/accuracy_round2",c("common.R","process_prepare_data.R")),
    file.path(repo,"validation/posterior_reference_4to8/reference.R"),
    file.path(root,"results/control_source_v1/bounded_run.py"),
    list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE),runtime$design_binary,runtime$native,sampler_dlls,case_files)
  files<-unique(files[file.info(files)$isdir==FALSE]);hash<-round2_sha256(files)
  data.frame(file=names(hash),sha256=unname(hash))
}
