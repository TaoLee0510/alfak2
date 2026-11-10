#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
src<-Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
source(file.path(src,"model.R"));source(file.path(src,"data.R"))
source(file.path(repo,"validation/accuracy_round2/truth_and_data.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
output<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(repo,"validation/observation_design/results"))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
runtime<-design_load(repo)
k<-as.matrix(expand.grid(rep(list(1:3),3L)))
context<-design_context(k,.005,runtime)
context_file<-file.path(output,"context.rds")
if(file.exists(context_file)&&!identical(readRDS(context_file),context))stop("Frozen design context changed")
if(!file.exists(context_file))saveRDS(context,context_file)
data_dir<-file.path(output,"data");dir.create(data_dir,showWarnings=FALSE)
rows<-list()
for(family in c("smooth","rough","interaction"))for(land in 1:2)for(replicate in 1:2) {
  unit<-design_generate_unit(context,family,land,replicate,round2_truth)
  uid<-unit$cases[[1]]$unit_id
  path<-file.path(data_dir,paste0(uid,"_latent.rds"))
  if(file.exists(path)&&!identical(readRDS(path),unit))stop("Frozen design data changed")
  if(!file.exists(path))saveRDS(unit,path)
  for(case in unit$cases) {
    path<-file.path(data_dir,paste0(case$id,".rds"))
    if(file.exists(path)&&!identical(readRDS(path),case))stop("Frozen design case changed")
    if(!file.exists(path))saveRDS(case,path)
    rows[[case$id]]<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,
      family=family,biological_replicate=replicate,design=case$design,
      typed_cells=case$total_typed_cells,initial_cells=case$total_initial_cells,
      cultures=case$cultures,observations=case$sampling_occasions,
      actual_max_distance=max(case$actual_distance),baseline_max_distance=max(case$baseline_distance),
      file=normalizePath(path),sha256=unname(round2_sha256(path)))
  }
  cat("Prepared",uid,"with all four count designs\n")
}
write.csv(do.call(rbind,rows),file.path(output,"dataset_manifest.csv"),row.names=FALSE)
cat("48 development cases generated; reserved calibration/test seeds not used\n")
