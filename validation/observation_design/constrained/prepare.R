#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT"))
source(file.path(root,"results/computation_source_v4/data.R"));source(file.path(repo,"validation/accuracy_round2/common.R"))
selected<-readRDS(file.path(out,"selected_designs.rds"));planning<-readRDS(file.path(out,"planning_context.rds"))
stopifnot(length(selected)>0L,!file.exists(file.path(out,"dataset_manifest.csv")))
context<-planning$context;nominal<-planning$nominal
original<-read.csv(file.path(root,"results/dataset_manifest.csv"));original<-original[original$design=="depth600",]
stopifnot(nrow(original)==12L)
dir.create(file.path(out,"data"),showWarnings=FALSE);dir.create(file.path(out,"latent_paths"),showWarnings=FALSE)
files<-c(file.path(src,"prepare.R"),file.path(root,"results/computation_source_v4/data.R"),
  file.path(repo,"validation/accuracy_round2/common.R"),file.path(out,c("selected_designs.rds","selected_designs.csv","planning_context.rds")),original$file)
write.csv(data.frame(file=files,sha256=unname(round2_sha256(files))),file.path(out,"generation_source_sha256.csv"),row.names=FALSE)
rows<-list();ordered<-character();preflight<-character();seeds<-list()
for(i in seq_len(nrow(original))) {
  row<-original[i,];stopifnot(unname(round2_sha256(row$file))==row$sha256)
  base<-readRDS(row$file);dest<-file.path(out,"data",paste0(base$id,".rds"));stopifnot(file.copy(row$file,dest))
  cases<-list(depth600=base)
  fi<-match(base$family,c("smooth","rough","interaction"));landscape<-base$landscape_replicate;rep<-base$biological_replicate
  ps<-8100000L+fi*10000L+landscape*100L+rep*10L;ss<-8200000L+fi*10000L+landscape*100L+rep*10L
  set.seed(ss);uniform<-matrix(runif(600L*4L),600L,4L)
  for(class in names(selected)) {
    design<-selected[[class]];nr<-length(design$initial);paths<-list();pop<-list();depth<-integer();traj<-integer();times<-numeric()
    for(culture in seq_len(nr)) {
      path<-design_simulate_path(context$graph$karyotypes,base$truth,.005,design$initial[[culture]],
        design$initial_cells[culture],ps+culture,times=c(0,1.5))
      paths[[culture]]<-path;pop<-c(pop,list(path$populations[,1],path$populations[,2]))
      depth<-c(depth,design$depths[[culture]]);traj<-c(traj,rep(culture-1L,2));times<-c(times,0,1.5)
    }
    y<-vapply(seq_along(pop),function(o)design_count_uniform(pop[[o]],uniform[,o],depth[o]),numeric(context$n))
    case<-base;case$id<-paste0(base$unit_id,"_opt_",class);case$design<-paste0("opt_",class)
    case$y<-y;case$times<-times;case$trajectory<-as.integer(traj)
    case$actual_distance<-apply(nominal$distance[,which(rowSums(y)>0),drop=FALSE],1,min)
    case$total_typed_cells<-sum(depth);case$total_initial_cells<-sum(design$initial_cells)
    case$cultures<-nr;case$sampling_occasions<-length(times);case$process_seeds<-ps+seq_len(nr)
    case$sampling_seeds<-ss;case$design_specification<-design;case$pairing_scope<-"same landscape/replicate index; baseline and candidate cultures are separately simulated"
    cases[[case$design]]<-case
    saveRDS(case,file.path(out,"data",paste0(case$id,".rds")))
    saveRDS(paths,file.path(out,"latent_paths",paste0(case$id,".rds")))
    seeds[[length(seeds)+1L]]<-data.frame(case_id=case$id,process_seed_A=ps+1L,
      process_seed_B=if(nr==2)ps+2L else NA_integer_,sampling_seed=ss)
    if(base$family=="interaction"&&landscape==1L&&rep==1L)preflight<-c(preflight,case$id)
  }
  for(case in cases) {
    file<-normalizePath(file.path(out,"data",paste0(case$id,".rds")))
    stopifnot(sum(case$y)==600,case$total_initial_cells==1200,
      identical(case$reference_weight,base$reference_weight),identical(case$baseline_distance,base$baseline_distance))
    rows[[length(rows)+1L]]<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,
      family=case$family,biological_replicate=case$biological_replicate,design=case$design,
      typed_cells=sum(case$y),initial_cells=case$total_initial_cells,cultures=case$cultures,observations=ncol(case$y),
      actual_max_distance=max(case$actual_distance),baseline_max_distance=max(case$baseline_distance),
      file=file,sha256=unname(round2_sha256(file)))
    ordered<-c(ordered,case$id)
  }
}
manifest<-do.call(rbind,rows);stopifnot(nrow(manifest)==12L*(1L+length(selected)),!anyDuplicated(manifest$case_id))
write.csv(manifest,file.path(out,"dataset_manifest.csv"),row.names=FALSE)
write.csv(do.call(rbind,seeds),file.path(out,"development_seed_ledger.csv"),row.names=FALSE)
queue<-list(preflight=as.list(preflight),main=as.list(setdiff(ordered,preflight)))
jsonlite::write_json(queue,file.path(out,"queue.json"),auto_unbox=TRUE,pretty=TRUE)
cat(nrow(manifest),"cases generated; preflight",paste(preflight,collapse=", "),"\n")
