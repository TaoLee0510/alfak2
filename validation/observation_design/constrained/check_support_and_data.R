#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT"))
source(file.path(root,"results/constrained_source_v2/planning_helpers.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
p<-readRDS(file.path(out,"planning_context.rds"));selected<-readRDS(file.path(out,"selected_designs.rds"))
original<-read.csv(file.path(root,"results/dataset_manifest.csv"));controls<-original[original$design=="control300",]
observed<-Reduce(`|`,lapply(controls$file,function(f)readRDS(f)$baseline_distance==0))
rows<-list()
for(class in names(selected))for(floor in p$floors)for(point in names(p$kernels)) {
  smooth<-function(x){x<-pmax(x,floor);x/sum(x)}
  b<-constrained_risk(list(smooth(p$nominal$A)),list(c(300,300)),p$kernels[[point]],p$context,p$A)
  x<-constrained_risk(lapply(selected[[class]]$initial,smooth),selected[[class]]$depths,p$kernels[[point]],p$context,p$A)
  rows[[length(rows)+1L]]<-data.frame(class=class,floor=floor,point=point,
    max_all_states=max(x$variance/b$variance),max_actual_observed=max(x$variance[observed]/b$variance[observed]))
}
information<-do.call(rbind,rows)
write.csv(information,file.path(out,"actual_observed_information_diagnostic.csv"),row.names=FALSE)
stopifnot(all(information$max_actual_observed<=1+1e-6))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"));audit<-list()
for(i in seq_len(nrow(manifest))) {
  row<-manifest[i,];case<-readRDS(row$file)
  base_row<-original[original$unit_id==row$unit_id & original$design=="depth600",];base<-readRDS(base_row$file)
  stopifnot(nrow(base_row)==1L,identical(case$truth,base$truth),identical(case$truth_relative,base$truth_relative),
    identical(case$reference_weight,base$reference_weight),identical(case$baseline_distance,base$baseline_distance),
    identical(case$planned_distance,base$planned_distance),identical(case$target_labels,base$target_labels),
    identical(case$target_edges,base$target_edges),sum(case$y)==600,all(case$y>=0),all(case$y==round(case$y)),
    nrow(case$y)==27L,length(case$times)==ncol(case$y),length(case$trajectory)==ncol(case$y),
    case$total_initial_cells==1200L,unname(round2_sha256(row$file))==row$sha256)
  if(row$design=="depth600")stopifnot(row$sha256==base_row$sha256)
  else {
    paths<-readRDS(file.path(out,"latent_paths",paste0(case$id,".rds")))
    stopifnot(sum(vapply(paths,function(a)sum(a$initial),numeric(1)))==1200,
      all(vapply(paths,function(a)identical(a$times,c(0,1.5)),logical(1))),
      length(case$process_seeds)==length(paths),all(case$process_seeds==vapply(paths,`[[`,numeric(1),"seed")))
  }
  audit[[row$case_id]]<-data.frame(case_id=row$case_id,design=row$design,sha256_unchanged=TRUE,
    original_reference_and_distance=TRUE,same_development_truth=TRUE,typed_cells=sum(case$y),
    total_initial_cells=case$total_initial_cells,actual_max_distance=max(case$actual_distance),
    baseline_max_distance=max(case$baseline_distance))
}
write.csv(do.call(rbind,audit),file.path(out,"data_integrity_audit.csv"),row.names=FALSE)
cat(nrow(manifest),"data cases passed input-budget, seed, target and immutable-reference checks; observed union",sum(observed),"states\n")
