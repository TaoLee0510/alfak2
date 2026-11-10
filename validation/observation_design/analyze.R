#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
out<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(repo,"validation/observation_design/results"))
parts<-file.path(out,paste0("fits_",0:1))
status<-do.call(rbind,lapply(parts,function(p)read.csv(file.path(p,"status.csv"))))
metrics<-do.call(rbind,lapply(parts,function(p)read.csv(file.path(p,"case_metrics.csv"))))
manifest<-read.csv(file.path(out,"dataset_manifest.csv"))
if(nrow(status)!=60L||anyDuplicated(status$case_id)||!setequal(status$case_id,manifest$case_id))
  stop("Incomplete or duplicated primary batch")
write.csv(status,file.path(out,"all_fit_status.csv"),row.names=FALSE)
write.csv(metrics,file.path(out,"all_case_metrics.csv"),row.names=FALSE)
resource_pass<-if("resource_pass"%in%names(status))tolower(as.character(status$resource_pass))=="true" else status$elapsed_seconds<=300
resource_pass[is.na(resource_pass)]<-FALSE
write.csv(status[!resource_pass,],file.path(out,"resource_failures.csv"),row.names=FALSE)
availability<-do.call(rbind,lapply(split(status,status$design),function(x) {
  data.frame(design=x$design[1],attempted=nrow(x),numerically_available=sum(x$status=="passed_quadrature_refinement"),
    raw_numerical_failures=if("raw_status"%in%names(x))sum(x$raw_status!="passed_quadrature_refinement")else sum(x$status!="passed_quadrature_refinement"),
    recovered=if("recovery"%in%names(x))sum(x$recovery!="none")else 0L,
    budget_failures=sum(!resource_pass[match(x$case_id,status$case_id)]))
}))
write.csv(availability,file.path(out,"availability_by_design.csv"),row.names=FALSE)
if(any(status$status!="passed_quadrature_refinement")) {
  write.csv(status[status$status!="passed_quadrature_refinement",],file.path(out,"numerical_failures.csv"),row.names=FALSE)
  writeLines("NUMERICALLY_INCOMPLETE: all failures retained; no complete design-risk ranking is certified.",file.path(out,"stage_status.txt"))
  stop("Numerical failures remain; do not aggregate away unavailable predictions")
}
key<-paste(metrics$case_id,metrics$view,metrics$distance)
if(anyDuplicated(key))stop("Duplicated metric rows")
for(i in seq_len(nrow(manifest))) {
  case<-readRDS(manifest$file[i]);rows<-metrics[metrics$case_id==case$id,]
  expected<-unlist(lapply(c("baseline_distance","planned_distance","actual_distance"),function(view)
    paste(view,c(-1L,sort(unique(case[[view]]))))))
  if(!setequal(paste(rows$view,rows$distance),expected))stop("Missing metric strata")
}
risks<-c("rmse","mse","shape_rmse","shape_mse","offset_squared","edge_rmse",
  "score50","score80","score95","coverage50","coverage80","coverage95",
  "width50","width80","width95")
if(any(metrics$prediction_fraction!=1)||any(metrics$interval_fraction50!=1)||
   any(metrics$interval_fraction80!=1)||any(metrics$interval_fraction95!=1))stop("Incomplete predictions or intervals")
if(any(!is.finite(as.matrix(metrics[,setdiff(risks,"edge_rmse")]))))stop("Nonfinite primary risks")
if(any(!is.finite(metrics$edge_rmse[metrics$n_edges>0])))stop("Missing supported edge risks")
primary<-metrics[metrics$view=="baseline_distance",]
if(any(primary$n_edges==0))stop("Primary edge stratum unsupported; report explicitly before comparison")
land<-aggregate(primary[,risks],primary[,c("design","family","landscape_id","distance")],mean)
family<-aggregate(land[,risks],land[,c("design","family","distance")],mean)
overall<-aggregate(family[,risks],family[,c("design","distance")],mean)
write.csv(land,file.path(out,"baseline_landscape_metrics.csv"),row.names=FALSE)
write.csv(family,file.path(out,"baseline_family_metrics.csv"),row.names=FALSE)
write.csv(overall,file.path(out,"baseline_overall_metrics.csv"),row.names=FALSE)
comparisons<-differences<-list()
for(reference in c("control300","depth600")) {
  base<-family[family$design==reference,]
  for(design in setdiff(unique(family$design),reference)) {
    x<-family[family$design==design,]
    idx<-match(paste(x$family,x$distance),paste(base$family,base$distance))
    z<-x[,c("family","distance")];z$design<-design;z$reference<-reference
    for(metric in c("rmse","shape_rmse","edge_rmse","score50","score80","score95"))z[[metric]]<-x[[metric]]/base[[metric]][idx]
    comparisons[[paste(reference,design)]]<-z
    a<-land[land$design==reference,];b<-land[land$design==design,]
    id<-match(paste(b$landscape_id,b$distance),paste(a$landscape_id,a$distance))
    v<-b[,c("family","landscape_id","distance")];v$design<-design;v$reference<-reference
    for(metric in c("rmse","shape_rmse","edge_rmse","score95"))v[[paste0(metric,"_difference")]]<-b[[metric]]-a[[metric]][id]
    differences[[paste(reference,design)]]<-v
  }
}
ratios<-do.call(rbind,comparisons)
write.csv(ratios,file.path(out,"family_risk_ratios.csv"),row.names=FALSE)
write.csv(do.call(rbind,differences),file.path(out,"paired_landscape_differences.csv"),row.names=FALSE)
triage<-list()
reference<-overall[overall$design=="depth600"&overall$distance>=0,]
for(design in c("time600","composition_anchor600")) {
  x<-overall[overall$design==design&overall$distance>=0,]
  x<-x[match(reference$distance,x$distance),]
  r<-ratios[ratios$design==design&ratios$reference=="depth600"&ratios$distance>=0,]
  flags<-c(point_nonworse=all(x$rmse<=reference$rmse),
    shape_5percent_all_distances=all(x$shape_rmse<=.95*reference$shape_rmse),
    edge_nonworse=all(x$edge_rmse<=reference$edge_rmse),
    interval_scores_nonworse=all(as.matrix(x[,c("score50","score80","score95")])<=as.matrix(reference[,c("score50","score80","score95")])),
    family_safety=all(as.matrix(r[,c("rmse","shape_rmse","edge_rmse","score95")])<=1.05))
  within_budget<-all(resource_pass[status$design%in%c(design,"depth600")])
  triage[[design]]<-cbind(data.frame(design=design,reference="depth600",all_distance_promising=all(flags)&&within_budget,
    statistical_development_flags_pass=all(flags),resource_budget_pass=within_budget),as.data.frame(as.list(flags)))
}
triage<-do.call(rbind,triage)
write.csv(triage,file.path(out,"development_triage.csv"),row.names=FALSE)
denominators<-aggregate(cbind(rank_available=is.finite(metrics$spearman),amplitude_available=is.finite(metrics$amplitude_ratio)),
  metrics[,c("design","view","distance")],sum)
write.csv(denominators,file.path(out,"descriptive_metric_denominators.csv"),row.names=FALSE)
cost<-aggregate(status$elapsed_seconds,status[,"design",drop=FALSE],function(x)c(median=median(x),maximum=max(x)))
write.csv(cost,file.path(out,"runtime_by_design.csv"),row.names=FALSE)
writeLines(if(any(triage$all_distance_promising))"DEVELOPMENT_PROMISING_REQUIRES_LARGER_DOMAIN_VALIDATION"else
  "NO_ALL_DISTANCE_WINNER_IN_FIXED_BUDGET_PILOT",file.path(out,"stage_status.txt"))
print(overall[overall$distance>=0,c("design","distance","rmse","shape_rmse","edge_rmse","score95","coverage95")],row.names=FALSE)
print(triage,row.names=FALSE)
