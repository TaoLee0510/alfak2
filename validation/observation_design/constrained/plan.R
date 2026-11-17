#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_SOURCE",file.path(root,"constrained")))
out<-Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT",file.path(root,"results/constrained_v1"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(out,"selected_designs.rds")))
source(file.path(root,"results/computation_source_v4/data.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_profile_fisher.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"));source(file.path(src,"planning_helpers.R"))
context<-readRDS(file.path(root,"results/context.rds"));context$landscape_precision<-solve(context$covariance)
nominal<-design_nominal_initial(context);n<-context$n
P<-matrix(0,n,n);op<-context$operator;P[cbind(op$from0+1L,op$to0+1L)]<-op$weight
x<-context$graph$karyotypes-2
profiles<-list(prior=context$prior$offset,
  smooth=.08*x[,1]-.06*x[,2]+.03*x[,3]-.04*rowSums(x^2),
  interaction=.16*x[,1]*x[,2]-.12*x[,2]*x[,3]+.12*x[,1]*x[,2]*x[,3])
kernels<-list()
for(profile in names(profiles))for(scale in c(.8,1,1.2))
  kernels[[paste(profile,scale,sep="_")]]<-constrained_kernel(P,profiles[[profile]]+log(scale))
A<-diag(n)-matrix(nominal$A,n,n,byrow=TRUE);floors<-c(1e-4,1e-6)
grid<-constrained_grid(context,nominal);stopifnot(length(grid$spec)==nrow(grid$table),!anyDuplicated(grid$table$id))
saveRDS(list(context=context,nominal=nominal,profiles=profiles,kernels=kernels,P=P,A=A,floors=floors),file.path(out,"planning_context.rds"))
write.csv(grid$table,file.path(out,"candidate_grid.csv"),row.names=FALSE)
cat("Frozen grid",nrow(grid$table),"candidates\n")
result<-grid$table;result$worst_observed_ratio<-0;result$objective_worst_far_ratio<-0
result$objective_mean_far_ratio<-0;result$feasible<-TRUE
base_risk<-details<-list();num<-0L
for(floor in floors)for(point in names(kernels)) {
  ker<-kernels[[point]];smooth<-function(p){p<-pmax(p,floor);p/sum(p)}
  base<-constrained_risk(list(smooth(nominal$A)),list(c(300,300)),ker,context,A)
  id0<-nominal$planned_distance==0;far<-lapply(1:2,function(d)nominal$planned_distance==d)
  far_base<-vapply(far,function(ii)mean(base$variance[ii]),numeric(1))
  num<-num+1L;base_risk[[num]]<-data.frame(point=point,floor=floor,node=seq_len(n),variance=base$variance)
  for(i in seq_len(nrow(result))) {
    d<-grid$spec[[i]]
    risk<-constrained_risk(lapply(d$initial,smooth),d$depths,ker,context,A)
    near_ratio<-max(risk$variance[id0]/base$variance[id0])
    far_ratio<-vapply(far,function(ii)mean(risk$variance[ii]),numeric(1))/far_base
    result$worst_observed_ratio[i]<-max(result$worst_observed_ratio[i],near_ratio)
    result$objective_worst_far_ratio[i]<-max(result$objective_worst_far_ratio[i],far_ratio)
    result$objective_mean_far_ratio[i]<-result$objective_mean_far_ratio[i]+mean(far_ratio)/(length(floors)*length(kernels))
  }
  cat("Planning point",num,"/",length(floors)*length(kernels),point,"floor",floor,"\n");flush.console()
}
result$feasible<-result$worst_observed_ratio<=1+1e-6
write.csv(result,file.path(out,"fisher_all_candidates.csv"),row.names=FALSE)
write.csv(do.call(rbind,base_risk),file.path(out,"fisher_baseline_node_variances.csv"),row.names=FALSE)
selected<-list();selection<-list()
for(class in c("single","dual")) {
  z<-result[result$class==class & result$feasible,]
  cat(class,"feasible",nrow(z),"of",sum(result$class==class),"\n")
  if(!nrow(z))next
  best<-min(z$objective_worst_far_ratio);z<-z[z$objective_worst_far_ratio<=best+1e-10,]
  z<-z[order(z$objective_mean_far_ratio,z$id),,drop=FALSE];winner<-z[1,,drop=FALSE]
  selected[[class]]<-grid$spec[[winner$id]];selection[[class]]<-winner
}
if(length(selection)) {
  write.csv(do.call(rbind,selection),file.path(out,"selected_designs.csv"),row.names=FALSE)
} else write.csv(result[FALSE,],file.path(out,"selected_designs.csv"),row.names=FALSE)
saveRDS(selected,file.path(out,"selected_designs.rds"))
selected_details<-list();j<-0L
for(class in names(selected))for(floor in floors)for(point in names(kernels))for(prior in c(TRUE,FALSE)) {
  d<-selected[[class]];ker<-kernels[[point]];smooth<-function(p){p<-pmax(p,floor);p/sum(p)}
  for(name in c("depth600",class)) {
    init<-if(name=="depth600")list(nominal$A)else d$initial
    depth<-if(name=="depth600")list(c(300,300))else d$depths
    risk<-constrained_risk(lapply(init,smooth),depth,ker,context,A,prior)
    U<-diag(n)-matrix(1/n,n,n);Cshape<-U%*%risk$covariance%*%U
    e<-nominal$edges;edgevar<-diag(risk$covariance)[e[,1]]+diag(risk$covariance)[e[,2]]-2*risk$covariance[e]
    for(distance in 0:2) {
      ii<-nominal$planned_distance==distance;ee<-pmax(nominal$planned_distance[e[,1]],nominal$planned_distance[e[,2]])==distance
      j<-j+1L;selected_details[[j]]<-data.frame(class=class,design=name,point=point,floor=floor,
        eta_prior=prior,distance=distance,relative_variance=mean(risk$variance[ii]),shape_variance=mean(diag(Cshape)[ii]),
        edge_variance=mean(edgevar[ee]))
    }
  }
}
if(length(selected_details))write.csv(do.call(rbind,selected_details),file.path(out,"selected_fisher_diagnostics.csv"),row.names=FALSE)
consumed<-c(list.files(src,full.names=TRUE),file.path(root,"results/context.rds"),
  file.path(root,"results/computation_source_v4/data.R"),file.path(repo,"validation/accuracy_round2/oracle_profile_fisher.R"),
  file.path(repo,"validation/accuracy_round2/common.R"))
consumed<-consumed[file.info(consumed)$isdir==FALSE]
write.csv(data.frame(file=consumed,sha256=unname(round2_sha256(consumed))),file.path(out,"planning_consumed_source_sha256.csv"),row.names=FALSE)
cat("Selected",length(selected),"candidate(s); no generating truth or actual estimate used in selection\n")
