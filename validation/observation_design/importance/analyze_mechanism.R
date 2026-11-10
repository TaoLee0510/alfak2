#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-file.path(normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT")),"mechanism")
src<-normalizePath(Sys.getenv("ALFAK2_IS_ANALYSIS_SOURCE"));source(file.path(src,"score.R"))
manifest<-read.csv(file.path(out,"manifest.csv"));stopifnot(nrow(manifest)==48)
statuses<-metrics<-nodes<-audits<-list()
for(i in seq_len(nrow(manifest))) {
 r<-manifest[i,];case<-readRDS(r$scoring_file)
 stopifnot(digest::digest(file=r$file,algo="sha256")==r$sha256,digest::digest(file=r$scoring_file,algo="sha256")==r$scoring_sha256)
 resource<-read.csv(file.path(r$output,"resource.csv"));rp<-tolower(as.character(resource$resource_pass))=="true"
 result<-file.path(r$output,"result.rds");fit<-if(file.exists(result))readRDS(result)else list(status="unavailable")
 valid<-fit$status=="passed_joint_importance"&&rp&&resource$exit_code==0
 statuses[[i]]<-data.frame(case_id=r$case_id,arm=r$arm,unit_id=r$unit_id,status=fit$status,
  numerical_pass=valid,resource_pass=rp,elapsed_seconds=resource$child_seconds,
  maximum_rss_gib=resource$maximum_rss_bytes/1024^3,ess=if(is.null(fit$ess))NA else fit$ess,
  max_pareto_k=if(is.null(fit$pareto_k))NA else max(fit$pareto_k),
  max_mean_mcse_ratio=if(is.null(fit$max_mean_mcse_ratio))NA else fit$max_mean_mcse_ratio,
  max_sd_mcse_ratio=if(is.null(fit$max_sd_mcse_ratio))NA else fit$max_sd_mcse_ratio)
 if(!valid)fit$status<-"unavailable"
 s<-design_score(case,fit);metrics[[i]]<-s$metrics;nodes[[i]]<-s$nodes
 if(valid) {
  A<-diag(length(fit$mean))-matrix(case$reference_weight,length(fit$mean),length(fit$mean),byrow=TRUE)
  er<-max(abs(s$nodes$pred_relative-as.numeric(A%*%fit$mean)))
  stopifnot(er<1e-10,all(fit$flags),all(apply(fit$quantiles$relative,1,function(q)all(diff(q)>=0))))
  full<-s$metrics[s$metrics$distance==-1,]
  decomp<-max(abs(full$mse-full$shape_mse-full$offset_squared));stopifnot(decomp<1e-12)
  audits[[length(audits)+1L]]<-data.frame(case_id=r$case_id,reference_error=er,mse_decomposition_error=decomp,
   original_prior=TRUE,weighted_quantiles=TRUE)
 }
 rm(fit);invisible(gc())
}
status<-do.call(rbind,statuses);m<-do.call(rbind,metrics)
write.csv(status,file.path(out,"all_fit_status.csv"),row.names=FALSE)
write.csv(m,file.path(out,"all_case_metrics.csv"),row.names=FALSE)
write.csv(do.call(rbind,nodes),file.path(out,"all_node_predictions.csv"),row.names=FALSE)
write.csv(do.call(rbind,audits),file.path(out,"posterior_integrity_audit.csv"),row.names=FALSE)
if(!all(status$numerical_pass))stop("Numerical/resource gate incomplete: all 48 cases retained, no complete aggregate")
stopifnot(all(m$prediction_fraction==1),all(m$interval_fraction50==1),all(m$interval_fraction80==1),all(m$interval_fraction95==1))
primary<-m[m$view=="baseline_distance",]
risks<-c("rmse","mse","shape_rmse","shape_mse","offset_squared","edge_rmse","score50","score80","score95", "coverage50","coverage80","coverage95","width50","width80","width95")
stopifnot(all(is.finite(as.matrix(primary[,risks]))))
land<-aggregate(primary[,risks],primary[,c("design","family","landscape_id","distance")],mean)
family<-aggregate(land[,risks],land[,c("design","family","distance")],mean)
overall<-aggregate(family[,risks],family[,c("design","distance")],mean)
write.csv(land,file.path(out,"landscape_metrics.csv"),row.names=FALSE)
write.csv(family,file.path(out,"family_metrics.csv"),row.names=FALSE)
write.csv(overall,file.path(out,"overall_metrics.csv"),row.names=FALSE)
comparisons<-list(process_unknown=c("B_mean_unknown","A_finite_unknown"),initial_finite=c("C_finite_oracle","A_finite_unknown"),
 process_oracle=c("D_mean_oracle","C_finite_oracle"),initial_mean=c("D_mean_oracle","B_mean_unknown"))
paired<-family_ratios<-overall_changes<-list()
for(name in names(comparisons)) {
 pair<-comparisons[[name]];a<-land[land$design==pair[1],];b<-land[land$design==pair[2],]
 b<-b[match(paste(a$landscape_id,a$distance),paste(b$landscape_id,b$distance)),]
 x<-a[,c("family","landscape_id","distance")];x$contrast<-name;x$oracle_probe<-any(grepl("oracle",pair))
 for(k in risks)x[[k]]<-a[[k]]-b[[k]]
 paired[[name]]<-x
 a<-family[family$design==pair[1],];b<-family[family$design==pair[2],];b<-b[match(paste(a$family,a$distance),paste(b$family,b$distance)),]
 x<-a[,c("family","distance")];x$contrast<-name;x$oracle_probe<-any(grepl("oracle",pair))
 for(k in risks)x[[k]]<-a[[k]]/b[[k]]
 family_ratios[[name]]<-x
 a<-overall[overall$design==pair[1],];b<-overall[overall$design==pair[2],];b<-b[match(a$distance,b$distance),]
 x<-data.frame(contrast=name,distance=a$distance,oracle_probe=any(grepl("oracle",pair)))
 for(k in risks)x[[k]]<-100*(a[[k]]/b[[k]]-1)
 overall_changes[[name]]<-x
}
write.csv(do.call(rbind,paired),file.path(out,"paired_landscape_differences.csv"),row.names=FALSE)
write.csv(do.call(rbind,family_ratios),file.path(out,"family_risk_ratios.csv"),row.names=FALSE)
write.csv(do.call(rbind,overall_changes),file.path(out,"overall_risk_changes.csv"),row.names=FALSE)
interaction<-land[land$design=="A_finite_unknown",c("family","landscape_id","distance")]
for(k in risks) {
 get<-function(arm){z<-land[land$design==arm,];z[[k]][match(paste(interaction$landscape_id,interaction$distance),paste(z$landscape_id,z$distance))]}
 interaction[[k]]<-get("D_mean_oracle")-get("C_finite_oracle")-get("B_mean_unknown")+get("A_finite_unknown")
}
write.csv(interaction,file.path(out,"paired_interaction.csv"),row.names=FALSE)
old<-read.csv(file.path(root,"results/constrained_v2/all_case_metrics.csv"));old<-old[old$design=="opt_single"&old$view=="baseline_distance",]
current<-primary[primary$design=="A_finite_unknown",];old<-old[match(paste(current$unit_id,current$distance),paste(old$unit_id,old$distance)),]
x<-current[,c("unit_id","family","landscape_id","distance")]
for(k in risks)x[[k]]<-current[[k]]-old[[k]]
write.csv(x,file.path(out,"numerical_correction_paired_differences.csv"),row.names=FALSE)
writeLines("ALL_48_NUMERICAL_RESOURCE_GATES_PASSED; descriptive six-landscape mechanism study, oracle arms not deployable",file.path(out,"stage_status.txt"))
print(overall[,c("design","distance","rmse","shape_rmse","edge_rmse","score95","coverage95","width95")],row.names=FALSE)
