#!/usr/bin/env Rscript
# The protocol prespecifies control300 to quantify doubled depth. This separate
# descriptive contrast requires all twelve complete pairs; it cannot rank the
# incomplete five-design study or qualify another design for promotion.
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-file.path(root,"results/analysis_input")
status<-read.csv(file.path(out,"all_fit_status.csv"))
metrics<-read.csv(file.path(out,"all_case_metrics.csv"))
designs<-c("control300","depth600")
s<-status[status$design%in%designs,]
stopifnot(nrow(s)==24L,all(s$status=="passed_quadrature_refinement"),
  length(unique(s$unit_id))==12L,length(unique(s$landscape_id))==6L,
  all(table(s$unit_id,s$design)==1L))
x<-metrics[metrics$design%in%designs & metrics$view=="baseline_distance",]
stopifnot(all(x$prediction_fraction==1),all(x$interval_fraction50==1),
  all(x$interval_fraction80==1),all(x$interval_fraction95==1),
  !anyDuplicated(paste(x$case_id,x$distance)),all(table(x$distance,x$design)==12L))
risks<-c("rmse","shape_rmse","edge_rmse","score50","score80","score95",
  "coverage50","coverage80","coverage95","width50","width80","width95")
stopifnot(all(is.finite(as.matrix(x[,risks]))))
land<-aggregate(x[,risks],x[,c("design","family","landscape_id","distance")],mean)
family<-aggregate(land[,risks],land[,c("design","family","distance")],mean)
overall<-aggregate(family[,risks],family[,c("design","distance")],mean)
a<-overall[overall$design=="control300",];b<-overall[overall$design=="depth600",]
b<-b[match(a$distance,b$distance),]
comparison<-data.frame(distance=a$distance,paired_biological_units=12L,independent_landscapes=6L)
for(metric in risks) {
  comparison[[paste0("control_",metric)]]<-a[[metric]]
  comparison[[paste0("depth_",metric)]]<-b[[metric]]
  comparison[[paste0(metric,"_percent_change")]]<-100*(b[[metric]]/a[[metric]]-1)
}
write.csv(comparison,file.path(out,"complete_depth_control_diagnostic.csv"),row.names=FALSE)
write.csv(family,file.path(out,"complete_depth_control_family_metrics.csv"),row.names=FALSE)
write.csv(land,file.path(out,"complete_depth_control_landscape_metrics.csv"),row.names=FALSE)
print(comparison[comparison$distance>=0,c("distance","rmse_percent_change","shape_rmse_percent_change",
  "edge_rmse_percent_change","score50_percent_change","score80_percent_change","score95_percent_change")],row.names=FALSE)
cat("Descriptive prespecified depth control only; all 12 pairs included; no five-design ranking or formal claim\n")
