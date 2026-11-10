#!/usr/bin/env Rscript
directory <- "validation/baseline_quadrature"
out <- file.path(directory,"results")
parts <- file.path(out,paste0("pilot_part",0:1))
status <- do.call(rbind,lapply(parts,function(p)read.csv(file.path(p,"status.csv"))))
case <- do.call(rbind,lapply(parts,function(p)read.csv(file.path(p,"case_metrics.csv"))))
if(nrow(status)!=12L || anyDuplicated(status$dataset_id))stop("Incomplete/duplicated pilot cases")
expected <- expand.grid(dataset_id=status$dataset_id,model=c("v1_joint","baseline_quadrature"),distance=c(-1L,0L,1L,2L))
key <- function(x)paste(x$dataset_id,x$model,x$distance)
if(!setequal(key(case),key(expected))||anyDuplicated(key(case)))stop("Incomplete model/distance risk grid")
risks <- c("rmse","mse","spatially_centered_rmse","spatially_centered_mse","edge_gradient_rmse",
  "interval_score95","coverage95","width95","global_offset_squared")
if(any(!is.finite(as.matrix(case[,risks])))||any(case$prediction_fraction!=1)||any(case$interval_fraction!=1))
  stop("Failed or incomplete predictions are present; no finite-mean deletion allowed")
landscape <- aggregate(case[,risks],case[,c("model","family","landscape_id","distance")],mean)
family <- aggregate(landscape[,risks],landscape[,c("model","family","distance")],mean)
overall <- aggregate(family[,risks],family[,c("model","distance")],mean)
base <- family[family$model=="v1_joint",]
candidate <- family[family$model=="baseline_quadrature",]
idx <- match(paste(candidate$family,candidate$distance),paste(base$family,base$distance))
ratios <- candidate[,c("family","distance")]
for(r in c("rmse","spatially_centered_rmse","edge_gradient_rmse","interval_score95"))ratios[[r]]<-candidate[[r]]/base[[r]][idx]
write.csv(status,file.path(out,"pilot_status.csv"),row.names=FALSE)
write.csv(case,file.path(out,"pilot_case_metrics.csv"),row.names=FALSE)
write.csv(landscape,file.path(out,"pilot_landscape_metrics.csv"),row.names=FALSE)
write.csv(family,file.path(out,"pilot_family_metrics.csv"),row.names=FALSE)
write.csv(overall,file.path(out,"pilot_overall_metrics.csv"),row.names=FALSE)
write.csv(ratios,file.path(out,"pilot_family_ratios.csv"),row.names=FALSE)
base_whole<-overall[overall$model=="v1_joint"&overall$distance== -1,]
quad_whole<-overall[overall$model=="baseline_quadrature"&overall$distance== -1,]
gain<-base_whole$mse-quad_whole$mse
shape_gain<-base_whole$spatially_centered_mse-quad_whole$spatially_centered_mse
offset_gain<-base_whole$global_offset_squared-quad_whole$global_offset_squared
stopifnot(abs(gain-shape_gain-offset_gain)<1e-12)
decision<-data.frame(status="NO_ALL_DISTANCE_SHAPE_UNCERTAINTY_WINNER",
  independent_test_opened=FALSE,independent_calibration_opened=FALSE,
  total_mse_gain=gain,shape_mse_gain=shape_gain,offset_squared_gain=offset_gain,
  offset_fraction_of_mse_gain=offset_gain/gain)
write.csv(decision,file.path(out,"pilot_stage_decision.csv"),row.names=FALSE)
num <- read.csv(file.path(out,"case_comparison_summary.csv"))
lines <- c("# Growth-scale integration: numerical and truth-development checks","",
  "The numerical candidate preserves the original Gaussian prior and integrates its existing independent intercept coefficient. It uses no additional growth observation in count-only fits. This separates a computational correction from changing the biological model or claiming new information.","",
  "## Small-graph numerical reference","",
  sprintf("All %d/16 cases pass quadrature refinement with every positive-weight inner node trusted. Against the independently checked MCMC references, relative-fitness SD ratios range %.3f–%.3f, versus %.3f–%.3f for the original Gaussian/Laplace approximation. The largest standardized mean error falls from %.3f to %.3f. Two nonzero-beta sparse count-only cases retain material mean discrepancies under the predeclared criterion; the approximation is improved but not exact.",sum(num$status=="passed_quadrature_refinement"),min(num$candidate_sd_ratio_min),max(num$candidate_sd_ratio_max),min(num$original_sd_ratio_min),max(num$original_sd_ratio_max),max(num$original_max_abs_mean_error_sd),max(num$candidate_max_abs_mean_error_sd)),"",
  "These are fixed-model computational comparisons using an explicitly fixed SD-2 initial-logit prior. They are not generating-truth RMSE, repeated-sampling coverage, SBC or a blanket deployment guarantee. [Reference report](../posterior_reference_4to8/REVIEW.md); the complete per-case numerical tables are retained in `results/`.","",
  "## Deployment-prior truth pilot","",
  "Twelve unchanged, revealed development samples represent six finite-population landscapes, two per family. Both estimators retain deployment eta priors, the same V1 fitness prior, complete 27-node domain and native support distances. Every V1 mean/covariance is replayed and every conditional eta-prior vector is checked against the frozen fit. Both models receive 40,000 actual joint-distribution draws for reference-contrast intervals, with common random numbers. No calibration is learned here.","",
  sprintf("All %d/12 integration fits pass numerical refinement. Median integration time is %.2fs per dataset with at most two fit processes. Input hashes, exact consumed datasets/frozen fits and node-level outputs are archived in the two pilot-part directories.",sum(status$status=="passed_quadrature_refinement"),median(status$elapsed_seconds)),"",
  "Risks average sampling repeats within landscape, landscapes within family and families equally. Distance -1 covers the whole target domain. Shape removes one spatial error mean over the whole domain, never an individual shell.","",
  "| Distance | Model | RMSE | Shape RMSE | Edge RMSE | Raw IS95 | Coverage95 |",
  "|---:|---|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(overall))) {
  x<-overall[i,];lines<-c(lines,sprintf("| %d | %s | %.5f | %.5f | %.5f | %.5f | %.3f |",x$distance,x$model,x$rmse,x$spatially_centered_rmse,x$edge_gradient_rmse,x$interval_score95,x$coverage95))
}
lines<-c(lines,"","## Descriptive family/distance safeguards","",
  "Ratios above 1.05 flag observed harm relative to frozen V1. These six development landscapes cannot establish formal noninferiority; all cells remain visible.","",
  "| Family | Distance | RMSE ratio | Shape ratio | Edge ratio | Raw IS95 ratio |",
  "|---|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(ratios))) {
  x<-ratios[i,];lines<-c(lines,sprintf("| %s | %d | %.4f | %.4f | %.4f | %.4f |",x$family,x$distance,x$rmse,x$spatially_centered_rmse,x$edge_gradient_rmse,x$interval_score95))
}
lines<-c(lines,"","## Decision boundaries","",
  sprintf("The stage remains closed: no all-distance shape/uncertainty winner. Distances 1 and 2 have slightly higher centered shape/edge error. Raw interval-score ratios exceed 1.05 in smooth-d0 (%.4f) and rough-d2 (%.4f). The common-reference offset reduction explains %.2f%% of the whole-domain MSE gain; the shape contribution is positive but small. This is a development stopping decision, not a formal independent-test rejection.",ratios$interval_score95[ratios$family=="smooth"&ratios$distance==0],ratios$interval_score95[ratios$family=="rough"&ratios$distance==2],100*offset_gain/gain),"",
  "A more faithful posterior computation does not automatically improve estimation of truths outside the model's structural prior. The pilot is used only to decide whether an all-distance shape/uncertainty claim is plausible enough for further development. No production default is promoted, no new accuracy calibration/test is opened, and no final test N/power is certified. Numerical correction alone cannot establish the original all-distance accuracy objective.","",
  "The most important remaining distinctions are (i) residual inner-Laplace error for rare initial states, (ii) biological-prior mismatch and (iii) information absent from sparse two-timepoint observations. These must remain separate in any subsequent method or observation-design change.","",
  paste("Generated:",format(Sys.time(),tz="America/New_York",usetz=TRUE)))
writeLines(lines,file.path(directory,"REVIEW.md"))
print(overall,row.names=FALSE)
cat("Wrote numerical/truth development report; no promotion asserted.\n")
