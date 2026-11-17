#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "metrics.R"))
args <- commandArgs(TRUE)
results_dir <- if (length(args)) args[[1L]] else file.path(script_dir, "results")
report_file <- if (length(args) > 1L) args[[2L]] else file.path(results_dir, "review.md")
status_files <- list.files(results_dir, "_status\\.csv$", full.names = TRUE)
files <- sub("_status\\.csv$", "_metrics.csv", status_files)
files <- files[file.exists(files)]
if (!length(files)) stop("No model metrics found")
metrics <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
statuses <- do.call(rbind, lapply(status_files, utils::read.csv, stringsAsFactors = FALSE))
metric_key <- c("dataset_id", "model", "distance_group", "metric")
duplicate <- duplicated(metrics[, metric_key])
if (any(duplicate)) {
  groups <- split(seq_len(nrow(metrics)), interaction(metrics[, metric_key], drop = TRUE))
  for (ii in groups) if (length(unique(metrics$value[ii])) > 1L) stop("Conflicting duplicate metrics; use distinct model labels/output directories for different implementations")
  metrics <- metrics[!duplicate, , drop = FALSE]
}
statuses <- statuses[!duplicated(statuses[, c("dataset_id", "model")]), , drop = FALSE]
target_count <- unique(metrics$value[metrics$metric == "n_target" & metrics$distance_group == "all_node_micro"])
target_count <- if (length(target_count) == 1L) as.character(target_count) else "a fixed set of"
node_file <- sub("_status\\.csv$", "_nodes.csv", status_files[[1L]])
example_nodes <- utils::read.csv(node_file, stringsAsFactors = FALSE)
chromosome_count <- length(strsplit(example_nodes$karyotype[[1L]], ".", fixed = TRUE)[[1L]])
landscape_count <- length(unique(metrics$landscape_id))
dataset_count <- length(unique(metrics$dataset_id))
paired <- paired_landscape_summary(metrics)
utils::write.csv(paired, file.path(results_dir, "paired_landscape_summary.csv"), row.names = FALSE)
utils::write.csv(metrics, file.path(results_dir, "all_model_metrics.csv"), row.names = FALSE)
aggregate_metrics <- aggregate(value ~ model + distance_group + metric, metrics, finite_mean, na.action = na.pass)
utils::write.csv(aggregate_metrics, file.path(results_dir, "aggregate_metrics.csv"), row.names = FALSE)
family_metrics <- aggregate(value ~ model + family + distance_group + metric, metrics, finite_mean, na.action = na.pass)
utils::write.csv(family_metrics, file.path(results_dir, "family_metrics.csv"), row.names = FALSE)
lines <- c("# alfak2 accuracy development benchmark", "",
           sprintf("These are development results on %d landscapes and %d paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.", landscape_count, dataset_count), "",
           sprintf("The finite target domain contains %s %d-chromosome karyotypes. Every method is scored on the same frozen target nodes and reference observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.", target_count, chromosome_count), "",
           "Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.", "",
           "Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. This small development sample is insufficient for a general claim of improvement at every distance.", "",
           "## Completion", "", "| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |", "|---|---:|---:|---:|---:|")
for (model in unique(statuses$model)) {
  s <- statuses[statuses$model == model, ]
  lines <- c(lines, sprintf("| %s | %d | %d | %d | %.2f |", model, nrow(s), sum(s$fit_status == "ok"),
                           sum(s$covariance_status %in% c("TMB_sdreport", "joint_laplace", "component_Laplace_mixture")), stats::median(s$runtime_seconds)))
}
lines <- c(lines, "", "## Fitness and interval results", "", "| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |", "|---|---|---:|---:|---:|---:|---:|")
for (model in unique(metrics$model)) for (distance in unique(metrics$distance_group)) {
  if (distance %in% c("local_target", "heldout_counts")) next
  x <- aggregate_metrics[aggregate_metrics$model == model & aggregate_metrics$distance_group == distance, ]
  get <- function(metric) {z <- x$value[x$metric == metric]; if (length(z)) z[1] else NA_real_}
  lines <- c(lines, sprintf("| %s | %s | %.4f | %.4f | %.3f | %.4f | %.3f |", model, distance,
                           get("relative_rmse"), get("relative_mae"), get("coverage_95"),
                           get("interval_score_95"), get("prediction_coverage")))
}
lines <- c(lines, "", "## Independent held-out count predictions", "",
           "Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.", "",
           "| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |", "|---|---:|---:|---:|---:|")
for (model in unique(metrics$model)) {
  x <- aggregate_metrics[aggregate_metrics$model == model & aggregate_metrics$distance_group == "heldout_counts", ]
  get <- function(metric) {z <- x$value[x$metric == metric]; if (length(z)) z[1] else NA_real_}
  lines <- c(lines, sprintf("| %s | %.4f | %.4f | %.5f | %.5f |", model,
                           get("heldout_t0_nll_per_cell_floor_1e12"), get("heldout_t1_nll_per_cell_floor_1e12"),
                           get("pi0_rmse_full_target"), get("pi1_rmse_full_target")))
}
distance_groups <- sort(unique(metrics$distance_group[grepl("^d[0-9]+$", metrics$distance_group)]))
lines <- c(lines, "", "## Fitness RMSE by truth family", "",
           paste0("| Model | Family | ", paste(distance_groups, collapse = " | "), " | Equal-distance macro |"),
           paste0("|---|---|", paste(rep("---:", length(distance_groups) + 1L), collapse = "|"), "|"))
for (model in unique(metrics$model)) for (family in unique(metrics$family)) {
  x <- family_metrics[family_metrics$model == model & family_metrics$family == family & family_metrics$metric == "relative_rmse", ]
  get <- function(group) {z <- x$value[x$distance_group == group]; if (length(z)) z[1] else NA_real_}
  values <- vapply(c(distance_groups, "all_distance_macro"), get, numeric(1))
  lines <- c(lines, paste0("| ", model, " | ", family, " | ", paste(sprintf("%.4f", values), collapse = " | "), " |"))
}
if (nrow(paired)) {
  lines <- c(lines, "", "## Paired candidate minus baseline differences", "",
             "Negative RMSE and interval-score differences favor the candidate. Coverage differences have no universal favorable direction; assess distance from the nominal 0.95 and interval score jointly.", "",
             "| Candidate | Distance | Metric | Mean difference | 95% landscape bootstrap interval | Landscapes |", "|---|---|---|---:|---|---:|")
  p <- paired[paired$metric %in% c("relative_rmse", "interval_score_95", "coverage_95") & paired$distance_group != "all_node_micro", ]
  for (i in seq_len(nrow(p))) lines <- c(lines, sprintf("| %s | %s | %s | %.4f | [%.4f, %.4f] | %d |", p$model[i], p$distance_group[i], p$metric[i], p$paired_delta[i], p$delta_low[i], p$delta_high[i], p$n_landscapes[i]))
}
lines <- c(lines, "", "## Limits and next gate", "",
           "The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.", "",
           "Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.", "",
           paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, report_file)
artifact_dir <- file.path(dirname(report_file), "results")
if (basename(report_file) != "pilot_review.md") artifact_dir <- file.path(artifact_dir, tools::file_path_sans_ext(basename(report_file)))
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(aggregate_metrics, file.path(artifact_dir, "aggregate_metrics.csv"), row.names = FALSE)
utils::write.csv(family_metrics, file.path(artifact_dir, "family_metrics.csv"), row.names = FALSE)
utils::write.csv(paired, file.path(artifact_dir, "paired_landscape_summary.csv"), row.names = FALSE)
utils::write.csv(statuses, file.path(artifact_dir, "fit_status.csv"), row.names = FALSE)
cat(normalizePath(report_file), "\n")
