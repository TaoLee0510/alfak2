#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "metrics.R"))
args <- commandArgs(TRUE)
if (length(args) < 4L) stop("Usage: write_final_exact_review.R TEST_RESULTS TEST_DATA FROZEN_CALIBRATION OUTPUT_MD")
results_dir <- args[[1L]]; data_dir <- args[[2L]]; calibration_dir <- args[[3L]]; report_file <- args[[4L]]
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(all(manifest$split == "test"), all(manifest$generating_protocol == "v2_exact_continuous_division"))
status_files <- list.files(results_dir, "_status\\.csv$", full.names = TRUE)
statuses <- do.call(rbind, lapply(status_files, utils::read.csv, stringsAsFactors = FALSE))
statuses <- statuses[!duplicated(statuses[, c("dataset_id", "model")]), , drop = FALSE]
metrics <- do.call(rbind, lapply(sub("_status\\.csv$", "_metrics.csv", status_files), utils::read.csv, stringsAsFactors = FALSE))
metrics <- metrics[!duplicated(metrics[, c("dataset_id", "model", "distance_group", "metric")]), , drop = FALSE]
required <- c("baseline", "division_v1", "baseline_calibrated", "division_v1_calibrated")
for (model in required) if (!setequal(statuses$dataset_id[statuses$model == model], manifest$dataset_id)) stop(paste("Final test incomplete:", model))
freeze <- readRDS(file.path(calibration_dir, "freeze_record.rds"))
stopifnot(identical(freeze$test_generated, FALSE))
frozen_manifest <- utils::read.csv(file.path(calibration_dir, "frozen_calibration_manifest.csv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(frozen_manifest))) {
  file <- file.path(calibration_dir, basename(frozen_manifest$object_file[i]))
  stopifnot(identical(unname(tools::md5sum(file)), frozen_manifest$md5[i]))
}
landscape_metrics <- aggregate(value ~ model + family + landscape_id + distance_group + metric,
                                metrics, finite_mean, na.action = na.pass)
aggregate_metrics <- aggregate(value ~ model + distance_group + metric, landscape_metrics, finite_mean, na.action = na.pass)
family <- aggregate(value ~ model + family + distance_group + metric, landscape_metrics, finite_mean, na.action = na.pass)
paired_raw <- paired_landscape_summary(metrics)
paired_calibrated <- paired_landscape_summary(metrics, baseline_model = "baseline_calibrated")
artifact_dir <- file.path(dirname(report_file), "results", "final_exact")
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(aggregate_metrics, file.path(artifact_dir, "aggregate_metrics.csv"), row.names = FALSE)
utils::write.csv(family, file.path(artifact_dir, "family_metrics.csv"), row.names = FALSE)
utils::write.csv(paired_raw, file.path(artifact_dir, "paired_against_raw_baseline.csv"), row.names = FALSE)
utils::write.csv(paired_calibrated, file.path(artifact_dir, "paired_against_calibrated_baseline.csv"), row.names = FALSE)
utils::write.csv(statuses, file.path(artifact_dir, "fit_status.csv"), row.names = FALSE)
utils::write.csv(manifest, file.path(artifact_dir, "test_dataset_manifest.csv"), row.names = FALSE)
node_files <- sub("_status\\.csv$", "_nodes.csv", status_files)
node_tables <- lapply(node_files, utils::read.csv, stringsAsFactors = FALSE)
all_nodes <- do.call(rbind, lapply(node_tables, function(x) x[, c("dataset_id", "model", "karyotype", "distance", "truth_raw", "truth_relative", "pred_raw", "pred_relative", "sd_raw", "sd_relative")]))
all_nodes <- all_nodes[!duplicated(all_nodes[, c("dataset_id", "model", "karyotype")]), , drop = FALSE]
labels <- readRDS(manifest$dataset_file[1L])$target_labels
matrix_from <- function(x, value) {
  out <- matrix(NA_real_, length(labels), nrow(manifest), dimnames = list(labels, manifest$dataset_id))
  ii <- cbind(match(x$karyotype, labels), match(x$dataset_id, manifest$dataset_id))
  out[ii] <- x[[value]]
  out
}
reference <- all_nodes[all_nodes$model == "baseline", ]
node_matrices <- list(target_labels = labels, dataset_ids = manifest$dataset_id,
  truth_relative = matrix_from(reference, "truth_relative"), truth_birth_rate = matrix_from(reference, "truth_raw"),
  reference_support_distance = matrix_from(reference, "distance"), models = list(),
  estimand = "fixed t0-count-weighted observed-reference fitness contrast",
  scenario = freeze$scenario)
for (model in unique(statuses$model)) {
  tables <- node_tables[vapply(node_tables, function(x) identical(unique(x$model), model), logical(1))]
  x <- do.call(rbind, tables)
  x <- x[!duplicated(x[, c("dataset_id", "karyotype")]), , drop = FALSE]
  values <- intersect(c("pred_raw", "pred_relative", "sd_raw", "sd_relative",
    paste0("relative_low_", c("50", "80", "95")), paste0("relative_high_", c("50", "80", "95")),
    paste0("calibrated_low_", c("50", "80", "95")), paste0("calibrated_high_", c("50", "80", "95"))), names(x))
  node_matrices$models[[model]] <- stats::setNames(lapply(values, function(value) matrix_from(x, value)), values)
}
for (model in c("baseline", "division_v1", "endpoint_v1")) {
  stopifnot(identical(node_matrices$models[[model]]$pred_relative, node_matrices$models[[paste0(model, "_calibrated")]]$pred_relative),
            identical(node_matrices$models[[model]]$sd_relative, node_matrices$models[[paste0(model, "_calibrated")]]$sd_relative))
}
saveRDS(node_matrices, file.path(artifact_dir, "node_matrices.rds"))
groups <- list()
for (model in unique(statuses$model)) for (group_type in c("family", "landscape_id")) {
  for (group in unique(manifest[[group_type]])) {
    ids <- manifest$dataset_id[manifest[[group_type]] == group]
    s <- statuses[statuses$model == model & statuses$dataset_id %in% ids, ]
    groups[[length(groups) + 1L]] <- data.frame(model = model, group_type = group_type, group = group,
      expected_datasets = length(ids), returned_datasets = nrow(s), successful_fits = sum(s$fit_status == "ok"))
  }
}
utils::write.csv(do.call(rbind, groups), file.path(artifact_dir, "group_completeness.csv"), row.names = FALSE)
get <- function(model, distance, metric) {
  z <- aggregate_metrics$value[aggregate_metrics$model == model & aggregate_metrics$distance_group == distance & aggregate_metrics$metric == metric]
  if (length(z)) z[1] else NA_real_
}
distances <- sort(unique(metrics$distance_group[grepl("^d[0-9]+$", metrics$distance_group)]))
gain <- vapply(distances, function(d) get("division_v1", d, "relative_rmse") < get("baseline", d, "relative_rmse"), logical(1))
complete_predictions <- all(vapply(c("baseline", "division_v1"), function(model)
  all(metrics$value[metrics$model == model & metrics$metric == "prediction_coverage"] == 1), logical(1)))
outcome <- if (all(gain) && complete_predictions) {
  point_ci <- paired_raw[paired_raw$model == "division_v1" & paired_raw$metric == "relative_rmse" & paired_raw$distance_group %in% distances, ]
  includes_zero <- sum(point_ci$delta_low <= 0 & point_ci$delta_high >= 0)
  sprintf("The frozen joint division candidate lowers average relative-fitness RMSE at every evaluated distance, but %d/%d distance-specific paired 95%% intervals still include zero. Family-specific regressions and calibrated interval-score tradeoffs remain. These results do not establish uniformly better estimation for every karyotype or scenario.", includes_zero, length(distances))
} else {
  "The frozen candidate does not meet the criterion of lower average relative-fitness RMSE at every evaluated distance with complete prediction coverage. Distance/family regressions remain visible below; the candidate and calibration were not retuned after this test."
}
lines <- c("# alfak2: final independent exact-division validation", "", outcome, "",
           "The strict all-distance point-and-uncertainty promotion gate is not passed. The joint candidate remains an explicit opt-in; the default pipeline remains two_layer. No tuning or extra test repetitions were added after these results.", "",
           sprintf("The test contains %d independent landscapes and %d datasets: three truth families, four reserved landscapes per family, two biological trajectories and two independent count samples per landscape. The entire 81-node four-chromosome CN1–3 domain is scored. Beta=0.005, interval=1.5, initial population=1,200 and count depths=150/150 are fixed.", length(unique(manifest$landscape_id)), nrow(manifest)), "",
           "## Point estimation: raw candidate versus raw baseline", "",
           "All models use one t0-count-weighted observed-reference contrast. No per-distance centering, affine truth correction or calibrated change to a point prediction is used. Reported means first average biological/sampling repeats within each landscape, then weight landscapes equally. The confidence intervals resample those independent paired landscapes, with the prespecified 1,000-replicate bootstrap. They condition on the frozen candidate and calibration. One test dataset has no d3 target nodes; this is an absent stratum, not a missing prediction.", "",
           "| Distance | Original RMSE | Joint division RMSE | RMSE reduction | Candidate minus original [95% landscape interval] | Paired landscapes/datasets |", "|---|---:|---:|---:|---|---:|")
for (distance in distances) {
  base <- get("baseline", distance, "relative_rmse"); candidate <- get("division_v1", distance, "relative_rmse")
  ci <- paired_raw[paired_raw$model == "division_v1" & paired_raw$distance_group == distance & paired_raw$metric == "relative_rmse", ]
  if (nrow(ci)) lines <- c(lines, sprintf("| %s | %.5f | %.5f | %.1f%% | %.5f [%.5f, %.5f] | %d/%d |", distance, base, candidate, 100 * (1 - candidate / base), ci$paired_delta[1], ci$delta_low[1], ci$delta_high[1], ci$n_landscapes[1], ci$n_paired_datasets[1]))
}
lines <- c(lines, "", "## Uncertainty: original and candidate, before and after independent calibration", "",
           "Both methods receive their own empirical calibration on the same separate six-landscape/24-dataset exact-simulation calibration split. Calibration changes interval endpoints only. It leaves every fitness prediction and SD unchanged. Proper interval scores penalize both excessive width and missed truth; lower is better. Generic calibration gains are not attributed to the new estimator.", "",
           "| Distance | Model | 95% coverage | 95% interval score | Interval availability |", "|---|---|---:|---:|---:|")
for (distance in distances) for (model in required) {
  lines <- c(lines, sprintf("| %s | %s | %.3f | %.5f | %.3f |", distance, model,
                           get(model, distance, "coverage_95"), get(model, distance, "interval_score_95"), get(model, distance, "valid_interval_fraction")))
}
lines <- c(lines, "", "The matched-calibration comparison is candidate_calibrated minus baseline_calibrated:", "",
           "| Distance | Proper interval-score difference [95% landscape interval] | Paired landscapes |", "|---|---|---:|")
for (distance in distances) {
  ci <- paired_calibrated[paired_calibrated$model == "division_v1_calibrated" & paired_calibrated$distance_group == distance & paired_calibrated$metric == "interval_score_95", ]
  if (nrow(ci)) lines <- c(lines, sprintf("| %s | %.5f [%.5f, %.5f] | %d |", distance, ci$paired_delta[1], ci$delta_low[1], ci$delta_high[1], ci$n_landscapes[1]))
}
lines <- c(lines, "", "## Family-specific point errors", "",
           "| Truth family | Distance | Original RMSE | Joint division RMSE | RMSE reduction |", "|---|---|---:|---:|---:|")
regressions <- character()
for (name in unique(manifest$family)) for (distance in distances) {
  select <- function(model) {z <- family$value[family$model == model & family$family == name & family$distance_group == distance & family$metric == "relative_rmse"]; if (length(z)) z[1] else NA_real_}
  base <- select("baseline"); candidate <- select("division_v1")
  lines <- c(lines, sprintf("| %s | %s | %.5f | %.5f | %.1f%% |", name, distance, base, candidate, 100 * (1 - candidate / base)))
  if (is.finite(candidate) && is.finite(base) && candidate > base) regressions <- c(regressions, paste(name, distance))
}
lines <- c(lines, "", if (length(regressions)) paste0("Remaining family/distance RMSE regressions: ", paste(regressions, collapse = ", "), ".") else "No family/distance average RMSE regression is observed in this test; individual-node improvements are not guaranteed.", "",
           "## Count predictions and completion", "",
           "Held-out multinomial count samples use a separate fixed seed and never enter fitting. These are per-cell log scores of the returned point probability forecasts, with a declared 1e-12 floor; they are not a fully integrated joint count predictive likelihood. Zero-probability hits and frequency RMSE are retained in the CSV.", "",
           "| Model | t0 NLL/cell | t1 NLL/cell | Successful fits | Trusted local covariance |", "|---|---:|---:|---:|---:|")
for (model in c("baseline", "endpoint_v1", "division_v1")) {
  s <- statuses[statuses$model == model, ]
  lines <- c(lines, sprintf("| %s | %.5f | %.5f | %d/%d | %d/%d |", model,
                           get(model, "heldout_counts", "heldout_t0_nll_per_cell_floor_1e12"),
                           get(model, "heldout_counts", "heldout_t1_nll_per_cell_floor_1e12"),
                           sum(s$fit_status == "ok"), nrow(s), sum(s$covariance_status == "TMB_sdreport"), nrow(s)))
}
lines <- c(lines, "", "Failure and missing-interval cases remain in fit-status and availability denominators. Error averages and paired intervals use finite comparable predictions, so incomplete availability cannot be interpreted as successful universal accuracy.", "",
           "## Freeze record and limits", "",
           paste("Calibration frozen at:", freeze$frozen_at_utc, "; the saved freeze record has test_generated=FALSE. The test landscape and sampling seeds were reserved before development. Candidate-v1 source/DLL and calibration-object checksums are retained in [the protocol artifact directory](results/protocol_v2/)."), "",
           "The initial tau-leap simulator has measurable finite-step growth bias. Formal calibration/test therefore use the independently implemented exact Gillespie process, chosen before test generation. The historical tau outcomes are retained as robustness evidence rather than overwritten. See [protocol v2](protocol_v2.md) and [simulator validation](results/simulator_validation.csv).", "",
           "This result is confined to the specified simulation families, depth, beta, time interval and bounded target domain. It does not validate a complete 22-chromosome landscape or real biological measurements. Absolute birth-rate level remains prior-identified without external growth data. Calibration uses six independent landscapes and provides empirical intervals only, without a finite-sample or distribution-free coverage guarantee. Unknown distances/scenarios are explicitly unvalidated.", "",
           "## Retained negative results and audit trail", "",
           "- [Original pilot and Laplace ablation](pilot_review.md): repaired convergence did not automatically improve every fitness error.",
           "- [Shared-scale and evidence-mixture screen](scale_grid_review.md): point/interval and truth-family tradeoffs remain visible.",
           "- [Separated basis/residual screen](basis_residual_review.md) and [eight-dimensional Pareto table](results/basis_residual_pareto.csv).",
           "- [Frozen tau robustness comparison](frozen_train_review.md): retained separately from the exact independent test.",
           "- [Five-chromosome d0–d5 support-limit probe](farfield_review.md): remote states are not identified by one founder-only observation.",
           "- [Machine-readable final results](results/final_exact/): paired raw and matched-calibration comparisons, all distance/family metrics, fit statuses, group completeness, the reserved test manifest and complete node-wise prediction/SD/interval matrices. Matrix checks verify calibrated and raw point predictions/SDs are identical.", "",
           paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, report_file)
cat(normalizePath(report_file), "\n")
