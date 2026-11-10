#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
input <- round2_arg("results"); report <- round2_arg("report", file.path(script_dir, "candidate_pilot_review.md"))
if (is.null(input)) stop("Supply --results")
inputs <- strsplit(input, ",", fixed = TRUE)[[1L]]
case <- do.call(rbind, lapply(file.path(inputs, "case_metrics.csv"), utils::read.csv, stringsAsFactors = FALSE))
status <- do.call(rbind, lapply(file.path(inputs, "fit_status.csv"), utils::read.csv, stringsAsFactors = FALSE))
if (!"v1_joint" %in% case$model) stop("V1 comparator is required")
expected_ids <- sort(unique(case$dataset_id[case$model == "v1_joint"]))
for (model in unique(case$model)) {
  if (!identical(sort(unique(case$dataset_id[case$model == model])), expected_ids))
    stop(paste("INCOMPLETE_MODEL_DATASET_SET:", model, "-- do not compare partial means"))
  for (id in expected_ids) {
    reference_distances <- sort(unique(case$distance[case$model == "v1_joint" & case$dataset_id == id]))
    if (!identical(sort(unique(case$distance[case$model == model & case$dataset_id == id])), reference_distances))
      stop(paste("MISSING_NATIVE_DISTANCE_ROWS:", model, id))
  }
}
node_file <- file.path(inputs[1L], "node_tables.rds")
first_nodes <- readRDS(node_file)[[1L]]
n_target <- length(unique(first_nodes$karyotype))
n_chr <- length(strsplit(first_nodes$karyotype[1L], ".", fixed = TRUE)[[1L]])
landscapes_per_family <- table(unique(case[, c("family", "landscape_id")])$family)
measure <- setdiff(names(case), c("dataset_id", "landscape_id", "family", "split", "model", "distance"))
land <- aggregate(case[, measure], case[, c("model", "family", "landscape_id", "distance")], round2_mean)
family <- aggregate(land[, measure], land[, c("model", "family", "distance")], round2_mean)
overall <- aggregate(family[, measure], family[, c("model", "distance")], round2_mean)
output <- file.path(dirname(report), "results", if (n_target == 27L) "candidate_pilot" else "development81_comparison")
round2_write_csv(case, file.path(output, "case_metrics.csv")); round2_write_csv(status, file.path(output, "fit_status.csv"))
round2_write_csv(land, file.path(output, "landscape_metrics.csv")); round2_write_csv(family, file.path(output, "family_metrics.csv"))
round2_write_csv(overall, file.path(output, "overall_metrics.csv"))
descriptive_denominators <- list()
for (model in unique(case$model)) for (d in sort(unique(case$distance))) for (metric in c("spearman", "pearson", "amplitude_ratio", "variance_ratio")) {
  z <- case[case$model == model & case$distance == d, ]
  defined <- is.finite(z[[metric]])
  descriptive_denominators[[length(descriptive_denominators) + 1L]] <- data.frame(model = model,
    distance = d, metric = metric, expected_datasets = nrow(z), defined_datasets = sum(defined),
    expected_landscapes = length(unique(z$landscape_id)),
    defined_landscapes = length(unique(z$landscape_id[defined])),
    reason_for_undefined = "too few nodes, constant prediction, or zero truth spread; not a missing fitted prediction")
}
descriptive_denominators <- do.call(rbind, descriptive_denominators)
round2_write_csv(descriptive_denominators, file.path(output, "descriptive_metric_denominators.csv"))
get <- function(model, d, risk) {v <- overall[[risk]][overall$model == model & overall$distance == d]; if (length(v)) v[1L] else NA_real_}
lines <- c("# Round 2 finite-population candidate development comparison", "",
  sprintf("Development screening only: %d chromosomes, all %d target nodes, %d independent landscapes. Scores use each count sample's native observed-support distances, not common process-pair distances. Every old calibration/revealed-test input is explicitly relabeled development. No new calibration/test landscape is opened.", n_chr, n_target, length(unique(case$landscape_id))), "",
  "V1 is explicitly isotropic with the frozen small prior scales. A changes the interaction prior to its observed-design nullspace with unidentified SD 0.2. B adds the predeclared weak/rough residual components with prior weights 0.75/0.25. All use deployment eta priors, no population-growth/pi0 oracle and no generating-family model input. The current default two-layer API and an interval-free zero point predictor are separate controls.", "",
  "Intervals for B sample the actual component mixture and then apply the shared reference contrast. Family labels are evaluation metadata only. Means first average count samples within a landscape, then landscapes within a family, then families equally.", "",
  "## Completion", "", "| Model | Returned | Successful | Trusted covariance | Median seconds |", "|---|---:|---:|---:|---:|")
for (model in unique(status$model)) {
  z <- status[status$model == model, ]
  lines <- c(lines, sprintf("| %s | %d | %d | %d | %.3f |", model, nrow(z), sum(z$fit_status %in% c("ok", "analytic_zero", "reused_frozen_v1")), sum(z$covariance_status == "TMB_sdreport"), stats::median(z$elapsed_seconds)))
}
lines <- c(lines, "", "## Native-distance risk and shape", "",
  "Spatial centering removes one error mean over the entire target domain before stratification. It never recenters individual shells. All unordered one-copy graph edges enter the gradient metric, assigned to the largest endpoint support distance. Distance -1 denotes the complete target domain.", "",
  "| Distance | Model | RMSE | Centered shape RMSE | Edge gradient RMSE | Spearman | Amplitude ratio | 95% coverage | 95% interval score | Interval fraction |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for (d in sort(unique(overall$distance))) for (model in unique(status$model)) {
  lines <- c(lines, sprintf("| %d | %s | %.5f | %.5f | %.5f | %.3f | %.3f | %.3f | %.5f | %.3f |", d, model,
    get(model, d, "rmse"), get(model, d, "spatially_centered_rmse"), get(model, d, "edge_gradient_rmse"),
    get(model, d, "spearman"), get(model, d, "amplitude_ratio"), get(model, d, "coverage95"),
    get(model, d, "interval_score95"), get(model, d, "interval_fraction")))
}
lines <- c(lines, "", "Rank/correlation and amplitude/variance ratios average only their defined values; a small stratum or constant predictor can make them undefined while all point predictions remain available. Exact dataset/landscape denominators are retained in `descriptive_metric_denominators.csv`.")
if (n_target == 81L) lines <- c(lines,
  "At native d3, Spearman is defined for 15/18 datasets and amplitude for 16/18 in each fitted model. The zero predictor has undefined rank at every distance. RMSE, centered shape RMSE, edge RMSE and proper interval score in the four fitted models each retain all 18/18 datasets at d0–d3.")
lines <- c(lines, "", "## Candidate versus V1 family safety screen", "",
  "These are descriptive ratios, not a formal noninferiority test. Ratios above 1.05 flag the predeclared 5% family safety margin. All cells remain visible; this pilot does not certify the final test sample size.", "",
  "| Family | Distance | Candidate | RMSE ratio | Shape RMSE ratio | Edge RMSE ratio | Proper interval-score ratio |", "|---|---:|---|---:|---:|---:|---:|")
comparisons <- list()
for (model in intersect(c("A", "B"), unique(status$model))) for (f in unique(family$family)) for (d in sort(unique(family$distance))) {
  a <- family[family$model == model & family$family == f & family$distance == d, ]
  b <- family[family$model == "v1_joint" & family$family == f & family$distance == d, ]
  if (!nrow(a) || !nrow(b)) next
  risk <- c("rmse", "spatially_centered_rmse", "edge_gradient_rmse", "interval_score95")
  ratios <- unlist(a[, risk], use.names = FALSE) / unlist(b[, risk], use.names = FALSE)
  comparisons[[length(comparisons) + 1L]] <- data.frame(model = model, family = f, distance = d,
    risk = risk, ratio_to_v1 = ratios, exceeds_5percent_margin = ratios > 1.05)
  lines <- c(lines, sprintf("| %s | %d | %s | %.4f | %.4f | %.4f | %.4f |", f, d, model, ratios[1], ratios[2], ratios[3], ratios[4]))
}
if (length(comparisons)) round2_write_csv(do.call(rbind, comparisons), file.path(output, "ratios_to_v1.csv"))
availability <- do.call(rbind, lapply(unique(status$model), function(model) {
  z <- case[case$model == model, ]; s <- status[status$model == model, ]
  data.frame(model = model, expected_datasets = length(expected_ids), returned_datasets = nrow(s),
    all_fits_returned = all(s$fit_status %in% c("ok", "reused_frozen_v1", "analytic_zero")),
    all_predictions_available = all(is.finite(z$prediction_fraction) & z$prediction_fraction == 1),
    all_intervals_available = all(is.finite(z$interval_fraction) & z$interval_fraction == 1),
    formal_promotion_evaluated = FALSE)
}))
round2_write_csv(availability, file.path(output, "availability_gate.csv"))
weights <- list()
for (file in unlist(lapply(inputs, list.files, pattern = "^B_.*[.]rds$", full.names = TRUE))) {
  if (grepl("provenance[.]rds$", file)) next
  result <- readRDS(file)
  if (is.null(result$fit)) next
  w <- result$fit$global$hyperparameters$model_weights
  if (is.null(w)) next
  weights[[length(weights) + 1L]] <- data.frame(dataset_id = result$status$dataset_id,
    family = result$status$family, landscape_id = result$status$landscape_id,
    component = if (is.null(names(w))) seq_along(w) else names(w), weight = as.numeric(w))
}
if (length(weights)) round2_write_csv(do.call(rbind, weights), file.path(output, "B_posterior_component_weights.csv"))
lines <- c(lines, "", "This is a development comparison, not an independent test or SBC study. Reused V1 results retain their frozen source-cache SHA and unchanged counts/reference; zero elapsed time denotes reuse, not a newly timed fit. The current default is a separate actual API control, not an endpoint-ablation surrogate. No prior weight is adjusted from a sealed test. A single winner, if any, must be frozen after development and evaluated once on new calibration/test landscapes. Calibration hyperparameter heldout scores used to select gamma/kappa are not by themselves nested OOF loss vectors for joint power certification.", "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE); writeLines(lines, report)
cat(normalizePath(report), "\n")
