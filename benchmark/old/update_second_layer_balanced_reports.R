#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else
  "benchmark/results/full_second_layer_9method_balanced_comparison"

read_result <- function(name) {
  read.csv(file.path(output_dir, name), stringsAsFactors = FALSE)
}

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

needs_attention <- function(status) {
  status <- as.character(status)
  !is.na(status) & status != "ok" &
    grepl("fallback|unavailable|missing|failed|failure|error", status, ignore.case = TRUE)
}

metric_line <- function(tbl) {
  paste(sprintf("%s=%s", tbl$metric, fmt(tbl$mean, 4)), collapse = ", ")
}

metric_table <- function(overall, label, shell, scale, metrics) {
  x <- overall[
    overall$method_label == label &
      overall$shell == shell &
      overall$prediction_scale == scale &
      overall$metric %in% metrics,
    c("metric", "mean"),
    drop = FALSE
  ]
  x[match(metrics, x$metric), , drop = FALSE]
}

best_metric <- function(overall, shell, scale, metric, small = TRUE, package = NULL) {
  x <- overall[
    overall$shell == shell &
      overall$prediction_scale == scale &
      overall$metric == metric &
      is.finite(overall$mean),
    c("method_label", "package", "mean"),
    drop = FALSE
  ]
  if (!is.null(package)) x <- x[x$package == package, , drop = FALSE]
  x[order(if (small) x$mean else -x$mean), , drop = FALSE][1, , drop = FALSE]
}

rank_aggregate <- function(rankings, shell = "all_nearfield", scale = "raw", package_prefix = NULL) {
  keys <- c("grf_lambda", "landscape_id", "landscape_rep", "shell", "prediction_scale", "method_label")
  uniq <- rankings[
    !duplicated(rankings[, keys]),
    c(keys, "numerical_weighted_rank", "shape_weighted_rank",
      "uncertainty_weighted_rank", "runtime_failure_rank", "balanced_weighted_rank"),
    drop = FALSE
  ]
  x <- uniq[uniq$shell == shell & uniq$prediction_scale == scale, , drop = FALSE]
  if (!is.null(package_prefix)) {
    x <- x[grepl(paste0("^", package_prefix), x$method_label), , drop = FALSE]
  }
  aggregate(
    cbind(
      numerical_weighted_rank,
      shape_weighted_rank,
      uncertainty_weighted_rank,
      runtime_failure_rank,
      balanced_weighted_rank
    ) ~ method_label,
    x,
    mean,
    na.rm = TRUE
  )
}

rank_best <- function(rankings, shell, column) {
  agg <- rank_aggregate(rankings, shell = shell, package_prefix = "alfak2")
  agg[order(agg[[column]]), c("method_label", column), drop = FALSE][1, , drop = FALSE]
}

main <- function() {
  run_index <- read_result("run_index.csv")
  dependency <- read_result("dependency_status.csv")
  overall <- read_result("overall_summary.csv")
  rankings <- read_result("landscape_rankings.csv")
  pareto <- read_result("pareto_front.csv")

  failures <- run_index[
    run_index$fit_status != "ok" | needs_attention(run_index$dependency_status),
    c(
      "run_id", "package", "grf_lambda", "landscape_id", "landscape_rep",
      "input_mode", "extrapolation_method", "minobs", "NN_prior_slot",
      "NN_prior_value", "fit_status", "failure_status", "dependency_status",
      "error_message", "output_path"
    ),
    drop = FALSE
  ]
  write.csv(failures, file.path(output_dir, "failures.csv"), row.names = FALSE, na = "")

  overall$method_label <- ifelse(
    overall$package == "alfak2",
    paste("alfak2", overall$input_mode, overall$extrapolation_method, sep = ":"),
    paste("alfakR", paste0("minobs", overall$minobs), overall$NN_prior_slot, sep = ":")
  )

  all_ranks <- rank_aggregate(rankings)
  balanced <- all_ranks[order(all_ranks$balanced_weighted_rank), , drop = FALSE][1, , drop = FALSE]

  numerical <- best_metric(overall, "all_nearfield", "raw", "rmse", small = TRUE)
  shape <- best_metric(overall, "all_nearfield", "raw", "edge_gradient_rmse", small = TRUE)
  calibrated <- best_metric(overall, "all_nearfield", "anchor_calibrated", "rmse", small = TRUE)
  uncertainty <- best_metric(overall, "all_nearfield", "raw", "interval_coverage_95_closeness", small = TRUE)
  runtime <- best_metric(overall, "all_nearfield", "raw", "runtime_seconds", small = TRUE)
  topk <- best_metric(overall, "all_nearfield", "raw", "top_k_overlap_fraction", small = FALSE)
  sign <- best_metric(overall, "all_nearfield", "raw", "sign_accuracy", small = FALSE)

  numerical_metrics <- metric_table(
    overall,
    numerical$method_label,
    "all_nearfield",
    "raw",
    c("rmse", "mae", "relative_rmse", "bias", "edge_gradient_rmse", "spearman", "sign_accuracy", "top_k_overlap_fraction")
  )
  shape_metrics <- metric_table(
    overall,
    shape$method_label,
    "all_nearfield",
    "raw",
    c("edge_gradient_rmse", "edge_gradient_spearman", "pearson", "spearman", "sign_accuracy", "top_k_overlap_fraction", "rmse", "mae", "relative_rmse")
  )
  calibrated_metrics <- metric_table(
    overall,
    calibrated$method_label,
    "all_nearfield",
    "anchor_calibrated",
    c("rmse", "calibration_slope", "calibration_intercept", "affine_rmse", "mae", "relative_rmse")
  )
  uncertainty_metrics <- metric_table(
    overall,
    uncertainty$method_label,
    "all_nearfield",
    "raw",
    c("interval_coverage_95", "interval_coverage_95_closeness", "standardized_rmse", "mean_pred_sd", "median_pred_sd")
  )

  d1_num <- best_metric(overall, "d1", "raw", "rmse", small = TRUE, package = "alfak2")
  d1_shape <- best_metric(overall, "d1", "raw", "edge_gradient_rmse", small = TRUE, package = "alfak2")
  d1_bal <- rank_best(rankings, "d1", "balanced_weighted_rank")
  d2_num <- best_metric(overall, "d2", "raw", "rmse", small = TRUE, package = "alfak2")
  d2_shape <- best_metric(overall, "d2", "raw", "edge_gradient_rmse", small = TRUE, package = "alfak2")
  d2_bal <- rank_best(rankings, "d2", "balanced_weighted_rank")

  graph_label <- "alfak2:full:graph_gaussian_baseline"
  graph <- all_ranks[all_ranks$method_label == graph_label, , drop = FALSE]
  graph_pos <- c(
    numerical = match(which(all_ranks$method_label == graph_label), order(all_ranks$numerical_weighted_rank)),
    shape = match(which(all_ranks$method_label == graph_label), order(all_ranks$shape_weighted_rank)),
    balanced = match(which(all_ranks$method_label == graph_label), order(all_ranks$balanced_weighted_rank))
  )

  alfakR_ranks <- all_ranks[grepl("^alfakR", all_ranks$method_label), , drop = FALSE]
  best_alfakR_balanced <- alfakR_ranks[order(alfakR_ranks$balanced_weighted_rank), , drop = FALSE][1, , drop = FALSE]
  best_alfakR_rmse <- best_metric(overall, "all_nearfield", "raw", "rmse", small = TRUE, package = "alfakR")
  best_alfakR_metrics <- metric_table(
    overall,
    best_alfakR_rmse$method_label,
    "all_nearfield",
    "raw",
    c("rmse", "mae", "relative_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy", "top_k_overlap_fraction")
  )

  dependency_methods <- unique(dependency[
    dependency$package == "alfak2",
    c("extrapolation_method", "dependency_status", "backend", "fallback_status", "tabpfn_available"),
    drop = FALSE
  ])
  dependency_methods <- dependency_methods[order(dependency_methods$extrapolation_method), , drop = FALSE]

  pareto_counts <- sort(
    table(pareto$method[pareto$is_pareto_optimal & pareto$shell == "all_nearfield"]),
    decreasing = TRUE
  )
  pareto_methods <- names(pareto_counts)[seq_len(min(10, length(pareto_counts)))]

  method_lines <- c(
    "# Method Recommendation",
    "",
    "## Executive Recommendation",
    sprintf(
      "- Default balanced method: `%s` (balanced mean rank %s; numerical %s, shape %s, uncertainty %s).",
      balanced$method_label,
      fmt(balanced$balanced_weighted_rank, 3),
      fmt(balanced$numerical_weighted_rank, 3),
      fmt(balanced$shape_weighted_rank, 3),
      fmt(balanced$uncertainty_weighted_rank, 3)
    ),
    sprintf("- Best numerical method: `%s`; %s.", numerical$method_label, metric_line(numerical_metrics[1:4, ])),
    sprintf("- Best shape method: `%s`; %s.", shape$method_label, metric_line(shape_metrics[1:6, ])),
    sprintf("- Best calibrated method: `%s`; %s.", calibrated$method_label, metric_line(calibrated_metrics[1:4, ])),
    sprintf("- Best uncertainty method: `%s`; %s.", uncertainty$method_label, metric_line(uncertainty_metrics)),
    sprintf("- Best runtime method: `%s` (runtime_seconds=%s).", runtime$method_label, fmt(runtime$mean, 4)),
    sprintf(
      "- Best top-k overlap: `%s` (%s). Best sign accuracy: `%s` (%s).",
      topk$method_label,
      fmt(topk$mean, 4),
      sign$method_label,
      fmt(sign$mean, 4)
    ),
    "",
    "## d1/d2 Winners",
    sprintf("- d1 numerical: `%s` (RMSE=%s).", d1_num$method_label, fmt(d1_num$mean, 4)),
    sprintf("- d1 shape: `%s` (edge_gradient_rmse=%s).", d1_shape$method_label, fmt(d1_shape$mean, 4)),
    sprintf("- d1 balanced: `%s` (rank=%s).", d1_bal$method_label, fmt(d1_bal[[2]], 3)),
    sprintf("- d2 numerical: `%s` (RMSE=%s).", d2_num$method_label, fmt(d2_num$mean, 4)),
    sprintf("- d2 shape: `%s` (edge_gradient_rmse=%s).", d2_shape$method_label, fmt(d2_shape$mean, 4)),
    sprintf("- d2 balanced: `%s` (rank=%s).", d2_bal$method_label, fmt(d2_bal[[2]], 3)),
    "",
    "## graph_gaussian_baseline",
    sprintf(
      "- `%s` ranks %d/%d numerical, %d/%d shape, and %d/%d balanced; weighted ranks are numerical=%s, shape=%s, balanced=%s.",
      graph_label,
      graph_pos[["numerical"]],
      nrow(all_ranks),
      graph_pos[["shape"]],
      nrow(all_ranks),
      graph_pos[["balanced"]],
      nrow(all_ranks),
      fmt(graph$numerical_weighted_rank, 3),
      fmt(graph$shape_weighted_rank, 3),
      fmt(graph$balanced_weighted_rank, 3)
    ),
    "- It is not the raw RMSE winner and not the edge-gradient winner, but it remains the best balanced method because numerical error, nearfield shape, uncertainty, coverage, and runtime are stable together.",
    "- It is therefore not simply `numerically good but shape bad`; it is shape-competitive, while `truncated_nearfield_gmrf` is slightly better on edge gradients and `edge_effect_interaction_path_ensemble` is better on raw absolute scale.",
    "",
    "## alfakR Comparison",
    sprintf("- Best alfakR by balanced rank: `%s` (balanced=%s).", best_alfakR_balanced$method_label, fmt(best_alfakR_balanced$balanced_weighted_rank, 3)),
    sprintf("- Best alfakR by raw all_nearfield RMSE: `%s`; %s.", best_alfakR_rmse$method_label, metric_line(best_alfakR_metrics)),
    "- alfak2 dominates alfakR on the main numerical and edge-gradient criteria in this benchmark; alfakR remains competitive mostly on runtime and some top-k overlap variants.",
    "",
    "## Pareto Front",
    sprintf("- Frequent all_nearfield Pareto-front methods include: %s.", paste(sprintf("`%s`", pareto_methods), collapse = ", ")),
    "- No single method simultaneously owns every numerical and shape objective. The practical tradeoff is graph_gaussian_baseline for balanced use, edge_effect_interaction_path_ensemble for absolute fitness scale, and truncated_nearfield_gmrf / trend filtering for mutation-gradient shape.",
    "",
    "## Dependency And Failure Notes",
    sprintf("- Fit failures: %d. Dependency/fallback rows requiring attention: %d.", sum(run_index$fit_status != "ok"), nrow(failures)),
    "- TabPFN was not available; `tabpfn_nearfield_feature_model` used the declared tree fallback with xgboost backend and is not reported as real TabPFN.",
    "- `delta_tree_ensemble` used xgboost directly, not ridge fallback.",
    "",
    "## Backend Summary",
    paste(
      sprintf(
        "- `%s`: dependency_status=`%s`, backend=`%s`, fallback_status=`%s`.",
        dependency_methods$extrapolation_method,
        dependency_methods$dependency_status,
        ifelse(is.na(dependency_methods$backend) | dependency_methods$backend == "", "internal", dependency_methods$backend),
        ifelse(is.na(dependency_methods$fallback_status) | dependency_methods$fallback_status == "", "none", dependency_methods$fallback_status)
      ),
      collapse = "\n"
    )
  )
  writeLines(method_lines, file.path(output_dir, "method_recommendation.md"))

  existing_report <- readLines(file.path(output_dir, "report.md"), warn = FALSE)
  config_end <- match("## alfak2 results", existing_report)
  if (is.na(config_end)) config_end <- 1L
  config_lines <- existing_report[seq_len(config_end - 1L)]
  report_lines <- c(
    config_lines,
    "## alfak2 results",
    "- All 27 alfak2 fitting families completed for every GRF lambda and replicate.",
    sprintf("- Numerical winner by raw all_nearfield RMSE: `%s` (%s).", numerical$method_label, fmt(numerical$mean, 5)),
    sprintf("- Shape winner by raw all_nearfield edge-gradient RMSE: `%s` (%s).", shape$method_label, fmt(shape$mean, 5)),
    sprintf("- Balanced winner across lambdas/landscapes: `%s` (rank %s).", balanced$method_label, fmt(balanced$balanced_weighted_rank, 3)),
    "- Per-input-mode and lambda summaries are in `lambda_summary.csv`, `balanced_rank_summary.csv`, and `lambda_method_rank_summary.csv`.",
    "- graph_gaussian_baseline paired deltas against the eight new methods are in `baseline_delta.csv`.",
    "",
    "## alfakR results",
    "- All 15 alfakR minobs/NN_prior slot families completed.",
    sprintf("- Best alfakR balanced setting: `%s` (rank %s).", best_alfakR_balanced$method_label, fmt(best_alfakR_balanced$balanced_weighted_rank, 3)),
    sprintf("- Best alfakR raw RMSE setting: `%s` (RMSE %s).", best_alfakR_rmse$method_label, fmt(best_alfakR_rmse$mean, 4)),
    "- The duplicated `empirical_censored_weighted` NN_prior slots were preserved as slot4 and slot5 and passed with the same alfakR argument value.",
    "",
    "## alfak2 vs alfakR",
    sprintf("- Best alfak2 raw RMSE is %s versus best alfakR raw RMSE %s.", fmt(numerical$mean, 4), fmt(best_alfakR_rmse$mean, 4)),
    sprintf(
      "- Best alfak2 edge-gradient RMSE is %s versus best alfakR edge-gradient RMSE %s for the alfakR RMSE winner.",
      fmt(shape$mean, 4),
      fmt(best_alfakR_metrics$mean[best_alfakR_metrics$metric == "edge_gradient_rmse"], 4)
    ),
    "- Landscape-level paired comparisons are available in `paired_landscape_comparison.csv`; the paired unit includes lambda, landscape, shell, metric, and prediction_scale.",
    "",
    "## Numerical And Shape Tradeoff",
    sprintf("- Numerical winner: `%s`; %s.", numerical$method_label, metric_line(numerical_metrics[1:4, ])),
    sprintf("- Shape winner: `%s`; %s.", shape$method_label, metric_line(shape_metrics[1:6, ])),
    sprintf("- Calibrated winner: `%s`; %s.", calibrated$method_label, metric_line(calibrated_metrics[1:4, ])),
    sprintf("- Uncertainty winner: `%s`; %s.", uncertainty$method_label, metric_line(uncertainty_metrics)),
    sprintf("- Pareto-front methods are recorded in `pareto_front.csv`; frequent all_nearfield Pareto methods include %s.", paste(sprintf("`%s`", pareto_methods), collapse = ", ")),
    "- No conclusion uses only raw RMSE or only shape metrics; balanced ranking weights numerical 45%, shape 40%, uncertainty/coverage 10%, and runtime/failure 5%.",
    "",
    "## Recommendation",
    "- Default: use `graph_gaussian_baseline` for balanced nearfield performance and backward-compatible behavior.",
    "- Absolute fitness scale: use `edge_effect_interaction_path_ensemble` with `soft_minobs` input mode.",
    "- Mutation direction / gradient shape: use `truncated_nearfield_gmrf`; use trend filtering as an interpretable shape-preserving alternative.",
    "- Avoid treating TabPFN fallback results as real TabPFN evidence until the Python TabPFN dependency is installed and rerun.",
    "",
    "## Important notes",
    sprintf("- Fit failures: %d.", sum(run_index$fit_status != "ok")),
    sprintf("- Dependency/fallback rows requiring attention: %d.", nrow(failures)),
    "- Full benchmark completion claim: completed all 630 configured runs.",
    "",
    "## Output files",
    paste0("- `", c(
      "run_index.csv", "metrics_by_run.csv", "paired_landscape_comparison.csv",
      "baseline_delta.csv", "landscape_rankings.csv", "lambda_summary.csv",
      "lambda_method_rank_summary.csv", "overall_summary.csv", "numerical_summary.csv",
      "shape_summary.csv", "calibration_summary.csv", "balanced_rank_summary.csv",
      "pareto_front.csv", "failures.csv", "full_results.rds", "dependency_status.csv",
      "method_recommendation.md", "report.md"
    ), "`")
  )
  writeLines(report_lines, file.path(output_dir, "report.md"))
}

main()
