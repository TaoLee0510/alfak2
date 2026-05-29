#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Summarize HPC selected depth/lambda comparison metrics.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/summarize_hpc_selected_depth_lambda_comparison.R \\\n",
    "    --input-dir=benchmark/results/hpc_22chr_9method_depth_lambda_comparison \\\n",
    "    --threads=${SLURM_CPUS_PER_TASK:-1}\n",
    sep = ""
  )
}

cap_threads <- function(x, max_threads = 60L) {
  x <- suppressWarnings(as.integer(x))
  if (!is.finite(x) || x < 1L) x <- 1L
  min(as.integer(x), as.integer(max_threads))
}

default_threads <- function() {
  cap_threads(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
}

parallel_lapply <- function(x, fun, threads = default_threads(), ...) {
  threads <- cap_threads(threads)
  if (threads <= 1L || length(x) <= 1L || .Platform$OS.type == "windows") {
    return(lapply(x, fun, ...))
  }
  parallel::mclapply(
    x,
    fun,
    ...,
    mc.cores = min(threads, length(x)),
    mc.preschedule = FALSE
  )
}

parse_args <- function(args) {
  out <- list(
    input_dir = "benchmark/results/hpc_22chr_9method_depth_lambda_comparison",
    prediction_scale = "raw",
    threads = default_threads(),
    help = FALSE
  )
  for (arg in args) {
    if (arg %in% c("-h", "--help")) out$help <- TRUE
    else if (grepl("^--input-dir=", arg)) out$input_dir <- sub("^--input-dir=", "", arg)
    else if (grepl("^--prediction-scale=", arg)) out$prediction_scale <- sub("^--prediction-scale=", "", arg)
    else if (grepl("^--threads=", arg)) out$threads <- cap_threads(sub("^--threads=", "", arg))
    else stop("Unknown argument: ", arg, call. = FALSE)
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0)
}

repo_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_dir <- args$input_dir
if (!grepl("^/", input_dir)) input_dir <- file.path(repo_dir, input_dir)
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package `data.table` is required.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required.", call. = FALSE)
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package `cowplot` is required for Supplementary Fig. 2-style plots.", call. = FALSE)
}
if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  stop("Package `ggalluvial` is required for Supplementary Fig. 2-style alluvium plots.", call. = FALSE)
}

threads <- cap_threads(args$threads)
data.table::setDTthreads(threads)
message("Using ", threads, " thread(s) for HPC summary analysis.")

dt <- data.table::fread(file.path(input_dir, "metrics_by_run.csv"), na.strings = c("", "NA"), nThread = threads)
deps <- data.table::fread(file.path(input_dir, "dependency_status.csv"), na.strings = c("", "NA"), nThread = threads)

tables_dir <- file.path(input_dir, "tables")
figures_dir <- file.path(input_dir, "figures")
results_dir <- file.path(input_dir, "results")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", na = "NA")
}

rename_present <- function(x, old, new) {
  keep <- old %in% names(x)
  if (any(keep)) data.table::setnames(x, old[keep], new[keep])
  x
}

safe_slug <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

package_display <- function(package) {
  ifelse(package == "alfak2", "alfak_V2", ifelse(package == "alfakR", "alfak", as.character(package)))
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.finite(x), paste0(formatC(100 * x, format = "f", digits = digits), "%"), "NA")
}

as_markdown_table <- function(x) {
  x <- data.table::as.data.table(x)
  if (!nrow(x)) return("_No rows._")
  x <- x[, lapply(.SD, as.character)]
  x <- x[, lapply(.SD, function(v) gsub("\\|", "\\\\|", v))]
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

method_label <- function(package, input_mode, extrapolation_method, minobs, nn_prior_slot) {
  pkg <- package_display(package)
  ifelse(
    package == "alfak2",
    paste(pkg, input_mode, extrapolation_method, sep = ":"),
    paste0(pkg, ":minobs", minobs, ":", nn_prior_slot)
  )
}

short_method <- function(x) {
  x <- sub("^alfak_V2:", "alfak_V2:", x)
  x <- sub("^alfak:", "alfak:", x)
  x <- gsub("edge_effect_interaction_path_ensemble", "edge_path", x)
  x <- gsub("graph_gaussian_baseline", "graph_gaussian", x)
  x <- gsub("truncated_nearfield_gmrf", "trunc_gmrf", x)
  x <- gsub("tabpfn_nearfield_feature_model", "tabpfn_fallback", x)
  x <- gsub("delta_tree_ensemble", "delta_tree", x)
  x <- gsub("empirical_censored_weighted_slot", "ecw", x)
  x <- gsub("empirical_censored", "ec", x)
  x
}

slot_suffix_number <- function(x) {
  x <- as.character(x)
  has_suffix <- grepl("_slot[0-9]+$", x)
  out <- rep(NA_integer_, length(x))
  out[has_suffix] <- suppressWarnings(as.integer(sub("^.*_slot([0-9]+)$", "\\1", x[has_suffix])))
  out
}

label_slot_short <- function(x) {
  x <- as.character(x)
  x <- gsub("empirical_censored_weighted_slot", "ecw", x)
  x <- gsub("empirical_censored", "ec", x)
  x
}

metric_meta <- data.table::data.table(
  metric = c(
    "rmse", "mae", "relative_rmse", "median_absolute_error", "q90_absolute_error",
    "centered_rmse", "affine_rmse", "rescaled_r2", "pearson", "spearman",
    "edge_gradient_rmse", "edge_gradient_spearman", "sign_accuracy",
    "beneficial_sign_accuracy", "deleterious_sign_accuracy", "top_k_overlap_fraction"
  ),
  domain = c(
    rep("numerical", 5),
    rep("shape", 11)
  ),
  direction = c(
    rep("lower", 5),
    "lower", "lower", "higher", "higher", "higher",
    "lower", "higher", "higher", "higher", "higher", "higher"
  )
)

selected_metrics <- metric_meta$metric
key_cols <- c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat", "shell", "metric")

dt[, method_label := method_label(package, input_mode, extrapolation_method, minobs, NN_prior_slot)]
dt[, method_short := short_method(method_label)]
dt[, package_display := package_display(package)]
dt[, ruggedness_index := 1 / grf_lambda]
dt[, display_minobs := as.integer(minobs)]
if ("anchor_count_reference" %in% names(dt)) {
  dt[
    package == "alfak2" & input_mode == "soft_minobs" & is.finite(suppressWarnings(as.numeric(anchor_count_reference))),
    display_minobs := as.integer(as.numeric(anchor_count_reference))
  ]
}
dt[package == "alfak2" & input_mode == "soft_minobs" & !is.finite(display_minobs), display_minobs := 10L]
dt[
  package == "alfak2" & input_mode == "soft_minobs" & is.finite(display_minobs),
  method_label := paste(package_display, paste0(input_mode, display_minobs), extrapolation_method, sep = ":")
]
dt[, method_short := short_method(method_label)]
dt[, NN_prior_slot_number := slot_suffix_number(NN_prior_slot)]
dt[, display_slot_keep := TRUE]
weighted_slot_keep <- dt[
  package == "alfakR" &
    NN_prior == "empirical_censored_weighted" &
    is.finite(NN_prior_slot_number),
  .(display_keep_slot_number = min(NN_prior_slot_number, na.rm = TRUE)),
  by = .(NN_prior)
]
if (nrow(weighted_slot_keep)) {
  dt <- merge(dt, weighted_slot_keep, by = "NN_prior", all.x = TRUE)
  dt[
    package == "alfakR" &
      NN_prior == "empirical_censored_weighted" &
      is.finite(NN_prior_slot_number) &
      NN_prior_slot_number != display_keep_slot_number,
    display_slot_keep := FALSE
  ]
} else {
  dt[, display_keep_slot_number := NA_integer_]
}
analysis_source_dt <- dt[display_slot_keep == TRUE]
dt <- merge(dt, metric_meta, by = "metric", all.x = TRUE)
analysis_source_dt <- merge(analysis_source_dt, metric_meta, by = "metric", all.x = TRUE)
analysis_dt <- analysis_source_dt[
  prediction_scale == args$prediction_scale &
    metric %in% selected_metrics &
    fit_status == "ok" &
    failure_status %in% c("ok", "fallback_used") &
    is.finite(value)
]

run_status_summary <- deps[, .N, by = .(package, dependency_status, fit_status, failure_status)]
write_tsv(run_status_summary, file.path(tables_dir, "hpc_run_status_summary.tsv"))

summarise_values <- function(x, by) {
  x[, .(
    n = .N,
    mean = mean(value),
    sd = stats::sd(value),
    se = stats::sd(value) / sqrt(.N),
    median = stats::median(value),
    q25 = stats::quantile(value, 0.25, names = FALSE),
    q75 = stats::quantile(value, 0.75, names = FALSE),
    ci95_low = mean(value) - 1.96 * stats::sd(value) / sqrt(.N),
    ci95_high = mean(value) + 1.96 * stats::sd(value) / sqrt(.N)
  ), by = by]
}

method_summary <- summarise_values(
  analysis_dt,
  c("sample_depth", "grf_lambda", "shell", "package", "input_mode",
    "extrapolation_method", "minobs", "NN_prior_slot", "method_label",
    "method_short", "metric", "domain", "direction")
)
data.table::setorder(method_summary, sample_depth, grf_lambda, shell, metric, package, method_label)
write_tsv(method_summary, file.path(tables_dir, "hpc_method_metric_summary.tsv"))

landscape_method_summary <- analysis_dt[, .(
  n_fit_repeats = .N,
  fit_repeat_mean = mean(value),
  fit_repeat_median = stats::median(value),
  fit_repeat_sd = stats::sd(value),
  fit_repeat_q25 = stats::quantile(value, 0.25, names = FALSE),
  fit_repeat_q75 = stats::quantile(value, 0.75, names = FALSE)
), by = .(
  sample_depth, grf_lambda, landscape_id, landscape_rep, shell, package,
  input_mode, extrapolation_method, minobs, NN_prior_slot, method_label,
  method_short, metric, domain, direction
)]
data.table::setorder(
  landscape_method_summary,
  sample_depth, grf_lambda, landscape_rep, shell, metric, package, method_label
)
write_tsv(
  landscape_method_summary,
  file.path(tables_dir, "hpc_landscape_fit_repeat_summary.tsv")
)

landscape_block_method_summary <- landscape_method_summary[, .(
  n_landscapes = .N,
  total_fit_repeats = sum(n_fit_repeats),
  mean_of_fit_repeat_means = mean(fit_repeat_mean),
  median_of_fit_repeat_means = stats::median(fit_repeat_mean),
  sd_of_fit_repeat_means = stats::sd(fit_repeat_mean),
  se_of_fit_repeat_means = stats::sd(fit_repeat_mean) / sqrt(.N),
  ci95_low_of_fit_repeat_means = mean(fit_repeat_mean) - 1.96 * stats::sd(fit_repeat_mean) / sqrt(.N),
  ci95_high_of_fit_repeat_means = mean(fit_repeat_mean) + 1.96 * stats::sd(fit_repeat_mean) / sqrt(.N),
  mean_of_fit_repeat_medians = mean(fit_repeat_median),
  median_of_fit_repeat_medians = stats::median(fit_repeat_median),
  sd_of_fit_repeat_medians = stats::sd(fit_repeat_median),
  mean_of_fit_repeat_sds = mean(fit_repeat_sd, na.rm = TRUE),
  median_of_fit_repeat_sds = stats::median(fit_repeat_sd, na.rm = TRUE),
  sd_of_fit_repeat_sds = stats::sd(fit_repeat_sd, na.rm = TRUE)
), by = .(
  sample_depth, grf_lambda, shell, package, input_mode, extrapolation_method,
  minobs, NN_prior_slot, method_label, method_short, metric, domain, direction
)]
data.table::setorder(
  landscape_block_method_summary,
  sample_depth, grf_lambda, shell, metric, package, method_label
)
write_tsv(
  landscape_block_method_summary,
  file.path(tables_dir, "hpc_landscape_block_method_metric_summary.tsv")
)

overall_summary <- summarise_values(
  analysis_dt,
  c("shell", "package", "input_mode", "extrapolation_method", "minobs",
    "NN_prior_slot", "method_label", "method_short", "metric", "domain", "direction")
)
data.table::setorder(overall_summary, shell, metric, package, method_label)
write_tsv(overall_summary, file.path(tables_dir, "hpc_overall_method_metric_summary.tsv"))

best_by_package <- method_summary[
  , .SD[if (direction[1] == "lower") which.min(mean) else which.max(mean)],
  by = .(sample_depth, grf_lambda, shell, metric, package)
]
data.table::setorder(best_by_package, sample_depth, grf_lambda, shell, metric, package)
write_tsv(best_by_package, file.path(tables_dir, "hpc_best_method_by_package_metric.tsv"))

overall_best_by_package <- overall_summary[
  , .SD[if (direction[1] == "lower") which.min(mean) else which.max(mean)],
  by = .(shell, metric, package)
]
data.table::setorder(overall_best_by_package, shell, metric, package)
write_tsv(overall_best_by_package, file.path(tables_dir, "hpc_overall_best_method_by_package_metric.tsv"))

alfakR_best_condition <- analysis_dt[package == "alfakR", {
  if (direction[1] == "lower") .SD[which.min(value)] else .SD[which.max(value)]
}, by = key_cols]
alfakR_best_condition <- alfakR_best_condition[
  , c(key_cols, "method_label", "method_short", "value"), with = FALSE
]
data.table::setnames(
  alfakR_best_condition,
  c("method_label", "method_short", "value"),
  c("best_alfakR_method", "best_alfakR_short", "best_alfakR_value")
)

alfak2_condition <- analysis_dt[package == "alfak2"]
delta_best <- merge(alfak2_condition, alfakR_best_condition, by = key_cols, allow.cartesian = FALSE)
delta_best[, raw_delta := value - best_alfakR_value]
delta_best[, improvement := ifelse(direction == "lower", best_alfakR_value - value, value - best_alfakR_value)]
delta_best[, relative_improvement := improvement / pmax(abs(best_alfakR_value), .Machine$double.eps)]
delta_best[, improved := improvement > 0]
write_tsv(
  delta_best[
    , .(sample_depth, grf_lambda, landscape_id, landscape_rep, fit_repeat, shell, metric,
        domain, direction, alfak_V2_method = method_label, best_alfak_method = best_alfakR_method,
        alfak_V2_value = value, best_alfak_value = best_alfakR_value, raw_delta, improvement,
        relative_improvement, improved)
  ],
  file.path(tables_dir, "hpc_alfak_V2_delta_vs_best_alfak_by_condition.tsv")
)

delta_best_summary <- delta_best[, .(
  n = .N,
  alfak2_mean = mean(value),
  best_alfakR_mean = mean(best_alfakR_value),
  mean_raw_delta = mean(raw_delta),
  median_raw_delta = stats::median(raw_delta),
  mean_improvement = mean(improvement),
  median_improvement = stats::median(improvement),
  mean_relative_improvement = mean(relative_improvement),
  improvement_rate = mean(improved)
), by = .(sample_depth, grf_lambda, shell, metric, domain, direction,
          alfak2_method = method_label, alfak2_short = method_short)]
data.table::setorder(delta_best_summary, sample_depth, grf_lambda, shell, metric, -mean_improvement)
delta_best_summary_out <- data.table::copy(delta_best_summary)
delta_best_summary_out <- rename_present(
  delta_best_summary_out,
  c("alfak2_mean", "best_alfakR_mean", "alfak2_method", "alfak2_short"),
  c("alfak_V2_mean", "best_alfak_mean", "alfak_V2_method", "alfak_V2_short")
)
write_tsv(delta_best_summary_out, file.path(tables_dir, "hpc_alfak_V2_delta_vs_best_alfak_summary.tsv"))

metric_values <- sort(unique(analysis_dt$metric))
message("Building pairwise alfak_V2-vs-alfak delta summaries across ", length(metric_values), " metric(s) with up to ", min(threads, length(metric_values)), " worker(s).")
build_pairwise_delta_part <- function(mm) {
  old_dt_threads <- data.table::getDTthreads()
  data.table::setDTthreads(1L)
  on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
  a2 <- analysis_dt[package == "alfak2" & metric == mm]
  ar <- analysis_dt[package == "alfakR" & metric == mm]
  if (!nrow(a2) || !nrow(ar)) return(NULL)
  pair <- merge(
    a2,
    ar[, c(key_cols, "method_label", "method_short", "value"), with = FALSE],
    by = key_cols,
    allow.cartesian = TRUE,
    suffixes = c("_alfak2", "_alfakR")
  )
  pair[, raw_delta := value_alfak2 - value_alfakR]
  pair[, improvement := ifelse(direction == "lower", value_alfakR - value_alfak2, value_alfak2 - value_alfakR)]
  pair[, relative_improvement := improvement / pmax(abs(value_alfakR), .Machine$double.eps)]
  pair[, improved := improvement > 0]
  pair[, .(
    n = .N,
    alfak2_mean = mean(value_alfak2),
    alfakR_mean = mean(value_alfakR),
    mean_raw_delta = mean(raw_delta),
    median_raw_delta = stats::median(raw_delta),
    mean_improvement = mean(improvement),
    median_improvement = stats::median(improvement),
    mean_relative_improvement = mean(relative_improvement),
    improvement_rate = mean(improved)
  ), by = .(sample_depth, grf_lambda, shell, metric, domain, direction,
            alfak2_method = method_label_alfak2,
            alfak2_short = method_short_alfak2,
            alfakR_method = method_label_alfakR,
            alfakR_short = method_short_alfakR)]
}
pairwise_delta_parts <- parallel_lapply(
  metric_values,
  build_pairwise_delta_part,
  threads = min(threads, length(metric_values))
)
pairwise_delta_summary <- data.table::rbindlist(pairwise_delta_parts, use.names = TRUE, fill = TRUE)
data.table::setorder(pairwise_delta_summary, sample_depth, grf_lambda, shell, metric, -mean_improvement)
pairwise_delta_summary_out <- data.table::copy(pairwise_delta_summary)
pairwise_delta_summary_out <- rename_present(
  pairwise_delta_summary_out,
  c("alfak2_mean", "alfakR_mean", "alfak2_method", "alfak2_short", "alfakR_method", "alfakR_short"),
  c("alfak_V2_mean", "alfak_mean", "alfak_V2_method", "alfak_V2_short", "alfak_method", "alfak_short")
)
write_tsv(pairwise_delta_summary_out, file.path(tables_dir, "hpc_pairwise_alfak_V2_vs_alfak_delta_summary.tsv"))

rank_metrics <- c("rmse", "mae", "relative_rmse", "centered_rmse", "rescaled_r2", "pearson", "spearman",
                  "edge_gradient_rmse", "edge_gradient_spearman", "sign_accuracy")
rank_dt <- analysis_dt[metric %in% rank_metrics]
rank_dt[, metric_rank := ifelse(
  direction == "lower",
  data.table::frank(value, ties.method = "average"),
  data.table::frank(-value, ties.method = "average")
), by = c(key_cols)]
rank_summary <- rank_dt[, .(
  metric_mean_rank = mean(metric_rank),
  metric_median_rank = stats::median(metric_rank),
  metric_top1_rate = mean(metric_rank == 1),
  metric_top3_rate = mean(metric_rank <= 3)
), by = .(sample_depth, grf_lambda, shell, package, method_label, method_short, metric, domain)]
domain_rank_summary <- rank_summary[, .(
  domain_mean_rank = mean(metric_mean_rank),
  domain_median_rank = stats::median(metric_mean_rank),
  top3_rate = mean(metric_top3_rate)
), by = .(sample_depth, grf_lambda, shell, package, method_label, method_short, domain)]
rank_wide <- data.table::dcast(
  domain_rank_summary,
  sample_depth + grf_lambda + shell + package + method_label + method_short ~ domain,
  value.var = "domain_mean_rank"
)
if (!"numerical" %in% names(rank_wide)) rank_wide[, numerical := NA_real_]
if (!"shape" %in% names(rank_wide)) rank_wide[, shape := NA_real_]
rank_wide[, balanced_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = c("numerical", "shape")]
data.table::setorder(rank_wide, sample_depth, grf_lambda, shell, balanced_rank)
write_tsv(rank_wide, file.path(tables_dir, "hpc_balanced_rank_summary.tsv"))

overall_rank_dt <- analysis_dt[metric %in% rank_metrics]
overall_rank_dt[, metric_rank := ifelse(
  direction == "lower",
  data.table::frank(value, ties.method = "average"),
  data.table::frank(-value, ties.method = "average")
), by = .(sample_depth, grf_lambda, landscape_id, landscape_rep, fit_repeat, shell, metric)]
overall_domain_rank <- overall_rank_dt[, .(
  domain_mean_rank = mean(metric_rank)
), by = .(shell, package, method_label, method_short, domain)]
overall_rank_wide <- data.table::dcast(
  overall_domain_rank,
  shell + package + method_label + method_short ~ domain,
  value.var = "domain_mean_rank"
)
if (!"numerical" %in% names(overall_rank_wide)) overall_rank_wide[, numerical := NA_real_]
if (!"shape" %in% names(overall_rank_wide)) overall_rank_wide[, shape := NA_real_]
overall_rank_wide[, balanced_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = c("numerical", "shape")]
data.table::setorder(overall_rank_wide, shell, balanced_rank)
write_tsv(overall_rank_wide, file.path(tables_dir, "hpc_overall_balanced_rank_summary.tsv"))

trade_metrics <- c("rmse", "mae", "relative_rmse", "centered_rmse", "rescaled_r2", "edge_gradient_rmse",
                   "spearman", "edge_gradient_spearman", "sign_accuracy")
trade_all <- overall_summary[shell == "all_nearfield" & package == "alfak2" & metric %in% trade_metrics]
trade_wide <- data.table::dcast(
  trade_all,
  method_label + method_short ~ metric,
  value.var = "mean"
)
trade_rank <- overall_rank_wide[shell == "all_nearfield" & package == "alfak2",
                                .(method_label, numerical_rank = numerical,
                                  shape_rank = shape, balanced_rank)]
alfak2_tradeoff <- merge(trade_wide, trade_rank, by = "method_label", all.x = TRUE)
gap_lower <- function(x) x - min(x, na.rm = TRUE)
gap_higher <- function(x) max(x, na.rm = TRUE) - x
for (col in intersect(c("rmse", "mae", "relative_rmse", "centered_rmse", "edge_gradient_rmse"), names(alfak2_tradeoff))) {
  alfak2_tradeoff[, paste0(col, "_gap_to_best") := gap_lower(get(col))]
}
for (col in intersect(c("rescaled_r2", "spearman", "edge_gradient_spearman", "sign_accuracy"), names(alfak2_tradeoff))) {
  alfak2_tradeoff[, paste0(col, "_gap_to_best") := gap_higher(get(col))]
}
data.table::setorder(alfak2_tradeoff, balanced_rank)
write_tsv(alfak2_tradeoff, file.path(tables_dir, "hpc_alfak_V2_all_nearfield_tradeoff.tsv"))

lambda_trend <- method_summary[, .(
  n = sum(n),
  mean = stats::weighted.mean(mean, n),
  median = stats::weighted.mean(median, n)
), by = .(sample_depth, grf_lambda, ruggedness_index = 1 / grf_lambda, shell, package,
          method_label, method_short, metric, domain, direction)]
write_tsv(lambda_trend, file.path(tables_dir, "hpc_lambda_trend_summary.tsv"))

rugged_effect <- analysis_dt[, {
  ok <- is.finite(value) & is.finite(ruggedness_index)
  if (sum(ok) >= 5L && length(unique(ruggedness_index[ok])) >= 2L) {
    fit <- stats::lm(value[ok] ~ ruggedness_index[ok])
    slope <- unname(stats::coef(fit)[[2]])
    corr <- suppressWarnings(stats::cor(ruggedness_index[ok], value[ok]))
  } else {
    slope <- NA_real_
    corr <- NA_real_
  }
  .(
    n = sum(ok),
    mean_value = mean(value[ok]),
    slope_vs_ruggedness = slope,
    cor_vs_ruggedness = corr,
    worse_when_more_rugged = ifelse(direction[1] == "lower", slope > 0, slope < 0)
  )
}, by = .(sample_depth, shell, package, method_label, method_short, metric, domain, direction)]
data.table::setorder(rugged_effect, sample_depth, shell, metric, package, method_label)
write_tsv(rugged_effect, file.path(tables_dir, "hpc_ruggedness_effect_summary.tsv"))

alfak_style_metric_map <- data.table::data.table(
  shell = rep(c("full_lscape", "d0", "d1", "d2"), each = 3L),
  metric = rep(c("pearson", "spearman", "rescaled_r2"), times = 4L),
  alfak_metric = c(
    "r", "rho", "R",
    "rfq", "rhofq", "Rfq",
    "rnn", "rhonn", "Rnn",
    "rd2", "rhod2", "Rd2"
  ),
  alfak_subset = rep(c("full", "fq", "nn", "d2"), each = 3L),
  alfak_family = rep(c("Pearson r", "Spearman rho", "Rescaled R^2"), times = 4L)
)

build_alfak_style_ground_truth_summary <- function(source_dt) {
  source_dt <- source_dt[
    prediction_scale == args$prediction_scale &
      fit_status == "ok" &
      failure_status %in% c("ok", "fallback_used")
  ]
  id_cols <- intersect(
    c(
      "task_id", "run_id", "sample_depth", "grf_lambda", "landscape_id",
      "landscape_rep", "fit_repeat", "package", "package_display", "input_mode",
      "extrapolation_method", "minobs", "display_minobs", "NN_prior_slot", "method_label",
      "method_short"
    ),
    names(source_dt)
  )
  metric_long <- merge(
    source_dt[is.finite(value)],
    alfak_style_metric_map,
    by = c("shell", "metric"),
    allow.cartesian = FALSE
  )
  if (!nrow(metric_long)) {
    return(data.table::data.table())
  }
  style <- data.table::dcast(
    metric_long,
    stats::as.formula(paste(paste(id_cols, collapse = " + "), "~ alfak_metric")),
    value.var = "value"
  )

  count_specs <- data.table::data.table(
    shell = c("full_lscape", "d0", "d1", "d2"),
    n_col = c("n_full_lscape", "nfq", "nnn", "nd2"),
    n_valid_col = c("n_valid_full_lscape", "n_valid_fq", "n_valid_nn", "n_valid_d2")
  )
  coverage <- source_dt[metric == "coverage" & shell %in% count_specs$shell]
  for (i in seq_len(nrow(count_specs))) {
    spec <- count_specs[i]
    cc <- coverage[shell == spec$shell[[1]]]
    if (!nrow(cc)) next
    cc <- cc[, c(id_cols, "n_eval_nodes", "n_valid_predictions"), with = FALSE]
    cc <- unique(cc)
    data.table::setnames(cc, c("n_eval_nodes", "n_valid_predictions"), c(spec$n_col[[1]], spec$n_valid_col[[1]]))
    style <- merge(style, cc, by = id_cols, all.x = TRUE)
  }
  style[, Rxv := NA_real_]

  preferred <- c(
    id_cols,
    "r", "rfq", "rnn", "rd2",
    "rho", "rhofq", "rhonn", "rhod2",
    "R", "Rfq", "Rnn", "Rd2",
    "Rxv", "n_full_lscape", "nfq", "nnn", "nd2",
    "n_valid_full_lscape", "n_valid_fq", "n_valid_nn", "n_valid_d2"
  )
  data.table::setcolorder(style, intersect(preferred, names(style)))
  data.table::setorder(style, package_display, input_mode, minobs, extrapolation_method, NN_prior_slot, sample_depth, grf_lambda, landscape_rep, fit_repeat)
  style
}

alfak_style_summary <- build_alfak_style_ground_truth_summary(analysis_source_dt)
write_tsv(alfak_style_summary, file.path(tables_dir, "hpc_alfak_style_ground_truth_summary.tsv"))

saveRDS(
  list(
    run_status_summary = run_status_summary,
    method_summary = method_summary,
    landscape_method_summary = landscape_method_summary,
    landscape_block_method_summary = landscape_block_method_summary,
    overall_summary = overall_summary,
    best_by_package = best_by_package,
    overall_best_by_package = overall_best_by_package,
    delta_best_summary = delta_best_summary,
    pairwise_delta_summary = pairwise_delta_summary,
    rank_summary = rank_wide,
    overall_rank_summary = overall_rank_wide,
    alfak2_tradeoff = alfak2_tradeoff,
    alfak_style_summary = alfak_style_summary,
    lambda_trend = lambda_trend,
    ruggedness_effect = rugged_effect
  ),
  file.path(results_dir, "hpc_selected_depth_lambda_analysis.rds")
)

ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))

plot_save <- function(plot, filename, width = 12, height = 7) {
  ggplot2::ggsave(file.path(figures_dir, filename), plot = plot, width = width, height = height, dpi = 180)
}

metric_label <- function(metric) {
  labels <- c(
    rmse = "RMSE",
    mae = "MAE",
    relative_rmse = "Relative RMSE",
    median_absolute_error = "Median absolute error",
    q90_absolute_error = "Q90 absolute error",
    centered_rmse = "Centered RMSE",
    affine_rmse = "Affine RMSE",
    rescaled_r2 = "Rescaled R^2",
    pearson = "Pearson r",
    spearman = "Spearman rho",
    edge_gradient_rmse = "Edge-gradient RMSE",
    edge_gradient_spearman = "Edge-gradient Spearman rho",
    sign_accuracy = "Sign accuracy",
    beneficial_sign_accuracy = "Beneficial sign accuracy",
    deleterious_sign_accuracy = "Deleterious sign accuracy",
    top_k_overlap_fraction = "Top-k overlap fraction"
  )
  out <- labels[metric]
  out[is.na(out)] <- metric[is.na(out)]
  unname(out)
}

alfak_metric_label <- function(metric) {
  labels <- c(
    r = "Pearson r",
    rfq = "Pearson r (fq/d0)",
    rnn = "Pearson r (nn/d1)",
    rd2 = "Pearson r (d2)",
    rho = "Spearman rho",
    rhofq = "Spearman rho (fq/d0)",
    rhonn = "Spearman rho (nn/d1)",
    rhod2 = "Spearman rho (d2)",
    R = "Rescaled R^2",
    Rfq = "Rescaled R^2 (fq/d0)",
    Rnn = "Rescaled R^2 (nn/d1)",
    Rd2 = "Rescaled R^2 (d2)",
    Rxv = "Cross-validation R^2"
  )
  out <- labels[metric]
  out[is.na(out)] <- metric[is.na(out)]
  unname(out)
}

quantile_summary <- function(x, by) {
  x[, {
    qq <- stats::quantile(value, probs = c(0.1, 0.5, 0.9), names = FALSE, na.rm = TRUE)
    .(lo = qq[[1]], med = qq[[2]], hi = qq[[3]], n = .N)
  }, by = by]
}

shell_order <- c("full_lscape", "all_eval", "all_nearfield", "d0", "d1", "d2")

method_artifact_path <- function(row) {
  pkg <- as.character(row$package[[1]])
  if (identical(pkg, "alfak2")) {
    file.path(
      input_dir,
      "method_reports",
      "alfak_V2",
      safe_slug(row$input_mode[[1]]),
      safe_slug(row$extrapolation_method[[1]])
    )
  } else {
    file.path(
      input_dir,
      "method_reports",
      "alfak",
      paste0("minobs", safe_slug(row$minobs[[1]])),
      safe_slug(row$NN_prior_slot[[1]])
    )
  }
}

karyotype_matrix <- function(karyotypes) {
  parts <- strsplit(as.character(karyotypes), ".", fixed = TRUE)
  n_chr <- unique(lengths(parts))
  if (length(n_chr) != 1L) {
    stop("Karyotypes have inconsistent chromosome counts.", call. = FALSE)
  }
  mat <- do.call(rbind, lapply(parts, as.numeric))
  storage.mode(mat) <- "double"
  mat
}

predict_cached_grf_truth <- function(landscape, karyotypes) {
  karyotypes <- as.character(karyotypes)
  out <- rep(NA_real_, length(karyotypes))
  if (!is.null(landscape$cache) && is.environment(landscape$cache)) {
    cached <- mget(karyotypes, envir = landscape$cache, ifnotfound = list(NA_real_))
    out <- as.numeric(unlist(cached, use.names = FALSE))
  }
  miss <- which(!is.finite(out))
  if (!length(miss)) return(out)

  if (is.null(landscape$centroids) || is.null(landscape$lambda) ||
      is.null(landscape$scale) || is.null(landscape$founder_raw) ||
      is.null(landscape$founder_fitness)) {
    stop("Landscape object does not contain the GRF fields needed for truth prediction.", call. = FALSE)
  }
  kmat <- karyotype_matrix(karyotypes[miss])
  centroids <- as.matrix(landscape$centroids)
  if (ncol(kmat) != ncol(centroids)) {
    stop("Karyotype chromosome count does not match landscape centroid dimension.", call. = FALSE)
  }
  m <- nrow(centroids)
  raw <- apply(kmat, 1, function(k) {
    delta <- sweep(centroids, 2, k, "-")
    dist <- sqrt(rowSums(delta^2))
    sum(sin(dist / as.numeric(landscape$lambda))) / (pi * sqrt(m))
  })
  out[miss] <- as.numeric(landscape$founder_fitness) +
    as.numeric(landscape$scale) * (raw - as.numeric(landscape$founder_raw))
  out
}

landscape_cache_candidates <- function(row) {
  row <- as.list(row)
  sample_depth <- as.integer(row$sample_depth)
  fit_repeat <- as.integer(row$fit_repeat)
  count_seed <- as.integer(row$count_seed)
  landscape_id <- as.character(row$landscape_id)
  n_chr <- if ("n_chr" %in% names(row) && is.finite(suppressWarnings(as.numeric(row$n_chr)))) {
    as.integer(row$n_chr)
  } else {
    NA_integer_
  }
  labels <- character()
  if (is.finite(n_chr)) {
    labels <- c(labels, paste(
      landscape_id,
      paste0("chr", n_chr),
      paste0("depth", sample_depth),
      paste0("fitrep", fit_repeat),
      paste0("countseed", count_seed),
      sep = "_"
    ))
  }
  labels <- c(
    labels,
    paste(
      landscape_id,
      "chr4",
      paste0("depth", sample_depth),
      paste0("fitrep", fit_repeat),
      paste0("countseed", count_seed),
      sep = "_"
    ),
    paste(
      landscape_id,
      paste0("depth", sample_depth),
      paste0("fitrep", fit_repeat),
      paste0("countseed", count_seed),
      sep = "_"
    )
  )
  unique(file.path(input_dir, "landscape_cache", paste0(labels, ".rds")))
}

read_representative_scatter <- function(rows, max_points_per_panel = 400L) {
  if (!nrow(rows)) return(data.table::data.table())
  pieces <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    rr <- rows[i]
    cache_path <- file.path(input_dir, "run_cache", paste0(rr$run_id[[1]], ".rds"))
    if (!file.exists(cache_path)) next
    run_cache <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.null(run_cache) || is.null(run_cache$result$predictions) || !nrow(run_cache$result$predictions)) next
    landscape_path <- landscape_cache_candidates(run_cache$row[1, , drop = FALSE])
    landscape_path <- landscape_path[file.exists(landscape_path)][1]
    if (is.na(landscape_path) || !nzchar(landscape_path)) next
    shared <- tryCatch(readRDS(landscape_path), error = function(e) NULL)
    if (is.null(shared$landscape)) next

    pred <- data.table::as.data.table(run_cache$result$predictions)
    pred <- pred[is.finite(fitness_mean) & !is.na(karyotype) & nzchar(as.character(karyotype))]
    if (!nrow(pred)) next
    pred <- pred[!duplicated(karyotype)]
    pred[, truth := predict_cached_grf_truth(shared$landscape, karyotype)]
    pred <- pred[is.finite(truth)]
    if (!nrow(pred)) next
    if (nrow(pred) > max_points_per_panel) {
      set.seed(as.integer(rr$task_id[[1]]) %% .Machine$integer.max)
      pred <- pred[sample.int(nrow(pred), max_points_per_panel)]
    }
    pred[, `:=`(
      pred_centered = fitness_mean - mean(fitness_mean, na.rm = TRUE),
      truth_centered = truth - mean(truth, na.rm = TRUE),
      sample_depth = rr$sample_depth[[1]],
      grf_lambda = rr$grf_lambda[[1]],
      minobs = rr$minobs[[1]],
      display_minobs = if ("display_minobs" %in% names(rr)) rr$display_minobs[[1]] else rr$minobs[[1]],
      package = rr$package[[1]],
      method_label = rr$method_label[[1]]
    )]
    pieces[[i]] <- pred[, .(
      sample_depth, grf_lambda, minobs, display_minobs, package, method_label,
      karyotype, truth_centered, pred_centered
    )]
  }
  data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
}

select_lambda_values <- function(available, requested = NULL, extremes = FALSE) {
  available <- sort(unique(available[is.finite(available)]))
  if (!length(available)) return(numeric())
  if (isTRUE(extremes)) return(unique(c(available[[1]], available[[length(available)]])))
  if (is.null(requested)) return(available[[ceiling(length(available) / 2)]])
  unique(vapply(requested, function(x) available[which.min(abs(available - x))], numeric(1)))
}

pick_representative_rows <- function(method_rows, by_cols, lambda_values = NULL, lambda_extremes = FALSE) {
  if (!nrow(method_rows)) return(method_rows)
  target_lambdas <- select_lambda_values(method_rows$grf_lambda, requested = lambda_values, extremes = lambda_extremes)
  if (length(target_lambdas)) method_rows <- method_rows[grf_lambda %in% target_lambdas]
  by_cols <- unique(c("grf_lambda", by_cols))
  data.table::setorderv(
    method_rows,
    intersect(c("grf_lambda", "sample_depth", "display_minobs", "minobs", "landscape_rep", "fit_repeat", "task_id"), names(method_rows))
  )
  method_rows[, .SD[1], by = by_cols]
}

suppfig2_landscape_level <- function(method_rows) {
  if (!nrow(method_rows)) return(method_rows)
  metric_cols <- intersect(
    c("r", "rfq", "rnn", "rd2", "rho", "rhofq", "rhonn", "rhod2", "R", "Rfq", "Rnn", "Rd2", "Rxv"),
    names(method_rows)
  )
  count_cols <- intersect(
    c("n_full_lscape", "nfq", "nnn", "nd2", "n_valid_full_lscape", "n_valid_fq", "n_valid_nn", "n_valid_d2"),
    names(method_rows)
  )
  by_cols <- intersect(
    c(
      "package", "package_display", "input_mode", "extrapolation_method",
      "minobs", "display_minobs", "NN_prior_slot", "method_label", "method_short",
      "sample_depth", "grf_lambda", "landscape_id", "landscape_rep"
    ),
    names(method_rows)
  )
  method_rows[, c(
    lapply(.SD[, metric_cols, with = FALSE], function(v) mean(v, na.rm = TRUE)),
    lapply(.SD[, count_cols, with = FALSE], function(v) round(mean(v, na.rm = TRUE)))
  ), by = by_cols]
}

plot_suppfig2_panel_a <- function(method_rows, use_minobs_axis) {
  metrics <- intersect(c("r", "R", "rho"), names(method_rows))
  if (!nrow(method_rows) || !length(metrics)) return(ggplot2::ggplot() + ggplot2::theme_void())
  z <- data.table::melt(
    method_rows,
    id.vars = intersect(c("sample_depth", "grf_lambda", "minobs", "display_minobs"), names(method_rows)),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
  )
  z <- z[is.finite(value)]
  if (!nrow(z)) return(ggplot2::ggplot() + ggplot2::theme_void())
  z[, value := pmax(-1, value)]
  z[, metric := factor(alfak_metric_label(metric), levels = alfak_metric_label(c("r", "R", "rho")))]
  z[, grf_lambda := factor(grf_lambda, levels = sort(unique(grf_lambda)))]
  if (isTRUE(use_minobs_axis)) {
    x_minobs <- if ("display_minobs" %in% names(z)) z$display_minobs else z$minobs
    z[, x_group := sprintf("%02d", as.integer(x_minobs))]
    z[, color_group := factor(sample_depth)]
    x_lab <- "Minimum observations"
    color_lab <- "Sample depth"
    group_cols <- c("grf_lambda", "metric", "x_group", "color_group")
  } else {
    z[, x_group := factor(sample_depth)]
    z[, color_group := factor(sample_depth)]
    x_lab <- "Sample depth"
    color_lab <- "Sample depth"
    group_cols <- c("grf_lambda", "metric", "x_group", "color_group")
  }
  q <- z[, {
    qq <- stats::quantile(value, probs = c(0.1, 0.5, 0.9), names = FALSE, na.rm = TRUE)
    .(lo = qq[[1]], med = qq[[2]], hi = qq[[3]])
  }, by = group_cols]
  ggplot2::ggplot(q, ggplot2::aes(x = x_group, y = med, color = color_group, group = color_group)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.18, position = ggplot2::position_dodge(width = 0.65)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.65), size = 1.4) +
    ggplot2::facet_grid(metric ~ grf_lambda, scales = "free_y", labeller = ggplot2::labeller(grf_lambda = function(x) paste0("lambda=", x))) +
    ggplot2::scale_y_continuous("Metric value") +
    ggplot2::scale_x_discrete(x_lab, expand = ggplot2::expansion(add = c(0.2, 0.8))) +
    ggplot2::scale_color_viridis_d(color_lab) +
    ggplot2::theme(legend.position = "top")
}

plot_suppfig2_panel_b <- function(method_rows, nfq_threshold) {
  needed <- c("R", "Rfq", "nfq")
  if (!nrow(method_rows) || length(setdiff(needed, names(method_rows)))) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  z <- method_rows[is.finite(R) & is.finite(Rfq) & is.finite(nfq)]
  if (!nrow(z)) return(ggplot2::ggplot() + ggplot2::theme_void())
  z[, R_plot := pmax(-1, R)]
  nfq_levels <- c(paste0("<=", nfq_threshold), paste0(">", nfq_threshold))
  rfq_levels <- c("Rfq > 0", "Rfq < 0")
  z[, nfq_bins := factor(
    ifelse(nfq > nfq_threshold, paste0(">", nfq_threshold), paste0("<=", nfq_threshold)),
    levels = nfq_levels
  )]
  z[, rfq_bin := factor(
    ifelse(Rfq > 0, "Rfq > 0", "Rfq < 0"),
    levels = rfq_levels
  )]
  skeleton <- data.table::CJ(
    nfq_bins = factor(nfq_levels, levels = nfq_levels),
    rfq_bin = factor(rfq_levels, levels = rfq_levels)
  )
  skeleton[, R_plot := min(z$R_plot, na.rm = TRUE)]
  ggplot2::ggplot(z, ggplot2::aes(x = nfq_bins, y = R_plot)) +
    ggplot2::geom_blank(data = skeleton, ggplot2::aes(x = nfq_bins, y = R_plot)) +
    ggplot2::geom_violin(fill = "grey85", color = "grey35", linewidth = 0.25, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.3, alpha = 0.85) +
    ggplot2::facet_grid(rfq_bin ~ ., scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_discrete("num. frequent karyotypes", drop = FALSE) +
    ggplot2::scale_y_continuous("Rescaled R^2 (landscape)") +
    ggplot2::coord_flip()
}

plot_suppfig2_panel_c <- function(scatter_dt, use_minobs_axis, lambda_page = FALSE) {
  if (!nrow(scatter_dt)) return(ggplot2::ggplot() + ggplot2::theme_void())
  scatter_dt[, depth_label := paste0("depth ", sample_depth)]
  scatter_dt[, lambda_label := paste0("lambda=", grf_lambda)]
  if (isTRUE(use_minobs_axis)) {
    c_minobs <- if ("display_minobs" %in% names(scatter_dt)) scatter_dt$display_minobs else scatter_dt$minobs
    scatter_dt[, minobs_label := paste0("N=", sprintf("%02d", as.integer(c_minobs)))]
    if (isTRUE(lambda_page)) {
      facet_formula <- stats::as.formula("depth_label ~ minobs_label")
      cor_by <- c("depth_label", "minobs_label")
    } else {
      facet_formula <- stats::as.formula("depth_label ~ lambda_label + minobs_label")
      cor_by <- c("depth_label", "lambda_label", "minobs_label")
    }
  } else {
    facet_formula <- stats::as.formula("depth_label ~ lambda_label")
    cor_by <- c("depth_label", "lambda_label")
  }
  cor_dt <- scatter_dt[, {
    ok <- is.finite(truth_centered) & is.finite(pred_centered)
    r <- if (sum(ok) >= 2L && stats::sd(truth_centered[ok]) > 0 && stats::sd(pred_centered[ok]) > 0) {
      suppressWarnings(stats::cor(truth_centered[ok], pred_centered[ok], method = "pearson"))
    } else {
      NA_real_
    }
    .(pearson_label = ifelse(is.finite(r), sprintf("r = %.2f", r), "r = NA"))
  }, by = cor_by]
  ggplot2::ggplot(scatter_dt, ggplot2::aes(x = truth_centered, y = pred_centered)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.25, color = "grey55") +
    ggplot2::geom_point(size = 0.45, alpha = 0.45) +
    ggplot2::geom_text(
      data = cor_dt,
      ggplot2::aes(x = -Inf, y = Inf, label = pearson_label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.2,
      size = 2.4
    ) +
    ggplot2::facet_grid(facet_formula, scales = "free") +
    ggplot2::scale_x_continuous("true fitness", breaks = scales::pretty_breaks(n = 2)) +
    ggplot2::scale_y_continuous("estimated fitness", breaks = scales::pretty_breaks(n = 2))
}

plot_suppfig2_panel_f <- function(method_rows, nfq_threshold, use_minobs_axis) {
  needed <- c("R", "Rfq", "nfq", "sample_depth", "grf_lambda")
  if (!nrow(method_rows) || length(setdiff(needed, names(method_rows)))) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  z <- method_rows[is.finite(R) & is.finite(Rfq) & is.finite(nfq)]
  if (!nrow(z)) return(ggplot2::ggplot() + ggplot2::theme_void())
  z[, `:=`(
    lambda_axis = paste0("lambda ", grf_lambda),
    depth_axis = paste0("depth ", sample_depth),
    rfq_axis = ifelse(Rfq > 0, "Rfq > 0", "Rfq < 0"),
    nfq_axis = ifelse(nfq > nfq_threshold, paste0("Nfq > ", nfq_threshold), paste0("Nfq <= ", nfq_threshold)),
    r_axis = ifelse(R > 0, "R2 > 0", "R2 < 0")
  )]
  z[, fillvar := "v1"]
  z[Rfq <= 0, fillvar := "v2"]
  z[Rfq > 0 & nfq <= nfq_threshold, fillvar := "v3"]
  if (isTRUE(use_minobs_axis)) {
    z[, minobs_axis := paste0("minobs ", if ("display_minobs" %in% names(z)) display_minobs else minobs)]
    agg <- z[, .N, by = .(lambda_axis, depth_axis, minobs_axis, rfq_axis, nfq_axis, r_axis, fillvar)]
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(axis1 = lambda_axis, axis2 = depth_axis, axis3 = minobs_axis,
                   axis4 = rfq_axis, axis5 = nfq_axis, axis6 = r_axis, y = N)
    ) +
      ggalluvial::geom_alluvium(ggplot2::aes(fill = fillvar), alpha = 0.75) +
      ggalluvial::geom_stratum(width = 1 / 4, color = "grey25", fill = "grey92") +
      ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 2) +
      ggplot2::scale_x_discrete(
        limits = c("axis1", "axis2", "axis3", "axis4", "axis5", "axis6"),
        labels = c("lambda", "sample\ndepth", "minobs", "Rfq > 0", paste0("Nfq > ", nfq_threshold), "R2 > 0"),
        expand = c(0.08, 0.08)
      )
  } else {
    agg <- z[, .N, by = .(lambda_axis, depth_axis, rfq_axis, nfq_axis, r_axis, fillvar)]
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(axis1 = lambda_axis, axis2 = depth_axis,
                   axis3 = rfq_axis, axis4 = nfq_axis, axis5 = r_axis, y = N)
    ) +
      ggalluvial::geom_alluvium(ggplot2::aes(fill = fillvar), alpha = 0.75) +
      ggalluvial::geom_stratum(width = 1 / 4, color = "grey25", fill = "grey92") +
      ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 2) +
      ggplot2::scale_x_discrete(
        limits = c("axis1", "axis2", "axis3", "axis4", "axis5"),
        labels = c("lambda", "sample\ndepth", "Rfq > 0", paste0("Nfq > ", nfq_threshold), "R2 > 0"),
        expand = c(0.08, 0.08)
      )
  }
  p +
    ggplot2::scale_fill_manual(
      name = "",
      values = c(v1 = "grey70", v2 = "#21908CFF", v3 = "#440154FF"),
      breaks = c("v2", "v3"),
      labels = c("Rfq < 0", paste0("Nfq <= ", nfq_threshold))
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

write_suppfig2_style_method_plots <- function() {
  supp_dir <- file.path(figures_dir, "suppfig2_style")
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  old_png <- list.files(supp_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
  if (length(old_png)) unlink(old_png, force = TRUE)
  if (!nrow(alfak_style_summary)) {
    write_tsv(data.table::data.table(), file.path(tables_dir, "hpc_suppfig2_style_plot_index.tsv"))
    return(data.table::data.table())
  }

  keys <- data.table::rbindlist(list(
    unique(alfak_style_summary[package == "alfak2", .(
      package, package_display, input_mode, extrapolation_method, minobs = NA_integer_,
      display_minobs = NA_integer_, NN_prior_slot = "",
      method_label = paste(package_display, input_mode, extrapolation_method, sep = ":"),
      method_short = short_method(paste(package_display, input_mode, extrapolation_method, sep = ":"))
    )]),
    unique(alfak_style_summary[package == "alfakR", .(
      package, package_display, input_mode = "", extrapolation_method = "",
      minobs = NA_integer_, display_minobs = NA_integer_, NN_prior_slot,
      method_label = paste0("alfak:", NN_prior_slot),
      method_short = label_slot_short(NN_prior_slot)
    )])
  ), use.names = TRUE, fill = TRUE)
  data.table::setorder(keys, package_display, input_mode, extrapolation_method, NN_prior_slot)

  nfq_threshold <- 8L
  alfak_c_lambdas <- c(0.2, 0.8)
  rows <- list()
  row_idx <- 0L
  for (i in seq_len(nrow(keys))) {
    key <- keys[i]
    is_alfak <- identical(as.character(key$package[[1]]), "alfakR")
    use_minobs_axis <- is_alfak || identical(as.character(key$input_mode[[1]]), "soft_minobs")
    if (is_alfak) {
      d <- alfak_style_summary[package == "alfakR" & NN_prior_slot == key$NN_prior_slot[[1]]]
      by_cols <- c("sample_depth", "display_minobs")
      lambda_pages <- select_lambda_values(d$grf_lambda, requested = alfak_c_lambdas)
    } else {
      d <- alfak_style_summary[
        package == "alfak2" &
          input_mode == key$input_mode[[1]] &
          extrapolation_method == key$extrapolation_method[[1]]
      ]
      by_cols <- if (use_minobs_axis) c("sample_depth", "display_minobs") else c("sample_depth")
      lambda_pages <- NA_real_
    }
    if (!nrow(d)) next
    d_landscape <- suppfig2_landscape_level(d)
    panel_a <- plot_suppfig2_panel_a(d_landscape, use_minobs_axis = use_minobs_axis)
    panel_b <- plot_suppfig2_panel_b(d_landscape, nfq_threshold = nfq_threshold)
    panel_f <- plot_suppfig2_panel_f(d_landscape, nfq_threshold = nfq_threshold, use_minobs_axis = use_minobs_axis)

    for (lambda_page in lambda_pages) {
      if (is_alfak) {
        d_scatter <- d[grf_lambda == lambda_page]
        lambda_tag <- paste0("lambda", gsub("\\.", "p", as.character(lambda_page)))
        out_file <- file.path("suppfig2_style", "alfak", safe_slug(key$NN_prior_slot[[1]]), paste0(lambda_tag, ".png"))
        lambda_title <- paste0(" (C panel lambda=", lambda_page, ")")
        rep_rows <- pick_representative_rows(
          d_scatter,
          by_cols = by_cols,
          lambda_values = lambda_page,
          lambda_extremes = FALSE
        )
      } else {
        d_scatter <- d
        out_file <- file.path(
          "suppfig2_style", "alfak_V2",
          safe_slug(key$input_mode[[1]]),
          paste0(safe_slug(key$extrapolation_method[[1]]), ".png")
        )
        lambda_title <- " (C panel all lambdas)"
        rep_rows <- pick_representative_rows(
          d_scatter,
          by_cols = by_cols,
          lambda_values = sort(unique(d_scatter$grf_lambda[is.finite(d_scatter$grf_lambda)])),
          lambda_extremes = FALSE
        )
      }
      dir.create(dirname(file.path(figures_dir, out_file)), recursive = TRUE, showWarnings = FALSE)
      scatter <- read_representative_scatter(rep_rows)
      panel_c <- plot_suppfig2_panel_c(scatter, use_minobs_axis = use_minobs_axis, lambda_page = is_alfak)
      title <- cowplot::ggdraw() +
        cowplot::draw_label(
          paste0(key$method_short[[1]], " Supplementary Fig. 2-style ground-truth validation", lambda_title),
          x = 0, y = 0.70, hjust = 0, fontface = "bold", size = 11
        ) +
        cowplot::draw_label(
          "Direct truth-vs-estimated fitness reliability metrics; no forward-prediction or trajectory validation is added. Panels A/B/F use landscape-level summaries.",
          x = 0, y = 0.28, hjust = 0, size = 8
        )
      top <- cowplot::plot_grid(panel_a, panel_b, labels = c("A", "B"), label_size = 10, rel_widths = c(2.5, 1))
      mid <- cowplot::plot_grid(panel_c, labels = c("C"), label_size = 10)
      bottom <- cowplot::plot_grid(panel_f, labels = c("F"), label_size = 10)
      plt <- cowplot::plot_grid(title, top, mid, bottom, nrow = 4, rel_heights = c(0.35, 3.2, 2.4, 2.2))
      ggplot2::ggsave(file.path(figures_dir, out_file), plot = plt, width = 10, height = 12, dpi = 180, bg = "white")
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.table::data.table(
        package = key$package[[1]],
        package_display = key$package_display[[1]],
        input_mode = key$input_mode[[1]],
        extrapolation_method = key$extrapolation_method[[1]],
        NN_prior_slot = key$NN_prior_slot[[1]],
        method_label = key$method_label[[1]],
        method_short = key$method_short[[1]],
        c_panel_lambda = if (is_alfak) lambda_page else NA_real_,
        figure = file.path("figures", out_file)
      )
    }
  }
  idx <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  write_tsv(idx, file.path(tables_dir, "hpc_suppfig2_style_plot_index.tsv"))
  idx
}

cleanup_excluded_method_artifacts <- function() {
  excluded <- unique(dt[
    display_slot_keep == FALSE & package == "alfakR",
    .(package, package_display, input_mode, extrapolation_method, minobs, NN_prior_slot, method_label, method_short)
  ])
  if (!nrow(excluded)) return(invisible(character()))
  removed <- character()
  for (i in seq_len(nrow(excluded))) {
    out_dir <- method_artifact_path(excluded[i])
    if (dir.exists(out_dir)) {
      unlink(out_dir, recursive = TRUE, force = TRUE)
      removed <- c(removed, out_dir)
    }
  }
  invisible(removed)
}

write_method_artifacts <- function() {
  method_keys <- unique(analysis_dt[, .(
    package, package_display, input_mode, extrapolation_method, minobs, NN_prior_slot,
    method_label, method_short
  )])
  data.table::setorder(method_keys, package_display, input_mode, minobs, extrapolation_method, NN_prior_slot)
  index_rows <- vector("list", nrow(method_keys))
  for (i in seq_len(nrow(method_keys))) {
    row <- method_keys[i]
    out_dir <- method_artifact_path(row)
    out_tables <- file.path(out_dir, "tables")
    out_figures <- file.path(out_dir, "figures")
    dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)
    d <- analysis_dt[method_label == row$method_label[[1]]]
    ms <- method_summary[method_label == row$method_label[[1]]]
    os <- overall_summary[method_label == row$method_label[[1]]]
    ls <- landscape_method_summary[method_label == row$method_label[[1]]]
    ak <- alfak_style_summary[method_label == row$method_label[[1]]]
    write_tsv(d, file.path(out_tables, "metrics_by_run.tsv"))
    write_tsv(ms, file.path(out_tables, "method_metric_summary.tsv"))
    write_tsv(os, file.path(out_tables, "overall_method_metric_summary.tsv"))
    write_tsv(ls, file.path(out_tables, "landscape_fit_repeat_summary.tsv"))
    write_tsv(ak, file.path(out_tables, "alfak_style_ground_truth_summary.tsv"))

    ak_full_plot <- FALSE
    ak_full_metrics <- intersect(c("r", "R", "rho"), names(ak))
    if (nrow(ak) && length(ak_full_metrics)) {
      ak_full <- data.table::melt(
        ak,
        id.vars = intersect(c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat"), names(ak)),
        measure.vars = ak_full_metrics,
        variable.name = "alfak_metric",
        value.name = "value"
      )
      ak_full <- ak_full[is.finite(value)]
      if (nrow(ak_full)) {
        ak_full[, value := pmax(-1, value)]
        ak_full[, alfak_metric_label := factor(
          alfak_metric_label(alfak_metric),
          levels = alfak_metric_label(c("r", "R", "rho"))
        )]
        ak_full_q <- quantile_summary(
          ak_full,
          by = c("sample_depth", "grf_lambda", "alfak_metric", "alfak_metric_label")
        )
        ak_full_q[, grf_lambda := factor(grf_lambda, levels = sort(unique(grf_lambda)))]
        p_ak_full <- ggplot2::ggplot(
          ak_full_q,
          ggplot2::aes(x = factor(sample_depth), y = med, group = 1)
        ) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.18) +
          ggplot2::geom_point(size = 1.8) +
          ggplot2::facet_grid(alfak_metric_label ~ grf_lambda, scales = "free_y") +
          ggplot2::scale_y_continuous("Metric value") +
          ggplot2::scale_x_discrete("Sample depth") +
          ggplot2::labs(
            title = paste0(row$method_label[[1]], " ALFA-K original-style r/R/rho"),
            subtitle = "Points are medians; bars are 10th-90th percentiles across landscapes and fit repeats."
          )
        ggplot2::ggsave(
          file.path(out_figures, "alfak_original_style_r_R_rho.png"),
          plot = p_ak_full,
          width = 12,
          height = 7,
          dpi = 180
        )
        ak_full_plot <- TRUE
      }
    }

    ak_subset_plot <- FALSE
    ak_subset_metrics <- intersect(alfak_style_metric_map$alfak_metric, names(ak))
    if (nrow(ak) && length(ak_subset_metrics)) {
      ak_subset <- data.table::melt(
        ak,
        id.vars = intersect(c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat"), names(ak)),
        measure.vars = ak_subset_metrics,
        variable.name = "alfak_metric",
        value.name = "value"
      )
      ak_subset <- merge(
        ak_subset[is.finite(value)],
        unique(alfak_style_metric_map[, .(alfak_metric, alfak_subset, alfak_family)]),
        by = "alfak_metric",
        all.x = TRUE
      )
      if (nrow(ak_subset)) {
        ak_subset[, value := pmax(-1, value)]
        ak_subset[, alfak_family := factor(alfak_family, levels = c("Pearson r", "Spearman rho", "Rescaled R^2"))]
        ak_subset[, alfak_subset := factor(alfak_subset, levels = c("full", "fq", "nn", "d2"))]
        ak_subset_q <- quantile_summary(
          ak_subset,
          by = c("sample_depth", "grf_lambda", "alfak_metric", "alfak_family", "alfak_subset")
        )
        p_ak_subset <- ggplot2::ggplot(
          ak_subset_q,
          ggplot2::aes(x = factor(grf_lambda), y = med, color = factor(sample_depth), group = sample_depth)
        ) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.08, position = ggplot2::position_dodge(width = 0.45)) +
          ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45), size = 1.6) +
          ggplot2::facet_grid(alfak_family ~ alfak_subset, scales = "free_y") +
          ggplot2::scale_y_continuous("Metric value") +
          ggplot2::scale_x_discrete("GRF lambda") +
          ggplot2::labs(
            title = paste0(row$method_label[[1]], " ALFA-K subset metrics"),
            subtitle = "ALFA-K fq/nn/d2 columns are mapped to benchmark d0/d1/d2 shells.",
            color = "Sample depth"
          )
        ggplot2::ggsave(
          file.path(out_figures, "alfak_original_style_subset_metrics.png"),
          plot = p_ak_subset,
          width = 13,
          height = 8,
          dpi = 180
        )
        ak_subset_plot <- TRUE
      }
    }

    gt <- d[
      metric %in% c("pearson", "spearman", "rescaled_r2") &
        shell %in% shell_order &
        is.finite(value)
    ]
    if (nrow(gt)) {
      gt[, metric_label := metric_label(metric)]
      gt[, shell := factor(shell, levels = shell_order)]
      gt[, violin_n := .N, by = .(metric_label, shell, grf_lambda, sample_depth)]
      gt_violin <- gt[violin_n >= 2L]
      p <- ggplot2::ggplot(
        gt,
        ggplot2::aes(x = factor(grf_lambda), y = value, fill = factor(sample_depth))
      )
      if (nrow(gt_violin)) {
        p <- p +
          ggplot2::geom_violin(
            data = gt_violin,
            position = ggplot2::position_dodge(width = 0.75),
            width = 0.72,
            alpha = 0.35,
            color = NA,
            trim = TRUE
          )
      }
      p <- p +
        ggplot2::geom_boxplot(
          position = ggplot2::position_dodge(width = 0.75),
          width = 0.18,
          outlier.size = 0.35,
          alpha = 0.85
        ) +
        ggplot2::geom_point(
          position = ggplot2::position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
          size = 0.45,
          alpha = 0.35,
          shape = 21,
          stroke = 0.1
        ) +
        ggplot2::facet_grid(metric_label ~ shell, scales = "free_y") +
        ggplot2::labs(
          title = paste0(row$method_label[[1]], " ground-truth landscape metrics"),
          subtitle = "ALFA-K-style fitness ranking/shape metrics; full_lscape uses the method's complete prediction table.",
          x = "GRF lambda",
          y = "Metric value",
          fill = "Depth"
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
      ggplot2::ggsave(
        file.path(out_figures, "alfak_style_ground_truth_metrics.png"),
        plot = p,
        width = 13,
        height = 8,
        dpi = 180
      )
    }

    top_overall <- os[
      metric %in% c("pearson", "spearman", "rescaled_r2", "rmse", "centered_rmse", "edge_gradient_spearman", "sign_accuracy") &
        shell %in% shell_order
    ][order(shell, metric)]
    top_report <- top_overall[, .(
      shell,
      metric,
      direction,
      n,
      mean = fmt_num(mean),
      median = fmt_num(median),
      q25 = fmt_num(q25),
      q75 = fmt_num(q75)
    )]
    report <- c(
      paste0("# ", row$method_label[[1]]),
      "",
      paste0("- package: `", row$package_display[[1]], "`"),
      if (identical(as.character(row$package[[1]]), "alfak2")) paste0("- input_mode: `", row$input_mode[[1]], "`") else paste0("- minobs: `", row$minobs[[1]], "`"),
      if (identical(as.character(row$package[[1]]), "alfak2")) paste0("- extrapolation_method: `", row$extrapolation_method[[1]], "`") else paste0("- NN_prior_slot: `", row$NN_prior_slot[[1]], "`"),
      "",
      "## Tables",
      "- `tables/metrics_by_run.tsv`",
      "- `tables/method_metric_summary.tsv`",
      "- `tables/overall_method_metric_summary.tsv`",
      "- `tables/landscape_fit_repeat_summary.tsv`",
      "- `tables/alfak_style_ground_truth_summary.tsv`",
      "",
      "## Figures",
      if (ak_full_plot) "- `figures/alfak_original_style_r_R_rho.png`" else "- ALFA-K original-style r/R/rho figure was skipped because no matching rows were available.",
      if (ak_subset_plot) "- `figures/alfak_original_style_subset_metrics.png`" else "- ALFA-K subset metric figure was skipped because no matching rows were available.",
      if (nrow(gt)) "- `figures/alfak_style_ground_truth_metrics.png`" else "- Ground-truth metric figure was skipped because no matching rows were available.",
      "",
      "## ALFA-K-Style Ground-Truth Metric Summary",
      as_markdown_table(top_report)
    )
    writeLines(report, file.path(out_dir, "report.md"))
    index_rows[[i]] <- data.table::data.table(
      package = row$package[[1]],
      package_display = row$package_display[[1]],
      input_mode = row$input_mode[[1]],
      extrapolation_method = row$extrapolation_method[[1]],
      minobs = row$minobs[[1]],
      NN_prior_slot = row$NN_prior_slot[[1]],
      method_label = row$method_label[[1]],
      method_short = row$method_short[[1]],
      report_dir = out_dir
    )
  }
  idx <- data.table::rbindlist(index_rows, use.names = TRUE, fill = TRUE)
  write_tsv(idx, file.path(input_dir, "method_reports", "method_report_index.tsv"))
  idx
}

write_metric_comparison_plots <- function() {
  metric_dir <- file.path(figures_dir, "metric_comparison")
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
  metric_files <- character()
  plot_metrics <- intersect(selected_metrics, sort(unique(analysis_dt$metric)))
  for (mm in plot_metrics) {
    pd <- analysis_dt[
      metric == mm &
        shell %in% shell_order &
        (
          package == "alfakR" |
            (package == "alfak2" & input_mode %in% c("full", "soft_minobs"))
        ) &
        is.finite(value)
    ]
    if (!nrow(pd) || length(unique(pd$method_label)) < 2L) next
    pd[, package_display := factor(package_display, levels = c("alfak_V2", "alfak"))]
    pd[, shell := factor(shell, levels = shell_order)]
    pd[, facet_group := ifelse(
      package == "alfak2",
      ifelse(input_mode == "soft_minobs", paste0("soft_minobs", as.integer(display_minobs)), as.character(input_mode)),
      paste0("minobs", as.integer(display_minobs))
    )]
    pd[, method_in_group := ifelse(
      package == "alfak2",
      short_method(extrapolation_method),
      label_slot_short(NN_prior_slot)
    )]
    pd[, facet_group := factor(facet_group, levels = c(
      "full", "soft_minobs5", "soft_minobs10", "soft_minobs20",
      "minobs5", "minobs10", "minobs20"
    ))]
    method_levels <- unique(pd[order(package_display, facet_group, method_in_group), method_in_group])
    pd[, method_in_group := factor(method_in_group, levels = method_levels)]
    pd[, violin_n := .N, by = .(package_display, facet_group, shell, method_in_group)]
    pd_violin <- pd[violin_n >= 2L]
    p <- ggplot2::ggplot(
      pd,
      ggplot2::aes(x = method_in_group, y = value, fill = package_display)
    )
    if (nrow(pd_violin)) {
      p <- p +
        ggplot2::geom_violin(
          data = pd_violin,
          width = 0.88,
          alpha = 0.35,
          trim = TRUE,
          color = NA,
          scale = "width"
        )
    }
    p <- p +
      ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.25, alpha = 0.85) +
      ggplot2::geom_jitter(width = 0.08, height = 0, size = 0.35, alpha = 0.22) +
      ggplot2::facet_grid(shell ~ package_display + facet_group, scales = "free_x", space = "free_x") +
      ggplot2::labs(
        title = paste0(metric_label(mm), " by method"),
        subtitle = "alfak_V2 is grouped by input/minobs; alfak is grouped by minobs and then NN_prior. ecw uses slot4 only.",
        x = NULL,
        y = metric_label(mm),
        fill = "Package"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"
      )
    fname <- file.path("metric_comparison", paste0("fig_metric_", safe_slug(mm), "_violin_boxplot.png"))
    ggplot2::ggsave(file.path(figures_dir, fname), plot = p, width = 14, height = 12, dpi = 180)
    metric_files <- c(metric_files, fname)
  }
  metric_files
}

cleanup_excluded_method_artifacts()
suppfig2_style_index <- write_suppfig2_style_method_plots()
method_report_index <- write_method_artifacts()
metric_comparison_files <- write_metric_comparison_plots()

numeric_plot_dt <- delta_best_summary[
  metric %in% c("rmse", "mae", "relative_rmse") &
    shell %in% c("d0", "d1", "d2", "all_nearfield", "all_eval", "full_lscape")
]
numeric_plot_dt[, facet := paste(shell, metric, sep = " / ")]
p1 <- ggplot2::ggplot(
  numeric_plot_dt,
  ggplot2::aes(x = stats::reorder(alfak2_short, mean_improvement), y = mean_improvement, fill = alfak2_short)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(metric ~ shell, scales = "free_x") +
  ggplot2::labs(
    title = "alfak_V2 numerical improvement vs best alfak",
    subtitle = "Positive means alfak_V2 has lower error than the best alfak setting in the same condition.",
    x = NULL,
    y = "Mean improvement"
  )
plot_save(p1, "fig01_numerical_improvement_vs_best_alfak.png", width = 13, height = 8)

shape_plot_dt <- delta_best_summary[
  metric %in% c("centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "edge_gradient_spearman", "sign_accuracy") &
    shell %in% c("d0", "d1", "d2", "all_nearfield", "all_eval", "full_lscape")
]
p2 <- ggplot2::ggplot(
  shape_plot_dt,
  ggplot2::aes(x = stats::reorder(alfak2_short, mean_improvement), y = mean_improvement, fill = alfak2_short)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(metric ~ shell, scales = "free_x") +
  ggplot2::labs(
    title = "alfak_V2 shape improvement vs best alfak",
    subtitle = "Positive is better after respecting each metric's direction.",
    x = NULL,
    y = "Mean improvement"
  )
plot_save(p2, "fig02_shape_improvement_vs_best_alfak.png", width = 13, height = 10)

trade_plot <- overall_rank_wide[shell == "all_nearfield"]
trade_plot[, package_display := package_display(package)]
p3 <- ggplot2::ggplot(
  trade_plot,
  ggplot2::aes(x = numerical, y = shape, color = package_display, label = method_short)
) +
  ggplot2::geom_point(size = 2.4, alpha = 0.85) +
  ggplot2::geom_text(
    data = trade_plot[package == "alfak2"],
    ggplot2::aes(label = method_short),
    hjust = 0, nudge_x = 0.15, size = 3, show.legend = FALSE
  ) +
  ggplot2::scale_y_reverse() +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "All-nearfield numerical vs shape rank",
    subtitle = "Ranks are averaged across selected numerical and shape metrics; closer to the upper-right is better.",
    x = "Mean numerical rank (lower is better)",
    y = "Mean shape rank (lower is better)"
  )
plot_save(p3, "fig03_all_nearfield_tradeoff_rank.png", width = 11, height = 7)

top_methods <- overall_rank_wide[shell == "all_nearfield"][order(balanced_rank)][1:min(.N, 8), method_label]
trend_plot_dt <- lambda_trend[
  shell == "all_nearfield" &
    metric %in% c("rmse", "centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "sign_accuracy") &
    method_label %in% top_methods
]
trend_line_dt <- data.table::copy(trend_plot_dt)
trend_line_dt[, line_n := data.table::uniqueN(grf_lambda), by = .(metric, sample_depth, method_short)]
trend_line_dt <- trend_line_dt[line_n >= 2L]
p4 <- ggplot2::ggplot(
  trend_plot_dt,
  ggplot2::aes(x = grf_lambda, y = mean, color = method_short, group = method_short)
)
if (nrow(trend_line_dt)) {
  p4 <- p4 +
    ggplot2::geom_line(data = trend_line_dt, linewidth = 0.55)
}
p4 <- p4 +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::facet_grid(metric ~ sample_depth, scales = "free_y") +
  ggplot2::labs(
    title = "Accuracy trend across landscape smoothness",
    subtitle = "Larger GRF lambda means smoother landscapes; lower lambda is more rugged.",
    x = "GRF lambda",
    y = "Mean metric value",
    color = "Method"
  )
plot_save(p4, "fig04_ruggedness_lambda_trend_top_methods.png", width = 13, height = 10)

gap_cols <- grep("_gap_to_best$", names(alfak2_tradeoff), value = TRUE)
gap_dt <- data.table::melt(
  alfak2_tradeoff[, c("method_label", "method_short", gap_cols), with = FALSE],
  id.vars = c("method_label", "method_short"),
  variable.name = "gap_metric",
  value.name = "gap"
)
gap_dt[, gap_metric := sub("_gap_to_best$", "", gap_metric)]
gap_dt <- gap_dt[is.finite(gap)]
p5 <- ggplot2::ggplot(
  gap_dt,
  ggplot2::aes(x = stats::reorder(method_short, gap), y = gap, fill = method_short)
) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ gap_metric, scales = "free_x") +
  ggplot2::labs(
    title = "alfak_V2 gap to best alfak_V2 all-nearfield method",
    subtitle = "For error metrics lower gaps are better; for correlation/accuracy metrics the gap is best minus method.",
    x = NULL,
    y = "Gap to best"
  )
plot_save(p5, "fig05_alfak_V2_gap_to_best_all_nearfield.png", width = 13, height = 8)

heat_dt <- delta_best_summary[
  shell %in% c("d0", "d1", "d2", "all_nearfield", "all_eval", "full_lscape") &
    metric %in% c("rmse", "mae", "centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "sign_accuracy")
]
heat_dt[, label := sprintf("%.0f%%", 100 * improvement_rate)]
p6 <- ggplot2::ggplot(
  heat_dt,
  ggplot2::aes(x = metric, y = alfak2_short, fill = improvement_rate)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
  ggplot2::facet_grid(shell ~ sample_depth) +
  ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0.5, limits = c(0, 1)) +
  ggplot2::labs(
    title = "Win rate vs best alfak by condition",
    subtitle = "Each cell is the fraction of matched conditions where alfak_V2 beats the best alfak setting.",
    x = NULL,
    y = NULL,
    fill = "Win rate"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
plot_save(p6, "fig06_win_rate_vs_best_alfak_heatmap.png", width = 14, height = 10)

best_overall <- overall_best_by_package[
  shell %in% c("d0", "d1", "d2", "all_nearfield", "all_eval", "full_lscape") &
    metric %in% c("rmse", "mae", "centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "sign_accuracy")
]
best_overall[, value_text := fmt_num(mean)]

best_compare <- data.table::dcast(
  best_overall[, .(shell, metric, direction, package, method_short, mean)],
  shell + metric + direction ~ package,
  value.var = c("method_short", "mean")
)
if (!"method_short_alfak2" %in% names(best_compare)) best_compare[, method_short_alfak2 := NA_character_]
if (!"method_short_alfakR" %in% names(best_compare)) best_compare[, method_short_alfakR := NA_character_]
if (!"mean_alfak2" %in% names(best_compare)) best_compare[, mean_alfak2 := NA_real_]
if (!"mean_alfakR" %in% names(best_compare)) best_compare[, mean_alfakR := NA_real_]
best_compare[, improvement_vs_fixed_best := ifelse(
  direction == "lower",
  mean_alfakR - mean_alfak2,
  mean_alfak2 - mean_alfakR
)]
best_compare[, relative_improvement_vs_fixed_best := improvement_vs_fixed_best / pmax(abs(mean_alfakR), .Machine$double.eps)]
best_compare[, result := ifelse(improvement_vs_fixed_best > 0, "alfak_V2 better", "alfak_V2 worse")]
best_compare_report <- best_compare[
  order(shell, metric),
  .(
    shell,
    metric,
    direction,
    alfak_V2_method = method_short_alfak2,
    alfak_V2_mean = fmt_num(mean_alfak2),
    alfak_method = method_short_alfakR,
    alfak_mean = fmt_num(mean_alfakR),
    improvement = fmt_num(improvement_vs_fixed_best),
    relative_improvement = fmt_pct(relative_improvement_vs_fixed_best),
    result
  )
]

best_delta_oracle <- delta_best_summary[
  shell %in% c("d0", "d1", "d2", "all_nearfield", "all_eval", "full_lscape") &
    metric %in% c("rmse", "centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "sign_accuracy"),
  .SD[which.max(mean_improvement)],
  by = .(shell, metric)
]
best_delta_oracle_report <- best_delta_oracle[
  order(shell, metric),
  .(
    shell,
    metric,
    alfak_V2_method = alfak2_short,
    mean_improvement = fmt_num(mean_improvement),
    win_rate = fmt_pct(improvement_rate)
  )
]

best_a2_tradeoff <- alfak2_tradeoff[1]
best_num <- alfak2_tradeoff[which.min(rmse)]
best_shape_edge <- if ("edge_gradient_rmse" %in% names(alfak2_tradeoff)) {
  alfak2_tradeoff[which.min(edge_gradient_rmse)]
} else data.table::data.table()
best_shape_center <- alfak2_tradeoff[which.min(centered_rmse)]

lambda_effect_compact <- rugged_effect[
  shell == "all_nearfield" &
    metric %in% c("rmse", "centered_rmse", "rescaled_r2", "edge_gradient_rmse", "spearman", "sign_accuracy")
][
  , .(
    n_methods = .N,
    frac_worse_when_rugged = mean(worse_when_more_rugged, na.rm = TRUE),
    median_cor = stats::median(cor_vs_ruggedness, na.rm = TRUE)
  ),
  by = .(sample_depth, package, metric, domain, direction)
]
write_tsv(lambda_effect_compact, file.path(tables_dir, "hpc_ruggedness_effect_compact.tsv"))

lambda_effect_report <- lambda_effect_compact[
  order(sample_depth, package, metric),
  .(
    sample_depth,
    package = package_display(package),
    metric,
    direction,
    n_methods,
    frac_worse_when_rugged = fmt_pct(frac_worse_when_rugged),
    median_cor_vs_ruggedness = fmt_num(median_cor)
  )
]

report_lines <- c(
  "# HPC selected depth/lambda comparison summary",
  "",
  paste0("- source-dir: `", input_dir, "`"),
  paste0("- prediction_scale: `", args$prediction_scale, "`"),
  paste0("- generated: `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "`"),
  paste0("- successful/fallback metric rows used: ", format(nrow(analysis_dt), big.mark = ",")),
  paste0("- fit status summary table: `tables/hpc_run_status_summary.tsv`"),
  "- duplicate alfak `empirical_censored_weighted` slots: display keeps the lowest numeric slot and excludes higher duplicate slots.",
  "",
  "## Key outputs",
  "- `tables/hpc_method_metric_summary.tsv`: method means by sample depth, lambda, shell, and metric.",
  "- `tables/hpc_alfak_style_ground_truth_summary.tsv`: ALFA-K original ground-truth columns (`r/rho/R`, `fq/nn/d2` variants) mapped to benchmark shells.",
  "- `tables/hpc_landscape_fit_repeat_summary.tsv`: within-landscape summaries across fit repeats.",
  "- `tables/hpc_landscape_block_method_metric_summary.tsv`: landscape-block summaries after first aggregating fit repeats within each landscape.",
  "- `tables/hpc_alfak_V2_delta_vs_best_alfak_summary.tsv`: alfak_V2 improvement against the best alfak setting in matched conditions.",
  "- `tables/hpc_pairwise_alfak_V2_vs_alfak_delta_summary.tsv`: pairwise alfak_V2-vs-alfak method deltas.",
  "- `tables/hpc_alfak_V2_all_nearfield_tradeoff.tsv`: alfak_V2 numerical/shape tradeoff and gap to best method.",
  "- `method_reports/`: per-method ALFA-K-style ground-truth summaries and figures.",
  "- `tables/hpc_ruggedness_effect_summary.tsv`: effect of ruggedness index, `1 / grf_lambda`, on each method/metric.",
  "",
  "## Figures",
  "- ![Numerical improvement](figures/fig01_numerical_improvement_vs_best_alfak.png)",
  "- ![Shape improvement](figures/fig02_shape_improvement_vs_best_alfak.png)",
  "- ![Tradeoff rank](figures/fig03_all_nearfield_tradeoff_rank.png)",
  "- ![Ruggedness trend](figures/fig04_ruggedness_lambda_trend_top_methods.png)",
  "- ![Gap to best](figures/fig05_alfak_V2_gap_to_best_all_nearfield.png)",
  "- ![Win rate heatmap](figures/fig06_win_rate_vs_best_alfak_heatmap.png)",
  paste0("- Supplementary Fig. 2-style ground-truth validation plots: ", if (exists("suppfig2_style_index") && nrow(suppfig2_style_index)) "`figures/suppfig2_style/`" else "_not generated_"),
  paste0("- Metric comparison violin/boxplots: ", if (length(metric_comparison_files)) "`figures/metric_comparison/`" else "_not generated_"),
  "",
  "## Direct answers",
  "",
  "### 1. Does alfak_V2 change fitness estimation accuracy by karyotype type?",
  "Use `d0` for directly informed nodes, `d1`/`d2` for one-hop/two-hop extrapolated nodes, `all_nearfield` for d1+d2, `all_eval` for d0+d1+d2, and `full_lscape` for the method's complete prediction table when available.",
  "The table below compares the best fixed alfak_V2 method with the best fixed alfak method after averaging over the selected depth/lambda grid. For lower-is-better metrics, positive `improvement` means error reduction. For higher-is-better metrics, positive `improvement` means higher correlation/accuracy.",
  "",
  as_markdown_table(best_compare_report),
  "",
  "The stricter oracle comparison below lets alfak choose its best setting separately in each matched condition. This is useful for seeing whether alfak_V2 still wins when alfak gets condition-specific tuning.",
  "",
  as_markdown_table(best_delta_oracle_report),
  "",
  "Interpretation: numerical fitness error is summarized by RMSE/MAE/relative RMSE. Shape is summarized by centered RMSE, ALFA-K-style Pearson/Spearman/rescaled R^2, edge-gradient RMSE, edge-gradient Spearman, and sign accuracy. d0 and full_lscape have no edge-gradient metric because no parent-to-child extrapolation edge is defined for those summaries.",
  "",
  "### 2. Which alfak_V2 extrapolation method is best overall?",
  paste0(
    "- Best all-nearfield balanced alfak_V2 method: `", best_a2_tradeoff$method_label,
    "`; numerical_rank=", fmt_num(best_a2_tradeoff$numerical_rank),
    ", shape_rank=", fmt_num(best_a2_tradeoff$shape_rank),
    ", balanced_rank=", fmt_num(best_a2_tradeoff$balanced_rank), "."
  ),
  paste0(
    "- Best all-nearfield RMSE alfak_V2 method: `", best_num$method_label,
    "`; RMSE=", fmt_num(best_num$rmse), ". Gap from balanced method to best RMSE = ",
    fmt_num(best_a2_tradeoff$rmse_gap_to_best), "."
  ),
  if (nrow(best_shape_edge)) paste0(
    "- Best all-nearfield edge-gradient RMSE alfak_V2 method: `", best_shape_edge$method_label,
    "`; edge_gradient_rmse=", fmt_num(best_shape_edge$edge_gradient_rmse),
    ". Gap from balanced method = ", fmt_num(best_a2_tradeoff$edge_gradient_rmse_gap_to_best), "."
  ) else "- Edge-gradient RMSE was unavailable for all-nearfield.",
  paste0(
    "- Best all-nearfield centered RMSE alfak_V2 method: `", best_shape_center$method_label,
    "`; centered_rmse=", fmt_num(best_shape_center$centered_rmse),
    ". Gap from balanced method = ", fmt_num(best_a2_tradeoff$centered_rmse_gap_to_best), "."
  ),
  "",
  "### 3. How does landscape ruggedness affect accuracy?",
  "The GRF generator defines larger `grf_lambda` as smoother landscapes, so this report uses `ruggedness_index = 1 / grf_lambda`. In `hpc_ruggedness_effect_summary.tsv`, `worse_when_more_rugged=TRUE` means error increases for lower-is-better metrics, or correlation/accuracy decreases for higher-is-better metrics, as ruggedness increases.",
  "",
  as_markdown_table(lambda_effect_report),
  "",
  "## Notes",
  "- `tabpfn_nearfield_feature_model` rows are included but should be interpreted as fallback evidence when dependency_status reports `tabpfn_unavailable_used_tree_fallback`.",
  "- The full benchmark combine step was not required for this report; this script uses `metrics_by_run.csv` directly to avoid constructing a very large all-method paired table."
)
writeLines(report_lines, file.path(input_dir, "hpc_selected_depth_lambda_summary_report.md"))

message("Wrote HPC summary under: ", input_dir)
