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
  ifelse(
    package == "alfak2",
    paste(package, input_mode, extrapolation_method, sep = ":"),
    paste0("alfakR:minobs", minobs, ":", nn_prior_slot)
  )
}

short_method <- function(x) {
  x <- sub("^alfak2:", "a2:", x)
  x <- sub("^alfakR:", "aR:", x)
  x <- gsub("edge_effect_interaction_path_ensemble", "edge_path", x)
  x <- gsub("graph_gaussian_baseline", "graph_gaussian", x)
  x <- gsub("truncated_nearfield_gmrf", "trunc_gmrf", x)
  x <- gsub("tabpfn_nearfield_feature_model", "tabpfn_fallback", x)
  x <- gsub("delta_tree_ensemble", "delta_tree", x)
  x <- gsub("empirical_censored_weighted_slot", "ecw", x)
  x <- gsub("empirical_censored", "ec", x)
  x
}

metric_meta <- data.table::data.table(
  metric = c(
    "rmse", "mae", "relative_rmse", "median_absolute_error", "q90_absolute_error",
    "centered_rmse", "affine_rmse", "pearson", "spearman",
    "edge_gradient_rmse", "edge_gradient_spearman", "sign_accuracy",
    "beneficial_sign_accuracy", "deleterious_sign_accuracy", "top_k_overlap_fraction"
  ),
  domain = c(
    rep("numerical", 5),
    rep("shape", 10)
  ),
  direction = c(
    rep("lower", 5),
    "lower", "lower", "higher", "higher",
    "lower", "higher", "higher", "higher", "higher", "higher"
  )
)

selected_metrics <- metric_meta$metric
key_cols <- c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat", "shell", "metric")

dt[, method_label := method_label(package, input_mode, extrapolation_method, minobs, NN_prior_slot)]
dt[, method_short := short_method(method_label)]
dt[, ruggedness_index := 1 / grf_lambda]
dt <- merge(dt, metric_meta, by = "metric", all.x = TRUE)
analysis_dt <- dt[
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
        domain, direction, alfak2_method = method_label, best_alfakR_method,
        alfak2_value = value, best_alfakR_value, raw_delta, improvement,
        relative_improvement, improved)
  ],
  file.path(tables_dir, "hpc_alfak2_delta_vs_best_alfakR_by_condition.tsv")
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
write_tsv(delta_best_summary, file.path(tables_dir, "hpc_alfak2_delta_vs_best_alfakR_summary.tsv"))

metric_values <- sort(unique(analysis_dt$metric))
message("Building pairwise alfak2-vs-alfakR delta summaries across ", length(metric_values), " metric(s) with up to ", min(threads, length(metric_values)), " worker(s).")
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
write_tsv(pairwise_delta_summary, file.path(tables_dir, "hpc_pairwise_alfak2_vs_alfakR_delta_summary.tsv"))

rank_metrics <- c("rmse", "mae", "relative_rmse", "centered_rmse", "pearson", "spearman",
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

trade_metrics <- c("rmse", "mae", "relative_rmse", "centered_rmse", "edge_gradient_rmse",
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
for (col in intersect(c("spearman", "edge_gradient_spearman", "sign_accuracy"), names(alfak2_tradeoff))) {
  alfak2_tradeoff[, paste0(col, "_gap_to_best") := gap_higher(get(col))]
}
data.table::setorder(alfak2_tradeoff, balanced_rank)
write_tsv(alfak2_tradeoff, file.path(tables_dir, "hpc_alfak2_all_nearfield_tradeoff.tsv"))

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
    lambda_trend = lambda_trend,
    ruggedness_effect = rugged_effect
  ),
  file.path(results_dir, "hpc_selected_depth_lambda_analysis.rds")
)

ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))

plot_save <- function(plot, filename, width = 12, height = 7) {
  ggplot2::ggsave(file.path(figures_dir, filename), plot = plot, width = width, height = height, dpi = 180)
}

numeric_plot_dt <- delta_best_summary[
  metric %in% c("rmse", "mae", "relative_rmse") &
    shell %in% c("d0", "d1", "d2", "all_nearfield")
]
numeric_plot_dt[, facet := paste(shell, metric, sep = " / ")]
p1 <- ggplot2::ggplot(
  numeric_plot_dt,
  ggplot2::aes(x = stats::reorder(alfak2_short, mean_improvement), y = mean_improvement, fill = alfak2_short)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(metric ~ shell, scales = "free_x") +
  ggplot2::labs(
    title = "alfak2 numerical improvement vs best alfakR",
    subtitle = "Positive means alfak2 has lower error than the best alfakR setting in the same condition.",
    x = NULL,
    y = "Mean improvement"
  )
plot_save(p1, "fig01_numerical_improvement_vs_best_alfakR.png", width = 13, height = 8)

shape_plot_dt <- delta_best_summary[
  metric %in% c("centered_rmse", "edge_gradient_rmse", "spearman", "edge_gradient_spearman", "sign_accuracy") &
    shell %in% c("d0", "d1", "d2", "all_nearfield")
]
p2 <- ggplot2::ggplot(
  shape_plot_dt,
  ggplot2::aes(x = stats::reorder(alfak2_short, mean_improvement), y = mean_improvement, fill = alfak2_short)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(metric ~ shell, scales = "free_x") +
  ggplot2::labs(
    title = "alfak2 shape improvement vs best alfakR",
    subtitle = "Positive is better after respecting each metric's direction.",
    x = NULL,
    y = "Mean improvement"
  )
plot_save(p2, "fig02_shape_improvement_vs_best_alfakR.png", width = 13, height = 10)

trade_plot <- overall_rank_wide[shell == "all_nearfield"]
p3 <- ggplot2::ggplot(
  trade_plot,
  ggplot2::aes(x = numerical, y = shape, color = package, label = method_short)
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
    metric %in% c("rmse", "centered_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy") &
    method_label %in% top_methods
]
p4 <- ggplot2::ggplot(
  trend_plot_dt,
  ggplot2::aes(x = grf_lambda, y = mean, color = method_short, group = method_short)
) +
  ggplot2::geom_line(linewidth = 0.55) +
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
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ gap_metric, scales = "free_x") +
  ggplot2::labs(
    title = "alfak2 gap to best alfak2 all-nearfield method",
    subtitle = "For error metrics lower gaps are better; for correlation/accuracy metrics the gap is best minus method.",
    x = NULL,
    y = "Gap to best"
  )
plot_save(p5, "fig05_alfak2_gap_to_best_all_nearfield.png", width = 13, height = 8)

heat_dt <- delta_best_summary[
  shell %in% c("d0", "d1", "d2", "all_nearfield") &
    metric %in% c("rmse", "mae", "centered_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy")
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
    title = "Win rate vs best alfakR by condition",
    subtitle = "Each cell is the fraction of matched conditions where alfak2 beats the best alfakR setting.",
    x = NULL,
    y = NULL,
    fill = "Win rate"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
plot_save(p6, "fig06_win_rate_vs_best_alfakR_heatmap.png", width = 14, height = 10)

best_overall <- overall_best_by_package[
  shell %in% c("d0", "d1", "d2", "all_nearfield") &
    metric %in% c("rmse", "mae", "centered_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy")
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
best_compare[, result := ifelse(improvement_vs_fixed_best > 0, "alfak2 better", "alfak2 worse")]
best_compare_report <- best_compare[
  order(shell, metric),
  .(
    shell,
    metric,
    direction,
    alfak2_method = method_short_alfak2,
    alfak2_mean = fmt_num(mean_alfak2),
    alfakR_method = method_short_alfakR,
    alfakR_mean = fmt_num(mean_alfakR),
    improvement = fmt_num(improvement_vs_fixed_best),
    relative_improvement = fmt_pct(relative_improvement_vs_fixed_best),
    result
  )
]

best_delta_oracle <- delta_best_summary[
  shell %in% c("d0", "d1", "d2", "all_nearfield") &
    metric %in% c("rmse", "centered_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy"),
  .SD[which.max(mean_improvement)],
  by = .(shell, metric)
]
best_delta_oracle_report <- best_delta_oracle[
  order(shell, metric),
  .(
    shell,
    metric,
    alfak2_method = alfak2_short,
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
    metric %in% c("rmse", "centered_rmse", "edge_gradient_rmse", "spearman", "sign_accuracy")
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
    package,
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
  "",
  "## Key outputs",
  "- `tables/hpc_method_metric_summary.tsv`: method means by sample depth, lambda, shell, and metric.",
  "- `tables/hpc_landscape_fit_repeat_summary.tsv`: within-landscape summaries across fit repeats.",
  "- `tables/hpc_landscape_block_method_metric_summary.tsv`: landscape-block summaries after first aggregating fit repeats within each landscape.",
  "- `tables/hpc_alfak2_delta_vs_best_alfakR_summary.tsv`: alfak2 improvement against the best alfakR setting in matched conditions.",
  "- `tables/hpc_pairwise_alfak2_vs_alfakR_delta_summary.tsv`: pairwise alfak2-vs-alfakR method deltas.",
  "- `tables/hpc_alfak2_all_nearfield_tradeoff.tsv`: alfak2 numerical/shape tradeoff and gap to best method.",
  "- `tables/hpc_ruggedness_effect_summary.tsv`: effect of ruggedness index, `1 / grf_lambda`, on each method/metric.",
  "",
  "## Figures",
  "- ![Numerical improvement](figures/fig01_numerical_improvement_vs_best_alfakR.png)",
  "- ![Shape improvement](figures/fig02_shape_improvement_vs_best_alfakR.png)",
  "- ![Tradeoff rank](figures/fig03_all_nearfield_tradeoff_rank.png)",
  "- ![Ruggedness trend](figures/fig04_ruggedness_lambda_trend_top_methods.png)",
  "- ![Gap to best](figures/fig05_alfak2_gap_to_best_all_nearfield.png)",
  "- ![Win rate heatmap](figures/fig06_win_rate_vs_best_alfakR_heatmap.png)",
  "",
  "## Direct answers",
  "",
  "### 1. Does alfak2 change fitness estimation accuracy by karyotype type?",
  "Use `d0` for directly informed nodes, `d1`/`d2` for one-hop/two-hop extrapolated nodes, and `all_nearfield` for d1+d2.",
  "The table below compares the best fixed alfak2 method with the best fixed alfakR method after averaging over the selected depth/lambda grid. For lower-is-better metrics, positive `improvement` means error reduction. For higher-is-better metrics, positive `improvement` means higher correlation/accuracy.",
  "",
  as_markdown_table(best_compare_report),
  "",
  "The stricter oracle comparison below lets alfakR choose its best setting separately in each matched condition. This is useful for seeing whether alfak2 still wins when alfakR gets condition-specific tuning.",
  "",
  as_markdown_table(best_delta_oracle_report),
  "",
  "Interpretation: numerical fitness error is summarized by RMSE/MAE/relative RMSE. Shape is summarized by centered RMSE, Spearman/Pearson-style ordering, edge-gradient RMSE, edge-gradient Spearman, and sign accuracy. d0 has no edge-gradient metric because no parent-to-child extrapolation edge is defined for d0.",
  "",
  "### 2. Which alfak2 extrapolation method is best overall?",
  paste0(
    "- Best all-nearfield balanced alfak2 method: `", best_a2_tradeoff$method_label,
    "`; numerical_rank=", fmt_num(best_a2_tradeoff$numerical_rank),
    ", shape_rank=", fmt_num(best_a2_tradeoff$shape_rank),
    ", balanced_rank=", fmt_num(best_a2_tradeoff$balanced_rank), "."
  ),
  paste0(
    "- Best all-nearfield RMSE alfak2 method: `", best_num$method_label,
    "`; RMSE=", fmt_num(best_num$rmse), ". Gap from balanced method to best RMSE = ",
    fmt_num(best_a2_tradeoff$rmse_gap_to_best), "."
  ),
  if (nrow(best_shape_edge)) paste0(
    "- Best all-nearfield edge-gradient RMSE alfak2 method: `", best_shape_edge$method_label,
    "`; edge_gradient_rmse=", fmt_num(best_shape_edge$edge_gradient_rmse),
    ". Gap from balanced method = ", fmt_num(best_a2_tradeoff$edge_gradient_rmse_gap_to_best), "."
  ) else "- Edge-gradient RMSE was unavailable for all-nearfield.",
  paste0(
    "- Best all-nearfield centered RMSE alfak2 method: `", best_shape_center$method_label,
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
