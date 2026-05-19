#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(htmltools)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: render_grf_downsampled_accuracy_report.R <result_dir>", call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
if (basename(run_dir) == "benchmark") {
  benchmark_dir <- run_dir
  run_dir <- dirname(run_dir)
} else if (dir.exists(file.path(run_dir, "benchmark"))) {
  benchmark_dir <- file.path(run_dir, "benchmark")
} else {
  benchmark_dir <- run_dir
}
tables_dir <- file.path(benchmark_dir, "tables")
if (!dir.exists(tables_dir)) {
  stop("Missing tables directory: ", tables_dir, call. = FALSE)
}

cmd <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE) else NA_character_
script_dir <- if (is.na(script_path)) getwd() else dirname(script_path)

required_extended <- file.path(tables_dir, c(
  "extended_fitness_accuracy_input_coverage.tsv",
  "extended_fitness_accuracy_method_condition_metrics.tsv",
  "extended_fitness_accuracy_method_summary.tsv",
  "extended_fitness_accuracy_pair_condition_metrics_by_alfakR_scope.tsv",
  "extended_fitness_accuracy_pair_summary_by_alfakR_scope.tsv"
))

ensure_extended_tables <- function() {
  if (all(file.exists(required_extended))) return(invisible(TRUE))
  evaluator <- file.path(script_dir, "evaluate_grf_fitness_accuracy_by_alfakR_scope.R")
  if (!file.exists(evaluator)) {
    stop("Missing extended accuracy evaluator: ", evaluator, call. = FALSE)
  }
  message("Missing extended accuracy tables; running evaluator first.")
  status <- system2("Rscript", c(evaluator, run_dir))
  if (!identical(status, 0L) || !all(file.exists(required_extended))) {
    stop("Failed to generate extended accuracy tables for: ", run_dir, call. = FALSE)
  }
  invisible(TRUE)
}

ensure_extended_tables()

report_dir <- file.path(run_dir, "report_assets", "grf_downsampled_accuracy_report")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_html <- file.path(run_dir, "grf_downsampled_accuracy_report.html")

fread_required <- function(path) {
  if (!file.exists(path)) stop("Missing required table: ", path, call. = FALSE)
  data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE, data.table = TRUE)
}

input_cov <- fread_required(file.path(tables_dir, "extended_fitness_accuracy_input_coverage.tsv"))
method_cond <- fread_required(file.path(tables_dir, "extended_fitness_accuracy_method_condition_metrics.tsv"))
method_summary <- fread_required(file.path(tables_dir, "extended_fitness_accuracy_method_summary.tsv"))
pair_cond <- fread_required(file.path(tables_dir, "extended_fitness_accuracy_pair_condition_metrics_by_alfakR_scope.tsv"))
pair_summary <- fread_required(file.path(tables_dir, "extended_fitness_accuracy_pair_summary_by_alfakR_scope.tsv"))
fit_results <- if (file.exists(file.path(tables_dir, "fit_results.tsv"))) {
  data.table::fread(file.path(tables_dir, "fit_results.tsv"), sep = "\t", header = TRUE, showProgress = FALSE, data.table = TRUE)
} else {
  data.table()
}

cfg_path <- file.path(benchmark_dir, "benchmark_config.rds")
cfg <- if (file.exists(cfg_path)) readRDS(cfg_path) else list()

scopes <- c("fq", "NN", "other", "whole")
scope_labels <- c(fq = "fq", NN = "NN", other = "other", whole = "whole")
standard_methods <- c("alfak2 full", "alfak2 minobs-matched", "alfakR empirical")
alfak2_methods <- c("alfak2 full", "alfak2 minobs-matched")

method_cond <- method_cond[method_label %in% standard_methods & support_scope %in% scopes]
method_summary <- method_summary[method_label %in% standard_methods & support_scope %in% scopes]
pair_cond <- pair_cond[alfakR_label == "alfakR empirical" & support_scope %in% scopes & alfak2_label %in% alfak2_methods]
pair_summary <- pair_summary[alfakR_label == "alfakR empirical" & support_scope %in% scopes & alfak2_label %in% alfak2_methods]

method_cond[, `:=`(
  support_scope = factor(as.character(support_scope), levels = scopes),
  method_label = factor(as.character(method_label), levels = standard_methods),
  minobs_f = factor(minobs, levels = sort(unique(minobs)))
)]
pair_cond[, `:=`(
  support_scope = factor(as.character(support_scope), levels = scopes),
  alfak2_label = factor(as.character(alfak2_label), levels = alfak2_methods),
  minobs_f = factor(minobs, levels = sort(unique(minobs)))
)]
input_cov[, minobs_f := factor(minobs, levels = sort(unique(minobs)))]

q_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(list(n = 0L, min = NA_real_, q25 = NA_real_, median = NA_real_, mean = NA_real_, q75 = NA_real_, max = NA_real_))
  }
  list(
    n = length(x),
    min = min(x),
    q25 = as.numeric(quantile(x, 0.25, names = FALSE)),
    median = median(x),
    mean = mean(x),
    q75 = as.numeric(quantile(x, 0.75, names = FALSE)),
    max = max(x)
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "fg", flag = "#"))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 1e-4, "p<1e-4",
      ifelse(p < 0.001, paste0("p=", formatC(p, digits = 2, format = "e")),
        paste0("p=", formatC(p, digits = 3, format = "fg"))
      )
    )
  )
}

p_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns")
      )
    )
  )
}

wilcox_one <- function(x, alternative) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || length(unique(x)) < 2L) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

metric_tests <- data.table::rbindlist(list(
  pair_cond[, .(
    metric = "delta_mae",
    metric_label = "MAE delta",
    direction = "less_than_zero",
    n = sum(is.finite(delta_mae)),
    median_delta = median(delta_mae, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta_mae, 0.25, na.rm = TRUE, names = FALSE)),
    q75_delta = as.numeric(quantile(delta_mae, 0.75, na.rm = TRUE, names = FALSE)),
    p_value = wilcox_one(delta_mae, "less")
  ), by = .(support_scope, minobs, alfak2_label)],
  pair_cond[, .(
    metric = "delta_rmse",
    metric_label = "RMSE delta",
    direction = "less_than_zero",
    n = sum(is.finite(delta_rmse)),
    median_delta = median(delta_rmse, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta_rmse, 0.25, na.rm = TRUE, names = FALSE)),
    q75_delta = as.numeric(quantile(delta_rmse, 0.75, na.rm = TRUE, names = FALSE)),
    p_value = wilcox_one(delta_rmse, "less")
  ), by = .(support_scope, minobs, alfak2_label)],
  pair_cond[, .(
    metric = "delta_centered_rmse",
    metric_label = "centered RMSE delta",
    direction = "less_than_zero",
    n = sum(is.finite(delta_centered_rmse)),
    median_delta = median(delta_centered_rmse, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta_centered_rmse, 0.25, na.rm = TRUE, names = FALSE)),
    q75_delta = as.numeric(quantile(delta_centered_rmse, 0.75, na.rm = TRUE, names = FALSE)),
    p_value = wilcox_one(delta_centered_rmse, "less")
  ), by = .(support_scope, minobs, alfak2_label)],
  pair_cond[, .(
    metric = "delta_pearson",
    metric_label = "Pearson delta",
    direction = "greater_than_zero",
    n = sum(is.finite(delta_pearson)),
    median_delta = median(delta_pearson, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta_pearson, 0.25, na.rm = TRUE, names = FALSE)),
    q75_delta = as.numeric(quantile(delta_pearson, 0.75, na.rm = TRUE, names = FALSE)),
    p_value = wilcox_one(delta_pearson, "greater")
  ), by = .(support_scope, minobs, alfak2_label)],
  pair_cond[, .(
    metric = "delta_spearman",
    metric_label = "Spearman delta",
    direction = "greater_than_zero",
    n = sum(is.finite(delta_spearman)),
    median_delta = median(delta_spearman, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta_spearman, 0.25, na.rm = TRUE, names = FALSE)),
    q75_delta = as.numeric(quantile(delta_spearman, 0.75, na.rm = TRUE, names = FALSE)),
    p_value = wilcox_one(delta_spearman, "greater")
  ), by = .(support_scope, minobs, alfak2_label)]
), fill = TRUE)
metric_tests[, p_adj_bh := p.adjust(p_value, method = "BH"), by = metric]
metric_tests[, p_label := paste0(fmt_p(p_value), " ", p_stars(p_value))]
metric_tests[, support_scope := factor(as.character(support_scope), levels = scopes)]
metric_tests[, alfak2_label := factor(as.character(alfak2_label), levels = alfak2_methods)]

input_summary <- input_cov[, c(
  q_summary(raw_input_rows),
  setNames(q_summary(input_rows_after_drop), paste0("after_drop_", names(q_summary(input_rows_after_drop)))),
  setNames(q_summary(input_rows_minobs), paste0("minobs_", names(q_summary(input_rows_minobs))))
), by = .(minobs)]
data.table::setorder(input_summary, minobs)

fit_status <- data.table()
if (nrow(fit_results) && all(c("engine", "status") %in% names(fit_results))) {
  fit_status <- fit_results[engine %in% c("alfak2", "alfakR"), .N, by = .(engine, status)]
  fit_status[, total := sum(N), by = engine]
  fit_status[, rate := N / total]
  data.table::setorder(fit_status, engine, status)
}

own_report <- method_summary[
  support_scope %in% scopes & method_label %in% standard_methods,
  .(support_scope, minobs, method_label, n_conditions, n_estimated_median,
    mae_median, rmse_median, centered_rmse_median,
    pearson_median, spearman_median, centered_r2_median,
    calibration_slope_median)
]
data.table::setorder(own_report, support_scope, minobs, method_label)

strict_report <- pair_cond[, .(
  n_conditions = .N,
  n_common_median = median(n_common, na.rm = TRUE),
  alfak2_mae_median = median(alfak2_mae, na.rm = TRUE),
  alfakR_mae_median = median(alfakR_mae, na.rm = TRUE),
  delta_mae_median = median(delta_mae, na.rm = TRUE),
  alfak2_rmse_median = median(alfak2_rmse, na.rm = TRUE),
  alfakR_rmse_median = median(alfakR_rmse, na.rm = TRUE),
  delta_rmse_median = median(delta_rmse, na.rm = TRUE),
  alfak2_centered_rmse_median = median(alfak2_centered_rmse, na.rm = TRUE),
  alfakR_centered_rmse_median = median(alfakR_centered_rmse, na.rm = TRUE),
  delta_centered_rmse_median = median(delta_centered_rmse, na.rm = TRUE),
  alfak2_pearson_median = median(alfak2_pearson, na.rm = TRUE),
  alfakR_pearson_median = median(alfakR_pearson, na.rm = TRUE),
  delta_pearson_median = median(delta_pearson, na.rm = TRUE),
  alfak2_spearman_median = median(alfak2_spearman, na.rm = TRUE),
  alfakR_spearman_median = median(alfakR_spearman, na.rm = TRUE),
  delta_spearman_median = median(delta_spearman, na.rm = TRUE),
  mae_condition_win_rate = mean(delta_mae < 0, na.rm = TRUE),
  centered_rmse_condition_win_rate = mean(delta_centered_rmse < 0, na.rm = TRUE),
  pearson_condition_win_rate = mean(delta_pearson > 0, na.rm = TRUE),
  spearman_condition_win_rate = mean(delta_spearman > 0, na.rm = TRUE)
), by = .(support_scope, minobs, alfak2_label)]
data.table::setorder(strict_report, support_scope, minobs, alfak2_label)

shape_report <- strict_report[, .(
  support_scope, minobs, alfak2_label, n_common_median,
  alfak2_centered_rmse_median, alfakR_centered_rmse_median,
  delta_centered_rmse_median,
  alfak2_pearson_median, alfakR_pearson_median, delta_pearson_median,
  alfak2_spearman_median, alfakR_spearman_median, delta_spearman_median,
  centered_rmse_condition_win_rate,
  pearson_condition_win_rate,
  spearman_condition_win_rate
)]

summary_files <- c(
  input = file.path(tables_dir, "report_input_minobs_summary.tsv"),
  own = file.path(tables_dir, "report_own_method_summary.tsv"),
  strict = file.path(tables_dir, "report_strict_common_empirical_summary.tsv"),
  shape = file.path(tables_dir, "report_strict_common_shape_summary.tsv"),
  tests = file.path(tables_dir, "report_strict_common_wilcoxon_tests.tsv")
)
data.table::fwrite(input_summary, summary_files[["input"]], sep = "\t")
data.table::fwrite(own_report, summary_files[["own"]], sep = "\t")
data.table::fwrite(strict_report, summary_files[["strict"]], sep = "\t")
data.table::fwrite(shape_report, summary_files[["shape"]], sep = "\t")
data.table::fwrite(metric_tests, summary_files[["tests"]], sep = "\t")

theme_report <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#eef2f7", colour = "#cbd5e1"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}

save_plot <- function(plot, filename, width = 13, height = 8) {
  path <- file.path(report_dir, filename)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 180, bg = "white")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

plot_paths <- list()

plot_paths$input <- save_plot(
  ggplot(input_cov, aes(x = minobs_f, y = input_rows_minobs)) +
    geom_violin(fill = "#8ecae6", colour = "#2563eb", alpha = 0.55, trim = FALSE) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "#1e3a8a") +
    geom_jitter(width = 0.08, height = 0, alpha = 0.25, size = 0.8, colour = "#0f172a") +
    labs(
      title = "Figure 1. Input karyotype counts after minobs filtering",
      x = "minobs",
      y = "retained input karyotypes"
    ) +
    theme_report(),
  "fig01_input_karyotypes_by_minobs.png", width = 8.5, height = 5.5
)

own_accuracy_long <- data.table::melt(
  method_cond,
  id.vars = c("support_scope", "minobs_f", "method_label"),
  measure.vars = c("mae", "rmse"),
  variable.name = "metric",
  value.name = "value"
)
own_accuracy_long[, metric := factor(metric, levels = c("mae", "rmse"), labels = c("MAE", "RMSE"))]
plot_paths$own_accuracy <- save_plot(
  ggplot(own_accuracy_long, aes(x = minobs_f, y = value, fill = method_label)) +
    geom_violin(position = position_dodge(width = 0.82), alpha = 0.55, trim = FALSE, scale = "width") +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.shape = NA, alpha = 0.9) +
    facet_grid(metric ~ support_scope, scales = "free_y") +
    labs(
      title = "Figure 2. Own-landscape numerical accuracy by method",
      x = "minobs",
      y = "error vs reference GRF",
      fill = "method"
    ) +
    scale_fill_manual(values = c("alfak2 full" = "#1f77b4", "alfak2 minobs-matched" = "#2ca02c", "alfakR empirical" = "#d62728")) +
    theme_report(),
  "fig02_own_landscape_accuracy_violin_box.png"
)

accuracy_delta_long <- data.table::melt(
  pair_cond,
  id.vars = c("support_scope", "minobs", "minobs_f", "alfak2_label"),
  measure.vars = c("delta_mae", "delta_rmse"),
  variable.name = "metric",
  value.name = "value"
)
accuracy_delta_long[, metric_label := factor(metric, levels = c("delta_mae", "delta_rmse"), labels = c("MAE delta", "RMSE delta"))]
acc_test <- metric_tests[metric %in% c("delta_mae", "delta_rmse")]
acc_y <- accuracy_delta_long[, .(y = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE))), by = .(support_scope, minobs, alfak2_label, metric)]
acc_test <- merge(acc_test, acc_y, by = c("support_scope", "minobs", "alfak2_label", "metric"), all.x = TRUE)
acc_test[, `:=`(minobs_f = factor(minobs, levels = levels(input_cov$minobs_f)),
                metric_label = factor(metric, levels = c("delta_mae", "delta_rmse"), labels = c("MAE delta", "RMSE delta")))]
plot_paths$accuracy_delta <- save_plot(
  ggplot(accuracy_delta_long, aes(x = minobs_f, y = value, fill = alfak2_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#475569") +
    geom_violin(position = position_dodge(width = 0.82), alpha = 0.55, trim = FALSE, scale = "width") +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.shape = NA, alpha = 0.9) +
    geom_text(
      data = acc_test,
      aes(x = minobs_f, y = y, label = p_label, group = alfak2_label),
      position = position_dodge(width = 0.82),
      inherit.aes = FALSE,
      size = 2.7,
      angle = 0,
      vjust = 0
    ) +
    facet_grid(metric_label ~ support_scope, scales = "free_y") +
    labs(
      title = "Figure 3. Strict-common numerical delta vs alfakR empirical",
      subtitle = "Delta = alfak2 - alfakR empirical; negative values favor alfak2. Labels are one-sided Wilcoxon p-values for delta < 0.",
      x = "minobs",
      y = "delta",
      fill = "alfak2 policy"
    ) +
    scale_fill_manual(values = c("alfak2 full" = "#1f77b4", "alfak2 minobs-matched" = "#2ca02c")) +
    theme_report(),
  "fig03_strict_common_accuracy_delta_violin_box.png"
)

shape_abs <- data.table::rbindlist(list(
  pair_cond[, .(support_scope, minobs_f, method_label = as.character(alfak2_label),
                centered_rmse = alfak2_centered_rmse,
                pearson = alfak2_pearson,
                spearman = alfak2_spearman)],
  pair_cond[, .(support_scope, minobs_f, method_label = "alfakR empirical",
                centered_rmse = alfakR_centered_rmse,
                pearson = alfakR_pearson,
                spearman = alfakR_spearman)]
), fill = TRUE)
shape_abs[, method_label := factor(method_label, levels = standard_methods)]
shape_abs_long <- data.table::melt(
  shape_abs,
  id.vars = c("support_scope", "minobs_f", "method_label"),
  measure.vars = c("centered_rmse", "pearson", "spearman"),
  variable.name = "metric",
  value.name = "value"
)
shape_abs_long[, metric := factor(metric, levels = c("centered_rmse", "pearson", "spearman"),
                                  labels = c("centered RMSE", "Pearson", "Spearman"))]
plot_paths$shape_abs <- save_plot(
  ggplot(shape_abs_long, aes(x = minobs_f, y = value, fill = method_label)) +
    geom_violin(position = position_dodge(width = 0.82), alpha = 0.55, trim = FALSE, scale = "width") +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.shape = NA, alpha = 0.9) +
    facet_grid(metric ~ support_scope, scales = "free_y") +
    labs(
      title = "Figure 4. Strict-common shape metrics by method",
      x = "minobs",
      y = "shape metric",
      fill = "method"
    ) +
    scale_fill_manual(values = c("alfak2 full" = "#1f77b4", "alfak2 minobs-matched" = "#2ca02c", "alfakR empirical" = "#d62728")) +
    theme_report(),
  "fig04_strict_common_shape_metrics_violin_box.png"
)

shape_delta_long <- data.table::melt(
  pair_cond,
  id.vars = c("support_scope", "minobs", "minobs_f", "alfak2_label"),
  measure.vars = c("delta_centered_rmse", "delta_pearson", "delta_spearman"),
  variable.name = "metric",
  value.name = "value"
)
shape_delta_long[, metric_label := factor(
  metric,
  levels = c("delta_centered_rmse", "delta_pearson", "delta_spearman"),
  labels = c("centered RMSE delta", "Pearson delta", "Spearman delta")
)]
shape_test <- metric_tests[metric %in% c("delta_centered_rmse", "delta_pearson", "delta_spearman")]
shape_y <- shape_delta_long[, .(y = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE))), by = .(support_scope, minobs, alfak2_label, metric)]
shape_test <- merge(shape_test, shape_y, by = c("support_scope", "minobs", "alfak2_label", "metric"), all.x = TRUE)
shape_test[, `:=`(
  minobs_f = factor(minobs, levels = levels(input_cov$minobs_f)),
  metric_label = factor(metric, levels = c("delta_centered_rmse", "delta_pearson", "delta_spearman"),
                        labels = c("centered RMSE delta", "Pearson delta", "Spearman delta"))
)]
plot_paths$shape_delta <- save_plot(
  ggplot(shape_delta_long, aes(x = minobs_f, y = value, fill = alfak2_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#475569") +
    geom_violin(position = position_dodge(width = 0.82), alpha = 0.55, trim = FALSE, scale = "width") +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.shape = NA, alpha = 0.9) +
    geom_text(
      data = shape_test,
      aes(x = minobs_f, y = y, label = p_label, group = alfak2_label),
      position = position_dodge(width = 0.82),
      inherit.aes = FALSE,
      size = 2.6,
      vjust = 0
    ) +
    facet_grid(metric_label ~ support_scope, scales = "free_y") +
    labs(
      title = "Figure 5. Strict-common shape delta vs alfakR empirical",
      subtitle = "Centered RMSE delta < 0 favors alfak2; Pearson/Spearman delta > 0 favors alfak2. Labels are one-sided Wilcoxon p-values.",
      x = "minobs",
      y = "delta",
      fill = "alfak2 policy"
    ) +
    scale_fill_manual(values = c("alfak2 full" = "#1f77b4", "alfak2 minobs-matched" = "#2ca02c")) +
    theme_report(),
  "fig05_strict_common_shape_delta_violin_box.png"
)

write_html_table <- function(x, digits = 3) {
  y <- as.data.frame(x)
  for (nm in names(y)) {
    if (is.numeric(y[[nm]])) y[[nm]] <- fmt_num(y[[nm]], digits)
  }
  paste(capture.output(print(knitr::kable(y, format = "html", escape = TRUE, table.attr = "class=\"report-table\""))), collapse = "\n")
}

rel_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
  file.path("report_assets", "grf_downsampled_accuracy_report", basename(path))
}

img_tag <- function(path, alt) {
  paste0("<img class=\"figure-img\" src=\"", htmlEscape(rel_path(path)), "\" alt=\"", htmlEscape(alt), "\">")
}

table_link <- function(path) {
  rel <- file.path("benchmark", "tables", basename(path))
  paste0("<a href=\"", htmlEscape(rel), "\">", htmlEscape(basename(path)), "</a>")
}

sample_depth <- cfg$sample_depth %||% NA
n_sim <- cfg$n_sim %||% NA
lambdas <- if (!is.null(cfg$lambdas)) paste(cfg$lambdas, collapse = ", ") else ""
time_gaps <- if (!is.null(cfg$time_gaps)) paste(cfg$time_gaps, collapse = ", ") else ""
minobs_cfg <- if (!is.null(cfg$minobs)) paste(cfg$minobs, collapse = ", ") else paste(sort(unique(input_cov$minobs)), collapse = ", ")

input_interp <- {
  first <- input_summary[which.min(minobs)]
  last <- input_summary[which.max(minobs)]
  paste0(
    "As minobs increases from ", first$minobs, " to ", last$minobs,
    ", median retained input karyotypes changes from ", fmt_num(first$minobs_median, 3),
    " to ", fmt_num(last$minobs_median, 3),
    ". This is the input-size context for interpreting fq/NN/other accuracy."
  )
}

accuracy_interp <- {
  whole <- strict_report[as.character(support_scope) == "whole" & alfak2_label == "alfak2 minobs-matched"]
  if (nrow(whole)) {
    paste0(
      "On strict-common whole-scope karyotypes, alfak2 minobs-matched median MAE is ",
      paste(paste0("minobs ", whole$minobs, ": ", fmt_num(whole$alfak2_mae_median, 3),
                   " vs alfakR ", fmt_num(whole$alfakR_mae_median, 3)), collapse = "; "),
      ". Negative deltas in Figure 3 and Table 4 indicate numerical improvement over alfakR empirical."
    )
  } else {
    "Strict-common whole-scope rows were not available."
  }
}

shape_interp <- {
  other <- shape_report[as.character(support_scope) == "other" & alfak2_label == "alfak2 minobs-matched"]
  if (nrow(other)) {
    paste0(
      "For far-field other nodes, alfak2 minobs-matched centered RMSE is generally lower than alfakR empirical, while Pearson/Spearman deltas are often smaller or negative. ",
      "This separates amplitude/offset-corrected shape error from rank-order preservation."
    )
  } else {
    "Other-scope strict-common rows were not available."
  }
}

css <- "
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 28px; color: #111827; line-height: 1.45; }
h1, h2, h3 { color: #0f172a; }
h1 { border-bottom: 2px solid #cbd5e1; padding-bottom: 8px; }
.meta, .note { background: #f8fafc; border: 1px solid #e2e8f0; padding: 12px 14px; border-radius: 6px; }
.figure-block, .table-block { margin: 26px 0 34px 0; }
.figure-img { max-width: 100%; border: 1px solid #d1d5db; border-radius: 6px; }
.caption { font-weight: 600; margin-top: 8px; }
.interpretation { margin-top: 6px; color: #334155; }
table.report-table { border-collapse: collapse; font-size: 12px; width: 100%; margin-top: 8px; }
table.report-table th, table.report-table td { border: 1px solid #d1d5db; padding: 4px 6px; text-align: right; }
table.report-table th:first-child, table.report-table td:first-child { text-align: left; }
table.report-table th { background: #eef2f7; color: #0f172a; }
.small { font-size: 12px; color: #475569; }
"

html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>GRF downsampled accuracy report</title><style>", css, "</style></head><body>",
  "<h1>GRF Downsampled Accuracy Report</h1>",
  "<div class=\"meta\"><strong>Result directory:</strong> ", htmlEscape(run_dir), "<br>",
  "<strong>Generated:</strong> ", htmlEscape(as.character(Sys.time())), "<br>",
  "<strong>sample_depth:</strong> ", htmlEscape(as.character(sample_depth)), " &nbsp; ",
  "<strong>n_sim:</strong> ", htmlEscape(as.character(n_sim)), " &nbsp; ",
  "<strong>lambda:</strong> ", htmlEscape(lambdas), " &nbsp; ",
  "<strong>time_gap:</strong> ", htmlEscape(time_gaps), " &nbsp; ",
  "<strong>minobs:</strong> ", htmlEscape(minobs_cfg), "</div>",
  "<h2>Analysis Scope</h2>",
  "<p>This report compares <strong>alfak2 full</strong>, <strong>alfak2 minobs-matched</strong>, and <strong>alfakR empirical</strong>. Numerical accuracy is summarized by MAE/RMSE against the reference GRF. Shape is summarized by centered RMSE, Pearson, Spearman, and centered R2. Strict-common comparisons use only karyotypes estimated by both alfak2 and alfakR empirical within the same alfakR-defined fq/NN/other scope.</p>",
  "<p class=\"note\"><strong>Statistical tests:</strong> paired condition-level deltas are tested by one-sided Wilcoxon signed-rank tests. For MAE/RMSE/centered RMSE, delta = alfak2 - alfakR empirical and the improvement alternative is delta &lt; 0. For Pearson/Spearman, the improvement alternative is delta &gt; 0.</p>",
  "<h2>Figures</h2>",
  "<div class=\"figure-block\"><h3>Figure 1. Input karyotype counts after minobs filtering</h3>", img_tag(plot_paths$input, "Input karyotype count violin plot"),
  "<p class=\"caption\">Content description: distribution of retained input karyotypes across benchmark conditions for each minobs threshold.</p>",
  "<p class=\"interpretation\">Interpretation: ", htmlEscape(input_interp), "</p></div>",
  "<div class=\"figure-block\"><h3>Figure 2. Own-landscape numerical accuracy</h3>", img_tag(plot_paths$own_accuracy, "Own landscape MAE RMSE violin box plot"),
  "<p class=\"caption\">Content description: MAE and RMSE distributions for each method using each method's own estimated landscape and support scope.</p>",
  "<p class=\"interpretation\">Interpretation: this view shows the method-level error distribution before forcing both methods onto identical karyotype sets. It highlights coverage and support-scope differences between alfak2 full, alfak2 minobs-matched, and alfakR empirical.</p></div>",
  "<div class=\"figure-block\"><h3>Figure 3. Strict-common numerical delta vs alfakR empirical</h3>", img_tag(plot_paths$accuracy_delta, "Strict common accuracy delta violin box plot"),
  "<p class=\"caption\">Content description: paired condition-level MAE/RMSE deltas on exactly shared karyotypes within alfakR fq/NN/other/whole scopes.</p>",
  "<p class=\"interpretation\">Interpretation: ", htmlEscape(accuracy_interp), "</p></div>",
  "<div class=\"figure-block\"><h3>Figure 4. Strict-common shape metrics by method</h3>", img_tag(plot_paths$shape_abs, "Strict common shape metrics violin box plot"),
  "<p class=\"caption\">Content description: centered RMSE, Pearson, and Spearman distributions on shared karyotypes for each method.</p>",
  "<p class=\"interpretation\">Interpretation: centered RMSE captures offset-corrected shape error, while Pearson/Spearman capture relative trend and rank-order preservation.</p></div>",
  "<div class=\"figure-block\"><h3>Figure 5. Strict-common shape delta vs alfakR empirical</h3>", img_tag(plot_paths$shape_delta, "Strict common shape delta violin box plot"),
  "<p class=\"caption\">Content description: paired shape deltas and one-sided Wilcoxon p-values across conditions.</p>",
  "<p class=\"interpretation\">Interpretation: ", htmlEscape(shape_interp), "</p></div>",
  "<h2>Tables</h2>",
  "<div class=\"table-block\"><h3>Table 1. Input minobs summary</h3>",
  "<p class=\"caption\">Content description: retained input karyotype count distribution by minobs. Source TSV: ", table_link(summary_files[["input"]]), ".</p>",
  write_html_table(input_summary), "</div>",
  if (nrow(fit_status)) paste0(
    "<div class=\"table-block\"><h3>Table 2. Fit status summary</h3>",
    "<p class=\"caption\">Content description: number and rate of fit outcomes by engine. This identifies failed alfakR variants before accuracy interpretation.</p>",
    write_html_table(fit_status), "</div>"
  ) else "",
  "<div class=\"table-block\"><h3>Table 3. Own-landscape method summary</h3>",
  "<p class=\"caption\">Content description: median numerical and shape metrics for each method on its own estimated support. Source TSV: ", table_link(summary_files[["own"]]), ".</p>",
  "<p class=\"interpretation\">Interpretation: use this table to compare typical method error while remembering that each method may estimate a different number of nodes.</p>",
  write_html_table(own_report), "</div>",
  "<div class=\"table-block\"><h3>Table 4. Strict-common numerical and shape summary vs alfakR empirical</h3>",
  "<p class=\"caption\">Content description: paired metrics on identical karyotypes. Negative error deltas favor alfak2; positive correlation deltas favor alfak2. Source TSV: ", table_link(summary_files[["strict"]]), ".</p>",
  "<p class=\"interpretation\">Interpretation: this is the primary fair comparison table for precision and shape because it removes coverage differences.</p>",
  write_html_table(strict_report), "</div>",
  "<div class=\"table-block\"><h3>Table 5. Shape-focused strict-common summary</h3>",
  "<p class=\"caption\">Content description: centered RMSE, Pearson, and Spearman comparison by alfakR fq/NN/other/whole scope. Source TSV: ", table_link(summary_files[["shape"]]), ".</p>",
  "<p class=\"interpretation\">Interpretation: use this table to distinguish lower shape error from better rank-order preservation.</p>",
  write_html_table(shape_report), "</div>",
  "<div class=\"table-block\"><h3>Table 6. Wilcoxon signed-rank tests on paired deltas</h3>",
  "<p class=\"caption\">Content description: one-sided tests for alfak2 improvement over alfakR empirical. Source TSV: ", table_link(summary_files[["tests"]]), ".</p>",
  "<p class=\"interpretation\">Interpretation: p-values test whether condition-level deltas consistently favor alfak2, not merely whether the median shown in summary tables is lower or higher.</p>",
  write_html_table(metric_tests[, .(support_scope, minobs, alfak2_label, metric_label, direction, n, median_delta, q25_delta, q75_delta, p_value, p_adj_bh)]), "</div>",
  "<h2>Output Files</h2>",
  "<ul>",
  paste0("<li>", table_link(summary_files), "</li>", collapse = ""),
  "</ul>",
  "<p class=\"small\">Generated by benchmark/scr/render_grf_downsampled_accuracy_report.R.</p>",
  "</body></html>"
)

writeLines(html, report_html, useBytes = TRUE)
message("Wrote report: ", report_html)
