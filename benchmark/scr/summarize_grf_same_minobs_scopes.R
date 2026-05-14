#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "benchmark/results/grf_alfak2_vs_alfakR_calibrated"
}
if (!dir.exists(results_dir)) {
  stop("Results directory does not exist: ", results_dir, call. = FALSE)
}

tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

safe_q25 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, 0.25, names = FALSE))
}

safe_q75 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, 0.75, names = FALSE))
}

method_label <- function(x) {
  out <- x
  out[x == "alfak2_effective_minobs_matched"] <- "alfak2 minobs-matched"
  out[x == "alfakR_none"] <- "alfakR none"
  out[x == "alfakR_empirical"] <- "alfakR empirical"
  out[x == "alfakR_empirical_censored"] <- "alfakR censored"
  out[x == "alfakR_empirical_censored_weighted"] <- "alfakR censored weighted"
  out[x == "alfakR_empirical_two_step"] <- "alfakR two-step"
  out
}

scope_label <- function(x) {
  out <- x
  out[x == "direct"] <- "fq"
  out[x == "nn"] <- "NN"
  out[x == "all"] <- "whole"
  out
}

summarize_metric <- function(df, group_cols, metric) {
  split_key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(df, split_key), function(x) {
    data.frame(
      x[1L, group_cols, drop = FALSE],
      n_conditions = nrow(x),
      median = safe_median(x[[metric]]),
      q25 = safe_q25(x[[metric]]),
      q75 = safe_q75(x[[metric]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

summary_path <- file.path(tables_dir, "summary_by_lambda_time_minobs_method.tsv")
summary_tbl <- read_tsv(summary_path)

numeric_cols <- c(
  "lambda", "time_start", "time_gap", "time_delta", "minobs", "sim_pm", "pm", "n_nodes", "n_scored",
  "observed_node_fraction", "pearson", "spearman", "rmse", "mae", "signed_bias",
  "centered_rmse", "centered_mae", "centered_r2", "sign_accuracy", "false_high_rate",
  "mean_estimated_sd"
)
for (nm in intersect(numeric_cols, names(summary_tbl))) {
  summary_tbl[[nm]] <- to_num(summary_tbl[[nm]])
}

methods_keep <- c(
  "alfak2_effective_minobs_matched",
  "alfakR_none",
  "alfakR_empirical",
  "alfakR_empirical_censored",
  "alfakR_empirical_censored_weighted",
  "alfakR_empirical_two_step"
)
scopes_keep <- c("direct", "nn", "all")
beta_cols <- intersect(c("sim_pm", "pm", "fit_beta_label"), names(summary_tbl))
key_cols <- c("lambda", "lambda_label", "time_start", "time_gap", "time_delta", "minobs", beta_cols, "support_scope")
metric_cols <- c(
  "n_scored", "pearson", "spearman", "rmse", "mae", "centered_rmse", "centered_mae",
  "centered_r2", "sign_accuracy", "false_high_rate", "mean_estimated_sd"
)

scope_tbl <- summary_tbl[
  summary_tbl$method %in% methods_keep & summary_tbl$support_scope %in% scopes_keep,
  ,
  drop = FALSE
]
scope_tbl$scope_label <- scope_label(scope_tbl$support_scope)
scope_tbl$method_label <- method_label(scope_tbl$method)
scope_tbl <- scope_tbl[order(scope_tbl$support_scope, scope_tbl$method, scope_tbl$minobs, scope_tbl$lambda, scope_tbl$time_gap), ]

summary_rows <- list()
idx <- 1L
for (metric in metric_cols) {
  metric_summary <- summarize_metric(
    scope_tbl,
    c("minobs", beta_cols, "support_scope", "scope_label", "method", "method_label"),
    metric
  )
  metric_summary$metric <- metric
  summary_rows[[idx]] <- metric_summary
  idx <- idx + 1L
}
performance_summary <- do.call(rbind, summary_rows)
performance_summary <- performance_summary[
  order(performance_summary$scope_label, performance_summary$metric, performance_summary$method_label, performance_summary$minobs),
  ,
  drop = FALSE
]

a2 <- scope_tbl[scope_tbl$method == "alfak2_effective_minobs_matched", c(key_cols, metric_cols), drop = FALSE]
names(a2)[match(metric_cols, names(a2))] <- paste0(metric_cols, "_alfak2")
ar <- scope_tbl[scope_tbl$method != "alfak2_effective_minobs_matched", c(key_cols, "method", "method_label", "nn_prior", metric_cols), drop = FALSE]
names(ar)[match(metric_cols, names(ar))] <- paste0(metric_cols, "_alfakR")
delta <- merge(ar, a2, by = key_cols, all = FALSE, sort = FALSE)
delta$scope_label <- scope_label(delta$support_scope)
delta$alfakR_method <- delta$method
delta$alfakR_label <- delta$method_label
delta$method <- NULL
delta$method_label <- NULL
for (metric in setdiff(metric_cols, "n_scored")) {
  delta[[paste0("delta_", metric)]] <- delta[[paste0(metric, "_alfak2")]] - delta[[paste0(metric, "_alfakR")]]
}
delta$scored_ratio_alfak2_over_alfakR <- delta$n_scored_alfak2 / delta$n_scored_alfakR

delta_metric_cols <- grep("^delta_", names(delta), value = TRUE)
delta_group_cols <- c("minobs", beta_cols, "support_scope", "scope_label", "alfakR_method", "alfakR_label", "nn_prior")
delta_rows <- list()
idx <- 1L
for (metric in delta_metric_cols) {
  split_key <- interaction(delta[delta_group_cols], drop = TRUE, lex.order = TRUE)
  metric_summary <- do.call(rbind, lapply(split(delta, split_key), function(x) {
    values <- x[[metric]]
    if (metric %in% c("delta_centered_rmse", "delta_rmse", "delta_mae", "delta_centered_mae", "delta_false_high_rate")) {
      win_rate <- mean(values < 0, na.rm = TRUE)
    } else if (metric %in% c("delta_spearman", "delta_pearson", "delta_centered_r2", "delta_sign_accuracy")) {
      win_rate <- mean(values > 0, na.rm = TRUE)
    } else {
      win_rate <- NA_real_
    }
    data.frame(
      x[1L, delta_group_cols, drop = FALSE],
      n_conditions = nrow(x),
      median = safe_median(values),
      q25 = safe_q25(values),
      q75 = safe_q75(values),
      metric = metric,
      alfak2_win_rate = win_rate,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  delta_rows[[idx]] <- metric_summary
  idx <- idx + 1L
}
delta_summary <- do.call(rbind, delta_rows)
delta_summary <- delta_summary[
  order(delta_summary$scope_label, delta_summary$metric, delta_summary$alfakR_label, delta_summary$minobs),
  ,
  drop = FALSE
]

write_tsv(scope_tbl, file.path(tables_dir, "same_minobs_scope_condition_metrics.tsv"))
write_tsv(performance_summary, file.path(tables_dir, "same_minobs_scope_performance_summary.tsv"))
write_tsv(delta, file.path(tables_dir, "same_minobs_scope_alfak2_vs_alfakR_condition_delta.tsv"))
write_tsv(delta_summary, file.path(tables_dir, "same_minobs_scope_alfak2_vs_alfakR_delta_summary.tsv"))

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("scales", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))
  method_levels <- c(
    "alfak2 minobs-matched",
    "alfakR none",
    "alfakR empirical",
    "alfakR censored",
    "alfakR censored weighted",
    "alfakR two-step"
  )
  scope_levels <- c("fq", "NN", "whole")
  x <- performance_summary[
    performance_summary$metric == "centered_rmse",
    ,
    drop = FALSE
  ]
  x$method_label <- factor(x$method_label, levels = method_levels)
  x$scope_label <- factor(x$scope_label, levels = scope_levels)
  x$beta_label <- if ("fit_beta_label" %in% names(x)) paste0("beta=", x$fit_beta_label) else "beta=single"
  p1 <- ggplot2::ggplot(
    x,
    ggplot2::aes(x = factor(minobs), y = median, color = method_label, group = method_label)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75, fill = method_label), alpha = 0.12, color = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::facet_grid(scope_label ~ beta_label, scales = "free_y") +
    ggplot2::labs(
      x = "matched min_obs / minobs",
      y = "centered RMSE, median and IQR",
      color = NULL,
      fill = NULL,
      title = "alfak2 and alfakR performance at the same minobs by scope"
    ) +
    ggplot2::theme(legend.position = "bottom")
  fig_width <- if (length(unique(x$beta_label)) > 1L) 14 else 10
  ggplot2::ggsave(file.path(figures_dir, "same_minobs_scope_centered_rmse.png"), p1, width = fig_width, height = 7.2, dpi = 300)
  ggplot2::ggsave(file.path(figures_dir, "same_minobs_scope_centered_rmse.pdf"), p1, width = fig_width, height = 7.2)

  y <- delta_summary[delta_summary$metric == "delta_centered_rmse", , drop = FALSE]
  y$alfakR_label <- factor(y$alfakR_label, levels = method_levels[-1L])
  y$scope_label <- factor(y$scope_label, levels = scope_levels)
  y$beta_label <- if ("fit_beta_label" %in% names(y)) paste0("beta=", y$fit_beta_label) else "beta=single"
  p2 <- ggplot2::ggplot(
    y,
    ggplot2::aes(x = factor(minobs), y = alfakR_label, fill = median)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f\nwin %.0f%%", median, 100 * alfak2_win_rate)), size = 2.8) +
    ggplot2::facet_grid(scope_label ~ beta_label) +
    ggplot2::scale_fill_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      midpoint = 0,
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggplot2::labs(
      x = "matched min_obs / minobs",
      y = NULL,
      fill = "delta centered RMSE",
      title = "alfak2 minobs-matched minus alfakR at the same minobs",
      subtitle = "Negative values favor alfak2; tiles show median delta and condition win rate"
    ) +
    ggplot2::theme(legend.position = "bottom")
  delta_fig_width <- if (length(unique(y$beta_label)) > 1L) 14 else 11
  ggplot2::ggsave(file.path(figures_dir, "same_minobs_scope_delta_centered_rmse.png"), p2, width = delta_fig_width, height = 7.2, dpi = 300)
  ggplot2::ggsave(file.path(figures_dir, "same_minobs_scope_delta_centered_rmse.pdf"), p2, width = delta_fig_width, height = 7.2)
}

message("Wrote same-minobs scope summaries to: ", tables_dir)
message("Wrote figures to: ", figures_dir)
