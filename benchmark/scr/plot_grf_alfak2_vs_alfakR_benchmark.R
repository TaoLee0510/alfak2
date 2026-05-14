#!/usr/bin/env Rscript

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) stop("Unexpected positional argument: ", arg, call. = FALSE)
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      val <- "true"
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- val
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || !length(value) || !nzchar(as.character(value[[1L]]))) default else value[[1L]]
}

usage <- function() {
  cat(
    "Plot GRF benchmark summaries for alfak2 vs alfakR.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/plot_grf_alfak2_vs_alfakR_benchmark.R \\\n",
    "    --results-dir=benchmark/results/grf_alfak2_vs_alfakR\n",
    sep = ""
  )
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

method_label <- function(x) {
  map <- c(
    alfak2_effective_full = "alfak2 full",
    alfak2_effective_minobs_matched = "alfak2 minobs-matched",
    alfakR_none = "alfakR none",
    alfakR_empirical = "alfakR empirical",
    alfakR_empirical_censored = "alfakR censored",
    alfakR_empirical_censored_weighted = "alfakR censored weighted",
    alfakR_empirical_two_step = "alfakR two-step"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

policy_label <- function(x) {
  out <- as.character(x)
  out[out == "full"] <- "full"
  out[out == "minobs_matched"] <- "minobs-matched"
  out
}

safe_median <- function(x) median(x, na.rm = TRUE)
safe_q25 <- function(x) stats::quantile(x, probs = 0.25, na.rm = TRUE, names = FALSE)
safe_q75 <- function(x) stats::quantile(x, probs = 0.75, na.rm = TRUE, names = FALSE)

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

read_common_node_source <- function(node_path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required for common-node paired comparisons.", call. = FALSE)
  }
  required_cols <- c(
    "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs", "method", "engine", "input_policy", "nn_prior", "k",
    "estimated_fitness", "true_fitness", "support_scope", "status", "estimation_error"
  )
  header <- names(data.table::fread(node_path, sep = "\t", header = TRUE, nrows = 0, showProgress = FALSE))
  optional_cols <- intersect(c("sim_pm", "pm", "fit_beta_label"), header)
  out <- data.table::fread(
    node_path,
    sep = "\t",
    header = TRUE,
    select = c(required_cols, optional_cols),
    showProgress = FALSE
  )
  out <- out[
    status == "ok" &
      support_scope %in% c("direct", "nn") &
      (method == "alfak2_effective_minobs_matched" | engine == "alfakR")
  ]
  missing <- setdiff(required_cols, names(out))
  if (length(missing)) {
    stop("Missing expected columns in node table: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

make_common_node_vs_alfakR <- function(node_path) {
  dt <- read_common_node_source(node_path)
  dt <- data.table::as.data.table(dt)
  num_cols <- c(
    "simulation_id", "lambda", "time_start", "time_gap", "time_delta", "minobs",
    intersect(c("sim_pm", "pm"), names(dt)),
    "estimated_fitness", "true_fitness", "estimation_error"
  )
  dt[, (num_cols) := lapply(.SD, to_num), .SDcols = num_cols]
  beta_cols <- intersect(c("sim_pm", "pm", "fit_beta_label"), names(dt))
  key_cols <- c("simulation_id", "lambda_label", "time_start", "time_gap", "time_delta", "minobs", beta_cols, "k")
  a2_cols <- c("simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
               "minobs", beta_cols, "k", "support_scope", "estimated_fitness",
               "true_fitness", "estimation_error")
  ar_cols <- c(key_cols, "method", "nn_prior", "support_scope", "estimated_fitness",
               "true_fitness", "estimation_error")

  a2 <- dt[
    method == "alfak2_effective_minobs_matched" & input_policy == "minobs_matched",
    ..a2_cols
  ]
  data.table::setnames(
    a2,
    c("support_scope", "estimated_fitness", "true_fitness", "estimation_error"),
    c("support_scope_a2", "estimated_fitness_a2", "true_fitness_a2", "estimation_error_a2")
  )
  ar <- dt[
    engine == "alfakR",
    ..ar_cols
  ]
  data.table::setnames(
    ar,
    c("method", "nn_prior", "support_scope", "estimated_fitness", "true_fitness", "estimation_error"),
    c("alfakR_method", "alfakR_nn_prior", "support_scope_ar", "estimated_fitness_alfakR", "true_fitness_alfakR", "estimation_error_alfakR")
  )
  rm(dt)
  paired <- merge(ar, a2, by = key_cols, all = FALSE, sort = FALSE, allow.cartesian = TRUE)
  paired <- paired[support_scope_ar == support_scope_a2]
  paired[, support_scope := support_scope_ar]
  paired[, true_fitness := true_fitness_a2]
  paired <- paired[is.finite(estimation_error_a2) & is.finite(estimation_error_alfakR) & is.finite(true_fitness)]

  group_cols <- c(
    "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs", beta_cols, "alfakR_method", "alfakR_nn_prior", "support_scope"
  )
  metrics <- paired[, {
    abs_a2 <- abs(estimation_error_a2)
    abs_ar <- abs(estimation_error_alfakR)
    truth_c <- true_fitness - mean(true_fitness, na.rm = TRUE)
    a2_c <- estimated_fitness_a2 - mean(estimated_fitness_a2, na.rm = TRUE)
    ar_c <- estimated_fitness_alfakR - mean(estimated_fitness_alfakR, na.rm = TRUE)
    centered_err_a2 <- a2_c - truth_c
    centered_err_ar <- ar_c - truth_c
    list(
      n_common = .N,
      n_alfak2_abs_better = sum(abs_a2 < abs_ar, na.rm = TRUE),
      n_alfakR_abs_better = sum(abs_ar < abs_a2, na.rm = TRUE),
      n_abs_ties = sum(abs_a2 == abs_ar, na.rm = TRUE),
      alfak2_mae = mean(abs_a2, na.rm = TRUE),
      alfakR_mae = mean(abs_ar, na.rm = TRUE),
      delta_mae = mean(abs_a2, na.rm = TRUE) - mean(abs_ar, na.rm = TRUE),
      median_delta_abs_error = stats::median(abs_a2 - abs_ar, na.rm = TRUE),
      alfak2_rmse = sqrt(mean(estimation_error_a2^2, na.rm = TRUE)),
      alfakR_rmse = sqrt(mean(estimation_error_alfakR^2, na.rm = TRUE)),
      delta_rmse = sqrt(mean(estimation_error_a2^2, na.rm = TRUE)) - sqrt(mean(estimation_error_alfakR^2, na.rm = TRUE)),
      alfak2_centered_rmse = sqrt(mean(centered_err_a2^2, na.rm = TRUE)),
      alfakR_centered_rmse = sqrt(mean(centered_err_ar^2, na.rm = TRUE)),
      delta_centered_rmse = sqrt(mean(centered_err_a2^2, na.rm = TRUE)) - sqrt(mean(centered_err_ar^2, na.rm = TRUE)),
      alfak2_spearman = safe_cor(estimated_fitness_a2, true_fitness, method = "spearman"),
      alfakR_spearman = safe_cor(estimated_fitness_alfakR, true_fitness, method = "spearman"),
      delta_spearman = safe_cor(estimated_fitness_a2, true_fitness, method = "spearman") -
        safe_cor(estimated_fitness_alfakR, true_fitness, method = "spearman"),
      estimate_pearson_between_methods = safe_cor(estimated_fitness_a2, estimated_fitness_alfakR, method = "pearson"),
      median_abs_estimate_diff = stats::median(abs(estimated_fitness_a2 - estimated_fitness_alfakR), na.rm = TRUE)
    )
  }, by = group_cols]

  summary <- metrics[, .(
    n_conditions = as.numeric(.N),
    median_common_nodes = as.numeric(stats::median(n_common, na.rm = TRUE)),
    total_common_nodes = as.numeric(sum(n_common, na.rm = TRUE)),
    pooled_node_win_rate = as.numeric(sum(n_alfak2_abs_better, na.rm = TRUE) / sum(n_common, na.rm = TRUE)),
    pooled_alfakR_node_win_rate = as.numeric(sum(n_alfakR_abs_better, na.rm = TRUE) / sum(n_common, na.rm = TRUE)),
    pooled_tie_rate = as.numeric(sum(n_abs_ties, na.rm = TRUE) / sum(n_common, na.rm = TRUE)),
    median_node_win_rate_by_condition = as.numeric(stats::median(n_alfak2_abs_better / n_common, na.rm = TRUE)),
    alfak2_median_mae = as.numeric(stats::median(alfak2_mae, na.rm = TRUE)),
    alfakR_median_mae = as.numeric(stats::median(alfakR_mae, na.rm = TRUE)),
    paired_median_delta_mae = as.numeric(stats::median(delta_mae, na.rm = TRUE)),
    mae_condition_win_rate = mean(delta_mae < 0, na.rm = TRUE),
    median_delta_abs_error = as.numeric(stats::median(median_delta_abs_error, na.rm = TRUE)),
    alfak2_median_rmse = as.numeric(stats::median(alfak2_rmse, na.rm = TRUE)),
    alfakR_median_rmse = as.numeric(stats::median(alfakR_rmse, na.rm = TRUE)),
    paired_median_delta_rmse = as.numeric(stats::median(delta_rmse, na.rm = TRUE)),
    rmse_condition_win_rate = mean(delta_rmse < 0, na.rm = TRUE),
    alfak2_median_centered_rmse = as.numeric(stats::median(alfak2_centered_rmse, na.rm = TRUE)),
    alfakR_median_centered_rmse = as.numeric(stats::median(alfakR_centered_rmse, na.rm = TRUE)),
    paired_median_delta_centered_rmse = as.numeric(stats::median(delta_centered_rmse, na.rm = TRUE)),
    centered_rmse_condition_win_rate = mean(delta_centered_rmse < 0, na.rm = TRUE),
    paired_median_delta_spearman = as.numeric(stats::median(delta_spearman, na.rm = TRUE)),
    spearman_condition_win_rate = mean(delta_spearman > 0, na.rm = TRUE),
    median_abs_estimate_diff = as.numeric(stats::median(median_abs_estimate_diff, na.rm = TRUE))
  ), by = c("minobs", beta_cols, "support_scope", "alfakR_method", "alfakR_nn_prior")]

  data.table::setorder(metrics, support_scope, alfakR_method, minobs, simulation_id, lambda_label, time_gap)
  data.table::setorder(summary, support_scope, alfakR_method, minobs)
  list(metrics = as.data.frame(metrics), summary = as.data.frame(summary))
}

summarize_metric <- function(df, group_cols, metric) {
  split_key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(df, split_key), function(x) {
    data.frame(
      x[1L, group_cols, drop = FALSE],
      n = nrow(x),
      median = safe_median(x[[metric]]),
      q25 = safe_q25(x[[metric]]),
      q75 = safe_q75(x[[metric]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

save_plot <- function(plot, path, width = 8, height = 5) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 300)
  ggplot2::ggsave(sub("[.]png$", ".pdf", path), plot = plot, width = width, height = height)
  invisible(path)
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package scales is required for plotting.", call. = FALSE)
  }

  repo_dir <- normalizePath(arg_value(args, "repo_dir", getwd()), winslash = "/", mustWork = FALSE)
  results_dir <- arg_value(args, "results_dir", "benchmark/results/grf_alfak2_vs_alfakR")
  if (!grepl("^/", results_dir)) {
    results_dir <- file.path(repo_dir, results_dir)
  }
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
  tables_dir <- file.path(results_dir, "tables")
  figures_dir <- file.path(results_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  summary_tbl <- read_tsv(file.path(tables_dir, "summary_by_lambda_time_minobs_method.tsv"))
  fit_tbl <- read_tsv(file.path(tables_dir, "fit_results.tsv"))
  input_tbl <- read_tsv(file.path(tables_dir, "input_table.tsv"))
  delta_tbl <- read_tsv(file.path(tables_dir, "alfak2_delta_vs_alfakR.tsv"))

  numeric_summary <- c(
    "lambda", "time_start", "time_gap", "time_delta", "minobs", "sim_pm", "pm", "n_nodes", "n_scored",
    "observed_node_fraction", "pearson", "spearman", "rmse", "mae", "signed_bias",
    "centered_rmse", "centered_mae", "centered_r2", "sign_accuracy", "false_high_rate",
    "mean_estimated_sd"
  )
  for (nm in intersect(numeric_summary, names(summary_tbl))) summary_tbl[[nm]] <- to_num(summary_tbl[[nm]])
  for (nm in intersect(c("minobs", "sim_pm", "pm", "local_nodes", "global_nodes", "elapsed_sec", "local_gradient_norm"), names(fit_tbl))) {
    fit_tbl[[nm]] <- to_num(fit_tbl[[nm]])
  }
  for (nm in intersect(c("minobs", "sim_pm", "raw_input_rows", "input_rows_after_drop", "input_rows_minobs"), names(input_tbl))) {
    input_tbl[[nm]] <- to_num(input_tbl[[nm]])
  }
  for (nm in grep("^delta_", names(delta_tbl), value = TRUE)) delta_tbl[[nm]] <- to_num(delta_tbl[[nm]])
  delta_tbl$minobs <- to_num(delta_tbl$minobs)
  for (nm in intersect(c("sim_pm", "pm"), names(delta_tbl))) delta_tbl[[nm]] <- to_num(delta_tbl[[nm]])

  summary_tbl$method_label <- method_label(summary_tbl$method)
  fit_tbl$method_label <- method_label(fit_tbl$method)
  delta_tbl$alfak2_policy_label <- policy_label(delta_tbl$alfak2_input_policy)
  delta_tbl$alfakR_label <- method_label(delta_tbl$alfakR_method)
  summary_tbl$minobs_f <- factor(summary_tbl$minobs, levels = sort(unique(summary_tbl$minobs)))
  input_tbl$minobs_f <- factor(input_tbl$minobs, levels = sort(unique(input_tbl$minobs)))
  fit_tbl$minobs_f <- factor(fit_tbl$minobs, levels = sort(unique(fit_tbl$minobs)))

  ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))
  method_colors <- c(
    "alfak2 minobs-matched" = "#0072B2",
    "alfak2 full" = "#009E73",
    "alfakR none" = "#D55E00",
    "alfakR empirical" = "#CC79A7",
    "alfakR censored" = "#E69F00",
    "alfakR censored weighted" = "#F0E442",
    "alfakR two-step" = "#999999"
  )

  input_long <- rbind(
    data.frame(minobs_f = input_tbl$minobs_f, nodes = input_tbl$input_rows_after_drop, type = "after diploid drop"),
    data.frame(minobs_f = input_tbl$minobs_f, nodes = input_tbl$input_rows_minobs, type = "after minobs filter")
  )
  p1 <- ggplot2::ggplot(input_long, ggplot2::aes(x = minobs_f, y = nodes, fill = type)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.35) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    ggplot2::scale_fill_manual(values = c("after diploid drop" = "#BBBBBB", "after minobs filter" = "#0072B2")) +
    ggplot2::labs(x = "minobs", y = "input nodes (log10 scale)", fill = NULL, title = "Input support shrinks under minobs filtering") +
    ggplot2::theme(legend.position = "top")
  save_plot(p1, file.path(figures_dir, "01_input_nodes_by_minobs.png"), width = 7, height = 4.8)

  x2 <- summary_tbl[
    summary_tbl$method %in% c("alfak2_effective_minobs_matched", "alfak2_effective_full") &
      summary_tbl$support_scope %in% c("direct", "all"),
    , drop = FALSE
  ]
  x2_sum <- summarize_metric(x2, c("method_label", "minobs", "minobs_f", "support_scope"), "centered_rmse")
  p2 <- ggplot2::ggplot(
    x2_sum,
    ggplot2::aes(x = minobs, y = median, color = method_label, group = method_label)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75, fill = method_label), alpha = 0.14, color = NA) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~support_scope, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(x2_sum$minobs))) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::labs(x = "minobs", y = "median centered RMSE", color = NULL, fill = NULL, title = "Effect of minobs on alfak2 accuracy") +
    ggplot2::theme(legend.position = "top")
  save_plot(p2, file.path(figures_dir, "02_alfak2_minobs_accuracy.png"), width = 8, height = 4.8)

  x3 <- summary_tbl[
    summary_tbl$method == "alfak2_effective_minobs_matched" &
      summary_tbl$support_scope %in% c("direct", "all"),
    , drop = FALSE
  ]
  x3_sum <- summarize_metric(x3, c("lambda_label", "time_gap", "minobs", "support_scope"), "centered_rmse")
  p3 <- ggplot2::ggplot(
    x3_sum,
    ggplot2::aes(x = minobs, y = median, color = factor(time_gap), group = factor(time_gap))
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_grid(support_scope ~ lambda_label, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(x3_sum$minobs))) +
    ggplot2::labs(x = "minobs", y = "median centered RMSE", color = "time_gap", title = "alfak2 minobs-matched accuracy by lambda and time gap") +
    ggplot2::theme(legend.position = "top")
  save_plot(p3, file.path(figures_dir, "03_alfak2_minobs_by_lambda.png"), width = 10.5, height = 5.8)

  x4 <- fit_tbl[fit_tbl$method == "alfak2_effective_minobs_matched", , drop = FALSE]
  diag_counts <- as.data.frame(table(x4$minobs_f, x4$local_covariance_status), stringsAsFactors = FALSE)
  names(diag_counts) <- c("minobs_f", "status", "n")
  diag_counts <- diag_counts[diag_counts$n > 0, , drop = FALSE]
  p4 <- ggplot2::ggplot(diag_counts, ggplot2::aes(x = minobs_f, y = n, fill = status)) +
    ggplot2::geom_col(position = "fill", width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = "minobs", y = "fraction of fits", fill = "local covariance status", title = "alfak2 minobs-matched diagnostics improve as minobs increases") +
    ggplot2::theme(legend.position = "top")
  save_plot(p4, file.path(figures_dir, "04_alfak2_minobs_diagnostics.png"), width = 7.5, height = 4.8)

  methods_keep <- c(
    "alfak2_effective_minobs_matched", "alfak2_effective_full", "alfakR_none",
    "alfakR_empirical", "alfakR_empirical_censored",
    "alfakR_empirical_censored_weighted", "alfakR_empirical_two_step"
  )
  x5 <- summary_tbl[
    summary_tbl$method %in% methods_keep &
      summary_tbl$support_scope %in% c("direct", "all"),
    , drop = FALSE
  ]
  x5_sum <- summarize_metric(x5, c("method_label", "support_scope"), "centered_rmse")
  x5_sum$method_label <- factor(x5_sum$method_label, levels = c(
    "alfak2 minobs-matched", "alfak2 full", "alfakR none", "alfakR empirical",
    "alfakR censored", "alfakR censored weighted", "alfakR two-step"
  ))
  p5 <- ggplot2::ggplot(x5_sum, ggplot2::aes(x = method_label, y = median, color = method_label)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = q25, ymax = q75), linewidth = 0.55) +
    ggplot2::facet_wrap(~support_scope, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(x = NULL, y = "centered RMSE, median and IQR", color = NULL, title = "Overall method accuracy") +
    ggplot2::theme(legend.position = "none")
  save_plot(p5, file.path(figures_dir, "05_method_accuracy_overall.png"), width = 8.5, height = 5.2)

  x6 <- delta_tbl[delta_tbl$support_scope %in% c("direct", "all"), , drop = FALSE]
  split_key <- interaction(x6[, c("support_scope", "alfak2_policy_label", "alfakR_label")], drop = TRUE, lex.order = TRUE)
  x6_sum <- do.call(rbind, lapply(split(x6, split_key), function(x) {
    data.frame(
      support_scope = x$support_scope[[1]],
      alfak2_policy_label = x$alfak2_policy_label[[1]],
      alfakR_label = x$alfakR_label[[1]],
      median_delta_rmse = safe_median(x$delta_centered_rmse),
      win_rate = mean(x$delta_centered_rmse < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  x6_sum$alfakR_label <- factor(x6_sum$alfakR_label, levels = c(
    "alfakR none", "alfakR empirical", "alfakR censored",
    "alfakR censored weighted", "alfakR two-step"
  ))
  x6_sum$alfak2_policy_label <- factor(x6_sum$alfak2_policy_label, levels = c("minobs-matched", "full"))
  p6 <- ggplot2::ggplot(
    x6_sum,
    ggplot2::aes(x = alfakR_label, y = alfak2_policy_label, fill = median_delta_rmse)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f\n%.0f%%", median_delta_rmse, 100 * win_rate)), size = 3) +
    ggplot2::facet_wrap(~support_scope) +
    ggplot2::scale_fill_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      midpoint = 0,
      labels = scales::label_number(accuracy = 0.001),
      guide = ggplot2::guide_colorbar(title.position = "top", barwidth = 12, barheight = 0.8)
    ) +
    ggplot2::labs(x = NULL, y = "alfak2 policy", fill = "median delta RMSE", title = "alfak2 minus alfakR centered RMSE") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_plot(p6, file.path(figures_dir, "06_alfak2_delta_vs_alfakR_heatmap.png"), width = 9, height = 4.8)

  x7_key <- c("lambda", "lambda_label", "time_start", "time_gap", "time_delta", "minobs",
              intersect(c("sim_pm", "pm", "fit_beta_label"), names(summary_tbl)), "support_scope")
  x7_cols <- c(
    x7_key, "centered_rmse", "spearman", "pearson", "false_high_rate",
    "n_scored", "n_nodes", "mean_estimated_sd"
  )
  x7_matched <- summary_tbl[
    summary_tbl$method == "alfak2_effective_minobs_matched" &
      summary_tbl$support_scope %in% c("direct", "nn"),
    x7_cols,
    drop = FALSE
  ]
  x7_full <- summary_tbl[
    summary_tbl$method == "alfak2_effective_full" &
      summary_tbl$support_scope %in% c("direct", "nn"),
    x7_cols,
    drop = FALSE
  ]
  x7_pair <- merge(x7_matched, x7_full, by = x7_key, suffixes = c("_matched", "_full"))
  x7_pair$delta_centered_rmse <- x7_pair$centered_rmse_matched - x7_pair$centered_rmse_full
  x7_pair$delta_spearman <- x7_pair$spearman_matched - x7_pair$spearman_full
  x7_pair$delta_false_high_rate <- x7_pair$false_high_rate_matched - x7_pair$false_high_rate_full
  x7_pair$scored_ratio <- x7_pair$n_scored_matched / x7_pair$n_scored_full

  x7_split <- interaction(x7_pair[, c("minobs", "support_scope")], drop = TRUE, lex.order = TRUE)
  x7_sum <- do.call(rbind, lapply(split(x7_pair, x7_split), function(x) {
    data.frame(
      minobs = x$minobs[[1L]],
      support_scope = x$support_scope[[1L]],
      n_conditions = nrow(x),
      matched_median_centered_rmse = safe_median(x$centered_rmse_matched),
      full_median_centered_rmse = safe_median(x$centered_rmse_full),
      paired_median_delta_centered_rmse = safe_median(x$delta_centered_rmse),
      rmse_win_rate = mean(x$delta_centered_rmse < 0, na.rm = TRUE),
      matched_median_spearman = safe_median(x$spearman_matched),
      full_median_spearman = safe_median(x$spearman_full),
      paired_median_delta_spearman = safe_median(x$delta_spearman),
      matched_median_false_high_rate = safe_median(x$false_high_rate_matched),
      full_median_false_high_rate = safe_median(x$false_high_rate_full),
      paired_median_delta_false_high_rate = safe_median(x$delta_false_high_rate),
      matched_median_n_scored = safe_median(x$n_scored_matched),
      full_median_n_scored = safe_median(x$n_scored_full),
      median_scored_ratio = safe_median(x$scored_ratio),
      matched_total_scored = sum(x$n_scored_matched, na.rm = TRUE),
      full_total_scored = sum(x$n_scored_full, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  x7_sum$support_label <- ifelse(x7_sum$support_scope == "direct", "fq / direct informed", "NN")
  x7_sum$minobs_f <- factor(x7_sum$minobs, levels = sort(unique(x7_sum$minobs)))
  x7_sum$support_label <- factor(x7_sum$support_label, levels = c("fq / direct informed", "NN"))
  utils::write.table(
    x7_sum,
    file.path(tables_dir, "alfak2_full_vs_minobs_matched_by_minobs_support.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p7 <- ggplot2::ggplot(
    x7_sum,
    ggplot2::aes(x = minobs_f, y = paired_median_delta_centered_rmse, fill = paired_median_delta_centered_rmse < 0)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linewidth = 0.35) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.3f\nwin %.0f%%", paired_median_delta_centered_rmse, 100 * rmse_win_rate)),
      size = 3,
      vjust = ifelse(x7_sum$paired_median_delta_centered_rmse < 0, 1.15, -0.15),
      color = "white"
    ) +
    ggplot2::facet_wrap(~support_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"), guide = "none") +
    ggplot2::labs(
      x = "minobs",
      y = "paired median delta centered RMSE\n(minobs-matched - full)",
      title = "alfak2 minobs-matched vs full at the same minobs",
      subtitle = "Negative values indicate lower RMSE for minobs-matched"
    )
  save_plot(p7, file.path(figures_dir, "07_alfak2_full_vs_minobs_matched_direct_nn.png"), width = 8, height = 4.8)

  common_metrics_path <- file.path(tables_dir, "alfak2_minobs_matched_vs_alfakR_common_node_metrics.tsv")
  common_summary_path <- file.path(tables_dir, "alfak2_minobs_matched_vs_alfakR_common_node_summary.tsv")
  node_path <- file.path(tables_dir, "node_accuracy.tsv")
  common_needs_update <- !file.exists(common_metrics_path) || !file.exists(common_summary_path)
  if (!common_needs_update && file.exists(node_path)) {
    common_needs_update <- file.info(common_metrics_path)$mtime < file.info(node_path)$mtime ||
      file.info(common_summary_path)$mtime < file.info(node_path)$mtime
  }
  if (common_needs_update) {
    message("Building common-karyotype alfak2 vs alfakR paired metrics from: ", node_path)
    common <- make_common_node_vs_alfakR(node_path)
    write_tsv(common$metrics, common_metrics_path)
    write_tsv(common$summary, common_summary_path)
    x8_sum <- common$summary
  } else {
    x8_sum <- read_tsv(common_summary_path)
  }

  numeric_common <- c(
    "minobs", "n_conditions", "median_common_nodes", "total_common_nodes",
    "pooled_node_win_rate", "pooled_alfakR_node_win_rate", "pooled_tie_rate",
    "median_node_win_rate_by_condition", "alfak2_median_mae", "alfakR_median_mae",
    "paired_median_delta_mae", "mae_condition_win_rate", "median_delta_abs_error",
    "alfak2_median_rmse", "alfakR_median_rmse", "paired_median_delta_rmse",
    "rmse_condition_win_rate", "alfak2_median_centered_rmse", "alfakR_median_centered_rmse",
    "paired_median_delta_centered_rmse", "centered_rmse_condition_win_rate",
    "paired_median_delta_spearman", "spearman_condition_win_rate", "median_abs_estimate_diff"
  )
  for (nm in intersect(numeric_common, names(x8_sum))) x8_sum[[nm]] <- to_num(x8_sum[[nm]])
  x8_sum$support_label <- ifelse(x8_sum$support_scope == "direct", "fq / direct informed", "NN")
  x8_sum$alfakR_label <- method_label(x8_sum$alfakR_method)
  x8_sum$minobs_f <- factor(x8_sum$minobs, levels = sort(unique(x8_sum$minobs)))
  x8_sum$support_label <- factor(x8_sum$support_label, levels = c("fq / direct informed", "NN"))
  x8_sum$alfakR_label <- factor(x8_sum$alfakR_label, levels = c(
    "alfakR none", "alfakR empirical", "alfakR censored",
    "alfakR censored weighted", "alfakR two-step"
  ))
  x8_sum <- x8_sum[order(x8_sum$support_label, x8_sum$alfakR_label, x8_sum$minobs), , drop = FALSE]

  p8 <- ggplot2::ggplot(
    x8_sum,
    ggplot2::aes(x = minobs_f, y = alfakR_label, fill = paired_median_delta_mae)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.3f\nwin %.0f%%", paired_median_delta_mae, 100 * pooled_node_win_rate)),
      size = 3
    ) +
    ggplot2::facet_wrap(~support_label) +
    ggplot2::scale_fill_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      midpoint = 0,
      labels = scales::label_number(accuracy = 0.001),
      guide = ggplot2::guide_colorbar(title.position = "top", barwidth = 12, barheight = 0.8)
    ) +
    ggplot2::labs(
      x = "matched min_obs / minobs",
      y = NULL,
      fill = "delta MAE",
      title = "Common-karyotype error: alfak2 minobs-matched minus alfakR",
      subtitle = "Only karyotypes estimated by both methods in the same support class are compared; negative values favor alfak2"
    ) +
    ggplot2::theme(legend.position = "bottom")
  save_plot(p8, file.path(figures_dir, "08_alfak2_minobs_matched_vs_alfakR_same_minobs_direct_nn.png"), width = 9.5, height = 5.5)

  p9 <- ggplot2::ggplot(
    x8_sum,
    ggplot2::aes(x = minobs, y = paired_median_delta_centered_rmse, color = alfakR_label, group = alfakR_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~support_label, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = sort(unique(x8_sum$minobs))) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(
      x = "matched min_obs / minobs",
      y = "paired median delta centered RMSE\n(alfak2 - alfakR)",
      color = NULL,
      title = "Common-karyotype centered RMSE difference"
    ) +
    ggplot2::theme(legend.position = "top")
  save_plot(p9, file.path(figures_dir, "09_alfak2_minobs_matched_vs_alfakR_common_nodes_centered_rmse.png"), width = 8.5, height = 5.2)

  message("Wrote figures to: ", figures_dir)
  invisible(figures_dir)
}

if (sys.nframe() == 0L) {
  main()
}
