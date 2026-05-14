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

repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
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

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

safe_centered_r2 <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L) return(NA_real_)
  truth_c <- truth[ok] - mean(truth[ok])
  pred_c <- pred[ok] - mean(pred[ok])
  den <- sum(truth_c^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  1 - sum((pred_c - truth_c)^2) / den
}

localize_path <- function(path) {
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(path)
  if (file.exists(path)) return(path)
  anchor <- file.path("benchmark", "results", basename(results_dir))
  pos <- regexpr(anchor, path, fixed = TRUE)[[1L]]
  if (pos > 0L) {
    candidate <- file.path(repo_dir, substr(path, pos, nchar(path)))
    if (file.exists(candidate)) return(candidate)
  }
  path
}

parse_karyotype_ids_base <- function(ids) {
  ids <- as.character(ids)
  pieces <- strsplit(ids, ".", fixed = TRUE)
  n_chr <- lengths(pieces)
  if (!length(ids) || any(n_chr == 0L) || length(unique(n_chr)) != 1L) {
    stop("Malformed karyotype IDs.", call. = FALSE)
  }
  mat <- do.call(rbind, lapply(pieces, as.integer))
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

compute_grf_fitness_truth <- function(karyotypes, centroids, lambda) {
  k_mat <- parse_karyotype_ids_base(karyotypes)
  if (!is.matrix(centroids) || !is.numeric(centroids) || !nrow(centroids)) {
    stop("`centroids` must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("`lambda` must be a positive finite scalar.", call. = FALSE)
  }
  if (ncol(k_mat) != ncol(centroids)) {
    stop("Karyotype dimension does not match GRF centroid dimension.", call. = FALSE)
  }
  out <- vapply(seq_len(nrow(k_mat)), function(i) {
    diffs <- sweep(centroids, 2L, as.numeric(k_mat[i, ]), FUN = "-")
    distances <- sqrt(rowSums(diffs^2))
    sum(sin(distances / lambda)) / (pi * sqrt(nrow(centroids)))
  }, numeric(1))
  names(out) <- rownames(k_mat)
  out
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
  out
}

summarize_estimates <- function(est, est_sd, truth) {
  ok <- is.finite(est) & is.finite(truth)
  err <- est[ok] - truth[ok]
  pred_c <- est[ok] - mean(est[ok], na.rm = TRUE)
  truth_c <- truth[ok] - mean(truth[ok], na.rm = TRUE)
  centered_err <- pred_c - truth_c
  data.frame(
    n_nodes = length(est),
    n_scored = sum(ok),
    pearson = safe_cor(est, truth, "pearson"),
    spearman = safe_cor(est, truth, "spearman"),
    rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    mae = if (length(err)) mean(abs(err)) else NA_real_,
    signed_bias = if (length(err)) mean(err) else NA_real_,
    centered_rmse = if (length(centered_err)) sqrt(mean(centered_err^2)) else NA_real_,
    centered_mae = if (length(centered_err)) mean(abs(centered_err)) else NA_real_,
    centered_r2 = safe_centered_r2(est, truth),
    sign_accuracy = if (length(centered_err)) mean(sign(pred_c) == sign(truth_c), na.rm = TRUE) else NA_real_,
    false_high_rate = if (length(centered_err)) mean(pred_c > 0 & truth_c <= 0, na.rm = TRUE) else NA_real_,
    mean_estimated_sd = mean(est_sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
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

extract_bootstrap_scope <- function(boot, scope) {
  mat <- if (identical(scope, "direct")) boot$final_fitness else boot$nn_fitness
  if (is.null(mat) || !is.matrix(mat) || !ncol(mat)) return(NULL)
  labels <- colnames(mat)
  if (is.null(labels) || !length(labels)) return(NULL)
  est_mean <- colMeans(mat, na.rm = TRUE)
  est_median <- apply(mat, 2L, stats::median, na.rm = TRUE)
  est_sd <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  data.frame(
    k = labels,
    estimated_fitness = as.numeric(est_mean),
    estimated_median = as.numeric(est_median),
    estimated_sd = as.numeric(est_sd),
    stringsAsFactors = FALSE
  )
}

extract_alfak2_nodes <- function(fit) {
  s <- fit$global$summary
  if (is.null(s) || !nrow(s)) return(NULL)
  est_col <- if ("fitness_mean_alfakR_scale" %in% names(s)) "fitness_mean_alfakR_scale" else "fitness_mean"
  sd_col <- if ("fitness_sd_alfakR_scale" %in% names(s)) "fitness_sd_alfakR_scale" else "fitness_sd"
  if (!all(c("karyotype", est_col, sd_col) %in% names(s))) return(NULL)
  data.frame(
    k = as.character(s$karyotype),
    estimated_fitness = as.numeric(s[[est_col]]),
    estimated_sd = as.numeric(s[[sd_col]]),
    stringsAsFactors = FALSE
  )
}

row_field <- function(row, name, default = NA) {
  if (!name %in% names(row)) return(default)
  value <- row[[name]]
  if (is.null(value) || !length(value)) return(default)
  value[[1L]]
}

optional_key_col <- function(x, name, default = NA_character_) {
  if (name %in% names(x)) as.character(x[[name]]) else rep(default, nrow(x))
}

fit_tbl <- read_tsv(file.path(tables_dir, "fit_results.tsv"))
summary_tbl <- read_tsv(file.path(tables_dir, "summary_by_lambda_time_minobs_method.tsv"))

num_fit <- c("simulation_id", "lambda", "time_start", "time_gap", "time_delta", "minobs", "sim_pm", "pm")
for (nm in intersect(num_fit, names(fit_tbl))) fit_tbl[[nm]] <- to_num(fit_tbl[[nm]])
num_summary <- c(
  "lambda", "time_start", "time_gap", "time_delta", "minobs", "sim_pm", "pm", "n_nodes", "n_scored",
  "observed_node_fraction", "pearson", "spearman", "rmse", "mae", "signed_bias",
  "centered_rmse", "centered_mae", "centered_r2", "sign_accuracy", "false_high_rate",
  "mean_estimated_sd"
)
for (nm in intersect(num_summary, names(summary_tbl))) summary_tbl[[nm]] <- to_num(summary_tbl[[nm]])

grf_cache <- new.env(parent = emptyenv())
load_grf <- function(path) {
  path <- localize_path(path)
  if (!file.exists(path)) stop("Missing GRF RDS: ", path, call. = FALSE)
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(key, envir = grf_cache, inherits = FALSE)) {
    assign(key, readRDS(path), envir = grf_cache)
  }
  get(key, envir = grf_cache, inherits = FALSE)
}

alfakR_rows <- fit_tbl[fit_tbl$engine == "alfakR" & fit_tbl$status == "ok", , drop = FALSE]
alfak2_rows <- fit_tbl[
  fit_tbl$engine == "alfak2" &
    fit_tbl$method == "alfak2_effective_minobs_matched" &
    fit_tbl$status == "ok",
  ,
  drop = FALSE
]
condition_key <- function(x) {
  paste(
    x$simulation_id,
    x$lambda_label,
    x$time_start,
    x$time_gap,
    x$time_delta,
    x$minobs,
    optional_key_col(x, "sim_pm"),
    optional_key_col(x, "pm"),
    optional_key_col(x, "fit_beta_label"),
    sep = "||"
  )
}
alfak2_rows$condition_key <- condition_key(alfak2_rows)
alfak2_key_index <- split(seq_len(nrow(alfak2_rows)), alfak2_rows$condition_key)

alfak2_cache <- new.env(parent = emptyenv())
load_alfak2_nodes <- function(path) {
  path <- localize_path(path)
  if (!file.exists(path)) stop("Missing alfak2 fit RDS: ", path, call. = FALSE)
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(key, envir = alfak2_cache, inherits = FALSE)) {
    fit <- readRDS(path)
    nodes <- extract_alfak2_nodes(fit)
    assign(key, nodes, envir = alfak2_cache)
  }
  get(key, envir = alfak2_cache, inherits = FALSE)
}

condition_rows <- vector("list", nrow(alfakR_rows) * 4L)
idx <- 0L
missing_alfak2 <- 0L
for (i in seq_len(nrow(alfakR_rows))) {
  fr <- alfakR_rows[i, , drop = FALSE]
  key <- condition_key(fr)
  a2_idx <- alfak2_key_index[[key]]
  if (is.null(a2_idx) || length(a2_idx) < 1L) {
    missing_alfak2 <- missing_alfak2 + 1L
    next
  }
  a2_row <- alfak2_rows[a2_idx[[1L]], , drop = FALSE]
  a2_nodes <- load_alfak2_nodes(a2_row$fit_path[[1L]])
  if (is.null(a2_nodes) || !nrow(a2_nodes)) {
    missing_alfak2 <- missing_alfak2 + 1L
    next
  }
  boot_path <- localize_path(fr$bootstrap_path[[1L]])
  if (!file.exists(boot_path)) next
  boot <- readRDS(boot_path)
  grf <- load_grf(fr$grf_rds[[1L]])
  for (scope in c("direct", "nn")) {
    nodes <- extract_bootstrap_scope(boot, scope)
    if (is.null(nodes) || !nrow(nodes)) next
    truth <- compute_grf_fitness_truth(nodes$k, grf$centroids, as.numeric(fr$lambda[[1L]]))
    truth_vec <- as.numeric(truth[match(nodes$k, names(truth))])
    ar_metrics <- summarize_estimates(
      est = nodes$estimated_fitness,
      est_sd = nodes$estimated_sd,
      truth = truth_vec
    )
    idx <- idx + 1L
    condition_rows[[idx]] <- data.frame(
      simulation_id = fr$simulation_id,
      lambda = fr$lambda,
      lambda_label = fr$lambda_label,
      time_start = fr$time_start,
      time_gap = fr$time_gap,
      time_delta = fr$time_delta,
      minobs = fr$minobs,
      sim_pm = row_field(fr, "sim_pm", NA_real_),
      pm = row_field(fr, "pm", NA_real_),
      fit_beta_label = row_field(fr, "fit_beta_label", NA_character_),
      method = fr$method,
      engine = "alfakR",
      input_policy = fr$input_policy,
      nn_prior = fr$nn_prior,
      bootstrap_method = fr$method,
      bootstrap_nn_prior = fr$nn_prior,
      support_scope = scope,
      source = "bootstrap_res",
      ar_metrics,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    m <- match(nodes$k, a2_nodes$k)
    a2_metrics <- summarize_estimates(
      est = a2_nodes$estimated_fitness[m],
      est_sd = a2_nodes$estimated_sd[m],
      truth = truth_vec
    )
    idx <- idx + 1L
    condition_rows[[idx]] <- data.frame(
      simulation_id = a2_row$simulation_id,
      lambda = a2_row$lambda,
      lambda_label = a2_row$lambda_label,
      time_start = a2_row$time_start,
      time_gap = a2_row$time_gap,
      time_delta = a2_row$time_delta,
      minobs = a2_row$minobs,
      sim_pm = row_field(a2_row, "sim_pm", NA_real_),
      pm = row_field(a2_row, "pm", NA_real_),
      fit_beta_label = row_field(a2_row, "fit_beta_label", NA_character_),
      method = a2_row$method,
      engine = "alfak2",
      input_policy = a2_row$input_policy,
      nn_prior = "alfak2",
      bootstrap_method = fr$method,
      bootstrap_nn_prior = fr$nn_prior,
      support_scope = scope,
      source = "alfak2_global_on_bootstrap_nodes",
      a2_metrics,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}
raw_condition_metrics <- do.call(rbind, condition_rows[seq_len(idx)])
if (missing_alfak2 > 0L) {
  warning("Skipped ", missing_alfak2, " alfakR rows without a paired alfak2 fit.", call. = FALSE)
}
raw_condition_metrics$scope_label <- scope_label(raw_condition_metrics$support_scope)
raw_condition_metrics$method_label <- method_label(raw_condition_metrics$method)
raw_condition_metrics <- raw_condition_metrics[
  order(raw_condition_metrics$support_scope, raw_condition_metrics$method, raw_condition_metrics$minobs,
        raw_condition_metrics$lambda, raw_condition_metrics$time_gap),
  ,
  drop = FALSE
]

metric_cols <- c(
  "n_scored", "pearson", "spearman", "rmse", "mae", "centered_rmse", "centered_mae",
  "centered_r2", "sign_accuracy", "false_high_rate", "mean_estimated_sd"
)
summary_rows <- list()
idx <- 1L
for (metric in metric_cols) {
  s <- summarize_metric(
    raw_condition_metrics,
    c("minobs", intersect(c("sim_pm", "pm", "fit_beta_label"), names(raw_condition_metrics)),
      "support_scope", "scope_label", "method", "method_label", "nn_prior", "source"),
    metric
  )
  s$metric <- metric
  summary_rows[[idx]] <- s
  idx <- idx + 1L
}
performance_summary <- do.call(rbind, summary_rows)
performance_summary <- performance_summary[
  order(performance_summary$scope_label, performance_summary$metric, performance_summary$method_label,
        performance_summary$minobs),
  ,
  drop = FALSE
]

key_cols <- c(
  "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
  "minobs", intersect(c("sim_pm", "pm", "fit_beta_label"), names(raw_condition_metrics)),
  "support_scope", "bootstrap_method", "bootstrap_nn_prior"
)
a2 <- raw_condition_metrics[
  raw_condition_metrics$method == "alfak2_effective_minobs_matched",
  c(key_cols, metric_cols),
  drop = FALSE
]
names(a2)[match(metric_cols, names(a2))] <- paste0(metric_cols, "_alfak2")
ar <- raw_condition_metrics[
  raw_condition_metrics$engine == "alfakR",
  c(key_cols, "method", "method_label", "nn_prior", metric_cols),
  drop = FALSE
]
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

delta_metric_cols <- grep("^delta_", names(delta), value = TRUE)
delta_group_cols <- c("minobs", intersect(c("sim_pm", "pm", "fit_beta_label"), names(delta)),
                      "support_scope", "scope_label", "alfakR_method", "alfakR_label", "nn_prior")
delta_rows <- list()
idx <- 1L
for (metric in delta_metric_cols) {
  split_key <- interaction(delta[delta_group_cols], drop = TRUE, lex.order = TRUE)
  s <- do.call(rbind, lapply(split(delta, split_key), function(x) {
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
  delta_rows[[idx]] <- s
  idx <- idx + 1L
}
delta_summary <- do.call(rbind, delta_rows)
delta_summary <- delta_summary[
  order(delta_summary$scope_label, delta_summary$metric, delta_summary$alfakR_label,
        delta_summary$minobs),
  ,
  drop = FALSE
]

write_tsv(raw_condition_metrics, file.path(tables_dir, "raw_bootstrap_fq_nn_condition_metrics.tsv"))
write_tsv(performance_summary, file.path(tables_dir, "raw_bootstrap_fq_nn_performance_summary.tsv"))
write_tsv(delta, file.path(tables_dir, "raw_bootstrap_fq_nn_alfak2_vs_alfakR_condition_delta.tsv"))
write_tsv(delta_summary, file.path(tables_dir, "raw_bootstrap_fq_nn_alfak2_vs_alfakR_delta_summary.tsv"))

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
  scope_levels <- c("fq", "NN")
  x <- performance_summary[performance_summary$metric == "centered_rmse", , drop = FALSE]
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
      y = "raw fitness centered RMSE, median and IQR",
      color = NULL,
      fill = NULL,
      title = "Raw fq and NN fitness performance",
      subtitle = "Node sets are the fq/NN columns of each alfakR bootstrap_res.Rds; alfak2 is evaluated on the same nodes"
    ) +
    ggplot2::theme(legend.position = "bottom")
  fig_width <- if (length(unique(x$beta_label)) > 1L) 14 else 10
  ggplot2::ggsave(file.path(figures_dir, "raw_bootstrap_fq_nn_centered_rmse.png"), p1, width = fig_width, height = 6.2, dpi = 300)
  ggplot2::ggsave(file.path(figures_dir, "raw_bootstrap_fq_nn_centered_rmse.pdf"), p1, width = fig_width, height = 6.2)

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
      title = "alfak2 minus raw alfakR bootstrap fitness",
      subtitle = "Each paired delta uses the same bootstrap_res.Rds fq/NN node set; negative values favor alfak2"
    ) +
    ggplot2::theme(legend.position = "bottom")
  delta_fig_width <- if (length(unique(y$beta_label)) > 1L) 14 else 10
  ggplot2::ggsave(file.path(figures_dir, "raw_bootstrap_fq_nn_delta_centered_rmse.png"), p2, width = delta_fig_width, height = 6.2, dpi = 300)
  ggplot2::ggsave(file.path(figures_dir, "raw_bootstrap_fq_nn_delta_centered_rmse.pdf"), p2, width = delta_fig_width, height = 6.2)
}

message("Wrote raw bootstrap fq/NN summaries to: ", tables_dir)
message("Wrote figures to: ", figures_dir)
