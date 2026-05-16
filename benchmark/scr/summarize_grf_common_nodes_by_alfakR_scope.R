#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "benchmark/results/grf_downsampled_input_benchmark_pm_5e_05/benchmark"
}
if (!dir.exists(results_dir)) {
  stop("Results directory does not exist: ", results_dir, call. = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required.", call. = FALSE)
}

dt <- data.table::data.table
tables_dir <- file.path(results_dir, "tables")
fit_path <- file.path(tables_dir, "fit_results.tsv")
if (!file.exists(fit_path)) {
  stop("Missing fit_results.tsv: ", fit_path, call. = FALSE)
}

point_path <- file.path(tables_dir, "common_node_pointwise_by_alfakR_scope.tsv.gz")
condition_path <- file.path(tables_dir, "common_node_condition_metrics_by_alfakR_scope.tsv")
summary_path <- file.path(tables_dir, "common_node_performance_summary_by_alfakR_scope.tsv")

read_tsv <- function(path) {
  data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE, data.table = TRUE)
}

scope_from_alfakR <- function(fq, nn) {
  data.table::fifelse(fq %in% TRUE, "fq", data.table::fifelse(nn %in% TRUE, "NN", "other"))
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

compute_grf_fitness_truth <- function(karyotypes, grf) {
  k_mat <- parse_karyotype_ids_base(karyotypes)
  centroids <- grf$centroids
  lambda <- as.numeric(grf$lambda)
  if (!is.matrix(centroids) || !is.numeric(centroids) || !nrow(centroids)) {
    stop("GRF centroids must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (!is.finite(lambda) || lambda <= 0) {
    stop("GRF lambda must be positive and finite.", call. = FALSE)
  }
  if (ncol(k_mat) != ncol(centroids)) {
    stop("Karyotype dimension does not match GRF centroid dimension.", call. = FALSE)
  }
  k_mat <- matrix(as.numeric(k_mat), nrow = nrow(k_mat), dimnames = dimnames(k_mat))
  acc <- numeric(nrow(k_mat))
  for (j in seq_len(nrow(centroids))) {
    diffs <- sweep(k_mat, 2L, centroids[j, ], FUN = "-")
    acc <- acc + sin(sqrt(rowSums(diffs^2)) / lambda)
  }
  out <- acc / (pi * sqrt(nrow(centroids)))
  names(out) <- rownames(k_mat)
  out
}

grf_cache <- new.env(parent = emptyenv())
truth_cache <- new.env(parent = emptyenv())

load_grf <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(key, envir = grf_cache, inherits = FALSE)) {
    assign(key, readRDS(key), envir = grf_cache)
  }
  get(key, envir = grf_cache, inherits = FALSE)
}

truth_for_k <- function(k, grf_path) {
  grf_key <- normalizePath(grf_path, winslash = "/", mustWork = TRUE)
  if (!exists(grf_key, envir = truth_cache, inherits = FALSE)) {
    assign(grf_key, new.env(parent = emptyenv(), hash = TRUE), envir = truth_cache)
  }
  env <- get(grf_key, envir = truth_cache, inherits = FALSE)
  current <- unlist(mget(k, envir = env, ifnotfound = as.list(rep(NA_real_, length(k))), inherits = FALSE), use.names = FALSE)
  missing <- unique(k[!is.finite(current)])
  if (length(missing)) {
    grf <- load_grf(grf_key)
    vals <- compute_grf_fitness_truth(missing, grf)
    list2env(as.list(vals), envir = env)
    current <- unlist(mget(k, envir = env, ifnotfound = as.list(rep(NA_real_, length(k))), inherits = FALSE), use.names = FALSE)
  }
  as.numeric(current)
}

read_alfak2_nodes <- function(path) {
  x <- readRDS(path)
  x <- data.table::as.data.table(x)
  if (!"karyotype" %in% names(x)) return(dt())
  est_col <- if ("fitness_mean_alfakR_scale" %in% names(x)) "fitness_mean_alfakR_scale" else "fitness_mean"
  sd_col <- if ("fitness_sd_alfakR_scale" %in% names(x)) "fitness_sd_alfakR_scale" else "fitness_sd"
  out <- x[, .(
    k = as.character(karyotype),
    alfak2_estimated_fitness = as.numeric(get(est_col)),
    alfak2_estimated_sd = if (sd_col %in% names(x)) as.numeric(get(sd_col)) else NA_real_,
    alfak2_support_tier = if ("support_tier" %in% names(x)) as.character(support_tier) else NA_character_
  )]
  out[is.finite(alfak2_estimated_fitness)]
}

read_alfakR_nodes <- function(path) {
  x <- readRDS(path)
  x <- data.table::as.data.table(x)
  if (!"k" %in% names(x) && "karyotype" %in% names(x)) x[, k := as.character(karyotype)]
  if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x[, mean := as.numeric(fitness_mean)]
  if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x[, sd := as.numeric(fitness_sd)]
  if (!"fq" %in% names(x)) x[, fq := FALSE]
  if (!"nn" %in% names(x)) x[, nn := FALSE]
  out <- x[, .(
    k = as.character(k),
    alfakR_estimated_fitness = as.numeric(mean),
    alfakR_estimated_sd = if ("sd" %in% names(x)) as.numeric(sd) else NA_real_,
    alfakR_fq = as.logical(fq),
    alfakR_nn = as.logical(nn)
  )]
  out[, alfakR_scope := scope_from_alfakR(alfakR_fq, alfakR_nn)]
  out[is.finite(alfakR_estimated_fitness)]
}

metric_row <- function(x, scope) {
  if (identical(scope, "whole")) {
    y <- x
  } else {
    y <- x[alfakR_scope == scope]
  }
  if (!nrow(y)) return(NULL)
  ok <- is.finite(y$alfak2_estimated_fitness) &
    is.finite(y$alfakR_estimated_fitness) &
    is.finite(y$true_fitness)
  y <- y[ok]
  if (!nrow(y)) return(NULL)
  a2_err <- y$alfak2_estimated_fitness - y$true_fitness
  ar_err <- y$alfakR_estimated_fitness - y$true_fitness
  truth_c <- y$true_fitness - mean(y$true_fitness)
  a2_c <- y$alfak2_estimated_fitness - mean(y$alfak2_estimated_fitness)
  ar_c <- y$alfakR_estimated_fitness - mean(y$alfakR_estimated_fitness)
  a2_cerr <- a2_c - truth_c
  ar_cerr <- ar_c - truth_c
  dt(
    support_scope = scope,
    n_common = nrow(y),
    n_fq = sum(y$alfakR_scope == "fq", na.rm = TRUE),
    n_nn = sum(y$alfakR_scope == "NN", na.rm = TRUE),
    n_other = sum(y$alfakR_scope == "other", na.rm = TRUE),
    alfak2_mae = mean(abs(a2_err)),
    alfakR_mae = mean(abs(ar_err)),
    delta_mae = mean(abs(a2_err)) - mean(abs(ar_err)),
    alfak2_rmse = sqrt(mean(a2_err^2)),
    alfakR_rmse = sqrt(mean(ar_err^2)),
    delta_rmse = sqrt(mean(a2_err^2)) - sqrt(mean(ar_err^2)),
    alfak2_centered_mae = mean(abs(a2_cerr)),
    alfakR_centered_mae = mean(abs(ar_cerr)),
    delta_centered_mae = mean(abs(a2_cerr)) - mean(abs(ar_cerr)),
    alfak2_centered_rmse = sqrt(mean(a2_cerr^2)),
    alfakR_centered_rmse = sqrt(mean(ar_cerr^2)),
    delta_centered_rmse = sqrt(mean(a2_cerr^2)) - sqrt(mean(ar_cerr^2)),
    alfak2_signed_bias = mean(a2_err),
    alfakR_signed_bias = mean(ar_err),
    alfak2_abs_better_n = sum(abs(a2_err) < abs(ar_err), na.rm = TRUE),
    alfakR_abs_better_n = sum(abs(ar_err) < abs(a2_err), na.rm = TRUE),
    abs_tie_n = sum(abs(a2_err) == abs(ar_err), na.rm = TRUE),
    alfak2_abs_better_rate = mean(abs(a2_err) < abs(ar_err), na.rm = TRUE),
    alfakR_abs_better_rate = mean(abs(ar_err) < abs(a2_err), na.rm = TRUE),
    alfak2_centered_abs_better_rate = mean(abs(a2_cerr) < abs(ar_cerr), na.rm = TRUE),
    alfakR_centered_abs_better_rate = mean(abs(ar_cerr) < abs(a2_cerr), na.rm = TRUE),
    alfak2_pearson = safe_cor(y$alfak2_estimated_fitness, y$true_fitness, "pearson"),
    alfakR_pearson = safe_cor(y$alfakR_estimated_fitness, y$true_fitness, "pearson"),
    alfak2_spearman = safe_cor(y$alfak2_estimated_fitness, y$true_fitness, "spearman"),
    alfakR_spearman = safe_cor(y$alfakR_estimated_fitness, y$true_fitness, "spearman"),
    alfak2_centered_r2 = safe_centered_r2(y$alfak2_estimated_fitness, y$true_fitness),
    alfakR_centered_r2 = safe_centered_r2(y$alfakR_estimated_fitness, y$true_fitness)
  )
}

method_label <- function(x) {
  out <- as.character(x)
  out[out == "alfakR_none"] <- "alfakR none"
  out[out == "alfakR_empirical"] <- "alfakR empirical"
  out[out == "alfakR_empirical_censored"] <- "alfakR censored"
  out[out == "alfakR_empirical_censored_weighted"] <- "alfakR censored weighted"
  out[out == "alfakR_empirical_two_step"] <- "alfakR two-step"
  out
}

fit_tbl <- read_tsv(fit_path)
for (nm in intersect(c("simulation_id", "lambda", "time_start", "time_gap", "time_delta", "minobs", "sim_pm", "pm"), names(fit_tbl))) {
  fit_tbl[, (nm) := suppressWarnings(as.numeric(get(nm)))]
}

key_cols <- c(
  "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
  "minobs", "sim_pm", "pm", "fit_beta_label", "input_md5", "grf_key"
)
required <- c(key_cols, "engine", "method", "status", "landscape_path", "grf_rds")
missing <- setdiff(required, names(fit_tbl))
if (length(missing)) {
  stop("fit_results.tsv is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

a2 <- fit_tbl[
  engine == "alfak2" & method == "alfak2_effective_minobs_matched" & status == "ok",
  c(key_cols, "landscape_path", "grf_rds", "task_order"),
  with = FALSE
]
data.table::setnames(a2, c("landscape_path", "task_order"), c("alfak2_landscape_path", "alfak2_task_order"))

ar <- fit_tbl[
  engine == "alfakR" & status == "ok",
  c(key_cols, "method", "nn_prior", "landscape_path", "task_order"),
  with = FALSE
]
data.table::setnames(ar, c("method", "nn_prior", "landscape_path", "task_order"),
                     c("alfakR_method", "alfakR_nn_prior", "alfakR_landscape_path", "alfakR_task_order"))

pairs <- merge(ar, a2, by = key_cols, allow.cartesian = TRUE, all = FALSE, sort = FALSE)
pairs <- pairs[file.exists(alfak2_landscape_path) & file.exists(alfakR_landscape_path) & file.exists(grf_rds)]
if (!nrow(pairs)) {
  stop("No paired alfak2/alfakR landscape files were found.", call. = FALSE)
}

if (file.exists(point_path)) invisible(file.remove(point_path))
condition_rows <- vector("list", nrow(pairs) * 4L)
condition_idx <- 0L
point_header <- TRUE

message("Common-node paired landscape comparisons: ", nrow(pairs), " alfak2/alfakR fit pairs.")
for (i in seq_len(nrow(pairs))) {
  pr <- pairs[i]
  if (i %% 50L == 0L) {
    message("  processed pairs: ", i, " / ", nrow(pairs))
  }
  a2_nodes <- read_alfak2_nodes(pr$alfak2_landscape_path[[1L]])
  ar_nodes <- read_alfakR_nodes(pr$alfakR_landscape_path[[1L]])
  if (!nrow(a2_nodes) || !nrow(ar_nodes)) next
  common <- merge(ar_nodes, a2_nodes, by = "k", all = FALSE, sort = FALSE)
  if (!nrow(common)) next
  common[, true_fitness := truth_for_k(k, pr$grf_rds[[1L]])]
  common <- common[is.finite(true_fitness)]
  if (!nrow(common)) next
  common[, `:=`(
    alfak2_error = alfak2_estimated_fitness - true_fitness,
    alfakR_error = alfakR_estimated_fitness - true_fitness
  )]
  common[, `:=`(
    alfak2_abs_error = abs(alfak2_error),
    alfakR_abs_error = abs(alfakR_error),
    delta_abs_error = abs(alfak2_error) - abs(alfakR_error),
    pointwise_abs_winner = data.table::fifelse(
      abs(alfak2_error) < abs(alfakR_error),
      "alfak2",
      data.table::fifelse(abs(alfakR_error) < abs(alfak2_error), "alfakR", "tie")
    ),
    estimate_diff_alfak2_minus_alfakR = alfak2_estimated_fitness - alfakR_estimated_fitness
  )]

  meta <- pr[, c(key_cols, "alfakR_method", "alfakR_nn_prior", "alfak2_task_order", "alfakR_task_order"), with = FALSE]
  common_out <- cbind(meta[rep(1L, nrow(common))], common)
  data.table::fwrite(
    common_out,
    file = point_path,
    sep = "\t",
    append = !point_header,
    col.names = point_header
  )
  point_header <- FALSE

  for (scope in c("whole", "fq", "NN", "other")) {
    row <- metric_row(common, scope)
    if (is.null(row)) next
    condition_idx <- condition_idx + 1L
    condition_rows[[condition_idx]] <- cbind(meta, row)
  }
}

condition_metrics <- data.table::rbindlist(condition_rows[seq_len(condition_idx)], fill = TRUE)
if (!nrow(condition_metrics)) {
  stop("No common-node condition metrics were generated.", call. = FALSE)
}
condition_metrics[, alfakR_label := method_label(alfakR_method)]
data.table::setcolorder(
  condition_metrics,
  c(intersect(c(key_cols, "alfakR_method", "alfakR_label", "alfakR_nn_prior", "support_scope"), names(condition_metrics)),
    setdiff(names(condition_metrics), c(key_cols, "alfakR_method", "alfakR_label", "alfakR_nn_prior", "support_scope")))
)
data.table::fwrite(condition_metrics, condition_path, sep = "\t")

metric_cols <- c(
  "n_common", "alfak2_mae", "alfakR_mae", "delta_mae",
  "alfak2_rmse", "alfakR_rmse", "delta_rmse",
  "alfak2_centered_mae", "alfakR_centered_mae", "delta_centered_mae",
  "alfak2_centered_rmse", "alfakR_centered_rmse", "delta_centered_rmse",
  "alfak2_signed_bias", "alfakR_signed_bias",
  "alfak2_abs_better_rate", "alfakR_abs_better_rate",
  "alfak2_centered_abs_better_rate", "alfakR_centered_abs_better_rate",
  "alfak2_pearson", "alfakR_pearson", "alfak2_spearman", "alfakR_spearman",
  "alfak2_centered_r2", "alfakR_centered_r2"
)
summary <- condition_metrics[, {
  out <- list(
    n_conditions = as.numeric(.N),
    total_common_nodes = as.numeric(sum(n_common, na.rm = TRUE)),
    median_common_nodes = as.numeric(stats::median(n_common, na.rm = TRUE)),
    pooled_alfak2_abs_better_rate = as.numeric(sum(alfak2_abs_better_n, na.rm = TRUE) / sum(n_common, na.rm = TRUE)),
    pooled_alfakR_abs_better_rate = as.numeric(sum(alfakR_abs_better_n, na.rm = TRUE) / sum(n_common, na.rm = TRUE)),
    pooled_abs_tie_rate = as.numeric(sum(abs_tie_n, na.rm = TRUE) / sum(n_common, na.rm = TRUE))
  )
  for (metric in metric_cols) {
    values <- get(metric)
    out[[paste0(metric, "_median")]] <- as.numeric(stats::median(values, na.rm = TRUE))
    out[[paste0(metric, "_q25")]] <- as.numeric(stats::quantile(values, 0.25, na.rm = TRUE, names = FALSE))
    out[[paste0(metric, "_q75")]] <- as.numeric(stats::quantile(values, 0.75, na.rm = TRUE, names = FALSE))
  }
  out$mae_condition_win_rate <- as.numeric(mean(delta_mae < 0, na.rm = TRUE))
  out$rmse_condition_win_rate <- as.numeric(mean(delta_rmse < 0, na.rm = TRUE))
  out$centered_rmse_condition_win_rate <- as.numeric(mean(delta_centered_rmse < 0, na.rm = TRUE))
  out$centered_mae_condition_win_rate <- as.numeric(mean(delta_centered_mae < 0, na.rm = TRUE))
  out
}, by = .(minobs, sim_pm, pm, fit_beta_label, support_scope, alfakR_method, alfakR_label, alfakR_nn_prior)]
data.table::setorder(summary, support_scope, alfakR_label, minobs)
data.table::fwrite(summary, summary_path, sep = "\t")

message("Wrote pointwise common-node table: ", point_path)
message("Wrote condition metrics: ", condition_path)
message("Wrote performance summary: ", summary_path)
