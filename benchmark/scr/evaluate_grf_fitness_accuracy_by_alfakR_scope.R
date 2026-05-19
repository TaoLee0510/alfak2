#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) >= 1L) args[[1L]] else "benchmark/results/grf_downsampled_method_blind_full_vs_minobs_2hop_2acba0a_pm_5e_05"
if (basename(run_dir) == "benchmark") {
  benchmark_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  run_dir <- dirname(benchmark_dir)
} else {
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  benchmark_dir <- file.path(run_dir, "benchmark")
}
if (!dir.exists(benchmark_dir)) stop("Missing benchmark directory: ", benchmark_dir, call. = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("Package data.table is required.", call. = FALSE)

dt <- data.table::data.table
tables_dir <- file.path(benchmark_dir, "tables")
parts_dir <- file.path(tables_dir, "fit_results_parts")
out_prefix <- "extended_fitness_accuracy"

path_or_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

resolve_run_path <- function(x) {
  x <- path_or_na(x)
  out <- x
  run_name <- basename(run_dir)
  for (i in seq_along(out)) {
    xi <- out[[i]]
    if (is.na(xi) || file.exists(xi)) next
    pos <- regexpr(run_name, xi, fixed = TRUE)[[1L]]
    if (pos > 0L) {
      candidate <- file.path(dirname(run_dir), substring(xi, pos))
      if (file.exists(candidate)) out[[i]] <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
    }
  }
  out
}

read_fit_parts <- function() {
  files <- list.files(parts_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(files)) stop("No fit result part TSV files found in: ", parts_dir, call. = FALSE)
  rows <- lapply(files, function(path) {
    tryCatch(data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE, data.table = TRUE),
             error = function(e) NULL)
  })
  x <- data.table::rbindlist(rows, fill = TRUE)
  num_cols <- intersect(
    c("simulation_id", "lambda", "time_start", "time_gap", "time_delta", "minobs",
      "sim_pm", "pm", "grf_centroid_min_cn", "grf_centroid_max_cn",
      "anchor_count_reference", "anchor_count_power", "local_nodes", "global_nodes"),
    names(x)
  )
  for (nm in num_cols) x[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  for (nm in intersect(c("grf_rds", "input_rds", "fit_path", "landscape_path", "outdir"), names(x))) {
    x[, (nm) := resolve_run_path(get(nm))]
  }
  x[]
}

read_input_table <- function() {
  path <- file.path(tables_dir, "input_table.tsv")
  x <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE, data.table = TRUE)
  num_cols <- intersect(
    c("simulation_id", "lambda", "time_start", "time_gap", "time_delta", "minobs",
      "sim_pm", "raw_input_rows", "input_rows_after_drop", "input_rows_minobs"),
    names(x)
  )
  for (nm in num_cols) x[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  for (nm in intersect(c("grf_rds", "input_rds", "input_csv"), names(x))) {
    x[, (nm) := resolve_run_path(get(nm))]
  }
  x[]
}

parse_karyotype_ids <- function(ids) {
  pieces <- strsplit(as.character(ids), ".", fixed = TRUE)
  n_chr <- lengths(pieces)
  if (!length(ids) || any(n_chr == 0L) || length(unique(n_chr)) != 1L) {
    stop("Malformed karyotype IDs.", call. = FALSE)
  }
  mat <- do.call(rbind, lapply(pieces, as.integer))
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

is_diploid_label <- function(ids) {
  if (!length(ids)) return(logical())
  mat <- parse_karyotype_ids(ids)
  rowSums(mat != 2L, na.rm = TRUE) == 0L
}

compute_grf_fitness_truth <- function(karyotypes, grf) {
  ids <- as.character(karyotypes)
  if (!length(ids)) return(stats::setNames(numeric(), character()))
  k_mat <- parse_karyotype_ids(ids)
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
  stats::setNames(acc / (pi * sqrt(nrow(centroids))), rownames(k_mat))
}

grf_cache <- new.env(parent = emptyenv())
truth_cache <- new.env(parent = emptyenv())
true_node_cache <- new.env(parent = emptyenv())

load_grf <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(key, envir = grf_cache, inherits = FALSE)) assign(key, readRDS(key), envir = grf_cache)
  get(key, envir = grf_cache, inherits = FALSE)
}

truth_for_k <- function(k, grf_path) {
  k <- as.character(k)
  grf_key <- normalizePath(grf_path, winslash = "/", mustWork = TRUE)
  if (!exists(grf_key, envir = truth_cache, inherits = FALSE)) {
    assign(grf_key, new.env(parent = emptyenv(), hash = TRUE), envir = truth_cache)
  }
  env <- get(grf_key, envir = truth_cache, inherits = FALSE)
  current <- unlist(mget(k, envir = env, ifnotfound = as.list(rep(NA_real_, length(k))), inherits = FALSE), use.names = FALSE)
  missing <- unique(k[!is.finite(current)])
  if (length(missing)) {
    vals <- compute_grf_fitness_truth(missing, load_grf(grf_key))
    list2env(as.list(vals), envir = env)
    current <- unlist(mget(k, envir = env, ifnotfound = as.list(rep(NA_real_, length(k))), inherits = FALSE), use.names = FALSE)
  }
  as.numeric(current)
}

true_sim_nodes_for_grf <- function(grf_path, drop_diploid = TRUE) {
  key <- normalizePath(grf_path, winslash = "/", mustWork = TRUE)
  cache_key <- paste(key, drop_diploid, sep = "\r")
  if (!exists(cache_key, envir = true_node_cache, inherits = FALSE)) {
    grf <- load_grf(key)
    nodes <- setdiff(names(grf$sim_wide), "time")
    if (isTRUE(drop_diploid) && length(nodes)) nodes <- nodes[!is_diploid_label(nodes)]
    assign(cache_key, unique(nodes), envir = true_node_cache)
  }
  get(cache_key, envir = true_node_cache, inherits = FALSE)
}

scope_from_alfakR <- function(fq, nn) {
  data.table::fifelse(fq %in% TRUE, "fq", data.table::fifelse(nn %in% TRUE, "NN", "other"))
}

scope_from_alfak2 <- function(support_tier) {
  tier <- as.character(support_tier)
  data.table::fifelse(tier == "directly_informed", "fq",
    data.table::fifelse(tier %in% c("local_borrowed", "weakly_supported", "graph_borrowed", "prior_dominated"), "NN", "other")
  )
}

read_landscape_nodes <- function(fr) {
  x <- data.table::as.data.table(readRDS(fr$landscape_path[[1L]]))
  if (identical(fr$engine[[1L]], "alfak2")) {
    if (!"karyotype" %in% names(x)) return(dt())
    est_col <- if ("fitness_mean_alfakR_scale" %in% names(x)) "fitness_mean_alfakR_scale" else "fitness_mean"
    sd_col <- if ("fitness_sd_alfakR_scale" %in% names(x)) "fitness_sd_alfakR_scale" else "fitness_sd"
    out <- x[, .(
      k = as.character(karyotype),
      estimated_fitness = as.numeric(get(est_col)),
      estimated_sd = if (sd_col %in% names(x)) as.numeric(get(sd_col)) else NA_real_,
      own_scope = scope_from_alfak2(if ("support_tier" %in% names(x)) support_tier else NA_character_),
      support_tier = if ("support_tier" %in% names(x)) as.character(support_tier) else NA_character_
    )]
    out[, `:=`(fq = own_scope == "fq", nn = own_scope == "NN")]
  } else {
    if (!"k" %in% names(x) && "karyotype" %in% names(x)) x[, k := as.character(karyotype)]
    if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x[, mean := as.numeric(fitness_mean)]
    if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x[, sd := as.numeric(fitness_sd)]
    if (!"fq" %in% names(x)) x[, fq := FALSE]
    if (!"nn" %in% names(x)) x[, nn := FALSE]
    out <- x[, .(
      k = as.character(k),
      estimated_fitness = as.numeric(mean),
      estimated_sd = if ("sd" %in% names(x)) as.numeric(sd) else NA_real_,
      own_scope = scope_from_alfakR(as.logical(fq), as.logical(nn)),
      fq = as.logical(fq),
      nn = as.logical(nn)
    )]
    out[, support_tier := own_scope]
  }
  out <- out[is.finite(estimated_fitness)]
  out[, true_fitness := truth_for_k(k, fr$grf_rds[[1L]])]
  out[is.finite(true_fitness)]
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

metric_row <- function(x, scope_col = "own_scope", scope = "whole", true_nodes = character()) {
  y <- if (identical(scope, "whole")) x else x[get(scope_col) == scope]
  if (!nrow(y)) return(NULL)
  ok <- is.finite(y$estimated_fitness) & is.finite(y$true_fitness)
  y <- y[ok]
  if (!nrow(y)) return(NULL)
  err <- y$estimated_fitness - y$true_fitness
  pred_c <- y$estimated_fitness - mean(y$estimated_fitness)
  truth_c <- y$true_fitness - mean(y$true_fitness)
  centered_err <- pred_c - truth_c
  truth_sd <- stats::sd(y$true_fitness)
  pred_sd <- stats::sd(y$estimated_fitness)
  slope <- if (is.finite(truth_sd) && truth_sd > 0) {
    stats::cov(y$estimated_fitness, y$true_fitness) / stats::var(y$true_fitness)
  } else {
    NA_real_
  }
  sd_ok <- is.finite(y$estimated_sd) & y$estimated_sd >= 0
  ci95 <- if (any(sd_ok)) mean(abs(err[sd_ok]) <= stats::qnorm(0.975) * y$estimated_sd[sd_ok]) else NA_real_
  ci99 <- if (any(sd_ok)) mean(abs(err[sd_ok]) <= stats::qnorm(0.995) * y$estimated_sd[sd_ok]) else NA_real_
  dt(
    support_scope = scope,
    n_estimated = nrow(y),
    n_ci_scored = sum(sd_ok),
    n_in_true_sim = sum(y$k %in% true_nodes),
    estimated_nodes_in_true_sim_fraction = if (nrow(y)) sum(y$k %in% true_nodes) / nrow(y) else NA_real_,
    true_sim_coverage_fraction = if (length(true_nodes)) sum(unique(y$k) %in% true_nodes) / length(true_nodes) else NA_real_,
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2)),
    signed_bias = mean(err),
    median_error = stats::median(err),
    centered_mae = mean(abs(centered_err)),
    centered_rmse = sqrt(mean(centered_err^2)),
    pearson = safe_cor(y$estimated_fitness, y$true_fitness, "pearson"),
    spearman = safe_cor(y$estimated_fitness, y$true_fitness, "spearman"),
    centered_r2 = safe_centered_r2(y$estimated_fitness, y$true_fitness),
    calibration_slope = slope,
    calibration_intercept = if (is.finite(slope)) mean(y$estimated_fitness) - slope * mean(y$true_fitness) else NA_real_,
    estimate_sd_ratio = if (is.finite(truth_sd) && truth_sd > 0) pred_sd / truth_sd else NA_real_,
    ci95_coverage = ci95,
    ci99_coverage = ci99,
    mean_estimated_sd = mean(y$estimated_sd, na.rm = TRUE)
  )
}

paired_metric_row <- function(common, scope = "whole") {
  y <- if (identical(scope, "whole")) common else common[alfakR_scope == scope]
  if (!nrow(y)) return(NULL)
  a2 <- y[, .(k, estimated_fitness = alfak2_estimated_fitness, estimated_sd = alfak2_estimated_sd, true_fitness)]
  ar <- y[, .(k, estimated_fitness = alfakR_estimated_fitness, estimated_sd = alfakR_estimated_sd, true_fitness)]
  a2m <- metric_row(a2[, own_scope := scope], "own_scope", scope, true_nodes = character())
  arm <- metric_row(ar[, own_scope := scope], "own_scope", scope, true_nodes = character())
  if (is.null(a2m) || is.null(arm)) return(NULL)
  a2_err <- a2$estimated_fitness - a2$true_fitness
  ar_err <- ar$estimated_fitness - ar$true_fitness
  a2_cerr <- (a2$estimated_fitness - mean(a2$estimated_fitness)) - (a2$true_fitness - mean(a2$true_fitness))
  ar_cerr <- (ar$estimated_fitness - mean(ar$estimated_fitness)) - (ar$true_fitness - mean(ar$true_fitness))
  dt(
    support_scope = scope,
    n_common = nrow(y),
    n_fq = sum(y$alfakR_scope == "fq"),
    n_NN = sum(y$alfakR_scope == "NN"),
    n_other = sum(y$alfakR_scope == "other"),
    alfak2_mae = a2m$mae,
    alfakR_mae = arm$mae,
    delta_mae = a2m$mae - arm$mae,
    alfak2_rmse = a2m$rmse,
    alfakR_rmse = arm$rmse,
    delta_rmse = a2m$rmse - arm$rmse,
    alfak2_centered_rmse = a2m$centered_rmse,
    alfakR_centered_rmse = arm$centered_rmse,
    delta_centered_rmse = a2m$centered_rmse - arm$centered_rmse,
    alfak2_signed_bias = a2m$signed_bias,
    alfakR_signed_bias = arm$signed_bias,
    delta_abs_bias = abs(a2m$signed_bias) - abs(arm$signed_bias),
    alfak2_pearson = a2m$pearson,
    alfakR_pearson = arm$pearson,
    delta_pearson = a2m$pearson - arm$pearson,
    alfak2_spearman = a2m$spearman,
    alfakR_spearman = arm$spearman,
    delta_spearman = a2m$spearman - arm$spearman,
    alfak2_centered_r2 = a2m$centered_r2,
    alfakR_centered_r2 = arm$centered_r2,
    delta_centered_r2 = a2m$centered_r2 - arm$centered_r2,
    alfak2_calibration_slope = a2m$calibration_slope,
    alfakR_calibration_slope = arm$calibration_slope,
    alfak2_estimate_sd_ratio = a2m$estimate_sd_ratio,
    alfakR_estimate_sd_ratio = arm$estimate_sd_ratio,
    alfak2_ci95_coverage = a2m$ci95_coverage,
    alfakR_ci95_coverage = arm$ci95_coverage,
    delta_ci95_coverage = a2m$ci95_coverage - arm$ci95_coverage,
    alfak2_ci99_coverage = a2m$ci99_coverage,
    alfakR_ci99_coverage = arm$ci99_coverage,
    delta_ci99_coverage = a2m$ci99_coverage - arm$ci99_coverage,
    alfak2_abs_better_rate = mean(abs(a2_err) < abs(ar_err), na.rm = TRUE),
    alfak2_centered_abs_better_rate = mean(abs(a2_cerr) < abs(ar_cerr), na.rm = TRUE),
    estimate_diff_alfak2_minus_alfakR = mean(a2$estimated_fitness - ar$estimated_fitness, na.rm = TRUE)
  )
}

method_label <- function(x) {
  out <- as.character(x)
  out[out == "alfak2_effective_full"] <- "alfak2 full"
  out[out == "alfak2_effective_minobs_matched"] <- "alfak2 minobs-matched"
  out[out == "alfakR_none"] <- "alfakR none"
  out[out == "alfakR_empirical"] <- "alfakR empirical"
  out[out == "alfakR_empirical_censored"] <- "alfakR censored"
  out[out == "alfakR_empirical_censored_weighted"] <- "alfakR censored weighted"
  out[out == "alfakR_empirical_two_step"] <- "alfakR two-step"
  out
}

summarize_numeric <- function(x, cols, by) {
  x[, {
    out <- list(n_conditions = .N)
    for (nm in cols) {
      v <- get(nm)
      out[[paste0(nm, "_median")]] <- as.numeric(stats::median(v, na.rm = TRUE))
      out[[paste0(nm, "_q25")]] <- as.numeric(stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE))
      out[[paste0(nm, "_q75")]] <- as.numeric(stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE))
    }
    out
  }, by = by]
}

fit_tbl <- read_fit_parts()
input_tbl <- read_input_table()
fit_tbl <- fit_tbl[status == "ok" & file.exists(landscape_path) & file.exists(grf_rds)]

key_cols <- c("simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
              "minobs", "sim_pm", "pm", "fit_beta_label", "input_md5", "grf_key")
input_key_cols <- c("simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
                    "minobs", "sim_pm", "input_md5", "grf_key")

input_cov <- unique(input_tbl[, c(input_key_cols, "raw_input_rows", "input_rows_after_drop", "input_rows_minobs", "grf_rds"), with = FALSE])
input_cov[, true_sim_nodes := vapply(grf_rds, function(p) length(true_sim_nodes_for_grf(p, drop_diploid = TRUE)), numeric(1))]
input_cov[, `:=`(
  sampled_input_fraction_of_true_sim = input_rows_after_drop / true_sim_nodes,
  minobs_input_fraction_of_true_sim = input_rows_minobs / true_sim_nodes
)]

node_cache <- new.env(parent = emptyenv())
read_nodes_cached <- function(task_order, fr) {
  key <- as.character(task_order)
  if (!exists(key, envir = node_cache, inherits = FALSE)) {
    assign(key, read_landscape_nodes(fr), envir = node_cache)
  }
  get(key, envir = node_cache, inherits = FALSE)
}

method_rows <- list()
pair_rows <- list()
method_idx <- 0L
pair_idx <- 0L

group_key <- interaction(fit_tbl[, ..key_cols], drop = TRUE, lex.order = TRUE)
levels_key <- levels(group_key)
message("Processing ", length(levels_key), " benchmark conditions and ", nrow(fit_tbl), " successful fits.")
for (g in seq_along(levels_key)) {
  if (g %% 25L == 0L) message("  processed conditions: ", g, " / ", length(levels_key))
  rows <- fit_tbl[group_key == levels_key[[g]]]
  if (!nrow(rows)) next
  true_nodes <- true_sim_nodes_for_grf(rows$grf_rds[[1L]], drop_diploid = TRUE)
  cov_row <- input_cov[
    simulation_id == rows$simulation_id[[1L]] &
      lambda == rows$lambda[[1L]] &
      time_start == rows$time_start[[1L]] &
      time_gap == rows$time_gap[[1L]] &
      minobs == rows$minobs[[1L]] &
      input_md5 == rows$input_md5[[1L]]
  ][1L]

  for (i in seq_len(nrow(rows))) {
    fr <- rows[i]
    nodes <- read_nodes_cached(fr$task_order[[1L]], fr)
    if (!nrow(nodes)) next
    for (scope in c("whole", "fq", "NN", "other")) {
      mr <- metric_row(nodes, "own_scope", scope, true_nodes)
      if (is.null(mr)) next
      method_idx <- method_idx + 1L
      method_rows[[method_idx]] <- cbind(
        fr[, c(key_cols, "engine", "method", "input_policy", "nn_prior", "task_order", "local_nodes", "global_nodes"), with = FALSE],
        mr,
        cov_row[, .(true_sim_nodes, raw_input_rows, input_rows_after_drop, input_rows_minobs,
                    sampled_input_fraction_of_true_sim, minobs_input_fraction_of_true_sim)]
      )
    }
  }

  a2_rows <- rows[engine == "alfak2"]
  ar_rows <- rows[engine == "alfakR"]
  if (!nrow(a2_rows) || !nrow(ar_rows)) next
  for (ia in seq_len(nrow(a2_rows))) {
    a2_fr <- a2_rows[ia]
    a2_nodes <- data.table::copy(read_nodes_cached(a2_fr$task_order[[1L]], a2_fr))
    if (!nrow(a2_nodes)) next
    data.table::setnames(
      a2_nodes,
      c("estimated_fitness", "estimated_sd", "own_scope", "support_tier"),
      c("alfak2_estimated_fitness", "alfak2_estimated_sd", "alfak2_own_scope", "alfak2_support_tier")
    )
    for (ir in seq_len(nrow(ar_rows))) {
      ar_fr <- ar_rows[ir]
      ar_nodes <- read_nodes_cached(ar_fr$task_order[[1L]], ar_fr)
      if (!nrow(ar_nodes)) next
      ar_nodes <- data.table::copy(ar_nodes)
      ar_nodes[, alfakR_scope := own_scope]
      data.table::setnames(
        ar_nodes,
        c("estimated_fitness", "estimated_sd", "support_tier"),
        c("alfakR_estimated_fitness", "alfakR_estimated_sd", "alfakR_support_tier")
      )
      common <- merge(
        ar_nodes[, .(k, alfakR_estimated_fitness, alfakR_estimated_sd, alfakR_scope, true_fitness)],
        a2_nodes[, .(k, alfak2_estimated_fitness, alfak2_estimated_sd, alfak2_own_scope)],
        by = "k", all = FALSE, sort = FALSE
      )
      if (!nrow(common)) next
      for (scope in c("whole", "fq", "NN", "other")) {
        pr <- paired_metric_row(common, scope)
        if (is.null(pr)) next
        pair_idx <- pair_idx + 1L
        pair_rows[[pair_idx]] <- cbind(
          ar_fr[, ..key_cols],
          dt(
            alfak2_method = a2_fr$method[[1L]],
            alfak2_label = method_label(a2_fr$method[[1L]]),
            alfak2_input_policy = a2_fr$input_policy[[1L]],
            alfakR_method = ar_fr$method[[1L]],
            alfakR_label = method_label(ar_fr$method[[1L]]),
            alfakR_nn_prior = ar_fr$nn_prior[[1L]]
          ),
          pr
        )
      }
    }
  }
  rm(list = ls(envir = node_cache), envir = node_cache)
}

method_metrics <- data.table::rbindlist(method_rows, fill = TRUE)
pair_metrics <- data.table::rbindlist(pair_rows, fill = TRUE)
if (!nrow(method_metrics)) stop("No method metrics generated.", call. = FALSE)
if (!nrow(pair_metrics)) stop("No paired metrics generated.", call. = FALSE)

method_metrics[, method_label := method_label(method)]
data.table::setcolorder(method_metrics, c(
  intersect(c(key_cols, "engine", "method", "method_label", "input_policy", "nn_prior", "support_scope"), names(method_metrics)),
  setdiff(names(method_metrics), c(key_cols, "engine", "method", "method_label", "input_policy", "nn_prior", "support_scope"))
))

method_summary_cols <- c(
  "n_estimated", "true_sim_coverage_fraction", "estimated_nodes_in_true_sim_fraction",
  "sampled_input_fraction_of_true_sim", "minobs_input_fraction_of_true_sim",
  "mae", "rmse", "signed_bias", "centered_rmse", "pearson", "spearman",
  "centered_r2", "calibration_slope", "estimate_sd_ratio", "ci95_coverage", "ci99_coverage"
)
method_summary <- summarize_numeric(
  method_metrics,
  intersect(method_summary_cols, names(method_metrics)),
  by = c("engine", "method", "method_label", "input_policy", "nn_prior", "minobs", "support_scope")
)
data.table::setorder(method_summary, support_scope, engine, method, minobs)

pair_summary_cols <- c(
  "n_common", "delta_mae", "delta_rmse", "delta_centered_rmse", "delta_abs_bias",
  "delta_pearson", "delta_spearman", "delta_centered_r2", "delta_ci95_coverage",
  "delta_ci99_coverage", "alfak2_abs_better_rate", "alfak2_centered_abs_better_rate",
  "alfak2_mae", "alfakR_mae", "alfak2_centered_rmse", "alfakR_centered_rmse",
  "alfak2_signed_bias", "alfakR_signed_bias", "alfak2_pearson", "alfakR_pearson",
  "alfak2_spearman", "alfakR_spearman", "alfak2_ci95_coverage", "alfakR_ci95_coverage",
  "alfak2_ci99_coverage", "alfakR_ci99_coverage"
)
pair_by <- c("alfak2_method", "alfak2_label", "alfak2_input_policy",
             "alfakR_method", "alfakR_label", "alfakR_nn_prior", "minobs", "support_scope")
pair_summary <- summarize_numeric(
  pair_metrics,
  intersect(pair_summary_cols, names(pair_metrics)),
  by = pair_by
)
pair_win <- pair_metrics[, .(
  mae_condition_win_rate = mean(delta_mae < 0, na.rm = TRUE),
  centered_rmse_condition_win_rate = mean(delta_centered_rmse < 0, na.rm = TRUE)
), by = pair_by]
pair_summary <- merge(pair_summary, pair_win, by = pair_by, all.x = TRUE, sort = FALSE)
data.table::setorder(pair_summary, support_scope, alfakR_label, alfak2_label, minobs)

data.table::fwrite(method_metrics, file.path(tables_dir, paste0(out_prefix, "_method_condition_metrics.tsv")), sep = "\t")
data.table::fwrite(method_summary, file.path(tables_dir, paste0(out_prefix, "_method_summary.tsv")), sep = "\t")
data.table::fwrite(pair_metrics, file.path(tables_dir, paste0(out_prefix, "_pair_condition_metrics_by_alfakR_scope.tsv")), sep = "\t")
data.table::fwrite(pair_summary, file.path(tables_dir, paste0(out_prefix, "_pair_summary_by_alfakR_scope.tsv")), sep = "\t")
data.table::fwrite(input_cov, file.path(tables_dir, paste0(out_prefix, "_input_coverage.tsv")), sep = "\t")

message("Wrote method condition metrics: ", file.path(tables_dir, paste0(out_prefix, "_method_condition_metrics.tsv")))
message("Wrote method summary: ", file.path(tables_dir, paste0(out_prefix, "_method_summary.tsv")))
message("Wrote paired condition metrics: ", file.path(tables_dir, paste0(out_prefix, "_pair_condition_metrics_by_alfakR_scope.tsv")))
message("Wrote paired summary: ", file.path(tables_dir, paste0(out_prefix, "_pair_summary_by_alfakR_scope.tsv")))
message("Wrote input coverage: ", file.path(tables_dir, paste0(out_prefix, "_input_coverage.tsv")))
