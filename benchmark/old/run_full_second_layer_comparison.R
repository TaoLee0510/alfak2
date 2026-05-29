#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(full = FALSE, quick = FALSE, resume = FALSE, n_workers = 1L)
  for (arg in args) {
    if (identical(arg, "--full")) out$full <- TRUE
    else if (identical(arg, "--quick")) out$quick <- TRUE
    else if (identical(arg, "--resume")) out$resume <- TRUE
    else if (grepl("^--n-workers=", arg)) out$n_workers <- as.integer(sub("^--n-workers=", "", arg))
    else if (identical(arg, "--n-workers")) out$need_n_workers <- TRUE
    else if (isTRUE(out$need_n_workers)) {
      out$n_workers <- as.integer(arg)
      out$need_n_workers <- FALSE
    } else if (grepl("^--output-dir=", arg)) out$output_dir <- sub("^--output-dir=", "", arg)
    else if (identical(arg, "--output-dir")) out$need_output_dir <- TRUE
    else if (isTRUE(out$need_output_dir)) {
      out$output_dir <- arg
      out$need_output_dir <- FALSE
    }
    else if (grepl("^--alfakR-repo=", arg)) out$alfakR_repo <- sub("^--alfakR-repo=", "", arg)
    else if (identical(arg, "--help") || identical(arg, "-h")) out$help <- TRUE
    else stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (isTRUE(out$help)) {
    cat(
      "Usage:\n",
      "  Rscript benchmark/run_full_second_layer_9method_balanced_comparison.R --quick [--resume] [--n-workers 1]\n",
      "  Rscript benchmark/run_full_second_layer_9method_balanced_comparison.R --full [--resume] [--n-workers 1]\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
  if (!isTRUE(out$full) && !isTRUE(out$quick)) out$quick <- TRUE
  if (isTRUE(out$full) && isTRUE(out$quick)) stop("Use only one of --full or --quick.", call. = FALSE)
  if (!is.finite(out$n_workers) || out$n_workers < 1L) out$n_workers <- 1L
  out
}

repo_root <- function() {
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "benchmark/run_full_second_layer_comparison.R"), ".."), mustWork = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

load_repos <- function(repo_dir, alfakR_repo) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package `pkgload` is required to run the benchmark from source.", call. = FALSE)
  }
  pkgload::load_all(repo_dir, quiet = TRUE)
  alfakR_status <- tryCatch({
    pkgload::load_all(alfakR_repo, quiet = TRUE)
    "loaded"
  }, error = function(e) paste("failed:", conditionMessage(e)))
  alfakR_status
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

rbind_fill <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  xs <- xs[vapply(xs, nrow, integer(1), USE.NAMES = FALSE) > 0L]
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

git_hash <- function(path) {
  tryCatch(system2("git", c("-C", path, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
}

package_version_safe <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

benchmark_graph_controls <- function(row) {
  n_chr <- if ("n_chr" %in% names(row) && is.finite(as.numeric(row$n_chr))) as.integer(row$n_chr) else 4L
  if (n_chr >= 20L) {
    return(list(
      local_shell_depth = 1L,
      global_extra_shell = 1L,
      eval_max_nodes = 1000000L,
      fit_max_nodes = 1000000L
    ))
  }
  list(
    local_shell_depth = 2L,
    global_extra_shell = 1L,
    eval_max_nodes = 5000L,
    fit_max_nodes = 5000L
  )
}

simulate_shared_input <- function(row, output_dir) {
  n_chr <- if ("n_chr" %in% names(row) && is.finite(as.numeric(row$n_chr))) {
    as.integer(row$n_chr)
  } else {
    4L
  }
  sample_depth <- if ("sample_depth" %in% names(row) && is.finite(as.numeric(row$sample_depth))) {
    as.integer(row$sample_depth)
  } else {
    600L
  }
  fit_repeat <- if ("fit_repeat" %in% names(row) && is.finite(as.numeric(row$fit_repeat))) {
    as.integer(row$fit_repeat)
  } else {
    1L
  }
  cache_label <- paste(
    as.character(row$landscape_id),
    paste0("chr", n_chr),
    paste0("depth", sample_depth),
    paste0("fitrep", fit_repeat),
    paste0("countseed", as.integer(row$count_seed)),
    sep = "_"
  )
  path <- file.path(output_dir, "landscape_cache", paste0(cache_label, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  landscape <- alfak2::simulate_grf_landscape(
    n_chr = n_chr,
    n_centroids = 12,
    lambda = as.numeric(row$grf_lambda),
    scale = 0.5,
    seed = as.integer(row$landscape_seed),
    min_cn = 1,
    max_cn = 4,
    centroid_min = -4,
    centroid_max = 8
  )
  sim <- alfak2::simulate_sparse_counts(
    landscape,
    beta = 0.01,
    dt = 1,
    n0 = sample_depth,
    n1 = sample_depth,
    detection_threshold = 1,
    seed = as.integer(row$count_seed),
    initial_population = 1500,
    initial_shell_depth = 1,
    time_step = 0.25,
    carrying_capacity = 5000,
    max_unique = 5000
  )
  colnames(sim$counts) <- c("0", "1")
  out <- list(landscape = landscape, sim = sim, counts = sim$counts, sample_depth = sample_depth, fit_repeat = fit_repeat)
  saveRDS(out, path)
  out
}

alfak2_counts_for_mode <- function(counts, input_mode) {
  totals <- rowSums(counts)
  if (identical(input_mode, "full")) {
    return(list(counts = counts, anchor_count_reference = NULL, input_depth = "raw"))
  }
  if (identical(input_mode, "minobs_matched")) {
    keep <- totals >= 10
    if (sum(keep) < 3L) keep[order(totals, decreasing = TRUE)[seq_len(min(3L, length(totals)))]] <- TRUE
    return(list(counts = counts[keep, , drop = FALSE], anchor_count_reference = NULL, input_depth = "raw"))
  }
  if (identical(input_mode, "soft_minobs")) {
    return(list(counts = counts, anchor_count_reference = 10, input_depth = "raw"))
  }
  stop("Unknown alfak2 input mode: ", input_mode, call. = FALSE)
}

coerce_alfakR_predictions <- function(path) {
  x <- readRDS(path)
  if (!is.data.frame(x) || !"k" %in% names(x)) {
    stop("alfakR landscape output has an unexpected schema.", call. = FALSE)
  }
  data.frame(
    karyotype = as.character(x$k),
    fitness_mean = if ("mean" %in% names(x)) as.numeric(x$mean) else if ("median" %in% names(x)) as.numeric(x$median) else NA_real_,
    fitness_sd = if ("sd" %in% names(x)) as.numeric(x$sd) else NA_real_,
    stringsAsFactors = FALSE
  )
}

run_alfak2_one <- function(row, shared, run_dir) {
  cfg <- alfak2_counts_for_mode(shared$counts, row$input_mode)
  graph_controls <- benchmark_graph_controls(row)
  started <- proc.time()[["elapsed"]]
  if ("fit_seed" %in% names(row) && is.finite(as.numeric(row$fit_seed))) {
    set.seed(as.integer(row$fit_seed))
  } else if ("count_seed" %in% names(row) && is.finite(as.numeric(row$count_seed))) {
    set.seed(as.integer(row$count_seed))
  }
  fit <- alfak2::fit_alfak2(
    cfg$counts,
    dt = 1,
    beta = 0.01,
    transition_kernel = "exact",
    local_shell_depth = graph_controls$local_shell_depth,
    global_extra_shell = graph_controls$global_extra_shell,
    min_cn = shared$landscape$min_cn,
    max_cn = shared$landscape$max_cn,
    max_nodes = graph_controls$fit_max_nodes,
    lambda_l_grid = 1,
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05,
    graph_edge_weight = "unit",
    extrapolation_method = row$extrapolation_method,
    max_prediction_distance = 2,
    anchor_count_reference = cfg$anchor_count_reference,
    input_depth = cfg$input_depth,
    control = list(eval.max = 160, iter.max = 160),
    retry_control = list(eval.max = 400, iter.max = 400)
  )
  fit_path <- file.path(run_dir, "fit.rds")
  saveRDS(fit, fit_path)
  dependency_status <- fit$global$diagnostics$dependency_status %||% "ok"
  failure_status <- if (grepl("fallback|unavailable|missing", dependency_status)) "fallback_used" else "ok"
  list(
    status = "ok",
    failure_status = failure_status,
    error_message = NA_character_,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    predictions = alfak2::summarize_alfak2(fit, layer = "global"),
    fit_path = fit_path,
    dependency_status = dependency_status,
    diagnostics = fit$global$diagnostics
  )
}

run_alfakR_one <- function(row, shared, run_dir) {
  started <- proc.time()[["elapsed"]]
  if ("fit_seed" %in% names(row) && is.finite(as.numeric(row$fit_seed))) {
    set.seed(as.integer(row$fit_seed))
  } else if ("count_seed" %in% names(row) && is.finite(as.numeric(row$count_seed))) {
    set.seed(as.integer(row$count_seed))
  }
  counts <- shared$counts
  colnames(counts) <- c("0", "1")
  yi <- list(x = counts, dt = 1)
  suppressMessages(
    alfakR::alfak(
      yi = yi,
      outdir = run_dir,
      passage_times = NULL,
      minobs = as.integer(row$minobs),
      nboot = 3,
      n0 = 1e5,
      nb = 1e7,
      pm = 0.01,
      nn_prior = as.character(row$NN_prior),
      nn_prior_grid_n = 21L,
      krig_bootstrap_mode = "marginal"
    )
  )
  landscape_path <- file.path(run_dir, "landscape.Rds")
  list(
    status = "ok",
    failure_status = "ok",
    error_message = NA_character_,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    predictions = coerce_alfakR_predictions(landscape_path),
    fit_path = landscape_path,
    dependency_status = "ok",
    diagnostics = list()
  )
}

run_one <- function(row, output_dir, resume = FALSE, alfakR_loaded = TRUE) {
  cache_path <- file.path(output_dir, "run_cache", paste0(row$run_id, ".rds"))
  if (isTRUE(resume) && file.exists(cache_path)) return(readRDS(cache_path))
  run_dir <- file.path(output_dir, "runs", row$run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  shared <- simulate_shared_input(row, output_dir)
  graph_controls <- benchmark_graph_controls(row)
  eval_graph <- alfak2:::second_layer_canonical_eval_graph(
    counts = shared$counts,
    landscape = shared$landscape,
    dt = 1,
    beta = 0.01,
    min_cn = shared$landscape$min_cn,
    max_cn = shared$landscape$max_cn,
    max_nodes = graph_controls$eval_max_nodes
  )
  result <- tryCatch({
    if (identical(row$package, "alfak2")) {
      run_alfak2_one(row, shared, run_dir)
    } else {
      if (!isTRUE(alfakR_loaded)) stop("alfakR source tree was not loaded.", call. = FALSE)
      run_alfakR_one(row, shared, run_dir)
    }
  }, error = function(e) {
    list(
      status = "failed",
      failure_status = "failed",
      error_message = conditionMessage(e),
      runtime_seconds = NA_real_,
      predictions = data.frame(karyotype = character(), fitness_mean = numeric(), fitness_sd = numeric()),
      fit_path = NA_character_,
      dependency_status = "failed",
      diagnostics = list()
    )
  })
  attached <- alfak2:::second_layer_attach_predictions(eval_graph, result$predictions)
  metrics_eval <- alfak2:::second_layer_metric_table(
    attached$nodes,
    attached$edges,
    runtime_seconds = result$runtime_seconds,
    failure_status = result$failure_status
  )
  full_eval <- alfak2:::second_layer_full_lscape_eval(result$predictions, shared$landscape)
  metrics_full <- alfak2:::second_layer_metric_table(
    full_eval$nodes,
    full_eval$edges,
    shells = "full_lscape",
    runtime_seconds = result$runtime_seconds,
    failure_status = result$failure_status
  )
  metrics <- rbind_fill(list(metrics_eval, metrics_full))
  metrics$failure_status <- result$failure_status
  metrics$fit_status <- result$status
  metrics$error_message <- result$error_message
  metrics$output_path <- result$fit_path
  metrics$dependency_status <- result$dependency_status %||% NA_character_
  metrics <- cbind(as.data.frame(row, stringsAsFactors = FALSE), metrics)
  out <- list(row = row, result = result, metrics = metrics)
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_path)
  out
}

method_label <- function(x) {
  pkg <- ifelse(x$package == "alfak2", "alfak_V2", ifelse(x$package == "alfakR", "alfak", x$package))
  ifelse(
    x$package == "alfak2",
    paste(pkg, x$input_mode, x$extrapolation_method, sep = ":"),
    paste(pkg, paste0("minobs", x$minobs), x$NN_prior_slot, sep = ":")
  )
}

metric_higher_is_better <- function(metric) {
  metric %in% c(
    "pearson", "spearman", "edge_gradient_spearman", "sign_accuracy",
    "beneficial_sign_accuracy", "deleterious_sign_accuracy",
    "top_k_overlap_count", "top_k_overlap_fraction",
    "coverage", "uncalibrated_r2", "rescaled_r2"
  )
}

metric_direction <- function(metric) {
  if (metric_higher_is_better(metric)) "higher_is_better" else "lower_is_better"
}

present_cols <- function(x, cols) {
  cols[cols %in% names(x)]
}

benchmark_unit_cols <- function(x, include_metric = FALSE) {
  present_cols(
    x,
    c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat",
      "shell", "prediction_scale", if (include_metric) "metric")
  )
}

build_paired_comparison <- function(metrics) {
  if (!nrow(metrics)) return(data.frame())
  metrics$method_label <- method_label(metrics)
  unit_cols <- benchmark_unit_cols(metrics, include_metric = TRUE)
  key <- interaction(metrics[unit_cols], drop = TRUE)
  rows <- list()
  levels_key <- levels(key)
  for (gi in seq_along(levels_key)) {
    k <- levels_key[[gi]]
    df <- metrics[key == k & is.finite(metrics$value), , drop = FALSE]
    df <- df[!duplicated(df$method_label), , drop = FALSE]
    n <- nrow(df)
    if (n < 2L) next
    p <- utils::combn(seq_len(n), 2L)
    ia <- p[1L, ]
    ib <- p[2L, ]
    va <- df$value[ia]
    vb <- df$value[ib]
    delta <- va - vb
    higher <- metric_higher_is_better(df$metric[[1]])
    a_better <- if (higher) delta >= 0 else delta <= 0
    better <- ifelse(a_better, df$method_label[ia], df$method_label[ib])
    rows[[length(rows) + 1L]] <- data.frame(
      df[1, unit_cols, drop = FALSE],
      comparison_scope = ifelse(df$package[ia] == df$package[ib], paste0(df$package[ia], "_internal"), "alfak_V2_vs_alfak"),
      package_a = df$package[ia],
      package_b = df$package[ib],
      method_a = df$method_label[ia],
      method_b = df$method_label[ib],
      input_mode_a = df$input_mode[ia],
      input_mode_b = df$input_mode[ib],
      minobs_a = df$minobs[ia],
      minobs_b = df$minobs[ib],
      NN_prior_slot_a = df$NN_prior_slot[ia],
      NN_prior_slot_b = df$NN_prior_slot[ib],
      value_a = va,
      value_b = vb,
      delta = delta,
      metric_direction = metric_direction(df$metric[[1]]),
      better_method = better,
      better_package = vapply(strsplit(better, ":", fixed = TRUE), `[[`, character(1), 1L),
      stringsAsFactors = FALSE
    )
  }
  rbind_fill(rows)
}

build_baseline_delta <- function(metrics) {
  if (!nrow(metrics)) return(data.frame())
  metrics$method_label <- method_label(metrics)
  baseline_label <- "alfak_V2:full:graph_gaussian_baseline"
  key_cols <- benchmark_unit_cols(metrics, include_metric = TRUE)
  row_cols <- present_cols(
    metrics,
    c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat",
      "shell", "prediction_scale", "metric", "package", "input_mode",
      "extrapolation_method", "minobs", "NN_prior_slot")
  )
  rows <- list()
  ri <- 0L
  key <- interaction(metrics[key_cols], drop = TRUE)
  for (k in levels(key)) {
    df <- metrics[key == k & is.finite(metrics$value), , drop = FALSE]
    base <- df$value[df$method_label == baseline_label][1]
    if (!is.finite(base)) next
    for (i in seq_len(nrow(df))) {
      ri <- ri + 1L
      rows[[ri]] <- data.frame(
        df[i, row_cols, drop = FALSE],
        baseline_method = baseline_label,
        method_label = df$method_label[[i]],
        value = df$value[[i]],
        baseline_value = base,
        delta_vs_graph_gaussian = df$value[[i]] - base,
        stringsAsFactors = FALSE
      )
    }
  }
  rbind_fill(rows)
}

build_landscape_rankings <- function(metrics) {
  weights <- c(
    rmse = 0.18, mae = 0.10, relative_rmse = 0.07, bias_abs = 0.04,
    q90_absolute_error = 0.04, uncalibrated_r2 = 0.02,
    edge_gradient_rmse = 0.12, edge_gradient_spearman = 0.08,
    centered_rmse = 0.07, rescaled_r2 = 0.04, spearman = 0.06, sign_accuracy = 0.04,
    top_k_overlap_fraction = 0.03,
    interval_coverage_95_closeness = 0.04, standardized_rmse = 0.03,
    coverage = 0.03, runtime_seconds = 0.02, failure_rate = 0.03
  )
  numerical <- c("rmse", "mae", "relative_rmse", "bias_abs", "q90_absolute_error", "uncalibrated_r2")
  shape <- c("edge_gradient_rmse", "edge_gradient_spearman", "centered_rmse", "rescaled_r2", "spearman", "sign_accuracy", "top_k_overlap_fraction")
  uncertainty <- c("interval_coverage_95_closeness", "standardized_rmse", "coverage")
  runtime <- c("runtime_seconds", "failure_rate")
  key_metrics <- names(weights)
  df <- metrics[metrics$metric %in% key_metrics & is.finite(metrics$value), , drop = FALSE]
  if (!nrow(df)) return(data.frame())
  df$method_label <- method_label(df)
  rows <- list()
  ri <- 0L
  unit_metric_cols <- benchmark_unit_cols(df, include_metric = TRUE)
  unit_cols <- benchmark_unit_cols(df, include_metric = FALSE)
  key <- interaction(df[unit_metric_cols], drop = TRUE)
  for (k in levels(key)) {
    x <- df[key == k, , drop = FALSE]
    ranks <- if (metric_higher_is_better(x$metric[[1]])) rank(-x$value, ties.method = "average") else rank(x$value, ties.method = "average")
    for (i in seq_len(nrow(x))) {
      ri <- ri + 1L
      rows[[ri]] <- data.frame(x[i, c(unit_cols, "metric", "method_label", "value"), drop = FALSE], rank = ranks[[i]])
    }
  }
  out <- rbind_fill(rows)
  weighted_one <- function(x, metrics_keep) {
    x <- x[x$metric %in% metrics_keep, , drop = FALSE]
    if (!nrow(x)) return(NA_real_)
    w <- weights[x$metric]
    sum(x$rank * w, na.rm = TRUE) / sum(w[is.finite(x$rank)], na.rm = TRUE)
  }
  out_unit_cols <- benchmark_unit_cols(out, include_metric = FALSE)
  split_key <- interaction(out[c(out_unit_cols, "method_label")], drop = TRUE)
  agg <- lapply(levels(split_key), function(k) {
    x <- out[split_key == k, , drop = FALSE]
    data.frame(
      x[1, c(out_unit_cols, "method_label"), drop = FALSE],
      numerical_weighted_rank = weighted_one(x, numerical),
      shape_weighted_rank = weighted_one(x, shape),
      uncertainty_weighted_rank = weighted_one(x, uncertainty),
      runtime_failure_rank = weighted_one(x, runtime),
      balanced_weighted_rank = weighted_one(x, key_metrics),
      stringsAsFactors = FALSE
    )
  })
  agg <- rbind_fill(agg)
  merge(out, agg, by = c(out_unit_cols, "method_label"), all.x = TRUE)
}

summary_stats <- function(metrics, by_cols) {
  if (!nrow(metrics)) return(data.frame())
  key_df <- metrics[by_cols]
  key_df[] <- lapply(key_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  split_key <- interaction(key_df, drop = TRUE, lex.order = TRUE)
  rows <- list()
  ri <- 0L
  for (k in levels(split_key)) {
    df <- metrics[split_key == k, , drop = FALSE]
    vals <- df$value[is.finite(df$value)]
    n_runs <- length(unique(df$run_id))
    n_success <- length(unique(df$run_id[df$failure_status == "ok"]))
    ri <- ri + 1L
    rows[[ri]] <- data.frame(
      df[1, by_cols, drop = FALSE],
      n_runs = n_runs,
      n_success = n_success,
      failure_rate = if (n_runs) 1 - n_success / n_runs else NA_real_,
      mean = if (length(vals)) mean(vals) else NA_real_,
      sd = if (length(vals) >= 2L) stats::sd(vals) else NA_real_,
      se = if (length(vals) >= 2L) stats::sd(vals) / sqrt(length(vals)) else NA_real_,
      median = if (length(vals)) stats::median(vals) else NA_real_,
      q25 = if (length(vals)) stats::quantile(vals, 0.25, names = FALSE) else NA_real_,
      q75 = if (length(vals)) stats::quantile(vals, 0.75, names = FALSE) else NA_real_,
      ci95_low = if (length(vals) >= 2L) mean(vals) - 1.96 * stats::sd(vals) / sqrt(length(vals)) else NA_real_,
      ci95_high = if (length(vals) >= 2L) mean(vals) + 1.96 * stats::sd(vals) / sqrt(length(vals)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  rbind_fill(rows)
}

build_rank_summary <- function(rankings) {
  if (!nrow(rankings)) return(data.frame())
  group_cols <- present_cols(rankings, c("sample_depth", "grf_lambda", "method_label"))
  x_cols <- unique(c(
    present_cols(rankings, c("sample_depth", "grf_lambda", "landscape_id", "landscape_rep", "fit_repeat", "shell", "prediction_scale", "method_label")),
    "numerical_weighted_rank", "shape_weighted_rank", "uncertainty_weighted_rank", "runtime_failure_rank", "balanced_weighted_rank"
  ))
  x <- unique(rankings[, x_cols, drop = FALSE])
  split_key <- interaction(x[group_cols], drop = TRUE)
  rows <- list()
  ri <- 0L
  for (k in levels(split_key)) {
    df <- x[split_key == k, , drop = FALSE]
    ri <- ri + 1L
    rows[[ri]] <- data.frame(
      df[1, group_cols, drop = FALSE],
      numerical_mean_rank = mean(df$numerical_weighted_rank, na.rm = TRUE),
      shape_mean_rank = mean(df$shape_weighted_rank, na.rm = TRUE),
      uncertainty_mean_rank = mean(df$uncertainty_weighted_rank, na.rm = TRUE),
      runtime_failure_mean_rank = mean(df$runtime_failure_rank, na.rm = TRUE),
      balanced_mean_rank = mean(df$balanced_weighted_rank, na.rm = TRUE),
      median_rank = stats::median(df$balanced_weighted_rank, na.rm = TRUE),
      win_rate = mean(df$balanced_weighted_rank <= 1, na.rm = TRUE),
      top3_rate = mean(df$balanced_weighted_rank <= 3, na.rm = TRUE),
      rank_sd = stats::sd(df$balanced_weighted_rank, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  rbind_fill(rows)
}

build_pareto_front <- function(metrics) {
  objectives <- c("rmse", "edge_gradient_rmse", "spearman", "sign_accuracy", "coverage", "runtime_seconds", "failure_rate")
  df <- metrics[metrics$metric %in% objectives & metrics$prediction_scale == "raw" & is.finite(metrics$value), , drop = FALSE]
  if (!nrow(df)) return(data.frame())
  df$method_label <- method_label(df)
  df$input_mode[is.na(df$input_mode) | df$input_mode == ""] <- "alfakR"
  group_cols <- present_cols(df, c("sample_depth", "grf_lambda", "shell", "input_mode", "method_label", "metric"))
  df <- aggregate(
    stats::as.formula(paste("value ~", paste(group_cols, collapse = " + "))),
    df,
    mean,
    na.rm = TRUE
  )
  id_cols <- present_cols(df, c("sample_depth", "grf_lambda", "shell", "input_mode", "method_label"))
  wide <- reshape(
    df[, c(id_cols, "metric", "value"), drop = FALSE],
    idvar = id_cols,
    timevar = "metric",
    direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))
  rows <- list()
  ri <- 0L
  pareto_key_cols <- present_cols(wide, c("sample_depth", "grf_lambda", "shell", "input_mode"))
  key <- interaction(wide[pareto_key_cols], drop = TRUE)
  for (k in levels(key)) {
    x <- wide[key == k, , drop = FALSE]
    if (!nrow(x)) next
    for (i in seq_len(nrow(x))) {
      dominated_by <- character()
      for (j in seq_len(nrow(x))) {
        if (i == j) next
        smaller <- c("rmse", "edge_gradient_rmse", "runtime_seconds", "failure_rate")
        larger <- c("spearman", "sign_accuracy", "coverage")
        a <- x[i, objectives, drop = FALSE]
        b <- x[j, objectives, drop = FALSE]
        ok_small <- all(mapply(function(bb, aa) !is.finite(bb) || !is.finite(aa) || bb <= aa, b[smaller], a[smaller]))
        ok_large <- all(mapply(function(bb, aa) !is.finite(bb) || !is.finite(aa) || bb >= aa, b[larger], a[larger]))
        strict_small <- any(mapply(function(bb, aa) is.finite(bb) && is.finite(aa) && bb < aa, b[smaller], a[smaller]))
        strict_large <- any(mapply(function(bb, aa) is.finite(bb) && is.finite(aa) && bb > aa, b[larger], a[larger]))
        if (ok_small && ok_large && (strict_small || strict_large)) dominated_by <- c(dominated_by, x$method_label[[j]])
      }
      ri <- ri + 1L
      rows[[ri]] <- data.frame(
        x[i, pareto_key_cols, drop = FALSE],
        method = x$method_label[[i]],
        is_pareto_optimal = !length(dominated_by),
        dominated_by = paste(dominated_by, collapse = ";"),
        x[i, objectives, drop = FALSE],
        stringsAsFactors = FALSE
      )
    }
  }
  rbind_fill(rows)
}

best_by_metric <- function(metrics, shell, metric, prediction_scale = "raw", higher = metric_higher_is_better(metric)) {
  x <- metrics[metrics$shell == shell & metrics$metric == metric & metrics$prediction_scale == prediction_scale & is.finite(metrics$value), , drop = FALSE]
  if (!nrow(x)) return(c(method_label = NA_character_, value = NA_character_))
  x$method_label <- method_label(x)
  agg <- aggregate(value ~ method_label, x, mean, na.rm = TRUE)
  agg <- agg[is.finite(agg$value), , drop = FALSE]
  if (!nrow(agg)) return(c(method_label = NA_character_, value = NA_character_))
  i <- if (higher) which.max(agg$value) else which.min(agg$value)
  c(method_label = agg$method_label[[i]], value = sprintf("%.5g", agg$value[[i]]))
}

dependency_requires_attention <- function(status) {
  status <- as.character(status %||% NA_character_)
  !is.na(status) & status != "ok" &
    grepl("fallback|unavailable|missing|failed|failure|error", status, ignore.case = TRUE)
}

write_recommendation <- function(metrics, rankings = NULL, path) {
  df <- metrics[metrics$shell == "all_nearfield" & is.finite(metrics$value), , drop = FALSE]
  if (!nrow(df)) {
    writeLines("# Method Recommendation\n\nNo finite metrics were available.", path)
    return(invisible(path))
  }
  numerical <- best_by_metric(metrics, "all_nearfield", "rmse", "raw", FALSE)
  shape <- best_by_metric(metrics, "all_nearfield", "edge_gradient_rmse", "raw", FALSE)
  calibrated <- best_by_metric(metrics, "all_nearfield", "rmse", "anchor_calibrated", FALSE)
  uncertainty <- best_by_metric(metrics, "all_nearfield", "interval_coverage_95_closeness", "raw", FALSE)
  runtime <- best_by_metric(metrics, "all_nearfield", "runtime_seconds", "raw", FALSE)
  balanced <- c(method_label = NA_character_, value = NA_character_)
  if (!is.null(rankings) && nrow(rankings)) {
    x <- unique(rankings[rankings$shell == "all_nearfield" & rankings$prediction_scale == "raw", c("method_label", "balanced_weighted_rank")])
    agg <- aggregate(balanced_weighted_rank ~ method_label, x, mean, na.rm = TRUE)
    if (nrow(agg)) {
      i <- which.min(agg$balanced_weighted_rank)
      balanced <- c(method_label = agg$method_label[[i]], value = sprintf("%.5g", agg$balanced_weighted_rank[[i]]))
    }
  }
  lines <- c(
    "# Method Recommendation",
    "",
    sprintf("- Best numerical method (raw all_nearfield RMSE): `%s` (%s).", numerical[["method_label"]], numerical[["value"]]),
    sprintf("- Best shape method (raw all_nearfield edge-gradient RMSE): `%s` (%s).", shape[["method_label"]], shape[["value"]]),
    sprintf("- Best calibrated method (anchor_calibrated all_nearfield RMSE): `%s` (%s).", calibrated[["method_label"]], calibrated[["value"]]),
    sprintf("- Best uncertainty method (95%% coverage closeness): `%s` (%s).", uncertainty[["method_label"]], uncertainty[["value"]]),
    sprintf("- Best runtime method: `%s` (%s seconds).", runtime[["method_label"]], runtime[["value"]]),
    sprintf("- Default balanced recommendation: `%s` (balanced mean rank %s).", balanced[["method_label"]], balanced[["value"]]),
    "",
    "Balanced ranking uses numerical accuracy (45%), shape accuracy (40%), uncertainty/coverage (10%), and runtime/failure (5%). Raw RMSE and shape metrics are both primary; affine or anchor-calibrated RMSE is reported only as diagnostic/calibratability evidence."
  )
  writeLines(lines, path)
  invisible(path)
}

write_report <- function(paths, cfg, run_index, run_results, metrics, lambda_summary, baseline_delta, path) {
  success <- sum(vapply(run_results, function(x) identical(x$result$status, "ok"), logical(1)))
  failed <- length(run_results) - success
  expected_alfak2 <- sum(run_index$package == "alfak2")
  expected_alfakR <- sum(run_index$package == "alfakR")
  fallback_rows <- metrics[dependency_requires_attention(metrics$dependency_status), , drop = FALSE]
  best_num <- best_by_metric(metrics, "all_nearfield", "rmse", "raw", FALSE)
  best_shape <- best_by_metric(metrics, "all_nearfield", "edge_gradient_rmse", "raw", FALSE)
  best_cal <- best_by_metric(metrics, "all_nearfield", "rmse", "anchor_calibrated", FALSE)
  lines <- c(
    "# Full Second-Layer 9-Method Balanced Comparison Report",
    "",
    "## Run configuration",
    sprintf("- Date/time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("- R version: `%s`", R.version.string),
    sprintf("- alfak_V2 git commit: `%s`", cfg$alfak2_git),
    sprintf("- alfak git commit: `%s`", cfg$alfakR_git),
    sprintf("- alfak_V2 version: `%s`", cfg$alfak2_version),
    sprintf("- alfak version: `%s`", cfg$alfakR_version),
    sprintf("- Python/TabPFN status: `%s`", cfg$tabpfn_status),
    sprintf("- tree backend status: `%s`", cfg$tree_backend_status),
    sprintf("- grf_lambda_values: %s", paste(sort(unique(run_index$grf_lambda)), collapse = ", ")),
    sprintf("- number of landscapes: %d", length(unique(run_index$landscape_id))),
    sprintf("- alfak_V2 input modes: %s", paste(unique(stats::na.omit(run_index$input_mode)), collapse = ", ")),
    sprintf("- alfak_V2 extrapolation methods: %s", paste(alfak2:::second_layer_alfak2_methods(), collapse = ", ")),
    sprintf("- alfak minobs values: %s", paste(unique(stats::na.omit(run_index$minobs)), collapse = ", ")),
    sprintf("- alfak NN_prior slots: %s", paste(alfak2:::second_layer_alfakR_slots()$NN_prior_slot, collapse = ", ")),
    sprintf("- expected_alfak_V2_runs: %d", expected_alfak2),
    sprintf("- expected_alfak_runs: %d", expected_alfakR),
    sprintf("- expected_total_runs: %d", nrow(run_index)),
    sprintf("- actual_started_runs: %d", length(run_results)),
    sprintf("- actual_completed_runs: %d", length(run_results)),
    sprintf("- actual_successful_runs: %d", success),
    sprintf("- actual_failed_runs: %d", failed),
    "",
    "## alfak_V2 results",
    "- Lambda-level summaries are written to `lambda_summary.csv` for the 27 alfak_V2 fitting families.",
    "- Paired graph_gaussian_baseline deltas for the eight new methods are written to `baseline_delta.csv`.",
    "- d1/d2/all_nearfield numerical and shape results are available in `metrics_by_run.csv`, `numerical_summary.csv`, and `shape_summary.csv`.",
    "",
    "## alfak results",
    "- The 15 minobs/NN_prior slot families are included in `lambda_summary.csv`.",
    "- The duplicated weighted NN_prior slots are preserved as slot4 and slot5 while passing the same alfak argument value.",
    "",
    "## alfak_V2 vs alfak",
    "- All methods are evaluated on the same canonical support-distance <= 2 graph per truth landscape, so landscape-level paired metric comparison is available in `paired_landscape_comparison.csv`.",
    "",
    "## Recommendation",
    sprintf("- Numerical winner by raw all_nearfield RMSE: `%s` (%s).", best_num[["method_label"]], best_num[["value"]]),
    sprintf("- Shape winner by raw all_nearfield edge-gradient RMSE: `%s` (%s).", best_shape[["method_label"]], best_shape[["value"]]),
    sprintf("- Calibrated winner by anchor_calibrated all_nearfield RMSE: `%s` (%s).", best_cal[["method_label"]], best_cal[["value"]]),
    "- See `method_recommendation.md` and `balanced_rank_summary.csv` for the balanced recommendation.",
    "",
    "## Important notes",
    sprintf("- Dependency/fallback rows requiring attention: %d.", nrow(fallback_rows)),
    sprintf("- Full benchmark completion claim: %s.", if (length(run_results) == nrow(run_index)) "the configured run index was fully attempted" else "not fully attempted"),
    "",
    "## Output files",
    paste0("- `", basename(paths), "`")
  )
  writeLines(lines, path)
  invisible(path)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_dir <- normalizePath(file.path(getwd()), mustWork = TRUE)
  alfakR_repo <- normalizePath(args$alfakR_repo %||% file.path(dirname(repo_dir), "alfakR"), mustWork = FALSE)
  output_dir <- normalizePath(args$output_dir %||% file.path(repo_dir, "benchmark/results/full_second_layer_9method_balanced_comparison"), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (args$n_workers != 1L) message("--n-workers currently runs deterministically in-process; requested value: ", args$n_workers)
  alfakR_status <- load_repos(repo_dir, alfakR_repo)
  alfakR_loaded <- identical(alfakR_status, "loaded")
  mode <- if (isTRUE(args$full)) "full" else "quick"
  run_index <- alfak2:::second_layer_build_run_index(mode)
  write_csv(run_index, file.path(output_dir, "run_index.csv"))
  results <- vector("list", nrow(run_index))
  for (i in seq_len(nrow(run_index))) {
    row <- run_index[i, , drop = FALSE]
    message(sprintf("[%d/%d] %s %s", i, nrow(run_index), row$run_id, method_label(row)))
    results[[i]] <- run_one(row, output_dir, resume = args$resume, alfakR_loaded = alfakR_loaded)
  }
  run_index <- rbind_fill(lapply(results, function(x) {
    data.frame(
      as.data.frame(x$row, stringsAsFactors = FALSE),
      fit_status = x$result$status,
      failure_status = x$result$failure_status,
      runtime_seconds = x$result$diagnostics$runtime_seconds %||% NA_real_,
      dependency_status = x$result$dependency_status %||% NA_character_,
      error_message = x$result$error_message,
      output_path = x$result$fit_path,
      cache_key = x$result$fit_path,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(run_index, file.path(output_dir, "run_index.csv"))
  metrics <- rbind_fill(lapply(results, `[[`, "metrics"))
  write_csv(metrics, file.path(output_dir, "metrics_by_run.csv"))
  failures <- rbind_fill(lapply(results, function(x) {
    needs_attention <- dependency_requires_attention(x$result$dependency_status %||% "ok")
    if (identical(x$result$status, "ok") && !needs_attention) return(NULL)
    data.frame(
      as.data.frame(x$row, stringsAsFactors = FALSE),
      fit_status = x$result$status,
      failure_status = x$result$failure_status,
      dependency_status = x$result$dependency_status %||% NA_character_,
      error_message = x$result$error_message,
      output_path = x$result$fit_path,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(failures, file.path(output_dir, "failures.csv"))
  dependency_status <- rbind_fill(lapply(results, function(x) {
    data.frame(
      as.data.frame(x$row, stringsAsFactors = FALSE),
      dependency_status = x$result$dependency_status %||% NA_character_,
      backend = x$result$diagnostics$backend %||%
        x$result$diagnostics$trend_filtering_backend %||%
        x$result$diagnostics$kernel %||%
        "internal",
      fallback_status = x$result$diagnostics$fallback_status %||% NA_character_,
      tabpfn_available = x$result$diagnostics$tabpfn_available %||% NA,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(dependency_status, file.path(output_dir, "dependency_status.csv"))
  paired <- build_paired_comparison(metrics)
  write_csv(paired, file.path(output_dir, "paired_landscape_comparison.csv"))
  baseline <- build_baseline_delta(metrics)
  write_csv(baseline, file.path(output_dir, "baseline_delta.csv"))
  rankings <- build_landscape_rankings(metrics)
  write_csv(rankings, file.path(output_dir, "landscape_rankings.csv"))
  pareto <- build_pareto_front(metrics)
  write_csv(pareto, file.path(output_dir, "pareto_front.csv"))
  summary_by <- c("grf_lambda", "package", "input_mode", "extrapolation_method", "minobs", "NN_prior_slot", "shell", "prediction_scale", "metric")
  lambda_summary <- summary_stats(
    metrics,
    summary_by
  )
  write_csv(lambda_summary, file.path(output_dir, "lambda_summary.csv"))
  rank_summary <- build_rank_summary(rankings)
  write_csv(rank_summary, file.path(output_dir, "lambda_method_rank_summary.csv"))
  write_csv(rank_summary, file.path(output_dir, "balanced_rank_summary.csv"))
  overall_summary <- summary_stats(
    metrics,
    c("package", "input_mode", "extrapolation_method", "minobs", "NN_prior_slot", "shell", "prediction_scale", "metric")
  )
  write_csv(overall_summary, file.path(output_dir, "overall_summary.csv"))
  numerical_metrics <- c("rmse", "mae", "relative_rmse", "bias", "bias_abs", "median_absolute_error", "q90_absolute_error", "uncalibrated_r2", "count_weighted_rmse", "count_weighted_mae")
  shape_metrics <- c("centered_rmse", "pearson", "spearman", "edge_gradient_rmse", "edge_gradient_spearman", "sign_accuracy", "beneficial_sign_accuracy", "deleterious_sign_accuracy", "top_k_overlap_fraction")
  calibration_metrics <- c("calibration_intercept", "calibration_slope", "affine_rmse", "interval_coverage_95", "interval_coverage_95_closeness", "standardized_rmse", "mean_pred_sd", "median_pred_sd")
  numerical_summary <- summary_stats(metrics[metrics$metric %in% numerical_metrics, , drop = FALSE], summary_by)
  shape_summary <- summary_stats(metrics[metrics$metric %in% shape_metrics, , drop = FALSE], summary_by)
  calibration_summary <- summary_stats(metrics[metrics$metric %in% calibration_metrics, , drop = FALSE], summary_by)
  write_csv(numerical_summary, file.path(output_dir, "numerical_summary.csv"))
  write_csv(shape_summary, file.path(output_dir, "shape_summary.csv"))
  write_csv(calibration_summary, file.path(output_dir, "calibration_summary.csv"))
  full_results <- list(
    run_index = run_index,
    run_results = results,
    metrics_by_run = metrics,
    paired_landscape_comparison = paired,
    baseline_delta = baseline,
    landscape_rankings = rankings,
    pareto_front = pareto,
    lambda_summary = lambda_summary,
    lambda_method_rank_summary = rank_summary,
    overall_summary = overall_summary,
    numerical_summary = numerical_summary,
    shape_summary = shape_summary,
    calibration_summary = calibration_summary,
    balanced_rank_summary = rank_summary,
    failures = failures,
    dependency_status = dependency_status,
    alfakR_status = alfakR_status
  )
  saveRDS(full_results, file.path(output_dir, "full_results.rds"))
  write_recommendation(metrics, rankings, file.path(output_dir, "method_recommendation.md"))
  tabpfn_status <- if (any(grepl("tabpfn", dependency_status$dependency_status %||% character()))) {
    paste(unique(stats::na.omit(dependency_status$dependency_status[grepl("tabpfn", dependency_status$dependency_status)])), collapse = ", ")
  } else {
    "not_attempted_or_not_recorded"
  }
  tree_status <- paste(unique(stats::na.omit(dependency_status$dependency_status[grepl("tree|ridge|xgboost|lightgbm|catboost|ranger|randomForest", dependency_status$dependency_status)])), collapse = ", ")
  cfg <- list(
    alfak2_git = git_hash(repo_dir),
    alfakR_git = git_hash(alfakR_repo),
    alfak2_version = package_version_safe("alfak2"),
    alfakR_version = package_version_safe("alfakR"),
    tabpfn_status = tabpfn_status,
    tree_backend_status = if (nzchar(tree_status)) tree_status else "not_recorded"
  )
  required_paths <- file.path(output_dir, c(
    "run_index.csv", "metrics_by_run.csv", "paired_landscape_comparison.csv",
    "baseline_delta.csv", "landscape_rankings.csv", "lambda_summary.csv",
    "lambda_method_rank_summary.csv", "overall_summary.csv", "numerical_summary.csv",
    "shape_summary.csv", "calibration_summary.csv", "balanced_rank_summary.csv",
    "pareto_front.csv", "failures.csv", "full_results.rds", "dependency_status.csv",
    "method_recommendation.md", "report.md"
  ))
  write_report(
    paths = required_paths,
    cfg = cfg,
    run_index = run_index,
    run_results = results,
    metrics = metrics,
    lambda_summary = lambda_summary,
    baseline_delta = baseline,
    path = file.path(output_dir, "report.md")
  )
  best_print <- function(shell, metric, scale = "raw", higher = metric_higher_is_better(metric)) {
    best_by_metric(metrics, shell, metric, scale, higher)[["method_label"]]
  }
  rank_for <- function(method, col) {
    x <- unique(rankings[rankings$shell == "all_nearfield" & rankings$prediction_scale == "raw", c("method_label", col)])
    if (!nrow(x)) return(NA_real_)
    agg <- aggregate(x[[col]], list(method_label = x$method_label), mean, na.rm = TRUE)
    names(agg)[2] <- col
    v <- agg[[col]][agg$method_label == method]
    if (length(v)) v[[1]] else NA_real_
  }
  baseline <- "alfak_V2:full:graph_gaussian_baseline"
  best_balanced_for_shell <- function(shell) {
    x <- unique(rankings[rankings$shell == shell & rankings$prediction_scale == "raw", c("method_label", "balanced_weighted_rank")])
    if (nrow(x)) {
      agg <- aggregate(balanced_weighted_rank ~ method_label, x, mean, na.rm = TRUE)
      agg$method_label[[which.min(agg$balanced_weighted_rank)]]
    } else NA_character_
  }
  best_balanced <- best_balanced_for_shell("all_nearfield")
  alfakR_best <- {
    x <- metrics[metrics$package == "alfakR" & metrics$shell == "all_nearfield" & metrics$prediction_scale == "raw" & metrics$metric == "rmse" & is.finite(metrics$value), , drop = FALSE]
    if (nrow(x)) {
      x$method_label <- method_label(x)
      agg <- aggregate(value ~ method_label, x, mean, na.rm = TRUE)
      agg$method_label[[which.min(agg$value)]]
    } else NA_character_
  }
  cat(
    "expected runs:", nrow(run_index), "\n",
    "completed runs:", length(results), "\n",
    "successful runs:", sum(vapply(results, function(x) identical(x$result$status, "ok"), logical(1))), "\n",
    "failed runs:", sum(!vapply(results, function(x) identical(x$result$status, "ok"), logical(1))), "\n",
    "output directory:", output_dir, "\n",
    "best numerical method by all_nearfield raw rmse:", best_print("all_nearfield", "rmse", "raw", FALSE), "\n",
    "best shape method by all_nearfield edge_gradient_rmse:", best_print("all_nearfield", "edge_gradient_rmse", "raw", FALSE), "\n",
    "best calibrated method by anchor_calibrated rmse:", best_print("all_nearfield", "rmse", "anchor_calibrated", FALSE), "\n",
    "best balanced method by balanced_weighted_rank:", best_balanced, "\n",
    "best d1 numerical method:", best_print("d1", "rmse", "raw", FALSE), "\n",
    "best d1 shape method:", best_print("d1", "edge_gradient_rmse", "raw", FALSE), "\n",
    "best d1 balanced method:", best_balanced_for_shell("d1"), "\n",
    "best d2 numerical method:", best_print("d2", "rmse", "raw", FALSE), "\n",
    "best d2 shape method:", best_print("d2", "edge_gradient_rmse", "raw", FALSE), "\n",
    "best d2 balanced method:", best_balanced_for_shell("d2"), "\n",
    "graph_gaussian_baseline numerical rank:", rank_for(baseline, "numerical_weighted_rank"), "\n",
    "graph_gaussian_baseline shape rank:", rank_for(baseline, "shape_weighted_rank"), "\n",
    "graph_gaussian_baseline balanced rank:", rank_for(baseline, "balanced_weighted_rank"), "\n",
    "alfakR best setting:", alfakR_best, "\n",
    "report.md path:", file.path(output_dir, "report.md"), "\n",
    sep = ""
  )
}

if (identical(environment(), globalenv())) main()
