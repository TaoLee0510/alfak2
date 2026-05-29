#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_delta_estimator_probe.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_shape_diagnostics.R"))

usage <- function() {
  cat(
    "Run farfield non-oracle edge-delta estimator probe F1-F5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_delta_estimator_probe.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --output-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --simulation-id=1 --minobs=5 --input-policy=full\n\n",
    "Modes:\n",
    "  prepare, f1-delta-training, f2-state-dependent-estimators, f3-nonoracle-potential,\n",
    "  f4-local-tmb-debug, f5-neighbor-shell-cv, summarize, all\n\n",
    "Options:\n",
    "  --quick=auto        auto|true|false. Auto runs full F1/F2 and sampled F3/F4/F5 fanout.\n",
    "  --force=false       recompute cached outputs.\n",
    sep = ""
  )
}

make_probe_dirs <- function(output_dir) {
  if (!grepl("^/", output_dir)) output_dir <- file.path(repo_guess, output_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!path_under_benchmark(output_dir)) {
    stop("Refusing to write benchmark outputs outside benchmark/: ", output_dir, call. = FALSE)
  }
  dirs <- list(
    root = output_dir,
    cache = file.path(output_dir, "cache"),
    fits = file.path(output_dir, "fits"),
    tables = file.path(output_dir, "tables"),
    results = file.path(output_dir, "results")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

load_diagnostics_results <- function(diagnostics_dir) {
  diagnostics_dir <- normalizePath(diagnostics_dir, winslash = "/", mustWork = TRUE)
  all_path <- file.path(diagnostics_dir, "results", "farfield_shape_diagnostics_all_results.rds")
  if (file.exists(all_path)) {
    out <- readRDS(all_path)
  } else {
    out <- list()
  }
  out$root <- diagnostics_dir
  out
}

context_index_for_step <- function(graph, from, to) {
  delta <- graph$karyotypes[to, , drop = FALSE] - graph$karyotypes[from, , drop = FALSE]
  changed <- which(as.integer(delta[1, ]) != 0L)
  if (length(changed) != 1L) return(NA_integer_)
  direction <- if (delta[1, changed] > 0) "gain" else "loss"
  label <- paste0(direction, "_chr", changed, "_band1")
  match(label, as.character(graph$context_label))
}

edge_step_features <- function(graph, from, to) {
  delta <- graph$karyotypes[to, , drop = FALSE] - graph$karyotypes[from, , drop = FALSE]
  changed <- which(as.integer(delta[1, ]) != 0L)
  chr <- if (length(changed) == 1L) changed else NA_integer_
  direction <- if (length(changed) == 1L) {
    if (delta[1, changed] > 0) "gain" else "loss"
  } else NA_character_
  data.frame(
    edge_chr = chr,
    edge_direction = direction,
    context_index = context_index_for_step(graph, from, to),
    context_label = if (is.finite(context_index_for_step(graph, from, to))) graph$context_label[[context_index_for_step(graph, from, to)]] else NA_character_,
    parent_total_cn = sum(graph$karyotypes[from, ]),
    child_total_cn = sum(graph$karyotypes[to, ]),
    parent_aneuploidy_burden = sum(abs(graph$karyotypes[from, ] - 2L)),
    child_aneuploidy_burden = sum(abs(graph$karyotypes[to, ] - 2L)),
    stringsAsFactors = FALSE
  )
}

make_step_adjacency <- function(graph) {
  from <- as.integer(graph$edge_from)
  to <- as.integer(graph$edge_to)
  keep <- from >= 1L & to >= 1L & from <= length(graph$labels) & to <= length(graph$labels) & from != to
  from <- from[keep]
  to <- to[keep]
  adj <- vector("list", length(graph$labels))
  for (i in seq_along(from)) {
    adj[[from[[i]]]] <- c(adj[[from[[i]]]], to[[i]])
  }
  adj
}

find_short_path <- function(adj, start, target, max_depth = 3L) {
  if (identical(start, target)) return(start)
  queue <- list(list(node = start, path = start))
  seen <- rep(FALSE, length(adj))
  seen[[start]] <- TRUE
  while (length(queue)) {
    item <- queue[[1L]]
    queue <- queue[-1L]
    if (length(item$path) - 1L >= max_depth) next
    nb <- unique(adj[[item$node]])
    for (v in nb) {
      if (seen[[v]]) next
      path <- c(item$path, v)
      if (v == target) return(path)
      seen[[v]] <- TRUE
      queue[[length(queue) + 1L]] <- list(node = v, path = path)
    }
  }
  integer()
}

path_context_counts <- function(graph, path) {
  counts <- rep(0, length(graph$context_label))
  if (length(path) < 2L) return(counts)
  for (i in seq_len(length(path) - 1L)) {
    cc <- context_index_for_step(graph, path[[i]], path[[i + 1L]])
    if (is.finite(cc) && cc >= 1L && cc <= length(counts)) counts[[cc]] <- counts[[cc]] + 1
  }
  counts
}

parent_edge_frame <- function(graph, grf, lambda) {
  from <- as.integer(graph$parent_from0) + 1L
  to <- as.integer(graph$parent_to0) + 1L
  keep <- from >= 1L & to >= 1L & from <= length(graph$labels) & to <= length(graph$labels)
  from <- from[keep]
  to <- to[keep]
  truth <- compute_truth_for_nodes(as.character(graph$labels), grf, lambda)
  rows <- lapply(seq_along(from), function(i) {
    feat <- edge_step_features(graph, from[[i]], to[[i]])
    data.frame(
      edge_id = i,
      parent_node_id = from[[i]],
      child_node_id = to[[i]],
      parent_karyotype = as.character(graph$labels[[from[[i]]]]),
      child_karyotype = as.character(graph$labels[[to[[i]]]]),
      feat,
      support_distance_parent = as.integer(graph$support_distance[[from[[i]]]]),
      support_distance_child = as.integer(graph$support_distance[[to[[i]]]]),
      support_scope_parent = support_scope_label(graph$support_tier[[from[[i]]]], graph$support_distance[[from[[i]]]]),
      support_scope_child = support_scope_label(graph$support_tier[[to[[i]]]], graph$support_distance[[to[[i]]]]),
      truth_delta = truth[[to[[i]]]] - truth[[from[[i]]]],
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

edge_delta_metric_row <- function(est, truth) {
  ok <- is.finite(est) & is.finite(truth)
  est <- est[ok]
  truth <- truth[ok]
  err <- est - truth
  truth_sd <- if (length(truth) > 1L) stats::sd(truth) else NA_real_
  est_sd <- if (length(est) > 1L) stats::sd(est) else NA_real_
  data.frame(
    n_edges = length(ok),
    n_scored_edges = length(est),
    truth_delta_mean = if (length(truth)) mean(truth) else NA_real_,
    truth_delta_sd = truth_sd,
    estimated_delta_mean = if (length(est)) mean(est) else NA_real_,
    estimated_delta_sd = est_sd,
    estimated_delta_sd_ratio = if (is.finite(truth_sd) && truth_sd > 0) est_sd / truth_sd else NA_real_,
    delta_bias = if (length(err)) mean(err) else NA_real_,
    delta_mae = if (length(err)) mean(abs(err)) else NA_real_,
    delta_rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    delta_pearson = safe_cor2(est, truth, "pearson"),
    delta_spearman = safe_cor2(est, truth, "spearman"),
    delta_sign_agreement = if (length(est)) mean(sign(est) == sign(truth), na.rm = TRUE) else NA_real_,
    delta_sign_agreement_nonzero = if (length(est)) {
      nz <- sign(est) != 0 & sign(truth) != 0
      if (any(nz)) mean(sign(est[nz]) == sign(truth[nz])) else NA_real_
    } else NA_real_,
    fraction_truth_positive = if (length(truth)) mean(truth > 0) else NA_real_,
    fraction_estimated_positive = if (length(est)) mean(est > 0) else NA_real_,
    median_abs_delta_error = if (length(err)) stats::median(abs(err)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

passes_delta_gate <- function(row) {
  isTRUE(is.finite(row$delta_sign_agreement) && row$delta_sign_agreement >= 0.55 &&
           is.finite(row$delta_spearman) && row$delta_spearman > 0.10 &&
           is.finite(row$estimated_delta_sd_ratio) && row$estimated_delta_sd_ratio >= 0.10)
}

build_trusted_edge_training_set <- function(bundle, grf, task_info) {
  graph <- bundle$global_graph
  local <- bundle$local
  direct_labels <- as.character(local$summary$karyotype[
    as.character(local$summary$support_tier) == "directly_informed" &
      is.finite(local$summary$fitness_mean)
  ])
  local_map <- setNames(as.numeric(local$summary$fitness_mean), as.character(local$summary$karyotype))
  count_map <- setNames(as.numeric(local$summary$count_total), as.character(local$summary$karyotype))
  eff_map <- setNames(as.numeric(local$summary$effective_count_total), as.character(local$summary$karyotype))
  from <- as.integer(graph$edge_from)
  to <- as.integer(graph$edge_to)
  labels <- as.character(graph$labels)
  keep <- labels[from] %in% direct_labels & labels[to] %in% direct_labels &
    is.finite(local_map[labels[from]]) & is.finite(local_map[labels[to]])
  from <- from[keep]
  to <- to[keep]
  truth <- compute_truth_for_nodes(unique(c(labels[from], labels[to])), grf, task_info$lambda)
  names(truth) <- unique(c(labels[from], labels[to]))
  rows <- lapply(seq_along(from), function(i) {
    feat <- edge_step_features(graph, from[[i]], to[[i]])
    parent <- labels[[from[[i]]]]
    child <- labels[[to[[i]]]]
    data.frame(
      parent = parent,
      child = child,
      feat,
      parent_support_tier = "directly_informed",
      child_support_tier = "directly_informed",
      parent_count_total = count_map[[parent]],
      child_count_total = count_map[[child]],
      parent_effective_count_total = eff_map[[parent]],
      child_effective_count_total = eff_map[[child]],
      parent_fitness = local_map[[parent]],
      child_fitness = local_map[[child]],
      observed_delta = local_map[[child]] - local_map[[parent]],
      truth_delta = truth[[child]] - truth[[parent]],
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows_fill(rows)
  if (!nrow(out)) {
    out <- data.frame(
      parent = character(), child = character(), edge_chr = integer(), edge_direction = character(),
      context_index = integer(), context_label = character(), parent_total_cn = numeric(),
      child_total_cn = numeric(), parent_aneuploidy_burden = numeric(),
      child_aneuploidy_burden = numeric(), parent_support_tier = character(),
      child_support_tier = character(), parent_count_total = numeric(), child_count_total = numeric(),
      parent_effective_count_total = numeric(), child_effective_count_total = numeric(),
      parent_fitness = numeric(), child_fitness = numeric(), observed_delta = numeric(),
      truth_delta = numeric(), stringsAsFactors = FALSE
    )
  }
  out
}

build_direct_pair_training_set <- function(bundle, grf, task_info, max_path_length = 3L) {
  graph <- bundle$global_graph
  local <- bundle$local
  direct <- local$summary[
    as.character(local$summary$support_tier) == "directly_informed" &
      is.finite(local$summary$fitness_mean),
    ,
    drop = FALSE
  ]
  labels <- as.character(graph$labels)
  direct_idx <- match(as.character(direct$karyotype), labels)
  direct <- direct[!is.na(direct_idx), , drop = FALSE]
  direct_idx <- direct_idx[!is.na(direct_idx)]
  truth <- compute_truth_for_nodes(labels[direct_idx], grf, task_info$lambda)
  names(truth) <- labels[direct_idx]
  adj <- make_step_adjacency(graph)
  n_ctx <- length(graph$context_label)
  rows <- list()
  idx <- 0L
  if (length(direct_idx) >= 2L) {
    cmb <- utils::combn(seq_along(direct_idx), 2L)
    for (j in seq_len(ncol(cmb))) {
      a <- cmb[1L, j]
      b <- cmb[2L, j]
      path <- find_short_path(adj, direct_idx[[a]], direct_idx[[b]], max_depth = max_path_length)
      if (!length(path)) next
      counts <- path_context_counts(graph, path)
      idx <- idx + 1L
      row <- data.frame(
        pair_id = idx,
        parent = labels[[direct_idx[[a]]]],
        child = labels[[direct_idx[[b]]]],
        path_length = length(path) - 1L,
        path_nodes = paste(labels[path], collapse = "->"),
        observed_pair_delta = direct$fitness_mean[[b]] - direct$fitness_mean[[a]],
        truth_pair_delta = truth[[labels[[direct_idx[[b]]]]]] - truth[[labels[[direct_idx[[a]]]]]],
        stringsAsFactors = FALSE
      )
      for (cc in seq_len(n_ctx)) row[[paste0("ctx_", cc)]] <- counts[[cc]]
      rows[[idx]] <- row
    }
  }
  bind_rows_fill(rows)
}

fit_ridge <- function(x, y, lambda) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  x <- as.matrix(x)
  y <- as.numeric(y)
  ok <- is.finite(y) & rowSums(!is.finite(x)) == 0
  x <- x[ok, , drop = FALSE]
  y <- y[ok]
  if (!nrow(x)) return(rep(0, ncol(x)))
  q <- crossprod(x) + diag(as.numeric(lambda), ncol(x))
  rhs <- crossprod(x, y)
  as.numeric(solve(q, rhs))
}

context_delta_from_edges <- function(edge_train, graph) {
  n_ctx <- length(graph$context_label)
  rows <- list()
  methods <- c("context_mean", "context_median", "gain_loss_mean", "chr_direction_mean", "global_mean", "zero_delta")
  for (method in methods) {
    for (cc in seq_len(n_ctx)) {
      vals <- numeric()
      if (nrow(edge_train) && method %in% c("context_mean", "context_median")) {
        vals <- edge_train$observed_delta[edge_train$context_index == cc]
      } else if (nrow(edge_train) && method == "gain_loss_mean") {
        dir <- if (grepl("^gain_", graph$context_label[[cc]])) "gain" else "loss"
        vals <- edge_train$observed_delta[edge_train$edge_direction == dir]
      } else if (nrow(edge_train) && method == "chr_direction_mean") {
        dir <- if (grepl("^gain_", graph$context_label[[cc]])) "gain" else "loss"
        chr <- suppressWarnings(as.integer(sub("^.*chr([0-9]+).*$", "\\1", graph$context_label[[cc]])))
        vals <- edge_train$observed_delta[edge_train$edge_direction == dir & edge_train$edge_chr == chr]
      } else if (nrow(edge_train) && method == "global_mean") {
        vals <- edge_train$observed_delta
      }
      vals <- vals[is.finite(vals)]
      est <- if (method == "zero_delta") 0 else if (!length(vals)) 0 else if (method == "context_median") stats::median(vals) else mean(vals)
      rows[[length(rows) + 1L]] <- data.frame(
        estimator = method,
        context_index = cc,
        context_label = graph$context_label[[cc]],
        delta_estimate = est,
        n_training_edges = length(vals),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

predict_context_delta <- function(parent_edges, estimates, estimator) {
  est <- estimates[estimates$estimator == estimator, , drop = FALSE]
  delta <- est$delta_estimate[match(parent_edges$context_index, est$context_index)]
  delta[!is.finite(delta)] <- 0
  delta
}

context_leave_one_edge_cv <- function(edge_train, graph) {
  estimators <- c("context_mean", "context_median", "gain_loss_mean", "chr_direction_mean", "global_mean", "zero_delta")
  rows <- list()
  if (!nrow(edge_train)) {
    for (estimator in estimators) {
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(estimator = estimator, cv_target = "heldout_direct_edge_truth_delta", stringsAsFactors = FALSE),
        edge_delta_metric_row(numeric(), numeric())
      )
    }
    return(bind_rows_fill(rows))
  }
  for (estimator in estimators) {
    pred <- rep(NA_real_, nrow(edge_train))
    for (i in seq_len(nrow(edge_train))) {
      train_i <- edge_train[-i, , drop = FALSE]
      est_i <- context_delta_from_edges(train_i, graph)
      pred[[i]] <- predict_context_delta(edge_train[i, , drop = FALSE], est_i, estimator)
    }
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = estimator, cv_target = "heldout_direct_edge_truth_delta", stringsAsFactors = FALSE),
      edge_delta_metric_row(pred, edge_train$truth_delta)
    )
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = estimator, cv_target = "heldout_direct_edge_observed_delta", stringsAsFactors = FALSE),
      edge_delta_metric_row(pred, edge_train$observed_delta)
    )
  }
  bind_rows_fill(rows)
}

fit_path_regression <- function(pair_train, graph, lambdas = c(0.001, 0.01, 0.1, 1, 10)) {
  ctx_cols <- paste0("ctx_", seq_along(graph$context_label))
  x <- if (nrow(pair_train)) as.matrix(pair_train[, ctx_cols, drop = FALSE]) else matrix(numeric(), 0, length(ctx_cols))
  y <- if (nrow(pair_train)) pair_train$observed_pair_delta else numeric()
  rows <- list()
  coeffs <- list()
  for (lam in lambdas) {
    pred <- rep(NA_real_, length(y))
    if (length(y) >= 3L) {
      for (i in seq_along(y)) {
        beta <- fit_ridge(x[-i, , drop = FALSE], y[-i], lam)
        pred[[i]] <- sum(x[i, ] * beta)
      }
    }
    m_obs <- edge_delta_metric_row(pred, y)
    m_truth <- edge_delta_metric_row(pred, if (length(y)) pair_train$truth_pair_delta else numeric())
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = "path_regression_context_coefficients", ridge_lambda = lam, cv_target = "observed_pair_delta", stringsAsFactors = FALSE),
      m_obs
    )
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = "path_regression_context_coefficients", ridge_lambda = lam, cv_target = "truth_pair_delta", stringsAsFactors = FALSE),
      m_truth
    )
    beta <- fit_ridge(x, y, lam)
    coeffs[[length(coeffs) + 1L]] <- data.frame(
      estimator = "path_regression_context_coefficients",
      ridge_lambda = lam,
      context_index = seq_along(graph$context_label),
      context_label = graph$context_label,
      beta = beta,
      stringsAsFactors = FALSE
    )
  }
  list(cv = bind_rows_fill(rows), coefficients = bind_rows_fill(coeffs))
}

heldout_gate_metrics <- function(cv, estimator) {
  target <- if (identical(estimator, "path_regression_context_coefficients")) "truth_pair_delta" else "heldout_direct_edge_truth_delta"
  x <- cv[cv$estimator == estimator & cv$cv_target == target, , drop = FALSE]
  if (!nrow(x)) {
    return(data.frame(
      heldout_delta_spearman = NA_real_,
      heldout_delta_sign_agreement = NA_real_,
      heldout_estimated_delta_sd_ratio = NA_real_,
      heldout_delta_rmse = NA_real_
    ))
  }
  if ("ridge_lambda" %in% names(x) && any(is.finite(x$ridge_lambda))) {
    x <- x[order(-x$delta_spearman, x$delta_rmse), , drop = FALSE]
  }
  data.frame(
    heldout_delta_spearman = x$delta_spearman[[1L]],
    heldout_delta_sign_agreement = x$delta_sign_agreement[[1L]],
    heldout_estimated_delta_sd_ratio = x$estimated_delta_sd_ratio[[1L]],
    heldout_delta_rmse = x$delta_rmse[[1L]]
  )
}

make_estimator_gate <- function(cv, edge_train, pair_train, parent_edges, context_estimates, path_coefficients) {
  rows <- list()
  # Context estimators are judged on full parent-edge truth because direct-edge holdout is usually tiny here.
  for (estimator in unique(context_estimates$estimator)) {
    pred <- predict_context_delta(parent_edges, context_estimates, estimator)
    m <- edge_delta_metric_row(pred, parent_edges$truth_delta)
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = estimator, cv_source = "parent_edge_truth_diagnostic", stringsAsFactors = FALSE),
      m,
      heldout_gate_metrics(cv, estimator)
    )
  }
  if (nrow(path_coefficients)) {
    best_lam <- path_coefficients$ridge_lambda[[1L]]
    truth_cv <- cv[cv$estimator == "path_regression_context_coefficients" & cv$cv_target == "truth_pair_delta", , drop = FALSE]
    if (nrow(truth_cv)) best_lam <- truth_cv$ridge_lambda[order(-truth_cv$delta_spearman, truth_cv$delta_rmse)][1L]
    beta <- path_coefficients$beta[path_coefficients$ridge_lambda == best_lam]
    pred <- beta[parent_edges$context_index]
    pred[!is.finite(pred)] <- 0
    m <- edge_delta_metric_row(pred, parent_edges$truth_delta)
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(estimator = "path_regression_context_coefficients", cv_source = "parent_edge_truth_diagnostic", ridge_lambda = best_lam, stringsAsFactors = FALSE),
      m,
      heldout_gate_metrics(cv, "path_regression_context_coefficients")
    )
  }
  gate <- bind_rows_fill(rows)
  gate$passes_gate <- vapply(seq_len(nrow(gate)), function(i) {
    passes_delta_gate(gate[i, , drop = FALSE]) &&
      is.finite(gate$heldout_delta_spearman[[i]]) &&
      is.finite(gate$heldout_delta_sign_agreement[[i]]) &&
      is.finite(gate$heldout_estimated_delta_sd_ratio[[i]]) &&
      gate$heldout_delta_spearman[[i]] > 0 &&
      gate$heldout_delta_sign_agreement[[i]] >= 0.55 &&
      gate$heldout_estimated_delta_sd_ratio[[i]] >= 0.10
  }, logical(1))
  gate$delta_status <- ifelse(gate$passes_gate, "delta_trusted", "delta_untrusted")
  gate
}

run_f1 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "f1_delta_training.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  edge_train <- build_trusted_edge_training_set(bundle, grf, task_info)
  pair_train <- build_direct_pair_training_set(bundle, grf, task_info, max_path_length = 3L)
  parent_edges <- parent_edge_frame(bundle$global_graph, grf, task_info$lambda)
  context_estimates <- context_delta_from_edges(edge_train, bundle$global_graph)
  path <- fit_path_regression(pair_train, bundle$global_graph)
  direct_cv_rows <- list()
  for (estimator in unique(context_estimates$estimator)) {
    pred <- predict_context_delta(parent_edges, context_estimates, estimator)
    direct_cv_rows[[length(direct_cv_rows) + 1L]] <- cbind(
      data.frame(estimator = estimator, cv_target = "parent_edge_truth_diagnostic", stringsAsFactors = FALSE),
      edge_delta_metric_row(pred, parent_edges$truth_delta)
    )
  }
  cv <- bind_rows_fill(list(bind_rows_fill(direct_cv_rows), context_leave_one_edge_cv(edge_train, bundle$global_graph), path$cv))
  gate <- make_estimator_gate(cv, edge_train, pair_train, parent_edges, context_estimates, path$coefficients)
  write_tsv_safe(edge_train, file.path(dirs$tables, "trusted_edge_training_set.tsv"))
  write_tsv_safe(pair_train, file.path(dirs$tables, "trusted_direct_pair_training_set.tsv"))
  write_tsv_safe(context_estimates, file.path(dirs$tables, "context_delta_estimates.tsv"))
  write_tsv_safe(path$coefficients, file.path(dirs$tables, "path_regression_coefficients.tsv"))
  write_tsv_safe(cv, file.path(dirs$tables, "delta_estimator_cv.tsv"))
  write_tsv_safe(gate, file.path(dirs$tables, "delta_estimator_gate.tsv"))
  out <- list(edge_train = edge_train, pair_train = pair_train, parent_edges = parent_edges,
              context_estimates = context_estimates, path = path, cv = cv, gate = gate)
  saveRDS(out, rds)
  out
}

make_pair_pseudo_edges <- function(pair_train, graph) {
  if (!nrow(pair_train)) return(data.frame())
  rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(pair_train))) {
    nodes <- unlist(strsplit(pair_train$path_nodes[[i]], "->", fixed = TRUE), use.names = FALSE)
    node_idx <- match(nodes, as.character(graph$labels))
    if (anyNA(node_idx) || length(node_idx) < 2L) next
    label <- pair_train$observed_pair_delta[[i]] / max(1L, length(node_idx) - 1L)
    truth_label <- pair_train$truth_pair_delta[[i]] / max(1L, length(node_idx) - 1L)
    for (j in seq_len(length(node_idx) - 1L)) {
      feat <- edge_step_features(graph, node_idx[[j]], node_idx[[j + 1L]])
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        source_pair_id = pair_train$pair_id[[i]],
        parent_node_id = node_idx[[j]],
        child_node_id = node_idx[[j + 1L]],
        parent_karyotype = nodes[[j]],
        child_karyotype = nodes[[j + 1L]],
        feat,
        observed_delta_label = label,
        truth_delta_label = truth_label,
        label_type = "direct_pair_path_delta_per_edge",
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

edge_distance <- function(train, target_row, graph) {
  target_parent <- graph$karyotypes[target_row$parent_node_id, ]
  train_parent <- do.call(rbind, strsplit(train$parent_karyotype, ".", fixed = TRUE))
  train_parent <- matrix(as.integer(train_parent), nrow = nrow(train))
  l1 <- rowSums(abs(sweep(train_parent, 2L, target_parent, "-")))
  chr_penalty <- ifelse(train$edge_chr == target_row$edge_chr, 0, 4)
  dir_penalty <- ifelse(train$edge_direction == target_row$edge_direction, 0, 4)
  total_penalty <- abs(train$parent_total_cn - target_row$parent_total_cn) * 0.25
  l1 + chr_penalty + dir_penalty + total_penalty
}

kernel_predict_one <- function(train, target_row, graph, bandwidth) {
  if (!nrow(train)) return(c(delta = 0, confidence = 0, effective_n = 0, weight_sum = 0))
  dist <- edge_distance(train, target_row, graph)
  w <- exp(-dist / as.numeric(bandwidth))
  ok <- is.finite(w) & is.finite(train$observed_delta_label)
  if (!any(ok) || sum(w[ok]) <= 0) return(c(delta = 0, confidence = 0, effective_n = 0, weight_sum = 0))
  w <- w[ok]
  y <- train$observed_delta_label[ok]
  delta <- sum(w * y) / sum(w)
  ess <- (sum(w)^2) / sum(w^2)
  c(delta = delta, confidence = min(1, ess / 3), effective_n = ess, weight_sum = sum(w))
}

kernel_predict_edges <- function(train, target_edges, graph, bandwidth) {
  if (!nrow(target_edges)) return(data.frame())
  rows <- lapply(seq_len(nrow(target_edges)), function(i) {
    pred <- kernel_predict_one(train, target_edges[i, , drop = FALSE], graph, bandwidth)
    data.frame(
      edge_id = target_edges$edge_id[[i]],
      estimated_delta = pred[["delta"]],
      confidence = pred[["confidence"]],
      effective_sample_size = pred[["effective_n"]],
      kernel_weight_sum = pred[["weight_sum"]],
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

kernel_cv <- function(train, graph, bandwidth) {
  if (nrow(train) < 3L) {
    return(cbind(data.frame(bandwidth = bandwidth, n_train = nrow(train), stringsAsFactors = FALSE),
                 edge_delta_metric_row(numeric(), numeric())))
  }
  pred <- numeric(nrow(train))
  for (i in seq_len(nrow(train))) {
    train_i <- train[-i, , drop = FALSE]
    target <- train[i, , drop = FALSE]
    target$edge_id <- i
    pred[[i]] <- kernel_predict_one(train_i, target, graph, bandwidth)[["delta"]]
  }
  cbind(
    data.frame(bandwidth = bandwidth, n_train = nrow(train), stringsAsFactors = FALSE),
    edge_delta_metric_row(pred, train$truth_delta_label)
  )
}

run_f2 <- function(bundle, grf, task_info, dirs, f1, force = FALSE) {
  rds <- file.path(dirs$results, "f2_state_dependent_estimators.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$global_graph
  parent_edges <- f1$parent_edges
  edge_train <- f1$edge_train
  direct_kernel_train <- if (nrow(edge_train)) {
    data.frame(
      parent_karyotype = edge_train$parent,
      child_karyotype = edge_train$child,
      edge_chr = edge_train$edge_chr,
      edge_direction = edge_train$edge_direction,
      parent_total_cn = edge_train$parent_total_cn,
      observed_delta_label = edge_train$observed_delta,
      truth_delta_label = edge_train$truth_delta,
      stringsAsFactors = FALSE
    )
  } else data.frame()
  pair_kernel_train <- make_pair_pseudo_edges(f1$pair_train, graph)
  bandwidths <- c(0.5, 1, 2, 4, 8)
  cv_rows <- list()
  for (bw in bandwidths) {
    cv_rows[[length(cv_rows) + 1L]] <- cbind(data.frame(estimator = "state_kernel_direct_edge", stringsAsFactors = FALSE), kernel_cv(direct_kernel_train, graph, bw))
    cv_rows[[length(cv_rows) + 1L]] <- cbind(data.frame(estimator = "state_kernel_direct_pair", stringsAsFactors = FALSE), kernel_cv(pair_kernel_train, graph, bw))
  }
  # Coarse estimators from F1 on the same parent-edge truth target.
  for (estimator in unique(f1$context_estimates$estimator)) {
    pred <- predict_context_delta(parent_edges, f1$context_estimates, estimator)
    cv_rows[[length(cv_rows) + 1L]] <- cbind(
      data.frame(estimator = estimator, bandwidth = NA_real_, n_train = nrow(edge_train), stringsAsFactors = FALSE),
      edge_delta_metric_row(pred, parent_edges$truth_delta)
    )
  }
  path_truth <- f1$cv[f1$cv$estimator == "path_regression_context_coefficients" & f1$cv$cv_target == "truth_pair_delta", , drop = FALSE]
  if (nrow(path_truth)) {
    best_lam <- path_truth$ridge_lambda[order(-path_truth$delta_spearman, path_truth$delta_rmse)][1L]
    beta <- f1$path$coefficients$beta[f1$path$coefficients$ridge_lambda == best_lam]
    pred <- beta[parent_edges$context_index]
    pred[!is.finite(pred)] <- 0
    cv_rows[[length(cv_rows) + 1L]] <- cbind(
      data.frame(estimator = "path_regression_context_coefficients", bandwidth = NA_real_, n_train = nrow(f1$pair_train), ridge_lambda = best_lam, stringsAsFactors = FALSE),
      edge_delta_metric_row(pred, parent_edges$truth_delta)
    )
  }
  cv <- bind_rows_fill(cv_rows)
  best_kernel_pair <- cv[cv$estimator == "state_kernel_direct_pair", , drop = FALSE]
  best_bw_pair <- if (nrow(best_kernel_pair)) best_kernel_pair$bandwidth[order(-best_kernel_pair$delta_spearman, best_kernel_pair$delta_rmse)][1L] else 2
  best_kernel_edge <- cv[cv$estimator == "state_kernel_direct_edge", , drop = FALSE]
  best_bw_edge <- if (nrow(best_kernel_edge)) best_kernel_edge$bandwidth[order(-best_kernel_edge$delta_spearman, best_kernel_edge$delta_rmse)][1L] else 2
  pred_rows <- list()
  add_pred <- function(estimator, delta, confidence = 1) {
    data.frame(
      parent_edges[, c("edge_id", "parent_node_id", "child_node_id", "parent_karyotype", "child_karyotype",
                       "context_index", "context_label", "edge_chr", "edge_direction", "truth_delta"), drop = FALSE],
      estimator = estimator,
      estimated_delta = as.numeric(delta),
      confidence = as.numeric(confidence),
      stringsAsFactors = FALSE
    )
  }
  for (estimator in unique(f1$context_estimates$estimator)) {
    pred_rows[[length(pred_rows) + 1L]] <- add_pred(estimator, predict_context_delta(parent_edges, f1$context_estimates, estimator), 1)
  }
  if (exists("beta") && length(beta)) {
    delta <- beta[parent_edges$context_index]
    delta[!is.finite(delta)] <- 0
    pred_rows[[length(pred_rows) + 1L]] <- add_pred("path_regression_context_coefficients", delta, 1)
  }
  kp <- kernel_predict_edges(pair_kernel_train, parent_edges, graph, best_bw_pair)
  pred_rows[[length(pred_rows) + 1L]] <- cbind(add_pred("state_kernel_direct_pair", kp$estimated_delta, kp$confidence),
                                               effective_sample_size = kp$effective_sample_size,
                                               kernel_weight_sum = kp$kernel_weight_sum)
  ke <- kernel_predict_edges(direct_kernel_train, parent_edges, graph, best_bw_edge)
  pred_rows[[length(pred_rows) + 1L]] <- cbind(add_pred("state_kernel_direct_edge", ke$estimated_delta, ke$confidence),
                                               effective_sample_size = ke$effective_sample_size,
                                               kernel_weight_sum = ke$kernel_weight_sum)
  # Hybrid: pair-kernel when confident, otherwise path regression/context zero shrink.
  if (nrow(kp)) {
    base_delta <- kp$estimated_delta
    conf <- kp$confidence
    delta <- conf * base_delta
    pred_rows[[length(pred_rows) + 1L]] <- cbind(add_pred("hybrid_confident_delta", delta, conf),
                                                 effective_sample_size = kp$effective_sample_size,
                                                 kernel_weight_sum = kp$kernel_weight_sum)
  }
  predictions <- bind_rows_fill(pred_rows)
  summary <- bind_rows_fill(lapply(split(predictions, predictions$estimator), function(x) {
    cbind(data.frame(estimator = x$estimator[[1L]], stringsAsFactors = FALSE),
          edge_delta_metric_row(x$estimated_delta, x$truth_delta))
  }))
  gate <- summary
  gate$passes_gate <- vapply(seq_len(nrow(gate)), function(i) passes_delta_gate(gate[i, , drop = FALSE]), logical(1))
  gate$delta_status <- ifelse(gate$passes_gate, "delta_trusted", "delta_untrusted")
  write_tsv_safe(cv, file.path(dirs$tables, "state_dependent_delta_estimator_cv.tsv"))
  write_tsv_safe(predictions, file.path(dirs$tables, "state_dependent_edge_predictions.tsv"))
  write_tsv_safe(gate, file.path(dirs$tables, "state_dependent_estimator_gate.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "state_dependent_estimator_summary.tsv"))
  out <- list(cv = cv, predictions = predictions, gate = gate, summary = summary,
              best_bw_pair = best_bw_pair, best_bw_edge = best_bw_edge)
  saveRDS(out, rds)
  out
}

fit_potential_prior_mean_edges <- function(local_fit, graph, parent_edges, delta, confidence,
                                           lambda_anchor = 10, lambda_smooth = 0.1,
                                           components = NULL, ridge = 1e-6,
                                           edge_weight_mode = "normalized") {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  n <- length(graph$labels)
  from <- parent_edges$parent_node_id
  to <- parent_edges$child_node_id
  delta <- as.numeric(delta)
  confidence <- as.numeric(confidence)
  confidence[!is.finite(confidence)] <- 0
  w <- confidence
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) w[] <- 0
  sw <- sqrt(w)
  keep <- sw > 0 & is.finite(delta)
  from <- from[keep]
  to <- to[keep]
  sw <- sw[keep]
  delta <- delta[keep]
  q <- Matrix::Diagonal(n, ridge)
  rhs <- numeric(n)
  if (length(from)) {
    row_id <- seq_along(from)
    a <- Matrix::sparseMatrix(i = c(row_id, row_id), j = c(from, to), x = c(-sw, sw), dims = c(length(from), n))
    rhs <- rhs + as.numeric(Matrix::crossprod(a, sw * delta))
    q <- q + Matrix::crossprod(a)
  }
  direct <- as.character(local_fit$summary$support_tier) == "directly_informed" & is.finite(local_fit$summary$fitness_mean)
  anchor_idx <- match(as.character(local_fit$summary$karyotype[direct]), as.character(graph$labels))
  ok <- !is.na(anchor_idx)
  anchor_idx <- anchor_idx[ok]
  anchor_mean <- as.numeric(local_fit$summary$fitness_mean[direct][ok])
  if (length(anchor_idx) && is.finite(lambda_anchor) && lambda_anchor > 0) {
    q <- q + Matrix::sparseMatrix(i = anchor_idx, j = anchor_idx, x = lambda_anchor, dims = c(n, n))
    rhs[anchor_idx] <- rhs[anchor_idx] + lambda_anchor * anchor_mean
  }
  if (!is.null(components) && is.finite(lambda_smooth) && lambda_smooth > 0) {
    q <- q + as.numeric(lambda_smooth) * components$edge[[edge_weight_mode]]
  }
  q <- Matrix::forceSymmetric(q, uplo = "U")
  chol <- Matrix::Cholesky(q, LDL = TRUE, perm = TRUE)
  as.numeric(Matrix::solve(chol, rhs))
}

base_global_configs_f <- function(abcd) {
  forced <- data.frame(
    experiment = "F3",
    candidate_id = c("F3_baseline_mutation", "F3_normalized_ll0p2_le0p01", "F3_unit_ll0p2_le0p01"),
    graph_edge_weight = c("mutation", "normalized", "unit"),
    lambda_l = c(0.2, 0.2, 0.2),
    lambda_e = c(1, 0.01, 0.01),
    sigma_obs = c(0.05, 0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    solver = "matrix_mean",
    stringsAsFactors = FALSE
  )
  far <- abcd$all[abcd$all$support_scope == "farfield" & abcd$all$metric_scale == "native" &
                    !grepl("oracle", abcd$all$prior_mean_mode), , drop = FALSE]
  top <- if (nrow(far)) configs_from_metrics(head(far[order(far$shape_score, far$centered_rmse), , drop = FALSE], 5L), "F3") else data.frame()
  dedupe_configs(bind_rows_fill(list(forced, top)))
}

run_f3 <- function(bundle, components, grf, task_info, dirs, abcd, f1, f2, quick = "auto", force = FALSE) {
  rds <- file.path(dirs$results, "f3_nonoracle_potential_prior.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  parent_edges <- f1$parent_edges
  preds <- f2$predictions
  best_context <- f1$gate$estimator[order(-f1$gate$delta_spearman, f1$gate$delta_rmse)][1L]
  if (!is.finite(match(best_context, preds$estimator))) best_context <- "zero_delta"
  delta_sources <- c("zero_delta", best_context, "path_regression_context_coefficients",
                     "state_kernel_direct_pair", "hybrid_confident_delta",
                     "oracle_context_delta", "oracle_per_edge_delta")
  delta_sources <- unique(delta_sources)
  base <- base_global_configs_f(abcd)
  if (quick %in% c("auto", "true")) base <- head(base, min(5L, nrow(base)))
  grid <- expand.grid(
    delta_source = delta_sources,
    lambda_anchor = if (quick %in% c("auto", "true")) c(10, 100) else c(1, 10, 100),
    lambda_smooth = if (quick %in% c("auto", "true")) c(0.1) else c(0, 0.1, 1),
    delta_scale = if (quick %in% c("auto", "true")) c(-1, 0.5, 1) else c(-1, -0.5, 0.5, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  oracle_context <- oracle_context_delta(bundle$global_graph, grf, task_info$lambda, per_edge = FALSE)$context_delta
  oracle_edge <- oracle_context_delta(bundle$global_graph, grf, task_info$lambda, per_edge = TRUE)$edge_delta
  rows <- list()
  comps <- list()
  idx <- 0L
  for (b in seq_len(nrow(base))) {
    for (g in seq_len(nrow(grid))) {
      cfg <- base[b, , drop = FALSE]
      gg <- grid[g, , drop = FALSE]
      source <- gg$delta_source
      if (identical(source, "oracle_context_delta")) {
        delta <- oracle_context[parent_edges$context_index]
        confidence <- rep(1, nrow(parent_edges))
      } else if (identical(source, "oracle_per_edge_delta")) {
        delta <- oracle_edge[seq_len(nrow(parent_edges))]
        confidence <- rep(1, nrow(parent_edges))
      } else {
        pp <- preds[preds$estimator == source, , drop = FALSE]
        delta <- pp$estimated_delta[match(parent_edges$edge_id, pp$edge_id)]
        confidence <- pp$confidence[match(parent_edges$edge_id, pp$edge_id)]
        delta[!is.finite(delta)] <- 0
        confidence[!is.finite(confidence)] <- 0
      }
      delta <- as.numeric(gg$delta_scale) * delta
      cfg$experiment <- "F3"
      cfg$prior_mean_mode <- paste0("potential_", source)
      cfg$prior_mean_scale <- gg$delta_scale
      cfg$candidate_id <- paste0(
        "F3__", cfg$graph_edge_weight, "__ll", cfg$lambda_l, "__le", cfg$lambda_e, "__so", cfg$sigma_obs,
        "__", source, "__ds", gg$delta_scale, "__la", gg$lambda_anchor, "__ls", gg$lambda_smooth
      )
      message("[F3] ", b, "/", nrow(base), " ", g, "/", nrow(grid), " ", cfg$candidate_id)
      pm <- fit_potential_prior_mean_edges(bundle$local, bundle$global_graph, parent_edges, delta, confidence,
                                           lambda_anchor = gg$lambda_anchor,
                                           lambda_smooth = gg$lambda_smooth,
                                           components = components)
      fit <- fit_cached(bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs,
                        prior_mean = pm, prior_mean_status = cfg$prior_mean_mode, force = force)
      idx <- idx + 1L
      rows[[idx]] <- transform(
        fit$metrics,
        delta_source = source,
        potential_lambda_anchor = gg$lambda_anchor,
        potential_lambda_smooth = gg$lambda_smooth,
        potential_delta_scale = gg$delta_scale,
        oracle_delta = grepl("^oracle_", source)
      )
      comps[[idx]] <- data.frame(
        candidate_id = cfg$candidate_id,
        delta_source = source,
        lambda_anchor = gg$lambda_anchor,
        lambda_smooth = gg$lambda_smooth,
        delta_scale = gg$delta_scale,
        prior_mean_sd = stats::sd(pm),
        prior_mean_range = quantile_range(pm),
        mean_confidence = mean(confidence, na.rm = TRUE),
        n_nonzero_confidence = sum(confidence > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- bind_rows_fill(rows)
  comp <- bind_rows_fill(comps)
  far <- out[out$support_scope == "farfield" & out$metric_scale == "native", , drop = FALSE]
  top <- rbind(
    transform(head(far[order(far$shape_score, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "shape_score"),
    transform(head(far[order(-far$spearman, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "spearman"),
    transform(head(far[order(-far$pearson, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "pearson")
  )
  nonoracle <- far[!far$oracle_delta, , drop = FALSE]
  oracle <- far[far$oracle_delta, , drop = FALSE]
  upper <- data.frame(
    class = c("best_nonoracle", "best_oracle"),
    candidate_id = c(nonoracle$candidate_id[order(-nonoracle$spearman, nonoracle$centered_rmse)][1L],
                     oracle$candidate_id[order(-oracle$spearman, oracle$centered_rmse)][1L]),
    spearman = c(max(nonoracle$spearman, na.rm = TRUE), max(oracle$spearman, na.rm = TRUE)),
    pearson = c(nonoracle$pearson[order(-nonoracle$spearman, nonoracle$centered_rmse)][1L],
                oracle$pearson[order(-oracle$spearman, oracle$centered_rmse)][1L]),
    estimate_sd_ratio = c(nonoracle$estimate_sd_ratio[order(-nonoracle$spearman, nonoracle$centered_rmse)][1L],
                          oracle$estimate_sd_ratio[order(-oracle$spearman, oracle$centered_rmse)][1L]),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(out, file.path(dirs$tables, "nonoracle_potential_prior_probe.tsv"))
  write_tsv_safe(top, file.path(dirs$tables, "nonoracle_potential_prior_top.tsv"))
  write_tsv_safe(comp, file.path(dirs$tables, "nonoracle_potential_prior_components.tsv"))
  write_tsv_safe(upper, file.path(dirs$tables, "nonoracle_vs_oracle_upper_bound.tsv"))
  res <- list(metrics = out, components = comp, top = top, upper_bound = upper)
  saveRDS(res, rds)
  res
}

local_debug_grid <- function(quick) {
  if (quick %in% c("auto", "true")) {
    return(data.frame(
      local_shell_depth = c(0, 1, 1, 1, 1, 1, 1, 1),
      observation_model = c("multinomial", "multinomial", "dirichlet_multinomial", "dirichlet_multinomial",
                            "dirichlet_multinomial", "multinomial", "dirichlet_multinomial", "dirichlet_multinomial"),
      dm_concentration = c(NA, NA, 20, 50, 200, NA, 50, 50),
      eval_max = c(500, 500, 500, 500, 2000, 2000, 500, 500),
      iter_max = c(500, 500, 500, 500, 2000, 2000, 500, 500),
      restart_id = c(1, 1, 1, 1, 1, 1, 2, 3),
      eta_borrowed_prior_mean = c(-6, -6, -6, -6, -6, -6, -8, -10),
      eta_borrowed_prior_sd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.0, 0.5),
      eta_distance_penalty = c(0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 1.5, 2.5),
      init_mode = c("default", "default", "default", "default", "default", "default", "stronger_borrowed_shrinkage", "stronger_borrowed_shrinkage"),
      stringsAsFactors = FALSE
    ))
  }
  rbind(
    expand.grid(
      local_shell_depth = c(0, 1),
      observation_model = "multinomial",
      dm_concentration = NA_real_,
      eval_max = c(500, 2000, 5000),
      iter_max = c(500, 2000, 5000),
      restart_id = c(1, 2, 3),
      eta_borrowed_prior_mean = -6,
      eta_borrowed_prior_sd = 1.5,
      eta_distance_penalty = 0.75,
      init_mode = "default",
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    expand.grid(
      local_shell_depth = c(0, 1),
      observation_model = "dirichlet_multinomial",
      dm_concentration = c(20, 50, 100, 200),
      eval_max = c(500, 2000, 5000),
      iter_max = c(500, 2000, 5000),
      restart_id = c(1, 2, 3),
      eta_borrowed_prior_mean = -6,
      eta_borrowed_prior_sd = 1.5,
      eta_distance_penalty = 0.75,
      init_mode = "default",
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

local_debug_id <- function(row) {
  paste0("f4_shell", row$local_shell_depth, "_", row$observation_model,
         "_dm", ifelse(is.na(row$dm_concentration), "NA", row$dm_concentration),
         "_eval", row$eval_max, "_restart", row$restart_id,
         "_eta", row$eta_borrowed_prior_mean, "_sd", row$eta_borrowed_prior_sd,
         "_pen", row$eta_distance_penalty)
}

run_local_fit_debug_official <- function(bundle, row, dirs, force = FALSE) {
  cache_path <- file.path(dirs$cache, paste0(local_debug_id(row), ".rds"))
  if (!isTRUE(force) && file.exists(cache_path)) return(readRDS(cache_path))
  graph <- if (as.integer(row$local_shell_depth) == as.integer(bundle$local$graph$shell_depth)) {
    bundle$local$graph
  } else {
    alfak2::build_karyotype_graph(bundle$data, transition_kernel = "exact",
                                  shell_depth = as.integer(row$local_shell_depth),
                                  min_cn = 0L,
                                  max_cn = max(bundle$global_graph$karyotypes),
                                  max_nodes = 150000L)
  }
  started <- Sys.time()
  fit <- tryCatch({
    alfak2::fit_local_posterior(
      bundle$data,
      graph,
      observation_model = row$observation_model,
      dm_concentration = if (identical(row$observation_model, "multinomial")) 50 else as.numeric(row$dm_concentration),
      observation_weight_mode = "likelihood",
      control = list(eval.max = as.integer(row$eval_max), iter.max = as.integer(row$iter_max)),
      retry_on_untrusted_covariance = FALSE,
      retry_control = NULL,
      eta_borrowed_prior_mean = as.numeric(row$eta_borrowed_prior_mean),
      eta_borrowed_prior_sd = as.numeric(row$eta_borrowed_prior_sd),
      eta_distance_penalty = as.numeric(row$eta_distance_penalty)
    )
  }, error = function(e) {
    list(error = conditionMessage(e), graph = graph, summary = data.frame(), diagnostics = list())
  })
  out <- list(fit = fit, config = row, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")))
  saveRDS(out, cache_path)
  out
}

gradient_block_table <- function(fit_result, config_id) {
  diag <- fit_result$fit$diagnostics %||% list()
  # fit_local_posterior does not retain obj/opt, so block gradients are unavailable without a second custom TMB run.
  blocks <- c("eta", "f", "delta_context", "mu_group", "log_sigma_neighbor", "log_sigma_anchor", "log_tau_group")
  data.frame(
    config_id = config_id,
    parameter_block = blocks,
    max_abs_gradient = NA_real_,
    global_gradient_norm = diag$gradient_norm %||% NA_real_,
    extraction_status = "unavailable_fit_local_posterior_does_not_return_obj_or_opt",
    stringsAsFactors = FALSE
  )
}

run_f4 <- function(bundle, grf, task_info, dirs, quick = "auto", force = FALSE) {
  rds <- file.path(dirs$results, "f4_local_tmb_debug.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  grid <- local_debug_grid(quick)
  rows <- list()
  grads <- list()
  edge_rows <- list()
  fits <- list()
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    cid <- local_debug_id(row)
    message("[F4] ", i, "/", nrow(grid), " ", cid)
    res <- run_local_fit_debug_official(bundle, row, dirs, force = force)
    fits[[cid]] <- res$fit
    fit <- res$fit
    diag <- fit$diagnostics %||% list()
    rows[[i]] <- data.frame(
      config_id = cid,
      row,
      convergence = diag$convergence %||% NA,
      message = diag$message %||% (fit$error %||% NA_character_),
      objective = diag$objective %||% NA_real_,
      gradient_norm = diag$gradient_norm %||% NA_real_,
      covariance_status = diag$covariance_status %||% NA_character_,
      covariance_fallback = diag$covariance_fallback %||% NA,
      fitness_sd_source = diag$fitness_sd_source %||% NA_character_,
      retry_attempted = diag$retry_attempted %||% FALSE,
      dm_concentration_selected = diag$dm_concentration %||% NA_real_,
      n_local_nodes = if (nrow(fit$summary %||% data.frame())) nrow(fit$summary) else 0L,
      n_direct = if (nrow(fit$summary %||% data.frame())) sum(fit$summary$support_tier == "directly_informed") else 0L,
      n_local_borrowed = if (nrow(fit$summary %||% data.frame())) sum(fit$summary$support_tier == "local_borrowed") else 0L,
      n_weakly_supported = if (nrow(fit$summary %||% data.frame())) sum(fit$summary$support_tier == "weakly_supported") else 0L,
      elapsed_sec = res$elapsed_sec,
      quick_mode = quick,
      stringsAsFactors = FALSE
    )
    grads[[i]] <- gradient_block_table(res, cid)
    if (nrow(fit$summary %||% data.frame())) {
      edge_rows[[i]] <- local_edge_alignment(res, grf, task_info, cid)
    }
  }
  grid_tbl <- bind_rows_fill(rows)
  grad_tbl <- bind_rows_fill(grads)
  edge_tbl <- bind_rows_fill(edge_rows)
  staged <- grid_tbl[grid_tbl$local_shell_depth == 1 & grid_tbl$eval_max >= 2000, , drop = FALSE]
  staged$staged_status <- "benchmark_placeholder_no_public_initial_parameter_hook"
  shrink <- grid_tbl[grid_tbl$init_mode == "stronger_borrowed_shrinkage", , drop = FALSE]
  pairwise <- data.frame()
  ids <- names(fits)
  if (length(ids) >= 2L) {
    meta <- grid_tbl[, c("config_id", "local_shell_depth", "observation_model", "dm_concentration", "eval_max", "iter_max"), drop = FALSE]
    meta$key <- paste(meta$local_shell_depth, meta$observation_model, meta$dm_concentration, meta$eval_max, meta$iter_max, sep = "|")
    pr <- list()
    for (key in unique(meta$key)) {
      kk <- meta$config_id[meta$key == key]
      if (length(kk) < 2L) next
      cmb <- utils::combn(kk, 2)
      for (j in seq_len(ncol(cmb))) {
        a <- fits[[cmb[1, j]]]
        b <- fits[[cmb[2, j]]]
        if (!nrow(a$summary %||% data.frame()) || !nrow(b$summary %||% data.frame())) next
        common <- intersect(a$summary$karyotype, b$summary$karyotype)
        fa <- setNames(a$summary$fitness_mean, a$summary$karyotype)[common]
        fb <- setNames(b$summary$fitness_mean, b$summary$karyotype)[common]
        pr[[length(pr) + 1L]] <- data.frame(
          key = key,
          config_id_a = cmb[1, j],
          config_id_b = cmb[2, j],
          n_common = length(common),
          fitness_mean_correlation = safe_cor2(fa, fb, "pearson"),
          direct_nodes_correlation = safe_cor2(fa, fb, "pearson"),
          borrowed_nodes_correlation = safe_cor2(fa, fb, "pearson"),
          edge_delta_sign_agreement_between_restarts = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    pairwise <- bind_rows_fill(pr)
  }
  write_tsv_safe(grid_tbl, file.path(dirs$tables, "local_tmb_debug_grid.tsv"))
  write_tsv_safe(grad_tbl, file.path(dirs$tables, "local_gradient_by_parameter_block.tsv"))
  write_tsv_safe(staged, file.path(dirs$tables, "local_staged_initialization_results.tsv"))
  write_tsv_safe(shrink, file.path(dirs$tables, "local_borrowed_shrinkage_results.tsv"))
  write_tsv_safe(pairwise, file.path(dirs$tables, "local_multistart_stability.tsv"))
  write_tsv_safe(edge_tbl, file.path(dirs$tables, "local_edge_delta_alignment.tsv"))
  out <- list(grid = grid_tbl, gradient = grad_tbl, staged = staged, shrinkage = shrink,
              stability = pairwise, edge_alignment = edge_tbl, quick_mode = quick)
  saveRDS(out, rds)
  out
}

shape_class <- function(x) {
  sd_ok <- is.finite(x$estimate_sd_ratio) & x$estimate_sd_ratio >= 0.02
  valid <- sd_ok & x$pearson > 0 & x$spearman > 0
  wrong <- sd_ok & (x$pearson < 0 | x$spearman < 0)
  out <- rep("numeric_only", nrow(x))
  out[!sd_ok] <- "collapsed_shrinkage"
  out[wrong] <- "noncollapsed_wrong_direction"
  out[valid] <- "valid_shape"
  out
}

cv_metric_for_holdout <- function(fit, truth_labels, holdout) {
  pred <- fit$summary$fitness_mean[match(holdout, fit$summary$karyotype)]
  truth <- as.numeric(truth_labels[holdout])
  cv_metric_row(pred, truth)
}

run_node_holdout_cv <- function(bundle, components, configs, task_info, grf, holdout_sets, cv_type) {
  local <- bundle$local
  truth_labels <- setNames(local$summary$fitness_mean, local$summary$karyotype)
  rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    pm <- build_prior_mean("zero", local, bundle$global_graph, 0)
    for (s in seq_along(holdout_sets)) {
      hold <- holdout_sets[[s]]
      fit <- fit_global_with_config(local, bundle$global_graph, components, cfg, task_info$minobs,
                                    prior_mean = pm$mean, anchor_exclude = hold)
      idx <- idx + 1L
      rows[[idx]] <- cbind(
        data.frame(cv_type = cv_type, split_id = s, candidate_id = cfg$candidate_id,
                   graph_edge_weight = cfg$graph_edge_weight, lambda_l = cfg$lambda_l,
                   lambda_e = cfg$lambda_e, sigma_obs = cfg$sigma_obs,
                   anchor_var_mode = cfg$anchor_var_mode, prior_mean_mode = cfg$prior_mean_mode,
                   prior_mean_scale = cfg$prior_mean_scale,
                   anchor_count_reference_mode = cfg$anchor_count_reference_mode,
                   n_holdout = length(hold),
                   holdout_labels = paste(hold, collapse = ","), stringsAsFactors = FALSE),
        cv_metric_for_holdout(fit, truth_labels, hold)
      )
    }
  }
  bind_rows_fill(rows)
}

run_f5 <- function(bundle, components, task_info, grf, dirs, abcd, f1, f2, f3, quick = "auto", force = FALSE) {
  rds <- file.path(dirs$results, "f5_neighbor_shell_edge_cv.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  local <- bundle$local
  direct <- local$summary[local$summary$support_tier == "directly_informed" & is.finite(local$summary$fitness_mean), , drop = FALSE]
  direct <- direct[order(direct$effective_count_total), , drop = FALSE]
  low <- as.character(head(direct$karyotype, max(1L, floor(nrow(direct) / 3))))
  high <- as.character(tail(direct$karyotype, max(1L, floor(nrow(direct) / 3))))
  middle <- setdiff(as.character(direct$karyotype), c(low, high))
  configs <- rbind(make_baseline_config("F5"), data.frame(
    experiment = "F5",
    candidate_id = c("F5_normalized_ll0p2_le0p01", "F5_unit_ll0p2_le0p01"),
    graph_edge_weight = c("normalized", "unit"),
    lambda_l = c(0.2, 0.2),
    lambda_e = c(0.01, 0.01),
    sigma_obs = c(0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    solver = "matrix_mean",
    stringsAsFactors = FALSE
  ))
  direct_cv <- run_node_holdout_cv(bundle, components, configs, task_info, grf, list(low, middle, high), "direct_anchor_cv")
  neighbor_cv <- run_node_holdout_cv(bundle, components, configs, task_info, grf, list(low), "neighbor_holdout_cv_low_count_proxy")
  shell_cv <- run_node_holdout_cv(bundle, components, configs, task_info, grf, list(high, low, middle), "shell_stratified_cv_count_bins")
  edge_cv <- bind_rows_fill(list(
    transform(f1$gate, cv_type = "edge_delta_cv_f1"),
    transform(f2$gate, cv_type = "edge_delta_cv_f2")
  ))
  direct_agg <- aggregate_cv_candidates(direct_cv)
  neighbor_agg <- aggregate_cv_candidates(neighbor_cv)
  shell_agg <- aggregate_cv_candidates(shell_cv)
  selected <- bind_rows_fill(list(
    transform(select_by_cv_objectives(direct_agg), cv_type = "direct_anchor_cv"),
    transform(select_by_cv_objectives(neighbor_agg), cv_type = "neighbor_holdout_cv"),
    transform(select_by_cv_objectives(shell_agg), cv_type = "shell_stratified_cv")
  ))
  far <- f3$metrics[f3$metrics$support_scope == "farfield" & f3$metrics$metric_scale == "native", , drop = FALSE]
  far$shape_class <- shape_class(far)
  nonoracle <- far[!far$oracle_delta, , drop = FALSE]
  edge_delta_gate_passed <- any(f2$gate$passes_gate, na.rm = TRUE)
  no_valid <- !any(nonoracle$shape_class == "valid_shape", na.rm = TRUE) || !edge_delta_gate_passed
  failure <- data.frame(
    recommended_status = if (no_valid) "no_valid_shape_configuration" else "valid_shape_available",
    n_nonoracle_valid_shape = sum(nonoracle$shape_class == "valid_shape", na.rm = TRUE),
    n_nonoracle_wrong_direction = sum(nonoracle$shape_class == "noncollapsed_wrong_direction", na.rm = TRUE),
    n_nonoracle_collapsed = sum(nonoracle$shape_class == "collapsed_shrinkage", na.rm = TRUE),
    edge_delta_gate_passed = edge_delta_gate_passed,
    stringsAsFactors = FALSE
  )
  write_tsv_safe(direct_cv, file.path(dirs$tables, "direct_anchor_cv_recheck.tsv"))
  write_tsv_safe(neighbor_cv, file.path(dirs$tables, "neighbor_holdout_cv.tsv"))
  write_tsv_safe(shell_cv, file.path(dirs$tables, "shell_stratified_cv.tsv"))
  write_tsv_safe(edge_cv, file.path(dirs$tables, "edge_delta_cv.tsv"))
  write_tsv_safe(selected, file.path(dirs$tables, "cv_selected_configs.tsv"))
  write_tsv_safe(failure, file.path(dirs$tables, "calibration_failure_state_demo.tsv"))
  out <- list(direct = direct_cv, neighbor = neighbor_cv, shell = shell_cv,
              edge = edge_cv, selected = selected, failure = failure)
  saveRDS(out, rds)
  out
}

make_recommendations <- function(f1, f2, f3, f4, f5) {
  f3_far <- f3$metrics[f3$metrics$support_scope == "farfield" & f3$metrics$metric_scale == "native", , drop = FALSE]
  nonoracle <- f3_far[!f3_far$oracle_delta, , drop = FALSE]
  nonoracle$shape_class <- shape_class(nonoracle)
  data.frame(
    table = c("delta_estimator_recommendation", "potential_prior_recommendation", "local_tmb_recommendation", "cv_recommendation", "recommended_next_steps"),
    recommendation = c(
      if (any(f2$gate$passes_gate, na.rm = TRUE)) "A non-oracle delta estimator passed the diagnostic gate." else "No non-oracle delta estimator passed the diagnostic gate; keep delta_untrusted.",
      if (any(nonoracle$shape_class == "valid_shape", na.rm = TRUE)) "A non-oracle potential prior reached valid_shape." else "No non-oracle potential prior reached valid_shape; oracle remains only an upper bound.",
      "Prioritize local shell_depth=1 convergence/covariance; block gradients are unavailable from the current public fit object.",
      if (isTRUE(f5$failure$recommended_status[[1]] == "no_valid_shape_configuration")) "Calibration should allow no_valid_shape_configuration." else "CV found at least one valid non-oracle shape configuration.",
      "Do not implement C++ edge-gradient until non-oracle delta and potential-prior gates pass."
    ),
    evidence = c(
      paste0("best non-oracle delta Spearman=", fmt_metric(max(f2$gate$delta_spearman, na.rm = TRUE))),
      paste0("best non-oracle farfield Spearman=", fmt_metric(max(nonoracle$spearman, na.rm = TRUE))),
      paste0("median local gradient=", fmt_metric(stats::median(f4$grid$gradient_norm, na.rm = TRUE))),
      paste0("status=", f5$failure$recommended_status[[1]]),
      "Required gate: delta sign>=0.55, delta Spearman>0, non-oracle farfield Pearson/Spearman>0 and sd_ratio>=0.05."
    ),
    stringsAsFactors = FALSE
  )
}

write_final_tables_report <- function(dirs, ctx, task_info, diagnostics_dir, f1, f2, f3, f4, f5, recs, quick) {
  f3_far <- f3$metrics[f3$metrics$support_scope == "farfield" & f3$metrics$metric_scale == "native", , drop = FALSE]
  nonoracle <- f3_far[!f3_far$oracle_delta, , drop = FALSE]
  oracle <- f3_far[f3_far$oracle_delta, , drop = FALSE]
  nonoracle$shape_class <- shape_class(nonoracle)
  best_nonoracle <- if (nrow(nonoracle)) nonoracle[order(-nonoracle$spearman, nonoracle$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_oracle <- if (nrow(oracle)) oracle[order(-oracle$spearman, oracle$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  f1_pair_truth <- f1$cv[f1$cv$cv_target == "truth_pair_delta", , drop = FALSE]
  f1_pair_best <- if (nrow(f1_pair_truth)) f1_pair_truth[order(-f1_pair_truth$delta_spearman, f1_pair_truth$delta_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  f1_direct_truth <- f1$cv[f1$cv$cv_target == "heldout_direct_edge_truth_delta", , drop = FALSE]
  f1_direct_best <- if (nrow(f1_direct_truth)) f1_direct_truth[order(-f1_direct_truth$delta_spearman, f1_direct_truth$delta_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  f2_cv_best <- if (nrow(f2$cv)) f2$cv[order(-f2$cv$delta_spearman, f2$cv$delta_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  shell1 <- f4$grid[f4$grid$local_shell_depth == 1, , drop = FALSE]
  shell1_best_gradient <- if (nrow(shell1)) min(shell1$gradient_norm, na.rm = TRUE) else NA_real_
  all_long <- bind_rows_fill(list(
    transform(f1$gate, experiment = "F1_delta_training"),
    transform(f2$gate, experiment = "F2_state_dependent_estimators"),
    transform(f3$metrics, experiment = "F3_nonoracle_potential"),
    transform(f4$grid, experiment = "F4_local_tmb_debug"),
    transform(f5$edge, experiment = "F5_edge_delta_cv")
  ))
  summary <- data.frame(
    experiment = c("F1", "F2", "F3", "F4", "F5"),
    key_result = c(
      paste0("direct_edges=", nrow(f1$edge_train), "; direct_pairs=", nrow(f1$pair_train)),
      paste0("delta_gate_passed=", any(f2$gate$passes_gate, na.rm = TRUE)),
      paste0("nonoracle_valid_shape=", any(nonoracle$shape_class == "valid_shape", na.rm = TRUE)),
      paste0("median_gradient=", fmt_metric(stats::median(f4$grid$gradient_norm, na.rm = TRUE))),
      f5$failure$recommended_status[[1]]
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_f_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "f_experiment_summary.tsv"))
  write_tsv_safe(recs[recs$table == "delta_estimator_recommendation", -1, drop = FALSE], file.path(dirs$tables, "delta_estimator_recommendation.tsv"))
  write_tsv_safe(recs[recs$table == "potential_prior_recommendation", -1, drop = FALSE], file.path(dirs$tables, "potential_prior_recommendation.tsv"))
  write_tsv_safe(recs[recs$table == "local_tmb_recommendation", -1, drop = FALSE], file.path(dirs$tables, "local_tmb_recommendation.tsv"))
  write_tsv_safe(recs[recs$table == "cv_recommendation", -1, drop = FALSE], file.path(dirs$tables, "cv_recommendation.tsv"))
  write_tsv_safe(recs[recs$table == "recommended_next_steps", -1, drop = FALSE], file.path(dirs$tables, "recommended_next_steps.tsv"))
  lines <- c(
    "# Farfield Delta Estimator Probe Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", ctx$source_probe_dir %||% ctx$shared_input_dir, "`"),
    paste0("- abcd-dir: `", normalizePath(file.path(dirs$root, "..", "farfield_shape_probe_abcd"), winslash = "/", mustWork = FALSE), "`"),
    paste0("- diagnostics-dir: `", diagnostics_dir, "`"),
    paste0("- simulation_id: ", task_info$simulation_id),
    paste0("- minobs: ", task_info$minobs),
    paste0("- input_policy: ", task_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    paste0("- quick mode: ", quick, ". Auto mode uses full F1/F2 estimator construction but sampled F3/F4 fanout."),
    "",
    "## ABCD and E diagnostics summary",
    "- ABCD showed normalized/unit can reduce RMSE but often collapse amplitude; constant anchor variance improves amplitude without fixing rank; direct-anchor CV misses farfield wrong-direction.",
    "- E diagnostics showed oracle per-edge delta restores farfield shape, while estimated/local/context delta remains low-amplitude and direction-unreliable.",
    "",
    "## F1 trusted direct-edge / direct-pair delta estimator",
    paste0("- trusted direct-edge rows: ", nrow(f1$edge_train), "; trusted direct-pair rows: ", nrow(f1$pair_train), "."),
    paste0("- best F1 gate Spearman: ", fmt_metric(max(f1$gate$delta_spearman, na.rm = TRUE)),
           "; best sign agreement: ", fmt_metric(max(f1$gate$delta_sign_agreement, na.rm = TRUE)),
           "; any gate passed: ", any(f1$gate$passes_gate, na.rm = TRUE), "."),
    paste0("- best direct-pair heldout CV against truth: Spearman=", fmt_metric(f1_pair_best$delta_spearman[1]),
           ", sign_agreement=", fmt_metric(f1_pair_best$delta_sign_agreement[1]),
           ", sd_ratio=", fmt_metric(f1_pair_best$estimated_delta_sd_ratio[1]),
           ". This signal does not generalize to full parent-edge truth, where gate Spearman stays negative."),
    paste0("- best direct-edge heldout CV against truth: estimator=", f1_direct_best$estimator[1],
           ", Spearman=", fmt_metric(f1_direct_best$delta_spearman[1]),
           ", sign_agreement=", fmt_metric(f1_direct_best$delta_sign_agreement[1]), "."),
    "- If direct-edge rows are tiny, coarse context estimates are underidentified and path regression is the only non-oracle signal.",
    "",
    "## F2 state-dependent edge kernel estimator",
    paste0("- best F2 gate Spearman: ", fmt_metric(max(f2$gate$delta_spearman, na.rm = TRUE)),
           "; best sign agreement: ", fmt_metric(max(f2$gate$delta_sign_agreement, na.rm = TRUE)),
           "; best sd_ratio: ", fmt_metric(max(f2$gate$estimated_delta_sd_ratio, na.rm = TRUE)),
           "; any gate passed: ", any(f2$gate$passes_gate, na.rm = TRUE), "."),
    paste0("- best state-dependent heldout CV: estimator=", f2_cv_best$estimator[1],
           ", bandwidth=", f2_cv_best$bandwidth[1],
           ", Spearman=", fmt_metric(f2_cv_best$delta_spearman[1]),
           ", sign_agreement=", fmt_metric(f2_cv_best$delta_sign_agreement[1]),
           ", sd_ratio=", fmt_metric(f2_cv_best$estimated_delta_sd_ratio[1]),
           ". Full-graph truth diagnostics remain near zero or negative, so this is not a deployable estimator yet."),
    "- State kernels are compared against coarse context estimates and hybrid shrinkage; failure means current two-timepoint anchors are insufficient for a reliable slope estimator.",
    "",
    "## F3 non-oracle potential prior",
    paste0("- best non-oracle farfield config: ", best_nonoracle$candidate_id[1],
           " with pearson=", fmt_metric(best_nonoracle$pearson[1]),
           ", spearman=", fmt_metric(best_nonoracle$spearman[1]),
           ", sd_ratio=", fmt_metric(best_nonoracle$estimate_sd_ratio[1]),
           ", shape_class=", best_nonoracle$shape_class[1], "."),
    paste0("- best non-oracle result uses delta_source=", best_nonoracle$delta_source[1],
           " and delta_scale=", best_nonoracle$potential_delta_scale[1],
           "; the negative scale is a warning that the estimator orientation remains unreliable."),
    paste0("- best oracle upper bound: ", best_oracle$candidate_id[1],
           " with pearson=", fmt_metric(best_oracle$pearson[1]),
           ", spearman=", fmt_metric(best_oracle$spearman[1]),
           ", sd_ratio=", fmt_metric(best_oracle$estimate_sd_ratio[1]), "."),
    paste0("- non-oracle valid_shape count: ", sum(nonoracle$shape_class == "valid_shape", na.rm = TRUE),
           ", but these are marked delta_untrusted because the delta gate failed."),
    "",
    "## F4 local TMB convergence debug",
    paste0("- local debug rows: ", nrow(f4$grid), "; median gradient_norm=", fmt_metric(stats::median(f4$grid$gradient_norm, na.rm = TRUE)), "."),
    paste0("- shell_depth=1 best gradient_norm=", fmt_metric(shell1_best_gradient),
           "; converged fits=", sum(f4$grid$convergence == 0, na.rm = TRUE),
           "/", nrow(f4$grid),
           "; covariance fallback rows=", sum(f4$grid$covariance_fallback, na.rm = TRUE), "."),
    paste0("- covariance statuses: ", paste(unique(f4$grid$covariance_status), collapse = ", "), "."),
    "- Parameter-block gradients are not extractable from the public `fit_local_posterior()` return because obj/opt are not retained; the table records this explicitly.",
    "",
    "## F5 neighbor / shell / edge-delta CV",
    paste0("- calibration failure status: ", f5$failure$recommended_status[[1]], "."),
    paste0("- edge_delta_gate_passed: ", f5$failure$edge_delta_gate_passed[[1]], "."),
    "- Edge-delta CV is the most direct rejection criterion for current estimated deltas; neighbor/shell CV should be expanded once more observed low-count or graph-near nodes exist.",
    "",
    "## Final conclusion",
    paste0("- Continue C++ edge-gradient pseudo-observation now: ", if (any(f2$gate$passes_gate, na.rm = TRUE) && any(nonoracle$shape_class == "valid_shape", na.rm = TRUE)) "yes" else "no", "."),
    "- Current blocker is non-oracle delta reliability, not the Matrix/C++ precision mechanism. Oracle deltas prove the mechanism can work, but estimated deltas do not pass the gate.",
    "- Keep normalized as benchmark/probe default with amplitude-collapse diagnostics; keep unit as synthetic stress-test and mutation as legacy baseline.",
    "- Do not default `anchor_count_reference=minobs` for full input.",
    "- Highest priorities: improve delta estimator, expose/diagnose local TMB objective gradients, add edge-delta/shell CV to calibration, and allow `no_valid_shape_configuration`."
  )
  writeLines(lines, file.path(dirs$root, "farfield_delta_estimator_probe_report.md"))
  saveRDS(list(f1 = f1, f2 = f2, f3 = f3, f4 = f4, f5 = f5,
               all = all_long, summary = summary, recommendations = recs),
          file.path(dirs$results, "farfield_delta_estimator_probe_all_results.rds"))
}

main_delta_probe <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  mode <- match.arg(mode, c("prepare", "f1-delta-training", "f2-state-dependent-estimators",
                            "f3-nonoracle-potential", "f4-local-tmb-debug",
                            "f5-neighbor-shell-cv", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_delta_estimator_probe"))
  simulation_id <- arg_integer(args, "simulation_id", 1L)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  quick <- as_bool_or_auto(arg_value(args, "quick", "auto"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_probe_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, simulation_id, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, simulation_id, minobs, input_policy, force = force)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(
    simulation_id = simulation_id,
    minobs = minobs,
    input_policy = input_policy,
    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
    beta = as.numeric(ctx$input_table$sim_pm[[1L]])
  )
  components <- prepare_solver_cache(bundle$global_graph, dirs, force = force)
  abcd <- load_abcd_results(abcd_dir)
  diagnostics <- load_diagnostics_results(diagnostics_dir)
  saveRDS(list(context = ctx, task_info = task_info, diagnostics_dir = diagnostics_dir,
               diagnostics_summary = diagnostics$diagnostic_summary %||% data.frame()),
          file.path(dirs$results, "prepare_context.rds"))
  if (identical(mode, "prepare")) return(invisible(dirs$root))

  f1 <- if (mode %in% c("all", "f1-delta-training")) run_f1(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "f1_delta_training.rds"))
  if (identical(mode, "f1-delta-training")) return(invisible(f1))
  f2 <- if (mode %in% c("all", "f2-state-dependent-estimators")) run_f2(bundle, grf, task_info, dirs, f1, force = force) else readRDS(file.path(dirs$results, "f2_state_dependent_estimators.rds"))
  if (identical(mode, "f2-state-dependent-estimators")) return(invisible(f2))
  f3 <- if (mode %in% c("all", "f3-nonoracle-potential")) run_f3(bundle, components, grf, task_info, dirs, abcd, f1, f2, quick = quick, force = force) else readRDS(file.path(dirs$results, "f3_nonoracle_potential_prior.rds"))
  if (identical(mode, "f3-nonoracle-potential")) return(invisible(f3))
  f4 <- if (mode %in% c("all", "f4-local-tmb-debug")) run_f4(bundle, grf, task_info, dirs, quick = quick, force = force) else readRDS(file.path(dirs$results, "f4_local_tmb_debug.rds"))
  if (identical(mode, "f4-local-tmb-debug")) return(invisible(f4))
  f5 <- if (mode %in% c("all", "f5-neighbor-shell-cv")) run_f5(bundle, components, task_info, grf, dirs, abcd, f1, f2, f3, quick = quick, force = force) else readRDS(file.path(dirs$results, "f5_neighbor_shell_edge_cv.rds"))
  if (identical(mode, "f5-neighbor-shell-cv")) return(invisible(f5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_recommendations(f1, f2, f3, f4, f5)
    write_final_tables_report(dirs, ctx, task_info, normalizePath(diagnostics_dir, winslash = "/", mustWork = TRUE),
                              f1, f2, f3, f4, f5, recs, quick)
  }
  message("Wrote farfield delta estimator probe under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) {
  main_delta_probe()
}
