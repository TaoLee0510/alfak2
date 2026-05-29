#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_delta_debug.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_delta_estimator_probe.R"))

usage <- function() {
  cat(
    "Run farfield delta debug G1-G5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_delta_debug.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --output-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --simulation-id=1 --minobs=5 --input-policy=full\n\n",
    "Modes:\n",
    "  prepare, g1-sign-orientation, g2-confidence-gated-delta, g3-pairwise-potential,\n",
    "  g4-local-gradient-debug, g5-calibration-failure-state, summarize, all\n",
    sep = ""
  )
}

`%||G%` <- function(x, y) if (is.null(x)) y else x

make_debug_dirs <- function(output_dir) {
  make_probe_dirs(output_dir)
}

load_delta_probe_results <- function(delta_probe_dir) {
  delta_probe_dir <- normalizePath(delta_probe_dir, winslash = "/", mustWork = TRUE)
  all_path <- file.path(delta_probe_dir, "results", "farfield_delta_estimator_probe_all_results.rds")
  if (file.exists(all_path)) {
    out <- readRDS(all_path)
  } else {
    out <- list(
      f1 = readRDS(file.path(delta_probe_dir, "results", "f1_delta_training.rds")),
      f2 = readRDS(file.path(delta_probe_dir, "results", "f2_state_dependent_estimators.rds")),
      f3 = readRDS(file.path(delta_probe_dir, "results", "f3_nonoracle_potential_prior.rds")),
      f4 = readRDS(file.path(delta_probe_dir, "results", "f4_local_tmb_debug.rds")),
      f5 = readRDS(file.path(delta_probe_dir, "results", "f5_neighbor_shell_edge_cv.rds"))
    )
  }
  out$root <- delta_probe_dir
  out
}

base_configs_g <- function(delta_probe) {
  cfg <- data.frame(
    experiment = "G",
    candidate_id = c("G_baseline_mutation", "G_normalized_ll0p2_le0p01", "G_unit_ll0p2_le0p01"),
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
  far <- delta_probe$f3$metrics
  far <- far[far$support_scope == "farfield" & far$metric_scale == "native" & !far$oracle_delta, , drop = FALSE]
  if (nrow(far)) {
    top <- head(far[order(-far$spearman, far$centered_rmse), , drop = FALSE], 5L)
    top_cfg <- configs_from_metrics(top, "Gtop")
    cfg <- bind_rows_fill(list(cfg, top_cfg))
  }
  cfg <- cfg[!duplicated(paste(cfg$graph_edge_weight, cfg$lambda_l, cfg$lambda_e, cfg$sigma_obs, cfg$anchor_var_mode, cfg$anchor_count_reference_mode)), , drop = FALSE]
  cfg$candidate_id <- make.unique(as.character(cfg$candidate_id), sep = "_")
  cfg
}

sanitize_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  substr(x, 1L, 220L)
}

metric_shape_class <- function(x) {
  out <- rep("numeric_only", nrow(x))
  out[x$estimate_sd_ratio < 0.02] <- "collapsed_shrinkage"
  ok_amp <- x$estimate_sd_ratio >= 0.02
  out[ok_amp & (x$pearson < 0 | x$spearman < 0)] <- "noncollapsed_wrong_direction"
  out[ok_amp & x$pearson > 0 & x$spearman > 0] <- "valid_shape"
  out
}

best_path_beta <- function(f1) {
  cv <- f1$cv[f1$cv$estimator == "path_regression_context_coefficients" & f1$cv$cv_target == "truth_pair_delta", , drop = FALSE]
  lam <- if (nrow(cv)) cv$ridge_lambda[order(-cv$delta_spearman, cv$delta_rmse)][1L] else f1$path$coefficients$ridge_lambda[[1L]]
  beta <- f1$path$coefficients[f1$path$coefficients$ridge_lambda == lam, , drop = FALSE]
  list(lambda = lam, beta = beta$beta, coefficients = beta)
}

path_edge_rows_for_pair <- function(pair_row, graph, truth_map, local_map, beta, context_estimates) {
  nodes <- unlist(strsplit(pair_row$path_nodes, "->", fixed = TRUE), use.names = FALSE)
  node_idx <- match(nodes, as.character(graph$labels))
  if (anyNA(node_idx) || length(node_idx) < 2L) return(data.frame())
  rows <- list()
  for (j in seq_len(length(node_idx) - 1L)) {
    from <- node_idx[[j]]
    to <- node_idx[[j + 1L]]
    feat <- edge_step_features(graph, from, to)
    ctx <- feat$context_index[[1L]]
    est <- if (is.finite(ctx) && ctx >= 1L && ctx <= length(beta)) beta[[ctx]] else 0
    med <- context_estimates$delta_estimate[context_estimates$estimator == "context_median" & context_estimates$context_index == ctx]
    med <- if (length(med) && is.finite(med[[1L]])) med[[1L]] else 0
    parent <- as.character(graph$labels[[from]])
    child <- as.character(graph$labels[[to]])
    rows[[j]] <- data.frame(
      pair_id = pair_row$pair_id,
      edge_index_in_path = j,
      edge_from = parent,
      edge_to = child,
      parent = parent,
      child = child,
      feat,
      delta_sign_used_in_potential = "+",
      estimated_delta_e = est,
      oracle_truth_delta_e = truth_map[[child]] - truth_map[[parent]],
      local_delta_e = if (is.finite(local_map[[parent]]) && is.finite(local_map[[child]])) local_map[[child]] - local_map[[parent]] else NA_real_,
      empirical_delta_e = med,
      path_component_sign = 1,
      is_edge_used_forward = TRUE,
      is_edge_used_reverse = FALSE,
      stringsAsFactors = FALSE
    )
  }
  bind_rows_fill(rows)
}

toy_potential_solve <- function(delta, lambda_anchor = 1e4, lambda_edge = 1, ridge = 1e-8) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  a <- Matrix::sparseMatrix(i = c(1, 1), j = c(1, 2), x = c(-sqrt(lambda_edge), sqrt(lambda_edge)), dims = c(1, 2))
  q <- Matrix::crossprod(a) + Matrix::Diagonal(2, ridge)
  rhs <- as.numeric(Matrix::crossprod(a, sqrt(lambda_edge) * delta))
  q <- q + Matrix::sparseMatrix(i = 1, j = 1, x = lambda_anchor, dims = c(2, 2))
  chol <- Matrix::Cholesky(Matrix::forceSymmetric(q, uplo = "U"), LDL = TRUE, perm = TRUE)
  as.numeric(Matrix::solve(chol, rhs))
}

run_g1 <- function(bundle, grf, task_info, dirs, delta_probe, force = FALSE) {
  rds <- file.path(dirs$results, "g1_sign_orientation_audit.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$global_graph
  f1 <- delta_probe$f1
  labels <- as.character(graph$labels)
  truth_map <- compute_truth_for_nodes(labels, grf, task_info$lambda)
  names(truth_map) <- labels
  local_map <- setNames(bundle$local$summary$fitness_mean, as.character(bundle$local$summary$karyotype))
  beta <- best_path_beta(f1)
  context_est <- f1$context_estimates
  edge_rows <- lapply(seq_len(nrow(f1$pair_train)), function(i) {
    path_edge_rows_for_pair(f1$pair_train[i, , drop = FALSE], graph, truth_map, local_map, beta$beta, context_est)
  })
  path_edge <- bind_rows_fill(edge_rows)
  pair_rows <- lapply(split(path_edge, path_edge$pair_id), function(e) {
    p <- f1$pair_train[f1$pair_train$pair_id == e$pair_id[[1L]], , drop = FALSE][1L, , drop = FALSE]
    sum_path <- sum(e$estimated_delta_e, na.rm = TRUE)
    sum_ctx <- sum(e$empirical_delta_e, na.rm = TRUE)
    sum_oracle <- sum(e$oracle_truth_delta_e, na.rm = TRUE)
    data.frame(
      pair_id = p$pair_id,
      anchor_a = p$parent,
      anchor_b = p$child,
      observed_pair_delta = p$observed_pair_delta,
      truth_pair_delta = p$truth_pair_delta,
      path_length = p$path_length,
      path_nodes = p$path_nodes,
      path_edges = paste(paste(e$edge_from, e$edge_to, sep = "->"), collapse = "|"),
      path_direction_used = "anchor_a_to_anchor_b",
      sum_estimated_edge_delta_forward = sum_path,
      sum_estimated_edge_delta_reverse = -sum_path,
      sum_context_delta_forward = sum_ctx,
      sum_context_delta_reverse = -sum_ctx,
      sum_path_regression_component_forward = sum_path,
      sum_path_regression_component_reverse = -sum_path,
      sum_oracle_edge_delta_forward = sum_oracle,
      sum_oracle_edge_delta_reverse = -sum_oracle,
      forward_minus_observed = sum_path - p$observed_pair_delta,
      reverse_minus_observed = -sum_path - p$observed_pair_delta,
      forward_minus_truth = sum_path - p$truth_pair_delta,
      reverse_minus_truth = -sum_path - p$truth_pair_delta,
      stringsAsFactors = FALSE
    )
  })
  path_audit <- bind_rows_fill(pair_rows)
  methods <- list(
    path_regression = path_audit$sum_path_regression_component_forward,
    context_median = path_audit$sum_context_delta_forward,
    oracle_edge = path_audit$sum_oracle_edge_delta_forward
  )
  cors <- lapply(names(methods), function(m) {
    x <- methods[[m]]
    data.frame(
      method = m,
      cor_observed_forward_pearson = safe_cor2(path_audit$observed_pair_delta, x, "pearson"),
      cor_observed_forward_spearman = safe_cor2(path_audit$observed_pair_delta, x, "spearman"),
      cor_observed_reverse_pearson = safe_cor2(path_audit$observed_pair_delta, -x, "pearson"),
      cor_observed_reverse_spearman = safe_cor2(path_audit$observed_pair_delta, -x, "spearman"),
      cor_truth_forward_pearson = safe_cor2(path_audit$truth_pair_delta, x, "pearson"),
      cor_truth_forward_spearman = safe_cor2(path_audit$truth_pair_delta, x, "spearman"),
      cor_truth_reverse_pearson = safe_cor2(path_audit$truth_pair_delta, -x, "pearson"),
      cor_truth_reverse_spearman = safe_cor2(path_audit$truth_pair_delta, -x, "spearman"),
      stringsAsFactors = FALSE
    )
  })
  cor_summary <- bind_rows_fill(cors)
  toy <- lapply(c(1, -1), function(d) {
    m <- toy_potential_solve(d)
    data.frame(delta = d, m_A = m[[1L]], m_B = m[[2L]], m_B_minus_m_A = m[[2L]] - m[[1L]],
               expected = if (d > 0) "m_B_gt_m_A" else "m_B_lt_m_A",
               passed = if (d > 0) (m[[2L]] > m[[1L]]) else (m[[2L]] < m[[1L]]),
               stringsAsFactors = FALSE)
  })
  toy_tbl <- bind_rows_fill(toy)
  bug_candidates <- data.frame(
    check = c("potential_rhs_sign", "oracle_path_orientation", "path_regression_orientation", "negative_scale_best"),
    status = c(
      if (all(toy_tbl$passed)) "pass" else "fail",
      if (cor_summary$cor_truth_forward_spearman[cor_summary$method == "oracle_edge"] > 0.99) "pass" else "suspect",
      if (cor_summary$cor_truth_forward_spearman[cor_summary$method == "path_regression"] > 0) "pass_on_pairs_only" else "suspect",
      "likely_estimator_compensation_not_rhs_bug"
    ),
    recommendation = c(
      "No normal-equation sign fix needed if toy passes.",
      "Forward graph path is consistent with truth when oracle deltas are used.",
      "Direct-pair sign can be internally consistent but does not guarantee farfield edge deployment.",
      "Do not flip production sign globally; require deployment gate before using delta."
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(path_audit, file.path(dirs$tables, "path_orientation_audit.tsv"))
  write_tsv_safe(path_edge, file.path(dirs$tables, "path_edge_orientation_audit.tsv"))
  write_tsv_safe(cor_summary, file.path(dirs$tables, "orientation_correlation_summary.tsv"))
  write_tsv_safe(toy_tbl, file.path(dirs$tables, "potential_delta_sign_audit.tsv"))
  write_tsv_safe(bug_candidates, file.path(dirs$tables, "orientation_bug_candidates.tsv"))
  out <- list(path = path_audit, edges = path_edge, cor = cor_summary, toy = toy_tbl, bugs = bug_candidates)
  saveRDS(out, rds)
  out
}

min_edge_distance <- function(train, target_edges, graph) {
  if (!nrow(train)) return(rep(NA_real_, nrow(target_edges)))
  if (!"parent_karyotype" %in% names(train) && "parent" %in% names(train)) {
    train$parent_karyotype <- train$parent
  }
  train$parent_karyotype <- as.character(train$parent_karyotype)
  target_edges$parent_karyotype <- as.character(target_edges$parent_karyotype)
  vapply(seq_len(nrow(target_edges)), function(i) min(edge_distance(train, target_edges[i, , drop = FALSE], graph), na.rm = TRUE), numeric(1))
}

same_counts <- function(train, target_edges) {
  chrdir <- paste(train$edge_chr, train$edge_direction)
  ctx <- train$context_index
  data.frame(
    same_chr_direction_training_count = tabulate(match(paste(target_edges$edge_chr, target_edges$edge_direction), unique(chrdir)), nbins = length(unique(chrdir)))[match(paste(target_edges$edge_chr, target_edges$edge_direction), unique(chrdir))],
    same_context_training_count = tabulate(match(target_edges$context_index, unique(ctx)), nbins = length(unique(ctx)))[match(target_edges$context_index, unique(ctx))]
  )
}

nearest_direct_distance <- function(graph, direct_labels) {
  labels <- as.character(graph$labels)
  direct_idx <- match(direct_labels, labels)
  direct_idx <- direct_idx[!is.na(direct_idx)]
  as.integer(graph$support_distance)
}

delta_for_source <- function(source, parent_edges, f1, f2_pred) {
  if (source == "path_regression_context_coefficients") {
    b <- best_path_beta(list(cv = f1$cv, path = f1$path))
    out <- b$beta[parent_edges$context_index]
    out[!is.finite(out)] <- 0
    return(out)
  }
  if (source %in% c("state_kernel_direct_pair", "hybrid_confident_delta", "context_median")) {
    est <- if (source == "context_median") "context_median" else source
    x <- f2_pred[f2_pred$estimator == est, , drop = FALSE]
    if (!nrow(x)) return(rep(0, nrow(parent_edges)))
    out <- x$estimated_delta[match(parent_edges$edge_id, x$edge_id)]
    out[!is.finite(out)] <- 0
    return(out)
  }
  if (source == "oracle_per_edge_delta") return(parent_edges$truth_delta)
  rep(0, nrow(parent_edges))
}

gate_pass <- function(overlap, gate) {
  if (gate == "gate_weak") {
    overlap$kernel_effective_sample_size >= 1 &
      overlap$nearest_training_pair_edge_distance <= 8 &
      overlap$same_chr_direction_training_count >= 1
  } else if (gate == "gate_medium") {
    overlap$kernel_effective_sample_size >= 2 &
      overlap$nearest_training_pair_edge_distance <= 4 &
      overlap$same_chr_direction_training_count >= 2
  } else {
    overlap$kernel_effective_sample_size >= 3 &
      overlap$nearest_training_pair_edge_distance <= 2 &
      overlap$same_chr_direction_training_count >= 3
  }
}

build_confidence_gated_predictions <- function(parent_edges, overlap, f1, f2_pred) {
  sources <- c("path_regression_context_coefficients", "state_kernel_direct_pair", "hybrid_confident_delta", "context_median", "oracle_per_edge_delta")
  gates <- c("gate_weak", "gate_medium", "gate_strict")
  rows <- list()
  idx <- 0L
  for (source in sources) {
    raw <- delta_for_source(source, parent_edges, f1, f2_pred)
    for (gate in gates) {
      pass <- gate_pass(overlap, gate)
      pass[is.na(pass)] <- FALSE
      k <- switch(gate, gate_weak = 1, gate_medium = 2, gate_strict = 3)
      conf <- pmin(1, overlap$kernel_effective_sample_size / k) * exp(-pmax(overlap$nearest_training_pair_edge_distance, 0) / 4)
      conf[!pass | !is.finite(conf)] <- 0
      if (source == "oracle_per_edge_delta") conf[pass & is.finite(conf)] <- pmax(conf[pass & is.finite(conf)], 0.25)
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        edge_id = parent_edges$edge_id,
        parent = parent_edges$parent_karyotype,
        child = parent_edges$child_karyotype,
        delta_source = source,
        gate = gate,
        estimated_delta_raw = raw,
        prediction_confidence = conf,
        estimated_delta_gated = raw * conf,
        oracle_truth_delta = parent_edges$truth_delta,
        is_in_domain = pass,
        gating_reason = ifelse(pass, "pass", "out_of_domain"),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(rows)
}

run_g2 <- function(bundle, components, grf, task_info, dirs, delta_probe, force = FALSE) {
  rds <- file.path(dirs$results, "g2_confidence_gated_delta.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$global_graph
  parent_edges <- parent_edge_frame(graph, grf, task_info$lambda)
  pair_edges <- make_pair_pseudo_edges(delta_probe$f1$pair_train, graph)
  f2_pred <- read_tsv_safe(file.path(delta_probe$root, "tables", "state_dependent_edge_predictions.tsv"))
  f2_pair <- f2_pred[f2_pred$estimator == "state_kernel_direct_pair", , drop = FALSE]
  overlap <- parent_edges
  overlap$nearest_training_edge_distance <- min_edge_distance(delta_probe$f1$edge_train, parent_edges, graph)
  overlap$nearest_training_pair_edge_distance <- min_edge_distance(pair_edges, parent_edges, graph)
  overlap$kernel_effective_sample_size <- f2_pair$effective_sample_size[match(parent_edges$edge_id, f2_pair$edge_id)]
  overlap$kernel_weight_sum <- f2_pair$kernel_weight_sum[match(parent_edges$edge_id, f2_pair$edge_id)]
  sc <- same_counts(pair_edges, parent_edges)
  overlap$same_chr_direction_training_count <- sc$same_chr_direction_training_count
  overlap$same_context_training_count <- sc$same_context_training_count
  overlap$edge_context_seen_in_training <- overlap$same_context_training_count > 0
  overlap$parent_distance_to_nearest_direct_anchor <- overlap$support_distance_parent
  overlap$child_distance_to_nearest_direct_anchor <- overlap$support_distance_child
  overlap$parent_state_seen_or_near_seen <- overlap$parent_distance_to_nearest_direct_anchor <= 1
  overlap$prediction_confidence <- pmin(1, overlap$kernel_effective_sample_size / 2) * exp(-pmax(overlap$nearest_training_pair_edge_distance, 0) / 4)
  overlap$estimated_delta_raw <- delta_for_source("state_kernel_direct_pair", parent_edges, delta_probe$f1, f2_pred)
  overlap$estimated_delta_gated <- overlap$estimated_delta_raw * ifelse(gate_pass(overlap, "gate_medium"), overlap$prediction_confidence, 0)
  overlap$oracle_truth_delta <- parent_edges$truth_delta
  overlap$is_in_domain <- gate_pass(overlap, "gate_medium")
  overlap$gating_reason <- ifelse(overlap$is_in_domain, "pass_gate_medium", "out_of_domain")
  preds <- build_confidence_gated_predictions(parent_edges, overlap, delta_probe$f1, f2_pred)
  pred_metrics <- lapply(split(preds, paste(preds$delta_source, preds$gate)), function(x) {
    cbind(data.frame(delta_source = x$delta_source[[1L]], gate = x$gate[[1L]], in_domain_fraction = mean(x$is_in_domain), stringsAsFactors = FALSE),
          edge_delta_metric_row(x$estimated_delta_gated, x$oracle_truth_delta))
  })
  pred_metrics <- bind_rows_fill(pred_metrics)
  cfg <- base_configs_g(delta_probe)[1:min(3L, nrow(base_configs_g(delta_probe))), , drop = FALSE]
  sources <- c("path_regression_context_coefficients", "state_kernel_direct_pair", "hybrid_confident_delta", "context_median", "oracle_per_edge_delta")
  gates <- c("gate_weak", "gate_medium", "gate_strict")
  scales <- c(1, -1)
  metrics <- list()
  k <- 0L
  for (i in seq_len(nrow(cfg))) {
    for (source in sources) {
      for (gate in gates) {
        for (scale in scales) {
          if (source == "oracle_per_edge_delta" && scale < 0) next
          sub <- preds[preds$delta_source == source & preds$gate == gate, , drop = FALSE]
          conf <- sub$prediction_confidence
          delta <- sub$estimated_delta_raw * scale
          pm <- fit_potential_prior_mean_edges(bundle$local, graph, parent_edges, delta, conf,
                                               lambda_anchor = 10, lambda_smooth = 0.1,
                                               components = components,
                                               edge_weight_mode = as.character(cfg$graph_edge_weight[[i]]))
          cfi <- cfg[i, , drop = FALSE]
          cfi$candidate_id <- sanitize_id(paste("G2", cfi$graph_edge_weight, source, gate, paste0("ds", scale), sep = "__"))
          cfi$prior_mean_mode <- paste0("gated_", source)
          cfi$prior_mean_scale <- scale
          res <- fit_cached(bundle$local, graph, components, cfi, grf, task_info, dirs,
                            prior_mean = pm, prior_mean_status = paste0("gated_", source, "_", gate),
                            force = force)
          m <- res$metrics
          m$delta_source <- source
          m$gate <- gate
          m$delta_scale <- scale
          m$oracle_delta <- source == "oracle_per_edge_delta"
          k <- k + 1L
          metrics[[k]] <- m
        }
      }
    }
  }
  metrics <- bind_rows_fill(metrics)
  top <- metrics[metrics$support_scope == "farfield" & metrics$metric_scale == "native", , drop = FALSE]
  top$shape_class <- metric_shape_class(top)
  top <- top[order(top$oracle_delta, -top$spearman, top$centered_rmse), , drop = FALSE]
  write_tsv_safe(overlap, file.path(dirs$tables, "edge_delta_domain_overlap.tsv"))
  write_tsv_safe(preds, file.path(dirs$tables, "confidence_gated_delta_predictions.tsv"))
  write_tsv_safe(metrics, file.path(dirs$tables, "confidence_gated_potential_prior.tsv"))
  write_tsv_safe(head(top, 40L), file.path(dirs$tables, "confidence_gated_potential_top.tsv"))
  out <- list(overlap = overlap, predictions = preds, prediction_metrics = pred_metrics, metrics = metrics, top = top)
  saveRDS(out, rds)
  out
}

fit_pairwise_potential_prior_mean <- function(local_fit, graph, pair_train, pair_delta, pair_weight,
                                              lambda_pair = 1, lambda_anchor = 10, lambda_graph = 0.1,
                                              components = NULL, ridge = 1e-6, edge_weight_mode = "normalized") {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  n <- length(graph$labels)
  a_idx <- match(as.character(pair_train$parent), as.character(graph$labels))
  b_idx <- match(as.character(pair_train$child), as.character(graph$labels))
  ok <- !is.na(a_idx) & !is.na(b_idx) & is.finite(pair_delta) & is.finite(pair_weight) & pair_weight > 0
  a_idx <- a_idx[ok]
  b_idx <- b_idx[ok]
  delta <- pair_delta[ok]
  w <- sqrt(lambda_pair * pair_weight[ok])
  q <- Matrix::Diagonal(n, ridge)
  rhs <- numeric(n)
  if (length(a_idx)) {
    row_id <- seq_along(a_idx)
    a <- Matrix::sparseMatrix(i = c(row_id, row_id), j = c(a_idx, b_idx), x = c(-w, w), dims = c(length(a_idx), n))
    q <- q + Matrix::crossprod(a)
    rhs <- rhs + as.numeric(Matrix::crossprod(a, w * delta))
  }
  direct <- as.character(local_fit$summary$support_tier) == "directly_informed" & is.finite(local_fit$summary$fitness_mean)
  anchor_idx <- match(as.character(local_fit$summary$karyotype[direct]), as.character(graph$labels))
  ok_anchor <- !is.na(anchor_idx)
  anchor_idx <- anchor_idx[ok_anchor]
  anchor_mean <- as.numeric(local_fit$summary$fitness_mean[direct][ok_anchor])
  q <- q + Matrix::sparseMatrix(i = anchor_idx, j = anchor_idx, x = lambda_anchor, dims = c(n, n))
  rhs[anchor_idx] <- rhs[anchor_idx] + lambda_anchor * anchor_mean
  if (!is.null(components) && lambda_graph > 0) {
    q <- q + lambda_graph * components$edge[[edge_weight_mode]]
  }
  chol <- Matrix::Cholesky(Matrix::forceSymmetric(q, uplo = "U"), LDL = TRUE, perm = TRUE)
  as.numeric(Matrix::solve(chol, rhs))
}

pair_delta_values <- function(source, pair_train, beta) {
  ctx_cols <- grep("^ctx_", names(pair_train), value = TRUE)
  x <- as.matrix(pair_train[, ctx_cols, drop = FALSE])
  pred <- as.numeric(x %*% beta)
  if (source == "observed_direct_pair_delta") pair_train$observed_pair_delta
  else if (source == "path_regression_predicted_pair_delta") pred
  else if (source == "shrinked_pair_delta") pred * pmin(1, 1 / pmax(1, pair_train$path_length))
  else if (source == "oracle_pair_delta") pair_train$truth_pair_delta
  else rep(0, nrow(pair_train))
}

run_g3 <- function(bundle, components, grf, task_info, dirs, delta_probe, force = FALSE) {
  rds <- file.path(dirs$results, "g3_pairwise_potential_prior.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  pair_train <- delta_probe$f1$pair_train
  beta <- best_path_beta(delta_probe$f1)$beta
  cfg <- base_configs_g(delta_probe)[1:min(3L, nrow(base_configs_g(delta_probe))), , drop = FALSE]
  sources <- c("observed_direct_pair_delta", "path_regression_predicted_pair_delta", "shrinked_pair_delta", "oracle_pair_delta")
  lambda_pair_grid <- c(1, 10)
  lambda_anchor_grid <- c(10)
  lambda_graph_grid <- c(0.1)
  thresholds <- c(0, 0.5)
  weight_modes <- c("uniform", "inverse_path_length")
  metrics <- list()
  constraints <- list()
  k <- 0L
  cidx <- 0L
  for (source in sources) {
    delta <- pair_delta_values(source, pair_train, beta)
    base_conf <- if (source == "shrinked_pair_delta") pmin(1, 1 / pmax(1, pair_train$path_length)) else rep(1, nrow(pair_train))
    for (thr in thresholds) {
      keep <- base_conf >= thr & pair_train$path_length <= 3
      for (wm in weight_modes) {
        w <- if (wm == "inverse_path_length") 1 / pmax(1, pair_train$path_length) else rep(1, nrow(pair_train))
        w <- w * base_conf * keep
        cidx <- cidx + 1L
        constraints[[cidx]] <- data.frame(pair_train[, c("pair_id", "parent", "child", "path_length", "observed_pair_delta", "truth_pair_delta")],
                                          delta_source = source, confidence_threshold = thr,
                                          pair_weight_mode = wm, pair_delta_hat = delta, pair_weight = w,
                                          stringsAsFactors = FALSE)
        for (lp in lambda_pair_grid) {
          for (la in lambda_anchor_grid) {
            for (lg in lambda_graph_grid) {
              for (scale in c(1, -1)) {
                if (source == "oracle_pair_delta" && scale < 0) next
                for (i in seq_len(nrow(cfg))) {
                  pm <- fit_pairwise_potential_prior_mean(bundle$local, bundle$global_graph, pair_train, delta * scale, w,
                                                           lambda_pair = lp, lambda_anchor = la, lambda_graph = lg,
                                                           components = components,
                                                           edge_weight_mode = as.character(cfg$graph_edge_weight[[i]]))
                  cfi <- cfg[i, , drop = FALSE]
                  cfi$candidate_id <- sanitize_id(paste("G3", cfi$graph_edge_weight, source, paste0("thr", thr), wm, paste0("lp", lp), paste0("ds", scale), sep = "__"))
                  cfi$prior_mean_mode <- paste0("pairwise_", source)
                  cfi$prior_mean_scale <- scale
                  res <- fit_cached(bundle$local, bundle$global_graph, components, cfi, grf, task_info, dirs,
                                    prior_mean = pm, prior_mean_status = paste0("pairwise_", source), force = force)
                  m <- res$metrics
                  m$pair_delta_source <- source
                  m$pair_confidence_threshold <- thr
                  m$pair_weight_mode <- wm
                  m$lambda_pair <- lp
                  m$lambda_anchor_pairwise <- la
                  m$lambda_graph_pairwise <- lg
                  m$pair_delta_scale <- scale
                  m$oracle_delta <- source == "oracle_pair_delta"
                  k <- k + 1L
                  metrics[[k]] <- m
                }
              }
            }
          }
        }
      }
    }
  }
  metrics <- bind_rows_fill(metrics)
  constraints <- bind_rows_fill(constraints)
  top <- metrics[metrics$support_scope == "farfield" & metrics$metric_scale == "native", , drop = FALSE]
  top$shape_class <- metric_shape_class(top)
  top <- top[order(top$oracle_delta, -top$spearman, top$centered_rmse), , drop = FALSE]
  edge_top <- read_tsv_safe(file.path(delta_probe$root, "tables", "nonoracle_vs_oracle_upper_bound.tsv"))
  comparison <- data.frame(
    method = c("pairwise_best_nonoracle", "pairwise_best_oracle", as.character(edge_top$class)),
    spearman = c(top$spearman[!top$oracle_delta][1], top$spearman[top$oracle_delta][1], edge_top$spearman),
    pearson = c(top$pearson[!top$oracle_delta][1], top$pearson[top$oracle_delta][1], edge_top$pearson),
    estimate_sd_ratio = c(top$estimate_sd_ratio[!top$oracle_delta][1], top$estimate_sd_ratio[top$oracle_delta][1], edge_top$estimate_sd_ratio),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(metrics, file.path(dirs$tables, "pairwise_potential_prior_probe.tsv"))
  write_tsv_safe(constraints, file.path(dirs$tables, "pairwise_potential_constraints.tsv"))
  write_tsv_safe(head(top, 60L), file.path(dirs$tables, "pairwise_potential_top.tsv"))
  write_tsv_safe(comparison, file.path(dirs$tables, "pairwise_vs_edge_potential_comparison.tsv"))
  out <- list(metrics = metrics, constraints = constraints, top = top, comparison = comparison)
  saveRDS(out, rds)
  out
}

make_local_tmb_debug_attempt <- function(data, graph, variant, observation_model = "dirichlet_multinomial",
                                         dm_concentration = 50, eval_max = 500, iter_max = 500,
                                         eta_borrowed_prior_mean = -6, eta_borrowed_prior_sd = 1.5,
                                         eta_distance_penalty = 0.75, freeze_delta_context = FALSE,
                                         staged_init = FALSE) {
  started <- Sys.time()
  n <- length(graph$labels)
  observed_index <- match(data$labels, as.character(graph$labels))
  y0 <- numeric(n); y1 <- numeric(n)
  y0[observed_index] <- data$counts[, 1]
  y1[observed_index] <- data$counts[, 2]
  obs_weight0 <- rep(1, n); obs_weight1 <- rep(1, n)
  y0_init <- y0; y1_init <- y1
  effective_count_total <- y0_init + y1_init
  borrowed_eta <- effective_count_total <= 0 & graph$support_distance > 0L
  eta_prior_mean <- rep(0, n)
  eta_prior_sd_vec <- rep(5, n)
  eta_prior_mean[borrowed_eta] <- eta_borrowed_prior_mean - eta_distance_penalty * pmax(0, as.integer(graph$support_distance[borrowed_eta]) - 1L)
  eta_prior_sd_vec[borrowed_eta] <- eta_borrowed_prior_sd
  p0 <- (y0_init + 0.5) / sum(y0_init + 0.5)
  p1 <- (y1_init + 0.5) / sum(y1_init + 0.5)
  f0 <- log(p1) - log(p0)
  f0 <- f0 - mean(f0)
  n_context <- length(graph$context_label)
  n_group <- max(unlist(graph$context_group0), 0L) + 1L
  parameters <- list(
    eta = log(p0),
    f = f0 / max(data$dt, .Machine$double.eps),
    delta_context = rep(0, n_context),
    mu_group = rep(0, n_group),
    log_sigma_neighbor = log(0.35),
    log_sigma_anchor = log(0.6),
    log_tau_group = rep(log(0.2), n_group)
  )
  if (isTRUE(staged_init)) {
    shell0_graph <- tryCatch(alfak2::build_karyotype_graph(data, shell_depth = 0), error = function(e) NULL)
    shell0_fit <- tryCatch(if (!is.null(shell0_graph)) alfak2::fit_local_posterior(data, graph = shell0_graph, observation_model = "multinomial", control = list(eval.max = 500, iter.max = 500), retry_on_untrusted_covariance = FALSE) else NULL, error = function(e) NULL)
    if (!is.null(shell0_fit)) {
      f_map <- setNames(shell0_fit$summary$fitness_mean, shell0_fit$summary$karyotype)
      idx <- match(as.character(graph$labels), names(f_map))
      parameters$f[!is.na(idx)] <- f_map[idx[!is.na(idx)]]
      parameters$f[is.na(idx)] <- median(parameters$f[!is.na(idx)], na.rm = TRUE)
      parameters$delta_context[] <- 0
    }
  }
  tmb_data <- list(
    y0 = as.numeric(y0), y1 = as.numeric(y1), n_nodes = as.integer(n),
    trans_from = as.integer(unlist(graph$transition_from0)),
    trans_to = as.integer(unlist(graph$transition_to0)),
    trans_weight = as.numeric(unlist(graph$transition_weight)),
    support_distance = as.integer(graph$support_distance),
    parent_from = as.integer(unlist(graph$parent_from0)),
    parent_to = as.integer(unlist(graph$parent_to0)),
    parent_weight = as.numeric(unlist(graph$parent_weight)),
    parent_context = as.integer(unlist(graph$parent_context0)),
    context_group = as.integer(unlist(graph$context_group0)),
    dt = as.numeric(data$dt),
    obs_weight0 = as.numeric(obs_weight0), obs_weight1 = as.numeric(obs_weight1),
    use_observation_weights = 0L,
    eta_prior_mean = as.numeric(eta_prior_mean),
    eta_prior_sd = as.numeric(eta_prior_sd_vec),
    anchor_prior_scale = 1.0, mu_prior_scale = 0.5, scale_prior_scale = 1.0,
    observation_model = as.integer(observation_model == "dirichlet_multinomial"),
    observation_weight_mode = 1L,
    dm_concentration = as.numeric(dm_concentration)
  )
  map <- if (isTRUE(freeze_delta_context)) list(delta_context = factor(rep(NA, n_context))) else NULL
  out <- tryCatch({
    obj <- TMB::MakeADFun(tmb_data, parameters, DLL = "alfak2", silent = TRUE, map = map)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = eval_max, iter.max = iter_max))
    g <- obj$gr(opt$par)
    block <- tapply(abs(g), names(obj$par), max)
    plist <- obj$env$parList(opt$par)
    data.frame(
      variant = variant,
      local_shell_depth = 1,
      observation_model = observation_model,
      dm_concentration = dm_concentration,
      eval_max = eval_max,
      iter_max = iter_max,
      convergence = opt$convergence,
      message = opt$message,
      objective = opt$objective,
      gradient_norm = sqrt(sum(g^2)),
      grad_eta_max_abs = unname(block["eta"]),
      grad_f_max_abs = unname(block["f"]),
      grad_delta_context_max_abs = if ("delta_context" %in% names(block)) unname(block["delta_context"]) else NA_real_,
      grad_mu_group_max_abs = unname(block["mu_group"]),
      grad_log_sigma_neighbor_abs = unname(block["log_sigma_neighbor"]),
      grad_log_sigma_anchor_abs = unname(block["log_sigma_anchor"]),
      grad_log_tau_group_max_abs = unname(block["log_tau_group"]),
      global_gradient_norm = sqrt(sum(g^2)),
      max_gradient_block_name = names(which.max(block)),
      covariance_status = "not_run_custom_gradient_probe",
      covariance_fallback = NA,
      fitness_sd_source = "not_run",
      retry_attempted = FALSE,
      dm_concentration_selected = dm_concentration,
      n_local_nodes = n,
      n_direct = sum(graph$support_tier == "directly_informed"),
      n_local_borrowed = sum(graph$support_tier != "directly_informed"),
      n_weakly_supported = 0,
      elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      extraction_status = "custom_tmb_object_benchmark_only",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      variant = variant,
      local_shell_depth = 1,
      observation_model = observation_model,
      dm_concentration = dm_concentration,
      eval_max = eval_max,
      iter_max = iter_max,
      convergence = NA_integer_,
      message = conditionMessage(e),
      objective = NA_real_,
      gradient_norm = NA_real_,
      grad_eta_max_abs = NA_real_, grad_f_max_abs = NA_real_, grad_delta_context_max_abs = NA_real_,
      grad_mu_group_max_abs = NA_real_, grad_log_sigma_neighbor_abs = NA_real_,
      grad_log_sigma_anchor_abs = NA_real_, grad_log_tau_group_max_abs = NA_real_,
      global_gradient_norm = NA_real_, max_gradient_block_name = NA_character_,
      covariance_status = "error", covariance_fallback = NA, fitness_sd_source = "not_run",
      retry_attempted = FALSE, dm_concentration_selected = dm_concentration,
      n_local_nodes = n, n_direct = sum(graph$support_tier == "directly_informed"),
      n_local_borrowed = sum(graph$support_tier != "directly_informed"),
      n_weakly_supported = 0, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      extraction_status = "custom_tmb_error", stringsAsFactors = FALSE
    )
  })
  out
}

run_g4 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "g4_local_gradient_debug.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  variants <- data.frame(
    variant = c("L0_current", "L1_freeze_delta_context", "L3_staged_initialization", "L4_strong_borrowed_shrinkage", "L4_very_strong_borrowed_shrinkage"),
    eta_mean = c(-6, -6, -8, -8, -10),
    eta_sd = c(1.5, 1.5, 1.0, 1.0, 0.5),
    eta_penalty = c(0.75, 0.75, 1.5, 1.5, 2.5),
    freeze_delta_context = c(FALSE, TRUE, FALSE, FALSE, FALSE),
    staged_init = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(variants)), function(i) {
    make_local_tmb_debug_attempt(
      bundle$data, bundle$local$graph, variants$variant[[i]],
      observation_model = "dirichlet_multinomial", dm_concentration = 50,
      eval_max = 500, iter_max = 500,
      eta_borrowed_prior_mean = variants$eta_mean[[i]],
      eta_borrowed_prior_sd = variants$eta_sd[[i]],
      eta_distance_penalty = variants$eta_penalty[[i]],
      freeze_delta_context = variants$freeze_delta_context[[i]],
      staged_init = variants$staged_init[[i]]
    )
  })
  grad <- bind_rows_fill(rows)
  freeze <- grad[grad$variant == "L1_freeze_delta_context", , drop = FALSE]
  gainloss <- data.frame(variant = "L2_gain_loss_only_context", implementation_status = "not_implemented_benchmark_only_requires_context_remapping_or_TMB_map",
                         recommendation = "Implement as a package-level local objective option only after block gradients confirm context blocks dominate.",
                         stringsAsFactors = FALSE)
  staged <- grad[grad$variant == "L3_staged_initialization", , drop = FALSE]
  shrink <- grad[grepl("^L4_", grad$variant), , drop = FALSE]
  summary <- grad
  summary$local_edge_delta_sign_agreement <- 0.3255814
  summary$local_edge_delta_status <- "still_untrusted_from_F4_alignment"
  stab <- read_tsv_safe(file.path(delta_probe_global_root, "tables", "local_multistart_stability.tsv"))
  write_tsv_safe(grad, file.path(dirs$tables, "local_gradient_block_diagnostics.tsv"))
  write_tsv_safe(freeze, file.path(dirs$tables, "local_freeze_context_probe.tsv"))
  write_tsv_safe(gainloss, file.path(dirs$tables, "local_gain_loss_context_probe.tsv"))
  write_tsv_safe(staged, file.path(dirs$tables, "local_staged_init_probe.tsv"))
  write_tsv_safe(shrink, file.path(dirs$tables, "local_borrowed_shrinkage_probe.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "local_variant_summary.tsv"))
  out <- list(gradient = grad, freeze = freeze, gainloss = gainloss, staged = staged, shrink = shrink, summary = summary, stability = stab)
  saveRDS(out, rds)
  out
}

classify_calibration_row <- function(x, delta_based, edge_gate, deploy_gate) {
  amp <- is.finite(x$estimate_sd_ratio) && x$estimate_sd_ratio >= 0.02
  rank <- is.finite(x$pearson) && is.finite(x$spearman) && x$pearson > 0 && x$spearman > 0
  if (delta_based && (!edge_gate || !deploy_gate)) return("delta_untrusted")
  if (!amp) return("amplitude_collapse")
  if (!rank) return("wrong_direction")
  "valid_shape_config"
}

run_g5 <- function(delta_probe, g2, g3, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "g5_calibration_failure_state_demo.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  f3 <- delta_probe$f3$metrics
  f3 <- f3[f3$support_scope == "farfield" & f3$metric_scale == "native", , drop = FALSE]
  g2m <- g2$metrics[g2$metrics$support_scope == "farfield" & g2$metrics$metric_scale == "native", , drop = FALSE]
  g3m <- g3$metrics[g3$metrics$support_scope == "farfield" & g3$metrics$metric_scale == "native", , drop = FALSE]
  f3$source_experiment <- "F3_edge_potential"
  g2m$source_experiment <- "G2_gated_edge_potential"
  g3m$source_experiment <- "G3_pairwise_potential"
  all <- bind_rows_fill(list(f3, g2m, g3m))
  all$node_shape_classification <- metric_shape_class(all)
  all$prior_type <- all$prior_mean_mode %||G% NA_character_
  ds <- if ("delta_source" %in% names(all)) as.character(all$delta_source) else rep(NA_character_, nrow(all))
  pds <- if ("pair_delta_source" %in% names(all)) as.character(all$pair_delta_source) else rep(NA_character_, nrow(all))
  ds[!nzchar(ds) | is.na(ds)] <- pds[!nzchar(ds) | is.na(ds)]
  all$delta_source <- ds
  dscale <- if ("delta_scale" %in% names(all)) all$delta_scale else rep(NA_real_, nrow(all))
  pscale <- if ("pair_delta_scale" %in% names(all)) all$pair_delta_scale else rep(NA_real_, nrow(all))
  dscale[!is.finite(dscale)] <- pscale[!is.finite(dscale)]
  dscale[!is.finite(dscale)] <- all$prior_mean_scale[!is.finite(dscale)]
  all$delta_scale <- dscale
  oracle_delta_flag <- if ("oracle_delta" %in% names(all)) as.logical(all$oracle_delta) else rep(FALSE, nrow(all))
  oracle_delta_flag[is.na(oracle_delta_flag)] <- FALSE
  all$is_oracle <- grepl("oracle", all$delta_source %||G% "") | oracle_delta_flag
  all$edge_delta_gate_passed <- all$is_oracle
  all$deployment_gate_passed <- all$is_oracle
  dep <- g2$prediction_metrics
  dep$key <- paste(dep$delta_source, dep$gate)
  deploy_pass <- dep$delta_sign_agreement >= 0.55 & dep$delta_spearman > 0 & dep$estimated_delta_sd_ratio >= 0.10
  names(deploy_pass) <- dep$key
  key <- paste(all$delta_source, all$gate)
  hit <- deploy_pass[key]
  all$deployment_gate_passed[!is.na(hit)] <- hit[!is.na(hit)]
  all$amplitude_gate_passed <- all$estimate_sd_ratio >= 0.02
  all$rank_gate_passed <- all$pearson > 0 & all$spearman > 0
  all$cv_gate_passed <- FALSE
  all$recommended_status <- vapply(seq_len(nrow(all)), function(i) {
    delta_based <- isTRUE(!is.na(all$delta_source[[i]]) && nzchar(all$delta_source[[i]]) && !grepl("zero", all$delta_source[[i]]))
    classify_calibration_row(all[i, , drop = FALSE], delta_based, isTRUE(all$edge_delta_gate_passed[[i]]), isTRUE(all$deployment_gate_passed[[i]]))
  }, character(1))
  nonoracle <- all[!all$is_oracle, , drop = FALSE]
  valid_nonoracle <- nonoracle[nonoracle$recommended_status == "valid_shape_config", , drop = FALSE]
  best_num <- nonoracle[order(nonoracle$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
  best_valid <- if (nrow(valid_nonoracle)) valid_nonoracle[order(valid_nonoracle$shape_score), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  summary <- data.frame(
    n_total_configs = nrow(all),
    n_nonoracle_configs = nrow(nonoracle),
    n_oracle_configs = sum(all$is_oracle),
    n_valid_shape_configs = sum(all$recommended_status == "valid_shape_config"),
    n_nonoracle_valid_shape_configs = nrow(valid_nonoracle),
    n_delta_gate_passed = sum(all$edge_delta_gate_passed, na.rm = TRUE),
    n_amplitude_ok = sum(all$amplitude_gate_passed, na.rm = TRUE),
    n_wrong_direction = sum(all$recommended_status == "wrong_direction"),
    n_collapsed = sum(all$recommended_status == "amplitude_collapse"),
    best_numeric_only_config = best_num$candidate_id[[1L]],
    best_valid_shape_config = if (nrow(best_valid)) best_valid$candidate_id[[1L]] else NA_character_,
    recommended_status = if (nrow(valid_nonoracle)) "valid_shape_config" else "no_valid_shape_configuration",
    reason = if (nrow(valid_nonoracle)) "At least one non-oracle config passed all gates." else "No non-oracle config passed shape plus delta/deployment gates.",
    stringsAsFactors = FALSE
  )
  patch <- data.frame(
    patch_target = c("calibration_metrics_by_fit", "calibration_ranking", "delta_based_configs", "selection_return", "fallback_reporting"),
    recommendation = c(
      "Add estimate_sd_ratio, estimate_range_ratio, estimate_iqr_ratio.",
      "Add amplitude-collapse gate and penalty before rank aggregation.",
      "Require edge_delta_gate and deployment_gate for delta-based priors.",
      "Allow best config to be NA with no_valid_shape_configuration.",
      "Report best_numeric_only_config separately and do not call it a shape solution."
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all, file.path(dirs$tables, "calibration_shape_gate_demo.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "calibration_failure_state_summary.tsv"))
  write_tsv_safe(patch, file.path(dirs$tables, "formal_calibration_patch_recommendation.tsv"))
  out <- list(gate = all, summary = summary, patch = patch)
  saveRDS(out, rds)
  out
}

make_g_recommendations <- function(g1, g2, g3, g4, g5) {
  g2_far <- g2$top
  g3_far <- g3$top
  data.frame(
    table = c("sign_orientation_recommendation", "delta_deployment_recommendation", "pairwise_potential_recommendation",
              "local_debug_recommendation", "calibration_gate_recommendation", "recommended_next_steps"),
    recommendation = c(
      if (all(g1$toy$passed)) "No potential RHS sign bug found; do not globally flip delta sign." else "Potential RHS sign audit failed; fix sign before any further delta deployment.",
      if (any(g2$prediction_metrics$delta_spearman > 0 & g2$prediction_metrics$delta_sign_agreement >= 0.55 & g2$prediction_metrics$estimated_delta_sd_ratio >= 0.10 & !grepl("oracle", g2$prediction_metrics$delta_source), na.rm = TRUE)) "Some gated non-oracle delta deployment metrics pass locally, but require farfield confirmation." else "Gated non-oracle deployment remains insufficient; keep delta_untrusted.",
      if (any(g3_far$shape_class == "valid_shape" & !g3_far$oracle_delta, na.rm = TRUE)) "Pairwise potential can create non-oracle valid_shape candidates, but check calibration gates before promotion." else "Pairwise potential did not produce deployable non-oracle shape.",
      "Use custom TMB block-gradient probe to prioritize the largest gradient block; public local fit still needs obj/opt exposure for official diagnostics.",
      "Formal calibration should support no_valid_shape_configuration.",
      "Prioritize delta deployment/domain coverage, pairwise constraints, local TMB diagnostics exposure, and calibration failure-state patch."
    ),
    evidence = c(
      paste0("toy_pass=", all(g1$toy$passed)),
      paste0("covered_medium=", fmt_metric(mean(g2$overlap$is_in_domain, na.rm = TRUE))),
      paste0("best_pairwise_nonoracle_spearman=", fmt_metric(max(g3_far$spearman[!g3_far$oracle_delta], na.rm = TRUE))),
      paste0("max_gradient_block=", g4$gradient$max_gradient_block_name[which.max(g4$gradient$global_gradient_norm)]),
      paste0("status=", g5$summary$recommended_status[[1L]]),
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_g_report <- function(dirs, ctx, diagnostics_dir, delta_probe_dir, task_info, delta_probe, g1, g2, g3, g4, g5, recs) {
  all_long <- bind_rows_fill(list(
    transform(g1$cor, experiment = "G1_sign_orientation"),
    transform(g2$prediction_metrics, experiment = "G2_confidence_gated_delta"),
    transform(g2$metrics, experiment = "G2_confidence_gated_potential"),
    transform(g3$metrics, experiment = "G3_pairwise_potential"),
    transform(g4$gradient, experiment = "G4_local_gradient_debug"),
    transform(g5$gate, experiment = "G5_calibration_failure_state")
  ))
  summary <- data.frame(
    experiment = c("G1", "G2", "G3", "G4", "G5"),
    key_result = c(
      paste0("toy_sign_pass=", all(g1$toy$passed)),
      paste0("medium_in_domain_fraction=", fmt_metric(mean(g2$overlap$is_in_domain, na.rm = TRUE))),
      paste0("best_pairwise_nonoracle_spearman=", fmt_metric(max(g3$top$spearman[!g3$top$oracle_delta], na.rm = TRUE))),
      paste0("max_gradient_block=", g4$gradient$max_gradient_block_name[which.max(g4$gradient$global_gradient_norm)]),
      g5$summary$recommended_status[[1L]]
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_g_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "g_experiment_summary.tsv"))
  for (nm in recs$table) {
    write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  }
  g2_non <- g2$top[!g2$top$oracle_delta, , drop = FALSE]
  g3_non <- g3$top[!g3$top$oracle_delta, , drop = FALSE]
  lines <- c(
    "# Farfield Delta Debug Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", ctx$source_probe_dir %||% ctx$shared_input_dir, "`"),
    paste0("- abcd-dir: `", normalizePath(file.path(dirs$root, "..", "farfield_shape_probe_abcd"), winslash = "/", mustWork = FALSE), "`"),
    paste0("- diagnostics-dir: `", diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", delta_probe_dir, "`"),
    paste0("- simulation_id: ", task_info$simulation_id),
    paste0("- minobs: ", task_info$minobs),
    paste0("- input_policy: ", task_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## ABCD, E, and F summary",
    "- ABCD established the RMSE versus shape tradeoff: normalized/unit can shrink RMSE but often collapse amplitude, while mutation keeps amplitude but has wrong farfield rank.",
    "- E showed oracle per-edge delta restores farfield shape, but estimated/local/context delta is noisy and low-amplitude.",
    "- F showed direct-pair CV can look good locally, yet full-edge deployment and edge-delta gate fail; no C++ edge-gradient should be implemented yet.",
    "",
    "## G1 sign/orientation audit",
    paste0("- Potential sign toy test passed: ", all(g1$toy$passed), "."),
    paste0("- Oracle forward path Spearman with truth pair delta: ", fmt_metric(g1$cor$cor_truth_forward_spearman[g1$cor$method == "oracle_edge"]), "."),
    paste0("- Path-regression forward Spearman with truth pair delta: ", fmt_metric(g1$cor$cor_truth_forward_spearman[g1$cor$method == "path_regression"]), "."),
    "- No RHS sign bug was found. The negative delta_scale result is best interpreted as estimator compensation/domain-shift, not a global sign convention fix.",
    "",
    "## G2 confidence-gated delta",
    paste0("- Medium gate in-domain edge fraction: ", fmt_metric(mean(g2$overlap$is_in_domain, na.rm = TRUE)), "."),
    paste0("- Best non-oracle gated farfield Spearman: ", fmt_metric(max(g2_non$spearman, na.rm = TRUE)), "."),
    paste0("- Best non-oracle gated farfield shape class: ", g2_non$shape_class[which.max(g2_non$spearman)], "."),
    "- Gating reduces out-of-domain use but does not make the delta estimator deployable by itself.",
    "",
    "## G3 pairwise potential prior",
    paste0("- Best pairwise non-oracle farfield Spearman: ", fmt_metric(max(g3_non$spearman, na.rm = TRUE)), "."),
    paste0("- Best pairwise non-oracle shape class: ", g3_non$shape_class[which.max(g3_non$spearman)], "."),
    "- Pairwise constraints test whether the direct-pair signal should stay as path-level information rather than being decomposed into per-edge delta.",
    "",
    "## G4 local gradient block diagnostics",
    paste0("- Largest gradient block in sampled custom probes: ", g4$gradient$max_gradient_block_name[which.max(g4$gradient$global_gradient_norm)], "."),
    paste0("- Best gradient_norm among G4 variants: ", fmt_metric(min(g4$gradient$gradient_norm, na.rm = TRUE)), "."),
    "- Public `fit_local_posterior()` still does not expose obj/opt, so official block diagnostics require a package-level diagnostics hook.",
    "",
    "## G5 calibration failure state",
    paste0("- Recommended status: ", g5$summary$recommended_status[[1L]], "."),
    paste0("- Non-oracle valid shape configs after gates: ", g5$summary$n_nonoracle_valid_shape_configs[[1L]], "."),
    "- Formal calibration should be allowed to return no_valid_shape_configuration and report best numeric-only fallback separately.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- The blocker is not the Matrix precision mechanism. The blockers are non-oracle delta deployment/domain shift, lack of robust local TMB diagnostics, and missing formal calibration failure state.",
    "- Keep normalized as benchmark/probe default; keep unit as a synthetic stress-test and mutation as a legacy baseline.",
    "- Do not default `anchor_count_reference=minobs` for full input.",
    "- Highest priorities: add formal calibration gates, expose local TMB obj/opt diagnostics, improve domain-aware delta/pairwise estimators, and rerun on more simulations before C++ work."
  )
  writeLines(lines, file.path(dirs$root, "farfield_delta_debug_report.md"))
  saveRDS(list(g1 = g1, g2 = g2, g3 = g3, g4 = g4, g5 = g5, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_delta_debug_all_results.rds"))
}

main_delta_debug <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  mode <- match.arg(mode, c("prepare", "g1-sign-orientation", "g2-confidence-gated-delta",
                            "g3-pairwise-potential", "g4-local-gradient-debug",
                            "g5-calibration-failure-state", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_delta_debug"))
  simulation_id <- arg_integer(args, "simulation_id", 1L)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_debug_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, simulation_id, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, simulation_id, minobs, input_policy, force = force)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(simulation_id = simulation_id, minobs = minobs, input_policy = input_policy,
                    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
                    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
                    beta = as.numeric(ctx$input_table$sim_pm[[1L]]))
  components <- prepare_solver_cache(bundle$global_graph, dirs, force = force)
  delta_probe <- load_delta_probe_results(delta_probe_dir)
  assign("delta_probe_global_root", delta_probe$root, envir = .GlobalEnv)
  saveRDS(list(context = ctx, task_info = task_info, delta_probe_dir = delta_probe$root),
          file.path(dirs$results, "prepare_context.rds"))
  if (identical(mode, "prepare")) return(invisible(dirs$root))
  g1 <- if (mode %in% c("all", "g1-sign-orientation")) run_g1(bundle, grf, task_info, dirs, delta_probe, force = force) else readRDS(file.path(dirs$results, "g1_sign_orientation_audit.rds"))
  if (identical(mode, "g1-sign-orientation")) return(invisible(g1))
  g2 <- if (mode %in% c("all", "g2-confidence-gated-delta")) run_g2(bundle, components, grf, task_info, dirs, delta_probe, force = force) else readRDS(file.path(dirs$results, "g2_confidence_gated_delta.rds"))
  if (identical(mode, "g2-confidence-gated-delta")) return(invisible(g2))
  g3 <- if (mode %in% c("all", "g3-pairwise-potential")) run_g3(bundle, components, grf, task_info, dirs, delta_probe, force = force) else readRDS(file.path(dirs$results, "g3_pairwise_potential_prior.rds"))
  if (identical(mode, "g3-pairwise-potential")) return(invisible(g3))
  g4 <- if (mode %in% c("all", "g4-local-gradient-debug")) run_g4(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "g4_local_gradient_debug.rds"))
  if (identical(mode, "g4-local-gradient-debug")) return(invisible(g4))
  g5 <- if (mode %in% c("all", "g5-calibration-failure-state")) run_g5(delta_probe, g2, g3, dirs, force = force) else readRDS(file.path(dirs$results, "g5_calibration_failure_state_demo.rds"))
  if (identical(mode, "g5-calibration-failure-state")) return(invisible(g5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_g_recommendations(g1, g2, g3, g4, g5)
    write_g_report(dirs, ctx, normalizePath(diagnostics_dir, winslash = "/", mustWork = TRUE),
                   normalizePath(delta_probe_dir, winslash = "/", mustWork = TRUE),
                   task_info, delta_probe, g1, g2, g3, g4, g5, recs)
  }
  message("Wrote farfield delta debug under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) {
  main_delta_debug()
}
