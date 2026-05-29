#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_shell1_empirical_uncertainty.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_local_hessian_diagnostics.R"))

usage <- function() {
  cat(
    "Run shell1 empirical/bootstrap uncertainty diagnostics L1-L5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_shell1_empirical_uncertainty.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --delta-debug-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --local-calibration-dir=benchmark/results/farfield_local_calibration_patch \\\n",
    "    --core-fix-dir=benchmark/results/farfield_core_fix_probe \\\n",
    "    --identifiability-dir=benchmark/results/farfield_local_identifiability_repair \\\n",
    "    --hessian-dir=benchmark/results/farfield_local_hessian_diagnostics \\\n",
    "    --output-dir=benchmark/results/farfield_shell1_empirical_uncertainty \\\n",
    "    --simulation-ids=1,2,3,4,5,6,7,8,9,10 --minobs=5 --input-policy=full --bootstrap-reps=50\n",
    sep = ""
  )
}

make_l_dirs <- function(output_dir) make_probe_dirs(output_dir)

l_config_table <- function() {
  core <- k_config_table("core")
  k4 <- k_config_table("k4")
  out <- rbind(
    core[core$config_id == "C0_shell0_control", , drop = FALSE],
    core[core$config_id == "C4_J4_M4_g_fixed_scale_borrowed", , drop = FALSE],
    k4[k4$config_id == "E7_no_context_direct_borrowed_residual", names(core), drop = FALSE],
    core[core$config_id == "C1_shell1_baseline", , drop = FALSE],
    core[core$config_id == "C2_all_scale_fixed", , drop = FALSE]
  )
  out$l_label <- c("safe_shell0", "J4_M4", "E7_no_context", "shell1_baseline", "fixed_scale_best")
  out
}

bootstrap_counts_multinomial_l <- function(counts, seed) {
  set.seed(seed)
  out <- counts
  for (j in seq_len(ncol(counts))) {
    depth <- sum(counts[, j])
    prob <- counts[, j] / max(depth, 1)
    if (!is.finite(sum(prob)) || sum(prob) <= 0) prob[] <- 1 / length(prob)
    out[, j] <- as.vector(stats::rmultinom(1, size = depth, prob = prob))
  }
  rownames(out) <- rownames(counts)
  colnames(out) <- colnames(counts)
  out
}

local_fit_from_counts_l <- function(counts, dt, cfg, return_tmb_objects = FALSE, eval_override = NULL) {
  data <- alfak2::prepare_alfak2_data(counts, dt = dt)
  cfg2 <- cfg
  if (!is.null(eval_override)) cfg2$eval_max <- as.integer(eval_override)
  fit_local_k(data, cfg2, return_tmb_objects = return_tmb_objects)
}

truth_for_fit_l <- function(fit, grf, lambda) {
  tm <- compute_grf_truth(fit$summary$karyotype, grf$centroids, lambda)
  as.numeric(tm[as.character(fit$summary$karyotype)])
}

edge_table_l <- function(fit, grf = NULL, lambda = NA_real_, config_id = NA_character_) {
  graph <- fit$graph
  parent_idx <- as.integer(unlist(graph$parent_from0)) + 1L
  child_idx <- as.integer(unlist(graph$parent_to0)) + 1L
  if (!length(parent_idx)) {
    return(data.frame(
      config_id = character(), edge_id = integer(), parent = character(), child = character(),
      parent_node_id = integer(), child_node_id = integer(), context_label = character(),
      parent_support_tier = character(), child_support_tier = character(),
      parent_support_distance = integer(), child_support_distance = integer(),
      estimated_delta = numeric(), truth_delta = numeric(), stringsAsFactors = FALSE
    ))
  }
  parent <- as.character(graph$labels)[parent_idx]
  child <- as.character(graph$labels)[child_idx]
  m <- setNames(as.numeric(fit$summary$fitness_mean), as.character(fit$summary$karyotype))
  truth <- rep(NA_real_, length(graph$labels))
  if (!is.null(grf)) {
    tm <- compute_grf_truth(as.character(graph$labels), grf$centroids, lambda)
    truth <- as.numeric(tm[as.character(graph$labels)])
  }
  names(truth) <- as.character(graph$labels)
  context_id <- as.integer(unlist(graph$parent_context0)) + 1L
  data.frame(
    config_id = config_id,
    edge_id = seq_along(parent),
    parent = parent,
    child = child,
    parent_node_id = parent_idx,
    child_node_id = child_idx,
    context_label = as.character(graph$context_label)[context_id],
    parent_support_tier = as.character(graph$support_tier)[parent_idx],
    child_support_tier = as.character(graph$support_tier)[child_idx],
    parent_support_distance = as.integer(graph$support_distance)[parent_idx],
    child_support_distance = as.integer(graph$support_distance)[child_idx],
    estimated_delta = as.numeric(m[child] - m[parent]),
    truth_delta = as.numeric(truth[child] - truth[parent]),
    stringsAsFactors = FALSE
  )
}

rank_stability_l <- function(fit, original_map, config_id, bootstrap_id = NA_integer_) {
  est <- setNames(as.numeric(fit$summary$fitness_mean), as.character(fit$summary$karyotype))
  labels <- intersect(names(est), names(original_map))
  scopes <- list(
    direct = as.character(fit$summary$support_tier) == "directly_informed",
    local_borrowed = as.character(fit$summary$support_tier) == "local_borrowed",
    weakly_supported = as.integer(fit$summary$support_distance) >= 2L,
    all = rep(TRUE, nrow(fit$summary))
  )
  rows <- list()
  for (nm in names(scopes)) {
    lab <- as.character(fit$summary$karyotype)[scopes[[nm]]]
    lab <- intersect(lab, labels)
    rows[[nm]] <- data.frame(
      config_id = config_id,
      bootstrap_id = bootstrap_id,
      support_scope = nm,
      n_nodes = length(lab),
      bootstrap_pairwise_correlation = safe_cor2(est[lab], original_map[lab], "pearson"),
      fitness_rank_stability = safe_cor2(est[lab], original_map[lab], "spearman"),
      stringsAsFactors = FALSE
    )
  }
  bind_rows_fill(rows)
}

summarize_node_boot_l <- function(fits, original_fit, grf, lambda, config_id) {
  labels <- as.character(original_fit$summary$karyotype)
  mat <- do.call(rbind, lapply(fits, function(f) setNames(as.numeric(f$summary$fitness_mean), as.character(f$summary$karyotype))[labels]))
  truth <- truth_for_fit_l(original_fit, grf, lambda)
  data.frame(
    config_id = config_id,
    node_id = as.integer(original_fit$summary$node_id),
    karyotype = labels,
    support_tier = original_fit$summary$support_tier,
    support_distance = original_fit$summary$support_distance,
    truth_fitness = unname(truth),
    original_map_fitness = unname(original_fit$summary$fitness_mean),
    bootstrap_mean = unname(colMeans(mat, na.rm = TRUE)),
    bootstrap_sd = unname(apply(mat, 2, stats::sd, na.rm = TRUE)),
    bootstrap_iqr = unname(apply(mat, 2, stats::IQR, na.rm = TRUE)),
    bootstrap_q05 = unname(apply(mat, 2, stats::quantile, probs = 0.05, na.rm = TRUE)),
    bootstrap_q50 = unname(apply(mat, 2, stats::quantile, probs = 0.50, na.rm = TRUE)),
    bootstrap_q95 = unname(apply(mat, 2, stats::quantile, probs = 0.95, na.rm = TRUE)),
    tmb_sd = unname(original_fit$summary$fitness_sd),
    count_total = unname(original_fit$summary$count_total),
    effective_count_total = unname(original_fit$summary$effective_count_total),
    stringsAsFactors = FALSE
  )
}

summarize_edge_boot_l <- function(fits, original_fit, grf, lambda, config_id) {
  orig <- edge_table_l(original_fit, grf, lambda, config_id)
  if (!nrow(orig)) {
    return(data.frame(
      config_id = character(), edge_id = integer(), parent = character(), child = character(),
      context_label = character(), parent_support_tier = character(), child_support_tier = character(),
      truth_delta = numeric(), original_map_delta = numeric(), bootstrap_delta_mean = numeric(),
      bootstrap_delta_sd = numeric(), bootstrap_sign_positive_rate = numeric(),
      bootstrap_sign_matches_truth_rate = numeric(), edge_delta_stable = logical(),
      stringsAsFactors = FALSE
    ))
  }
  key <- paste(orig$parent, orig$child, sep = "->")
  mat <- do.call(rbind, lapply(fits, function(f) {
    ed <- edge_table_l(f, grf, lambda, config_id)
    setNames(ed$estimated_delta, paste(ed$parent, ed$child, sep = "->"))[key]
  }))
  data.frame(
    config_id = config_id,
    edge_id = orig$edge_id,
    parent = orig$parent,
    child = orig$child,
    context_label = orig$context_label,
    parent_support_tier = orig$parent_support_tier,
    child_support_tier = orig$child_support_tier,
    truth_delta = orig$truth_delta,
    original_map_delta = orig$estimated_delta,
    bootstrap_delta_mean = colMeans(mat, na.rm = TRUE),
    bootstrap_delta_sd = apply(mat, 2, stats::sd, na.rm = TRUE),
    bootstrap_sign_positive_rate = colMeans(sign(mat) > 0, na.rm = TRUE),
    bootstrap_sign_matches_truth_rate = colMeans(sign(mat) == rep(sign(orig$truth_delta), each = nrow(mat)), na.rm = TRUE),
    edge_delta_stable = colMeans(sign(mat) == rep(sign(orig$truth_delta), each = nrow(mat)), na.rm = TRUE) >= 0.55,
    stringsAsFactors = FALSE
  )
}

run_l1 <- function(bundle, counts, dt, grf, task_info, dirs, bootstrap_reps, force = FALSE) {
  rds <- file.path(dirs$results, "l1_bootstrap_empirical_uncertainty.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- l_config_table()
  actual_reps <- min(as.integer(bootstrap_reps), 8L)
  rows <- list(); node_rows <- list(); edge_rows <- list(); rank_rows <- list()
  success_fits <- list(); original_fits <- list()
  idx <- 0L
  for (i in seq_len(nrow(cfg))) {
    cid <- cfg$config_id[[i]]
    orig <- local_fit_from_counts_l(counts, dt, cfg[i, , drop = FALSE], return_tmb_objects = FALSE)
    if (!inherits(orig$fit, "error")) original_fits[[cid]] <- orig$fit
    original_map <- if (!inherits(orig$fit, "error")) setNames(orig$fit$summary$fitness_mean, orig$fit$summary$karyotype) else numeric()
    fits <- list()
    for (b in seq_len(actual_reps)) {
      boot_counts <- bootstrap_counts_multinomial_l(counts, seed = 10000 + i * 1000 + b)
      res <- local_fit_from_counts_l(boot_counts, dt, cfg[i, , drop = FALSE], return_tmb_objects = FALSE, eval_override = min(500L, as.integer(cfg$eval_max[[i]])))
      idx <- idx + 1L
      rows[[idx]] <- local_diag_row_k(res, cfg[i, , drop = FALSE], grf, task_info$lambda,
                                      extra = list(bootstrap_id = b, bootstrap_mode = "multinomial_timepoint",
                                                   requested_bootstrap_reps = bootstrap_reps, actual_bootstrap_reps = actual_reps))
      if (!inherits(res$fit, "error")) {
        fits[[length(fits) + 1L]] <- res$fit
        rank_rows[[length(rank_rows) + 1L]] <- rank_stability_l(res$fit, original_map, cid, b)
      }
    }
    success_fits[[cid]] <- fits
    if (length(fits) && !inherits(orig$fit, "error")) {
      node_rows[[length(node_rows) + 1L]] <- summarize_node_boot_l(fits, orig$fit, grf, task_info$lambda, cid)
      edge_rows[[length(edge_rows) + 1L]] <- summarize_edge_boot_l(fits, orig$fit, grf, task_info$lambda, cid)
    }
  }
  run_tbl <- bind_rows_fill(rows)
  node_tbl <- bind_rows_fill(node_rows)
  edge_tbl <- bind_rows_fill(edge_rows)
  rank_tbl <- bind_rows_fill(rank_rows)
  sign_tbl <- stats::aggregate(cbind(bootstrap_sign_matches_truth_rate, edge_delta_stable) ~ config_id, edge_tbl, mean, na.rm = TRUE)
  names(sign_tbl)[names(sign_tbl) == "bootstrap_sign_matches_truth_rate"] <- "edge_delta_sign_stability"
  names(sign_tbl)[names(sign_tbl) == "edge_delta_stable"] <- "bootstrap_edge_gate_pass_rate"
  success_counts <- stats::aggregate(status ~ config_id, run_tbl, function(x) sum(x == "ok", na.rm = TRUE))
  names(success_counts)[names(success_counts) == "status"] <- "bootstrap_success_reps"
  summary <- data.frame(config_id = unique(run_tbl$config_id), stringsAsFactors = FALSE)
  summary$bootstrap_reps_requested <- bootstrap_reps
  summary$bootstrap_reps <- actual_reps
  summary <- merge(summary, success_counts, by = "config_id", all.x = TRUE)
  summary$bootstrap_success_reps[is.na(summary$bootstrap_success_reps)] <- 0L
  summary$bootstrap_failure_reps <- actual_reps - summary$bootstrap_success_reps
  summary <- merge(summary, sign_tbl, by = "config_id", all.x = TRUE)
  write_tsv_safe(run_tbl, file.path(dirs$tables, "l1_local_bootstrap_runs.tsv"))
  write_tsv_safe(node_tbl, file.path(dirs$tables, "l1_local_bootstrap_node_variance.tsv"))
  write_tsv_safe(edge_tbl, file.path(dirs$tables, "l1_local_bootstrap_edge_delta_variance.tsv"))
  write_tsv_safe(sign_tbl, file.path(dirs$tables, "l1_local_bootstrap_sign_stability.tsv"))
  write_tsv_safe(rank_tbl, file.path(dirs$tables, "l1_local_bootstrap_rank_stability.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "l1_bootstrap_empirical_uncertainty_summary.tsv"))
  out <- list(runs = run_tbl, nodes = node_tbl, edges = edge_tbl, sign = sign_tbl, rank = rank_tbl,
              summary = summary, success_fits = success_fits, original_fits = original_fits,
              actual_reps = actual_reps, requested_reps = bootstrap_reps)
  saveRDS(out, rds)
  out
}

regularized_cov_diag_l <- function(fit, method, value) {
  H <- fit$optimizer$obj$he(fit$optimizer$opt$par)
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  vals <- ev$values
  if (method == "ridge") vals2 <- vals + value
  else if (method == "eigen_floor") vals2 <- pmax(vals, value)
  else vals2 <- ifelse(abs(vals) < value, Inf, vals)
  inv_vals <- ifelse(is.finite(vals2) & vals2 > 0, 1 / vals2, 0)
  cov <- ev$vectors %*% (inv_vals * t(ev$vectors))
  pmap <- parameter_map_k(fit)
  frows <- pmap[pmap$parameter_block == "f" & !is.na(pmap$node_id), , drop = FALSE]
  sd <- sqrt(pmax(diag(cov)[frows$par_index], 0))
  if (identical(fit$diagnostics$local_parameterization, "g_equivalent")) sd <- sd / fit$data$dt
  data.frame(node_id = frows$node_id, karyotype = frows$karyotype,
             support_tier = frows$support_tier, support_distance = frows$support_distance,
             support_scope = frows$support_scope, nodewise_regularized_hessian_sd = sd,
             stringsAsFactors = FALSE)
}

run_l2 <- function(l1, counts, dt, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "l2_regularized_hessian_covariance.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- l_config_table()
  cfg <- cfg[cfg$config_id %in% c("C0_shell0_control", "C4_J4_M4_g_fixed_scale_borrowed", "E7_no_context_direct_borrowed_residual", "C1_shell1_baseline"), , drop = FALSE]
  methods <- expand.grid(method = c("ridge", "eigen_floor", "pseudo_inverse"),
                         value = c(1e-8, 1e-6, 1e-4), stringsAsFactors = FALSE)
  rows <- list(); idx <- 0L
  for (i in seq_len(nrow(cfg))) {
    res <- local_fit_from_counts_l(counts, dt, cfg[i, , drop = FALSE], return_tmb_objects = TRUE)
    if (inherits(res$fit, "error")) next
    for (j in seq_len(nrow(methods))) {
      rg <- tryCatch(regularized_cov_diag_l(res$fit, methods$method[[j]], methods$value[[j]]), error = function(e) e)
      if (inherits(rg, "error")) next
      rg$config_id <- cfg$config_id[[i]]
      rg$regularization_method <- methods$method[[j]]
      rg$regularization_value <- methods$value[[j]]
      rows[[idx <- idx + 1L]] <- rg
    }
  }
  reg <- bind_rows_fill(rows)
  boot <- l1$nodes[, c("config_id", "karyotype", "bootstrap_sd"), drop = FALSE]
  cmp <- merge(reg, boot, by = c("config_id", "karyotype"), all.x = TRUE)
  cmp$sd_ratio_regularized_to_bootstrap <- cmp$nodewise_regularized_hessian_sd / cmp$bootstrap_sd
  scope <- stats::aggregate(cbind(nodewise_regularized_hessian_sd, bootstrap_sd, sd_ratio_regularized_to_bootstrap) ~
                              config_id + regularization_method + regularization_value + support_scope,
                            cmp, median, na.rm = TRUE)
  rec <- stats::aggregate(abs(log(pmax(sd_ratio_regularized_to_bootstrap, 1e-8))) ~
                            config_id + regularization_method + regularization_value,
                          cmp, median, na.rm = TRUE)
  names(rec)[names(rec) == "abs(log(pmax(sd_ratio_regularized_to_bootstrap, 1e-08)))"] <- "median_abs_log_sd_ratio"
  rec <- rec[order(rec$median_abs_log_sd_ratio), , drop = FALSE]
  rec$recommendation <- ifelse(rec$median_abs_log_sd_ratio < 1, "regularized_hessian_same_order", "bootstrap_preferred")
  write_tsv_safe(reg, file.path(dirs$tables, "l2_regularized_hessian_covariance.tsv"))
  write_tsv_safe(cmp, file.path(dirs$tables, "l2_regularized_vs_bootstrap_sd.tsv"))
  write_tsv_safe(scope, file.path(dirs$tables, "l2_hessian_covariance_by_scope.tsv"))
  write_tsv_safe(head(rec, 20L), file.path(dirs$tables, "l2_regularized_hessian_recommendation.tsv"))
  out <- list(regularized = reg, comparison = cmp, by_scope = scope, recommendation = rec)
  saveRDS(out, rds)
  out
}

potential_mean_l <- function(local_fit, graph, edge_tbl, delta_scale = 1, lambda_anchor = 10, ridge = 1e-6) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix required", call. = FALSE)
  n <- length(graph$labels)
  key <- paste(edge_tbl$parent, edge_tbl$child, sep = "->")
  gparent <- as.character(graph$labels)[as.integer(unlist(graph$parent_from0)) + 1L]
  gchild <- as.character(graph$labels)[as.integer(unlist(graph$parent_to0)) + 1L]
  gkey <- paste(gparent, gchild, sep = "->")
  m <- match(gkey, key)
  keep <- !is.na(m) & is.finite(edge_tbl$original_map_delta[m])
  from <- as.integer(unlist(graph$parent_from0))[keep] + 1L
  to <- as.integer(unlist(graph$parent_to0))[keep] + 1L
  delta <- delta_scale * edge_tbl$original_map_delta[m[keep]]
  conf <- pmin(1, pmax(0, edge_tbl$bootstrap_sign_matches_truth_rate[m[keep]]))
  q <- Matrix::Diagonal(n, ridge)
  rhs <- numeric(n)
  if (length(from)) {
    sw <- sqrt(conf)
    rr <- seq_along(from)
    A <- Matrix::sparseMatrix(i = c(rr, rr), j = c(from, to), x = c(-sw, sw), dims = c(length(from), n))
    q <- q + Matrix::crossprod(A)
    rhs <- rhs + as.numeric(Matrix::crossprod(A, sw * delta))
  }
  direct <- as.character(local_fit$summary$support_tier) == "directly_informed"
  ai <- match(as.character(local_fit$summary$karyotype[direct]), as.character(graph$labels))
  ok <- !is.na(ai)
  ai <- ai[ok]
  am <- local_fit$summary$fitness_mean[direct][ok]
  q <- q + Matrix::sparseMatrix(i = ai, j = ai, x = lambda_anchor, dims = c(n, n))
  rhs[ai] <- rhs[ai] + lambda_anchor * am
  as.numeric(Matrix::solve(Matrix::Cholesky(Matrix::forceSymmetric(q), LDL = TRUE, perm = TRUE), rhs))
}

fit_global_with_prior_l <- function(local_fit, graph, prior_mean, cfg) {
  local_resid <- local_fit
  idx <- match(as.character(local_resid$summary$karyotype), as.character(graph$labels))
  local_resid$summary$fitness_mean <- local_resid$summary$fitness_mean - prior_mean[idx]
  gp <- alfak2::fit_graph_posterior(local_resid, graph,
                                    lambda_l_grid = cfg$lambda_l,
                                    lambda_e_grid = cfg$lambda_e,
                                    sigma_obs_grid = cfg$sigma_obs,
                                    graph_edge_weight = cfg$graph_edge_weight,
                                    compute_sd = FALSE)
  gp$summary$fitness_mean <- gp$summary$fitness_mean + prior_mean
  gp$summary$conf_low <- NA_real_
  gp$summary$conf_high <- NA_real_
  gp
}

run_l3 <- function(l1, counts, dt, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "l3_deterministic_borrowed_followup.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  e7_fit <- l1$original_fits[["E7_no_context_direct_borrowed_residual"]]
  m4_fit <- l1$original_fits[["C4_J4_M4_g_fixed_scale_borrowed"]]
  edge <- l1$edges[l1$edges$config_id == "E7_no_context_direct_borrowed_residual", , drop = FALSE]
  deployment <- data.frame(
    config_id = "E7_no_context_direct_borrowed_residual",
    deployment_sign_agreement = mean(sign(edge$original_map_delta) == sign(edge$truth_delta), na.rm = TRUE),
    deployment_spearman = safe_cor2(edge$original_map_delta, edge$truth_delta, "spearman"),
    estimated_delta_sd_ratio = stats::sd(edge$original_map_delta, na.rm = TRUE) / stats::sd(edge$truth_delta, na.rm = TRUE),
    bootstrap_stable_edge_fraction = mean(edge$edge_delta_stable, na.rm = TRUE),
    gate_passed = FALSE,
    stringsAsFactors = FALSE
  )
  deployment$gate_passed <- deployment$deployment_sign_agreement >= 0.55 &&
    deployment$deployment_spearman > 0 && deployment$estimated_delta_sd_ratio >= 0.10 &&
    deployment$bootstrap_stable_edge_fraction >= 0.55
  global_cfg <- data.frame(
    experiment = "L3",
    candidate_id = c("mutation_baseline", "normalized_default", "unit_stress"),
    graph_edge_weight = c("mutation", "normalized", "unit"),
    lambda_l = c(0.2, 0.2, 0.2),
    lambda_e = c(1, 0.01, 0.01),
    sigma_obs = c(0.05, 0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "E7_potential",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    stringsAsFactors = FALSE
  )
  scales <- c(0.25, 0.5, 1.0, -0.5)
  rows <- list(); idx <- 0L
  if (!is.null(e7_fit)) {
    graph <- alfak2::build_karyotype_graph(e7_fit$data, shell_depth = 2, max_nodes = 30000)
    for (scale in scales) {
      pm <- potential_mean_l(e7_fit, graph, edge, delta_scale = scale)
      for (j in seq_len(nrow(global_cfg))) {
        cfg <- global_cfg[j, , drop = FALSE]
        cfg$prior_mean_scale <- scale
        fit <- tryCatch(fit_global_with_prior_l(e7_fit, graph, pm, cfg), error = function(e) e)
        idx <- idx + 1L
        if (inherits(fit, "error")) {
          rows[[idx]] <- data.frame(prior_mean_scale = scale, candidate_id = cfg$candidate_id, status = "error", error_message = conditionMessage(fit))
        } else {
          m <- score_summary_abcd(fit$summary, graph, grf, task_info$lambda, task_info, cfg, "E7_potential")
          far <- m[m$support_scope == "farfield" & m$metric_scale == "native", , drop = FALSE]
          far$prior_mean_scale <- scale
          far$shape_classification <- metric_shape_class(far)
          far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_i(far[k, , drop = FALSE]), character(1))
          rows[[idx]] <- far
        }
      }
    }
  }
  pot <- bind_rows_fill(rows)
  stability <- l1$sign[l1$sign$config_id %in% c("E7_no_context_direct_borrowed_residual", "C4_J4_M4_g_fixed_scale_borrowed"), , drop = FALSE]
  follow <- rbind(
    data.frame(config_id = "E7_no_context_direct_borrowed_residual", comparison = "E7", stringsAsFactors = FALSE),
    data.frame(config_id = "C4_J4_M4_g_fixed_scale_borrowed", comparison = "J4_M4", stringsAsFactors = FALSE)
  )
  rec <- data.frame(
    recommendation = if (isTRUE(deployment$gate_passed) && any(pot$recommended_status == "valid_shape_config", na.rm = TRUE)) "E7_experimental_delta_candidate" else "E7_delta_not_deployable",
    deployment_gate_passed = deployment$gate_passed,
    any_valid_shape = any(pot$recommended_status == "valid_shape_config", na.rm = TRUE),
    negative_scale_best = any(pot$prior_mean_scale < 0 & pot$recommended_status == "valid_shape_config", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(follow, file.path(dirs$tables, "l3_deterministic_borrowed_delta_followup.tsv"))
  write_tsv_safe(stability, file.path(dirs$tables, "l3_deterministic_borrowed_bootstrap_stability.tsv"))
  write_tsv_safe(deployment, file.path(dirs$tables, "l3_deterministic_borrowed_deployment_gate.tsv"))
  write_tsv_safe(pot, file.path(dirs$tables, "l3_deterministic_borrowed_potential_probe.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "l3_deterministic_borrowed_recommendation.tsv"))
  out <- list(followup = follow, stability = stability, deployment = deployment, potential = pot, recommendation = rec)
  saveRDS(out, rds)
  out
}

run_l4 <- function(source_input_dir, simulation_ids, minobs, input_policy, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "l4_production_safe_calibration_dryrun.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  requested_sim_ids <- simulation_ids
  sim_ids <- requested_sim_ids[seq_len(min(length(requested_sim_ids), 5L))]
  global_sim_ids <- sim_ids
  cfgs <- data.frame(
    experiment = "L4",
    candidate_id = c("mutation_baseline", "normalized_default", "unit_stress"),
    graph_edge_weight = c("mutation", "normalized", "unit"),
    lambda_l = c(0.2, 0.2, 0.2),
    lambda_e = c(1, 0.01, 0.01),
    sigma_obs = c(0.05, 0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    stringsAsFactors = FALSE
  )
  rows <- list(); ridx <- 0L
  status_rows <- list()
  for (sim in sim_ids) {
    row <- resolve_shared_input_row(source_input_dir, sim, minobs)
    if (!nrow(row)) next
    yi <- readRDS(row$input_rds[[1L]])
    counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = input_policy, drop_diploid = TRUE)
    dt <- suppressWarnings(diff(as.numeric(colnames(counts))))
    if (length(dt) != 1L || !is.finite(dt) || dt <= 0) dt <- as.numeric(row$time_delta[[1L]])
    grf <- readRDS(row$grf_rds[[1L]])
    local_modes <- list(
      safe_shell0 = k_config_table("core")[k_config_table("core")$config_id == "C0_shell0_control", , drop = FALSE],
      experimental_shell1_map = k_config_table("core")[k_config_table("core")$config_id == "C4_J4_M4_g_fixed_scale_borrowed", , drop = FALSE]
    )
    sim_far <- list()
    for (lm in names(local_modes)) {
      res <- local_fit_from_counts_l(counts, dt, local_modes[[lm]], return_tmb_objects = FALSE, eval_override = 500)
      if (inherits(res$fit, "error") || !(sim %in% global_sim_ids)) next
      graph <- alfak2::build_karyotype_graph(res$fit$data, shell_depth = 2, max_nodes = 30000)
      for (j in seq_len(nrow(cfgs))) {
        started <- Sys.time()
        fit <- tryCatch(alfak2::fit_graph_posterior(res$fit, graph,
                                                    lambda_l_grid = cfgs$lambda_l[[j]],
                                                    lambda_e_grid = cfgs$lambda_e[[j]],
                                                    sigma_obs_grid = cfgs$sigma_obs[[j]],
                                                    graph_edge_weight = cfgs$graph_edge_weight[[j]],
                                                    compute_sd = FALSE), error = function(e) e)
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        ridx <- ridx + 1L
        if (inherits(fit, "error")) next
        task <- list(simulation_id = sim, minobs = minobs, input_policy = input_policy,
                     lambda = as.numeric(row$lambda[[1L]]), dt = dt,
                     beta = if ("sim_beta" %in% names(row)) as.numeric(row$sim_beta[[1L]]) else 0.00005)
        m <- score_summary_abcd(fit$summary, graph, grf, task$lambda, task, cfgs[j, , drop = FALSE], "production_safe_dryrun")
        far <- m[m$support_scope == "farfield" & m$metric_scale == "native", , drop = FALSE]
        far$simulation_id <- sim
        far$local_mode <- lm
        far$local_covariance_status <- res$fit$diagnostics$covariance_status
        far$runtime_sec <- elapsed
        far$shape_classification <- metric_shape_class(far)
        far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_i(far[k, , drop = FALSE]), character(1))
        far$failure_reason <- ifelse(far$recommended_status == "valid_shape_config", "passed", far$recommended_status)
        rows[[ridx]] <- far
        sim_far[[length(sim_far) + 1L]] <- far
      }
    }
    sf <- bind_rows_fill(sim_far)
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      simulation_id = sim,
      n_valid_shape_configs = sum(sf$recommended_status == "valid_shape_config", na.rm = TRUE),
      recommended_status = if (any(sf$recommended_status == "valid_shape_config", na.rm = TRUE)) "valid_shape_config" else "no_valid_shape_configuration",
      stringsAsFactors = FALSE
    )
  }
  tbl <- bind_rows_fill(rows)
  stat <- bind_rows_fill(status_rows)
  numeric <- tbl[order(tbl$centered_rmse), , drop = FALSE]
  numeric <- numeric[!duplicated(paste(numeric$simulation_id, numeric$local_mode)), , drop = FALSE]
  summary <- data.frame(
    normalized_collapse_rate = mean(tbl$amplitude_collapse[tbl$candidate_id == "normalized_default"], na.rm = TRUE),
    mutation_wrong_direction_rate = mean(tbl$recommended_status[tbl$candidate_id == "mutation_baseline"] == "wrong_direction", na.rm = TRUE),
    unit_collapse_rate = mean(tbl$amplitude_collapse[tbl$candidate_id == "unit_stress"], na.rm = TRUE),
    n_valid_shape_configs = sum(tbl$recommended_status == "valid_shape_config", na.rm = TRUE),
    n_no_valid_shape_configuration = sum(stat$recommended_status == "no_valid_shape_configuration", na.rm = TRUE),
    dryrun_note = "Runtime-bounded dry run on sim1-5 for both local and global scoring.",
    stringsAsFactors = FALSE
  )
  write_tsv_safe(tbl, file.path(dirs$tables, "l4_production_safe_calibration_dryrun.tsv"))
  write_tsv_safe(stat, file.path(dirs$tables, "l4_recommended_status_by_sim.tsv"))
  write_tsv_safe(numeric, file.path(dirs$tables, "l4_best_numeric_only_by_sim.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "l4_production_safe_summary.tsv"))
  out <- list(results = tbl, status = stat, numeric = numeric, summary = summary)
  saveRDS(out, rds)
  out
}

run_l5 <- function(l1, l2, l3, l4, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "l5_shape_uncertainty_separation.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  shape <- data.frame(
    protocol = "shape_ranking",
    compute_sd = FALSE,
    inputs = "all calibration candidates",
    outputs = "fitness_mean, centered_rmse, pearson, spearman, amplitude diagnostics, no_valid_shape_configuration",
    gate = "rank/amplitude/failure-state",
    stringsAsFactors = FALSE
  )
  uncert <- data.frame(
    protocol = "uncertainty_calibration",
    compute_sd = "selected_only",
    inputs = "safe_shell0 normalized, mutation baseline, unit stress-test, valid_shape if any, experimental shell1 bootstrap",
    outputs = "TMB sd for shell0 if finite, bootstrap/restart variance for shell1, coverage diagnostics",
    gate = "do not use nonfinite TMB covariance",
    stringsAsFactors = FALSE
  )
  selected <- data.frame(
    config_id = c("safe_shell0_normalized", "mutation_baseline", "unit_stress", "J4_M4_experimental_bootstrap", "E7_experimental_bootstrap"),
    reason = c("production default candidate", "legacy amplitude baseline", "synthetic stress-test", "stable MAP but TMB covariance nonfinite", "deterministic borrowed diagnostic"),
    uncertainty_method = c("TMB_if_finite_or_bootstrap", "TMB_if_finite", "TMB_if_finite", "bootstrap", "bootstrap"),
    stringsAsFactors = FALSE
  )
  runtime <- data.frame(
    stage = c("shape_ranking", "bootstrap_uncertainty", "regularized_hessian", "production_dryrun"),
    compute_sd = c(FALSE, FALSE, FALSE, FALSE),
    observed_or_expected_runtime = c("fast mean-only global solves", "bounded bootstrap local refits", "single-fit Hessian eigensolves", "sim1-5 mean-only scoring"),
    stringsAsFactors = FALSE
  )
  rec <- data.frame(
    recommendation = c("Use compute_sd=FALSE as default shape ranking mode.",
                       "Run uncertainty calibration only on selected configs.",
                       "Use bootstrap/restart variance for shell1; do not use nonfinite TMB sdreport.",
                       "Keep uncertainty failure separate from shape ranking failure.",
                       "Avoid all-config posterior SD loops in benchmark grids."),
    evidence = c(
      paste0("L4 valid shape count=", sum(l4$results$recommended_status == "valid_shape_config", na.rm = TRUE)),
      paste0("L1 actual bootstrap reps=", l1$actual_reps),
      paste0("L2 best method=", l2$recommendation$regularization_method[[1L]]),
      "K2 shell_depth=1 direct-only covariance was nonfinite.",
      "compute_sd=FALSE already validated as mean-equivalent."
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(shape, file.path(dirs$tables, "l5_shape_ranking_protocol.tsv"))
  write_tsv_safe(uncert, file.path(dirs$tables, "l5_uncertainty_calibration_protocol.tsv"))
  write_tsv_safe(selected, file.path(dirs$tables, "l5_selected_uncertainty_configs.tsv"))
  write_tsv_safe(runtime, file.path(dirs$tables, "l5_shape_vs_uncertainty_runtime.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "l5_shape_vs_uncertainty_recommendation.tsv"))
  out <- list(shape = shape, uncertainty = uncert, selected = selected, runtime = runtime, recommendation = rec)
  saveRDS(out, rds)
  out
}

make_l_recommendations <- function(l1, l2, l3, l4, l5) {
  data.frame(
    table = c("bootstrap_empirical_uncertainty_recommendation", "regularized_hessian_recommendation",
              "deterministic_borrowed_followup_recommendation", "production_safe_calibration_recommendation",
              "shape_uncertainty_separation_recommendation", "recommended_next_steps"),
    recommendation = c(
      "Shell1 MAP may be used only as an experimental mean provider; uncertainty must be empirical/bootstrap.",
      if (nrow(l2$recommendation) && l2$recommendation$recommendation[[1L]] == "regularized_hessian_same_order") "Regularized Hessian is a possible fast approximation but needs validation." else "Bootstrap variance is preferred over regularized Hessian.",
      if (isTRUE(l3$deployment$gate_passed)) "E7 passed deployment gate experimentally." else "E7/deterministic borrowed did not pass deployable delta gates.",
      "Production dry run should return no_valid_shape_configuration rather than force a numeric-only winner.",
      "Separate shape ranking from uncertainty calibration; use compute_sd=FALSE for ranking.",
      "Keep edge-gradient disabled; prioritize shell0 production fallback, bootstrap shell1 diagnostics, and delta deployment gates."
    ),
    evidence = c(
      paste0("actual_bootstrap_reps=", l1$actual_reps),
      paste0("best_regularized=", l2$recommendation$regularization_method[[1L]], "/", l2$recommendation$regularization_value[[1L]]),
      paste0("E7_gate=", l3$deployment$gate_passed, "; valid_shape=", any(l3$potential$recommended_status == "valid_shape_config", na.rm = TRUE)),
      paste0("L4 no_valid_count=", l4$summary$n_no_valid_shape_configuration[[1L]]),
      "L5 protocol tables written.",
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_l_report <- function(dirs, args_info, ctx, l1, l2, l3, l4, l5, recs) {
  all_long <- bind_rows_fill(list(
    transform(l1$runs, experiment = "L1_bootstrap_runs"),
    transform(l1$summary, experiment = "L1_bootstrap_summary"),
    transform(l2$recommendation, experiment = "L2_regularized_hessian"),
    transform(l3$deployment, experiment = "L3_E7_deployment"),
    transform(l3$potential, experiment = "L3_E7_potential"),
    transform(l4$results, experiment = "L4_production_dryrun"),
    transform(l5$recommendation, experiment = "L5_protocol")
  ))
  summary <- data.frame(
    experiment = c("L1", "L2", "L3", "L4", "L5"),
    key_result = c(
      paste0("actual_bootstrap_reps=", l1$actual_reps),
      paste0("best_regularized=", l2$recommendation$regularization_method[[1L]], "/", l2$recommendation$regularization_value[[1L]]),
      paste0("E7_gate=", l3$deployment$gate_passed[[1L]], "; valid_shape=", any(l3$potential$recommended_status == "valid_shape_config", na.rm = TRUE)),
      paste0("valid_shape=", sum(l4$results$recommended_status == "valid_shape_config", na.rm = TRUE)),
      "shape_ranking_compute_sd_false"
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_l_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "l_experiment_summary.tsv"))
  for (nm in recs$table) write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  m4sign <- l1$sign[l1$sign$config_id == "C4_J4_M4_g_fixed_scale_borrowed", , drop = FALSE]
  e7sign <- l1$sign[l1$sign$config_id == "E7_no_context_direct_borrowed_residual", , drop = FALSE]
  lines <- c(
    "# Farfield Shell1 Empirical Uncertainty Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", args_info$source_input_dir, "`"),
    paste0("- abcd-dir: `", args_info$abcd_dir, "`"),
    paste0("- diagnostics-dir: `", args_info$diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", args_info$delta_probe_dir, "`"),
    paste0("- delta-debug-dir: `", args_info$delta_debug_dir, "`"),
    paste0("- local-calibration-dir: `", args_info$local_calibration_dir, "`"),
    paste0("- core-fix-dir: `", args_info$core_fix_dir, "`"),
    paste0("- identifiability-dir: `", args_info$identifiability_dir, "`"),
    paste0("- hessian-dir: `", args_info$hessian_dir, "`"),
    paste0("- simulation_ids: ", paste(args_info$simulation_ids, collapse = ",")),
    paste0("- minobs: ", args_info$minobs),
    paste0("- input_policy: ", args_info$input_policy),
    paste0("- bootstrap_reps requested: ", args_info$bootstrap_reps, "; actual runtime-bounded reps: ", l1$actual_reps),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## Prior results summary",
    "- ABCD through K show no deployable non-oracle farfield shape. J4/M4 MAP is restart-stable but shell_depth=1 TMB covariance is nonfinite/near-singular.",
    "- K2 showed shell_depth=1 direct-only sdreport is also nonfinite, so shell1 uncertainty must be empirical until the local model is reparameterized.",
    "",
    "## L1 bootstrap empirical uncertainty",
    paste0("- Actual bootstrap reps per config: ", l1$actual_reps, " (requested ", l1$requested_reps, ")."),
    paste0("- J4/M4 bootstrap edge sign stability: ", if (nrow(m4sign)) fmt_metric(m4sign$edge_delta_sign_stability[[1L]]) else "NA", "."),
    paste0("- E7 bootstrap edge sign stability: ", if (nrow(e7sign)) fmt_metric(e7sign$edge_delta_sign_stability[[1L]]) else "NA", "."),
    "- Shell1 MAP should remain experimental: use it only as a mean-provider candidate and use bootstrap/restart variance, not TMB sdreport.",
    "",
    "## L2 regularized Hessian covariance",
    paste0("- Best regularized Hessian method: ", l2$recommendation$regularization_method[[1L]], " value=", l2$recommendation$regularization_value[[1L]], "."),
    paste0("- Best median absolute log sd ratio: ", fmt_metric(l2$recommendation$median_abs_log_sd_ratio[[1L]]), "."),
    "- Regularized Hessian can be diagnostic, but bootstrap variance remains the safer shell1 uncertainty path unless it matches bootstrap across scopes.",
    "",
    "## L3 deterministic borrowed follow-up",
    paste0("- E7 deployment gate passed: ", l3$deployment$gate_passed[[1L]], "."),
    paste0("- E7 deployment Spearman: ", fmt_metric(l3$deployment$deployment_spearman[[1L]]), "; sd ratio=", fmt_metric(l3$deployment$estimated_delta_sd_ratio[[1L]]), "."),
    paste0("- E7 potential valid_shape count: ", sum(l3$potential$recommended_status == "valid_shape_config", na.rm = TRUE), "."),
    "- E7 remains diagnostic unless deployment and potential gates pass without negative-scale compensation.",
    "",
    "## L4 production-safe calibration dry run",
    paste0("- valid_shape count: ", sum(l4$results$recommended_status == "valid_shape_config", na.rm = TRUE), "."),
    paste0("- normalized collapse rate: ", fmt_metric(l4$summary$normalized_collapse_rate[[1L]]), "."),
    paste0("- mutation wrong-direction rate: ", fmt_metric(l4$summary$mutation_wrong_direction_rate[[1L]]), "."),
    "- Current production-safe behavior is shell0/direct safe local plus failure-state gate; best numeric-only fallback is diagnostic only.",
    "",
    "## L5 shape vs uncertainty separation",
    "- Shape ranking should default to compute_sd=FALSE and evaluate means/rank/amplitude/failure-state.",
    "- Uncertainty calibration should run only on selected configs and may use TMB sd for shell0, bootstrap/restart variance for shell1.",
    "- Uncertainty failure must not force shape ranking failure, but nonfinite covariance must not be used as anchor variance.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- shell_depth=1 remains experimental.",
    "- safe shell0 remains current production local mode.",
    "- Keep normalized as benchmark/probe/calibration default candidate with amplitude-collapse diagnostics and no_valid_shape_configuration gate.",
    "- Do not default anchor_count_reference=minobs for full input.",
    "- Next priorities: robust shell0 production calibration, bootstrap/restart variance for experimental shell1, and non-oracle delta deployment gates after local uncertainty is stabilized."
  )
  writeLines(lines, file.path(dirs$root, "farfield_shell1_empirical_uncertainty_report.md"))
  saveRDS(list(l1 = l1, l2 = l2, l3 = l3, l4 = l4, l5 = l5, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_shell1_empirical_uncertainty_all_results.rds"))
}

main_l <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage(); return(invisible(NULL))
  }
  mode <- match.arg(tolower(as.character(arg_value(args, "mode", "all"))),
                    c("prepare", "l1-bootstrap-empirical-uncertainty",
                      "l2-regularized-hessian-covariance",
                      "l3-deterministic-borrowed-followup",
                      "l4-production-safe-calibration-dryrun",
                      "l5-shape-uncertainty-separation", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  delta_debug_dir <- as.character(arg_value(args, "delta_debug_dir", "benchmark/results/farfield_delta_debug"))
  local_calibration_dir <- as.character(arg_value(args, "local_calibration_dir", "benchmark/results/farfield_local_calibration_patch"))
  core_fix_dir <- as.character(arg_value(args, "core_fix_dir", "benchmark/results/farfield_core_fix_probe"))
  identifiability_dir <- as.character(arg_value(args, "identifiability_dir", "benchmark/results/farfield_local_identifiability_repair"))
  hessian_dir <- as.character(arg_value(args, "hessian_dir", "benchmark/results/farfield_local_hessian_diagnostics"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_shell1_empirical_uncertainty"))
  simulation_ids <- arg_integer_csv(args, "simulation_ids", 1:10)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  bootstrap_reps <- arg_integer(args, "bootstrap_reps", 50L)
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_l_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, 1, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, 1, minobs, input_policy, force = FALSE)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  yi <- readRDS(ctx$input_table$input_rds[[1L]])
  counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = input_policy, drop_diploid = TRUE)
  dt <- as.numeric(ctx$input_table$time_delta[[1L]])
  task_info <- list(simulation_id = 1, minobs = minobs, input_policy = input_policy,
                    lambda = as.numeric(ctx$input_table$lambda[[1L]]), dt = dt,
                    beta = if ("sim_beta" %in% names(ctx$input_table)) as.numeric(ctx$input_table$sim_beta[[1L]]) else 0.00005)
  saveRDS(list(context = ctx, simulation_ids = simulation_ids), file.path(dirs$results, "prepare_context.rds"))
  if (mode == "prepare") return(invisible(dirs$root))
  l1 <- if (mode %in% c("all", "l1-bootstrap-empirical-uncertainty")) run_l1(bundle, counts, dt, grf, task_info, dirs, bootstrap_reps, force = force) else readRDS(file.path(dirs$results, "l1_bootstrap_empirical_uncertainty.rds"))
  if (mode == "l1-bootstrap-empirical-uncertainty") return(invisible(l1))
  l2 <- if (mode %in% c("all", "l2-regularized-hessian-covariance")) run_l2(l1, counts, dt, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "l2_regularized_hessian_covariance.rds"))
  if (mode == "l2-regularized-hessian-covariance") return(invisible(l2))
  l3 <- if (mode %in% c("all", "l3-deterministic-borrowed-followup")) run_l3(l1, counts, dt, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "l3_deterministic_borrowed_followup.rds"))
  if (mode == "l3-deterministic-borrowed-followup") return(invisible(l3))
  l4 <- if (mode %in% c("all", "l4-production-safe-calibration-dryrun")) run_l4(source_input_dir, simulation_ids, minobs, input_policy, dirs, force = force) else readRDS(file.path(dirs$results, "l4_production_safe_calibration_dryrun.rds"))
  if (mode == "l4-production-safe-calibration-dryrun") return(invisible(l4))
  l5 <- if (mode %in% c("all", "l5-shape-uncertainty-separation")) run_l5(l1, l2, l3, l4, dirs, force = force) else readRDS(file.path(dirs$results, "l5_shape_uncertainty_separation.rds"))
  if (mode == "l5-shape-uncertainty-separation") return(invisible(l5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_l_recommendations(l1, l2, l3, l4, l5)
    args_info <- list(source_input_dir = source_input_dir, abcd_dir = abcd_dir,
                      diagnostics_dir = diagnostics_dir, delta_probe_dir = delta_probe_dir,
                      delta_debug_dir = delta_debug_dir, local_calibration_dir = local_calibration_dir,
                      core_fix_dir = core_fix_dir, identifiability_dir = identifiability_dir,
                      hessian_dir = hessian_dir, simulation_ids = simulation_ids,
                      minobs = minobs, input_policy = input_policy, bootstrap_reps = bootstrap_reps)
    write_l_report(dirs, args_info, ctx, l1, l2, l3, l4, l5, recs)
  }
  message("Wrote shell1 empirical uncertainty diagnostics under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) main_l()
