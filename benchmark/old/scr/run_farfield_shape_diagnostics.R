#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_shape_diagnostics.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_shape_abcd.R"))

fmt_metric <- function(x) {
  out <- vapply(as.numeric(x), function(xx) {
    if (!is.finite(xx)) "NA" else format(round(xx, 4), nsmall = 4)
  }, character(1))
  if (length(out) == 1L) out[[1L]] else out
}

usage <- function() {
  cat(
    "Run farfield shape diagnostics after the ABCD probe.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_shape_diagnostics.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --output-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --simulation-id=1 --minobs=5 --input-policy=full\n\n",
    "Modes:\n",
    "  prepare, edge-delta, oracle, local-stability, potential-prior, summarize, all\n\n",
    "Options:\n",
    "  --quick=auto        auto|true|false; auto limits the expensive local TMB grid.\n",
    "  --force=false       recompute cached diagnostic fits.\n",
    "  --cv-splits=8       CV splits for the follow-up diagnostic.\n",
    sep = ""
  )
}

make_diag_dirs <- function(output_dir) {
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

support_scope_label <- function(tier, distance) {
  tier <- as.character(tier)
  distance <- suppressWarnings(as.integer(distance))
  out <- rep("other", length(tier))
  out[tier == "directly_informed"] <- "direct"
  out[tier == "local_borrowed"] <- "local_borrowed"
  far <- distance >= 2L | tier %in% c("weakly_supported", "graph_borrowed", "prior_dominated")
  out[far] <- "farfield"
  out[tier == "weakly_supported"] <- "weakly_supported"
  out
}

compute_truth_for_nodes <- function(labels, grf, lambda) {
  truth_map <- compute_grf_truth(labels, grf$centroids, lambda)
  as.numeric(truth_map[as.character(labels)])
}

clone_local_with_truth_direct_anchors <- function(local_fit, graph, grf, lambda) {
  out <- local_fit
  direct <- as.character(out$summary$support_tier) == "directly_informed"
  truth <- compute_truth_for_nodes(out$summary$karyotype[direct], grf, lambda)
  out$summary$fitness_mean[direct] <- truth
  out$summary$conf_low[direct] <- truth - 1.959963984540054 * out$summary$fitness_sd[direct]
  out$summary$conf_high[direct] <- truth + 1.959963984540054 * out$summary$fitness_sd[direct]
  out$diagnostics$anchor_mean_source <- "truth_direct_anchors"
  out
}

edge_context_table <- function(graph, local_fit, grf, lambda, graph_layer,
                               empirical_delta_context = NULL,
                               prior_mean_local = NULL,
                               prior_mean_empirical = NULL) {
  from <- as.integer(graph$parent_from0) + 1L
  to <- as.integer(graph$parent_to0) + 1L
  ctx <- as.integer(graph$parent_context0) + 1L
  keep <- from >= 1L & to >= 1L & from <= length(graph$labels) & to <= length(graph$labels)
  from <- from[keep]
  to <- to[keep]
  ctx <- ctx[keep]
  if (!length(from)) {
    return(data.frame(
      graph_layer = character(),
      edge_id = integer(),
      parent_node_id = integer(),
      child_node_id = integer(),
      parent_karyotype = character(),
      child_karyotype = character(),
      context_index = integer(),
      context_label = character(),
      edge_direction = character(),
      edge_chr = integer(),
      ploidy_band = integer(),
      support_distance_parent = integer(),
      support_distance_child = integer(),
      support_scope_parent = character(),
      support_scope_child = character(),
      support_tier_parent = character(),
      support_tier_child = character(),
      edge_scope = character(),
      parent_observed = logical(),
      child_observed = logical(),
      truth_delta = numeric(),
      local_delta = numeric(),
      local_context_delta = numeric(),
      empirical_context_delta = numeric(),
      prior_mean_local_context_delta = numeric(),
      prior_mean_empirical_context_delta = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  labels <- as.character(graph$labels)
  truth <- compute_truth_for_nodes(unique(c(labels[from], labels[to])), grf, lambda)
  names(truth) <- unique(c(labels[from], labels[to]))
  local_map <- setNames(as.numeric(local_fit$summary$fitness_mean), as.character(local_fit$summary$karyotype))
  local_context <- local_fit$parameter_mode$delta_context
  if (is.null(local_context)) local_context <- rep(NA_real_, length(graph$context_label))
  if (is.null(empirical_delta_context)) empirical_delta_context <- rep(NA_real_, length(graph$context_label))
  context_label <- as.character(graph$context_label)
  ctx_label <- ifelse(ctx >= 1L & ctx <= length(context_label), context_label[ctx], NA_character_)
  chr <- suppressWarnings(as.integer(sub("^.*chr([0-9]+).*$", "\\1", ctx_label)))
  chr[!is.finite(chr)] <- NA_integer_
  direction <- ifelse(grepl("^gain_", ctx_label), "gain",
                      ifelse(grepl("^loss_", ctx_label), "loss", NA_character_))
  band <- suppressWarnings(as.integer(sub("^.*band([0-9]+).*$", "\\1", ctx_label)))
  band[!is.finite(band)] <- NA_integer_
  parent_tier <- as.character(graph$support_tier[from])
  child_tier <- as.character(graph$support_tier[to])
  parent_distance <- as.integer(graph$support_distance[from])
  child_distance <- as.integer(graph$support_distance[to])
  parent_scope <- support_scope_label(parent_tier, parent_distance)
  child_scope <- support_scope_label(child_tier, child_distance)
  edge_scope <- ifelse(parent_scope == "direct" & child_scope == "direct", "direct_direct",
                       ifelse(parent_scope == "direct" | child_scope == "direct", "direct_borrowed",
                              ifelse(parent_scope != "direct" & child_scope != "direct", "borrowed_borrowed", "other")))
  parent_label <- labels[from]
  child_label <- labels[to]
  truth_delta <- as.numeric(truth[child_label] - truth[parent_label])
  local_parent <- as.numeric(local_map[parent_label])
  local_child <- as.numeric(local_map[child_label])
  local_delta <- local_child - local_parent
  lctx <- rep(NA_real_, length(ctx))
  ectx <- rep(NA_real_, length(ctx))
  ok_ctx <- ctx >= 1L & ctx <= length(local_context)
  lctx[ok_ctx] <- as.numeric(local_context[ctx[ok_ctx]])
  ok_emp <- ctx >= 1L & ctx <= length(empirical_delta_context)
  ectx[ok_emp] <- as.numeric(empirical_delta_context[ctx[ok_emp]])
  pm_local_delta <- if (!is.null(prior_mean_local)) prior_mean_local[to] - prior_mean_local[from] else NA_real_
  pm_empirical_delta <- if (!is.null(prior_mean_empirical)) prior_mean_empirical[to] - prior_mean_empirical[from] else NA_real_
  data.frame(
    graph_layer = graph_layer,
    edge_id = seq_along(from),
    parent_node_id = from,
    child_node_id = to,
    parent_karyotype = parent_label,
    child_karyotype = child_label,
    context_index = ctx,
    context_label = ctx_label,
    edge_direction = direction,
    edge_chr = chr,
    ploidy_band = band,
    support_distance_parent = parent_distance,
    support_distance_child = child_distance,
    support_scope_parent = parent_scope,
    support_scope_child = child_scope,
    support_tier_parent = parent_tier,
    support_tier_child = child_tier,
    edge_scope = edge_scope,
    parent_observed = parent_tier == "directly_informed",
    child_observed = child_tier == "directly_informed",
    truth_delta = truth_delta,
    local_delta = local_delta,
    local_context_delta = lctx,
    empirical_context_delta = ectx,
    prior_mean_local_context_delta = pm_local_delta,
    prior_mean_empirical_context_delta = pm_empirical_delta,
    stringsAsFactors = FALSE
  )
}

empirical_context_delta_from_local <- function(local_fit, graph) {
  local_graph <- local_fit$graph
  ls <- local_fit$summary
  f <- setNames(as.numeric(ls$fitness_mean), as.character(ls$karyotype))
  labels <- as.character(local_graph$labels)
  from <- as.integer(local_graph$parent_from0) + 1L
  to <- as.integer(local_graph$parent_to0) + 1L
  ctx <- as.integer(local_graph$parent_context0) + 1L
  ok <- from >= 1L & to >= 1L & from <= length(labels) & to <= length(labels) &
    is.finite(f[labels[from]]) & is.finite(f[labels[to]]) &
    ctx >= 1L
  n_ctx <- length(graph$context_label)
  out <- rep(NA_real_, n_ctx)
  n_local_edges <- rep(0L, n_ctx)
  if (any(ok)) {
    delta <- as.numeric(f[labels[to[ok]]] - f[labels[from[ok]]])
    ctx_ok <- ctx[ok]
    global_med <- stats::median(delta[is.finite(delta)], na.rm = TRUE)
    if (!is.finite(global_med)) global_med <- 0
    for (cc in seq_len(n_ctx)) {
      vals <- delta[ctx_ok == cc]
      vals <- vals[is.finite(vals)]
      n_local_edges[[cc]] <- length(vals)
      if (length(vals) >= 2L) out[[cc]] <- stats::median(vals)
    }
    out[!is.finite(out)] <- global_med
  } else {
    out[] <- 0
  }
  data.frame(
    context_index = seq_len(n_ctx),
    context_label = as.character(graph$context_label),
    empirical_context_delta = out,
    n_local_edges = n_local_edges,
    stringsAsFactors = FALSE
  )
}

compute_edge_delta_metrics <- function(edge_df, estimator_col, group_cols = character()) {
  if (!nrow(edge_df) || !estimator_col %in% names(edge_df)) return(data.frame())
  x <- edge_df
  x$estimated_delta <- suppressWarnings(as.numeric(x[[estimator_col]]))
  x$truth_delta <- suppressWarnings(as.numeric(x$truth_delta))
  split_key <- if (length(group_cols)) interaction(x[group_cols], drop = TRUE, lex.order = TRUE) else factor(rep("all", nrow(x)))
  rows <- lapply(split(x, split_key), function(d) {
    ok <- is.finite(d$estimated_delta) & is.finite(d$truth_delta)
    est <- d$estimated_delta[ok]
    truth <- d$truth_delta[ok]
    err <- est - truth
    data.frame(
      d[1L, group_cols, drop = FALSE],
      estimator = estimator_col,
      n_edges = nrow(d),
      n_scored_edges = length(est),
      truth_delta_mean = if (length(truth)) mean(truth) else NA_real_,
      truth_delta_sd = if (length(truth) > 1L) stats::sd(truth) else NA_real_,
      estimated_delta_mean = if (length(est)) mean(est) else NA_real_,
      estimated_delta_sd = if (length(est) > 1L) stats::sd(est) else NA_real_,
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
  })
  bind_rows_fill(rows)
}

run_edge_delta_diagnostics <- function(bundle, grf, task_info, dirs) {
  emp <- empirical_context_delta_from_local(bundle$local, bundle$global_graph)
  pm_local <- build_prior_mean_local_context_delta(bundle$local, bundle$global_graph, 1)$mean
  pm_emp <- build_prior_mean_empirical_edge_delta(bundle$local, bundle$global_graph, 1)$mean
  local_edges <- edge_context_table(
    bundle$local$graph, bundle$local, grf, task_info$lambda, "local",
    empirical_delta_context = emp$empirical_context_delta,
    prior_mean_local = pm_local[match(bundle$local$graph$labels, bundle$global_graph$labels)],
    prior_mean_empirical = pm_emp[match(bundle$local$graph$labels, bundle$global_graph$labels)]
  )
  global_edges <- edge_context_table(
    bundle$global_graph, bundle$local, grf, task_info$lambda, "global",
    empirical_delta_context = emp$empirical_context_delta,
    prior_mean_local = pm_local,
    prior_mean_empirical = pm_emp
  )
  by_edge <- rbind(local_edges, global_edges)
  write_tsv_safe(by_edge, file.path(dirs$tables, "edge_delta_diagnostics_by_edge.tsv"))
  estimators <- c("local_delta", "local_context_delta", "empirical_context_delta",
                  "prior_mean_local_context_delta", "prior_mean_empirical_context_delta")
  by_context <- bind_rows_fill(lapply(estimators, function(est) {
    compute_edge_delta_metrics(by_edge, est, c("graph_layer", "context_label", "edge_direction", "edge_chr", "ploidy_band"))
  }))
  by_scope <- bind_rows_fill(lapply(estimators, function(est) {
    compute_edge_delta_metrics(by_edge, est, c("graph_layer", "edge_scope", "support_scope_parent", "support_scope_child"))
  }))
  summary <- bind_rows_fill(lapply(estimators, function(est) {
    compute_edge_delta_metrics(by_edge, est, c("graph_layer"))
  }))
  write_tsv_safe(by_context, file.path(dirs$tables, "edge_delta_diagnostics_by_context.tsv"))
  write_tsv_safe(by_scope, file.path(dirs$tables, "edge_delta_diagnostics_by_scope.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "edge_delta_direction_summary.tsv"))
  saveRDS(list(by_edge = by_edge, by_context = by_context, by_scope = by_scope,
               summary = summary, empirical_context_delta = emp),
          file.path(dirs$results, "edge_delta_diagnostics.rds"))
  list(by_edge = by_edge, by_context = by_context, by_scope = by_scope, summary = summary)
}

load_abcd_results <- function(abcd_dir) {
  abcd_dir <- normalizePath(abcd_dir, winslash = "/", mustWork = TRUE)
  list(
    root = abcd_dir,
    all = read_tsv_safe(file.path(abcd_dir, "tables", "all_experiments_long.tsv")),
    a = read_tsv_safe(file.path(abcd_dir, "tables", "experiment_A_global_scale_grid.tsv")),
    c = read_tsv_safe(file.path(abcd_dir, "tables", "experiment_C_prior_mean_edge_slope.tsv")),
    d_candidates = read_tsv_safe(file.path(abcd_dir, "tables", "experiment_D_shape_cv_candidates.tsv")),
    d_selected = read_tsv_safe(file.path(abcd_dir, "tables", "experiment_D_selected_farfield_evaluation.tsv"))
  )
}

base_config_table <- function(abcd) {
  base <- data.frame(
    experiment = "diag",
    candidate_id = c("diag_baseline_mutation", "diag_normalized_ll0p2_le0p01", "diag_unit_ll0p2_le0p01"),
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
  far <- abcd$all[abcd$all$support_scope == "farfield" &
                    abcd$all$metric_scale == "native", , drop = FALSE]
  top <- if (nrow(far)) {
    configs_from_metrics(head(far[order(far$shape_score, far$centered_rmse), , drop = FALSE], 5L), "diag")
  } else data.frame()
  cfg <- bind_rows_fill(list(base, top))
  dedupe_configs(cfg)
}

run_negative_scale_probe <- function(bundle, components, grf, task_info, dirs, abcd, force = FALSE) {
  base <- base_config_table(abcd)
  priors <- rbind(
    data.frame(prior_mean_mode = "zero", prior_mean_scale = 0, stringsAsFactors = FALSE),
    expand.grid(
      prior_mean_mode = c("local_context_delta", "empirical_edge_delta"),
      prior_mean_scale = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  configs <- do.call(rbind, lapply(seq_len(nrow(priors)), function(i) {
    x <- base
    x$experiment <- "E1b"
    x$prior_mean_mode <- priors$prior_mean_mode[[i]]
    x$prior_mean_scale <- priors$prior_mean_scale[[i]]
    x$candidate_id <- mapply(
      config_id,
      prefix = "E1b",
      graph_edge_weight = x$graph_edge_weight,
      lambda_l = x$lambda_l,
      lambda_e = x$lambda_e,
      sigma_obs = x$sigma_obs,
      anchor_var_mode = x$anchor_var_mode,
      prior_mean_mode = x$prior_mean_mode,
      prior_mean_scale = x$prior_mean_scale,
      anchor_count_reference_mode = x$anchor_count_reference_mode,
      USE.NAMES = FALSE
    )
    x
  }))
  configs <- dedupe_configs(configs)
  rows <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    message("[E1b] ", i, "/", nrow(configs), " ", cfg$candidate_id)
    pm <- build_prior_mean(cfg$prior_mean_mode, bundle$local, bundle$global_graph, cfg$prior_mean_scale)
    rows[[i]] <- fit_cached(
      bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs,
      prior_mean = pm$mean,
      prior_mean_status = pm$status,
      force = force
    )$metrics
  }
  out <- bind_rows_fill(rows)
  write_tsv_safe(out, file.path(dirs$tables, "negative_scale_prior_mean_probe.tsv"))
  far <- out[out$support_scope == "farfield" & out$metric_scale == "native", , drop = FALSE]
  top <- rbind(
    transform(head(far[order(far$shape_score, far$centered_rmse), , drop = FALSE], 20L), ranking_mode = "shape_score"),
    transform(head(far[order(-far$spearman, far$centered_rmse), , drop = FALSE], 20L), ranking_mode = "spearman"),
    transform(head(far[order(-far$pearson, far$centered_rmse), , drop = FALSE], 20L), ranking_mode = "pearson")
  )
  write_tsv_safe(top, file.path(dirs$tables, "negative_scale_prior_mean_top.tsv"))
  saveRDS(list(configs = configs, metrics = out, top = top),
          file.path(dirs$results, "negative_scale_prior_mean_probe.rds"))
  out
}

oracle_context_delta <- function(graph, grf, lambda, per_edge = FALSE) {
  from <- as.integer(graph$parent_from0) + 1L
  to <- as.integer(graph$parent_to0) + 1L
  ctx <- as.integer(graph$parent_context0) + 1L
  labels <- as.character(graph$labels)
  truth <- compute_truth_for_nodes(labels, grf, lambda)
  delta <- truth[to] - truth[from]
  if (isTRUE(per_edge)) {
    return(list(edge_delta = delta, context_delta = NULL))
  }
  n_ctx <- length(graph$context_label)
  context_delta <- rep(0, n_ctx)
  for (cc in seq_len(n_ctx)) {
    vals <- delta[ctx == cc]
    vals <- vals[is.finite(vals)]
    context_delta[[cc]] <- if (length(vals)) stats::median(vals) else 0
  }
  list(edge_delta = delta, context_delta = context_delta)
}

fit_potential_prior_mean <- function(local_fit, graph, delta_mode, grf, lambda,
                                     delta_scale = 1, lambda_anchor = 10,
                                     lambda_smooth = 0.1, components = NULL,
                                     ridge = 1e-6) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  n <- length(graph$labels)
  from <- as.integer(graph$parent_from0) + 1L
  to <- as.integer(graph$parent_to0) + 1L
  ctx <- as.integer(graph$parent_context0) + 1L
  keep <- from >= 1L & to >= 1L & from <= n & to <= n
  from <- from[keep]
  to <- to[keep]
  ctx <- ctx[keep]
  w <- as.numeric(graph$parent_weight)[keep]
  w[!is.finite(w) | w <= 0] <- 1
  med_w <- stats::median(w[w > 0], na.rm = TRUE)
  if (!is.finite(med_w) || med_w <= 0) med_w <- 1
  w <- w / med_w
  delta_mode <- as.character(delta_mode)
  delta <- rep(0, length(from))
  component_context <- NULL
  if (identical(delta_mode, "empirical_context_delta")) {
    component_context <- empirical_context_delta_from_local(local_fit, graph)
    delta <- component_context$empirical_context_delta[ctx]
  } else if (identical(delta_mode, "local_context_delta")) {
    dc <- local_fit$parameter_mode$delta_context
    if (is.null(dc)) dc <- rep(0, length(graph$context_label))
    delta <- as.numeric(dc[ctx])
  } else if (identical(delta_mode, "oracle_context_delta")) {
    od <- oracle_context_delta(graph, grf, lambda, per_edge = FALSE)
    delta <- od$context_delta[ctx]
  } else if (identical(delta_mode, "oracle_per_edge_delta")) {
    od <- oracle_context_delta(graph, grf, lambda, per_edge = TRUE)
    delta <- od$edge_delta[keep]
  } else if (identical(delta_mode, "zero")) {
    delta[] <- 0
  } else {
    stop("Unsupported potential delta mode: ", delta_mode, call. = FALSE)
  }
  delta[!is.finite(delta)] <- 0
  delta <- as.numeric(delta_scale) * delta
  sw <- sqrt(w)
  row_id <- seq_along(from)
  a <- Matrix::sparseMatrix(
    i = c(row_id, row_id),
    j = c(from, to),
    x = c(-sw, sw),
    dims = c(length(from), n),
    giveCsparse = TRUE
  )
  rhs_edges <- sw * delta
  q <- Matrix::crossprod(a)
  rhs <- as.numeric(Matrix::crossprod(a, rhs_edges))
  direct <- as.character(local_fit$summary$support_tier) == "directly_informed" &
    is.finite(local_fit$summary$fitness_mean)
  anchor_idx <- match(as.character(local_fit$summary$karyotype[direct]), as.character(graph$labels))
  anchor_ok <- !is.na(anchor_idx)
  anchor_idx <- anchor_idx[anchor_ok]
  anchor_mean <- as.numeric(local_fit$summary$fitness_mean[direct][anchor_ok])
  if (length(anchor_idx) && is.finite(lambda_anchor) && lambda_anchor > 0) {
    q <- q + Matrix::sparseMatrix(i = anchor_idx, j = anchor_idx, x = lambda_anchor, dims = c(n, n))
    rhs[anchor_idx] <- rhs[anchor_idx] + lambda_anchor * anchor_mean
  }
  if (!is.null(components) && is.finite(lambda_smooth) && lambda_smooth > 0) {
    q <- q + as.numeric(lambda_smooth) * components$edge$normalized
  }
  q <- q + Matrix::Diagonal(n, ridge)
  q <- Matrix::forceSymmetric(q, uplo = "U")
  started <- Sys.time()
  mean <- tryCatch({
    chol <- Matrix::Cholesky(q, LDL = TRUE, perm = TRUE)
    as.numeric(Matrix::solve(chol, rhs))
  }, error = function(e) {
    warning("Potential prior solve failed: ", conditionMessage(e), call. = FALSE)
    rep(0, n)
  })
  list(
    mean = mean,
    status = paste0("potential_", delta_mode),
    components = data.frame(
      delta_mode = delta_mode,
      delta_scale = delta_scale,
      lambda_anchor = lambda_anchor,
      lambda_smooth = lambda_smooth,
      n_edges = length(from),
      n_anchors = length(anchor_idx),
      mean_sd = stats::sd(mean),
      mean_range = quantile_range(mean),
      elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      stringsAsFactors = FALSE
    ),
    context_components = component_context
  )
}

build_potential_configs <- function(abcd) {
  base <- base_config_table(abcd)
  # Keep the potential grid diagnostic but bounded enough for repeated Matrix fits.
  grid <- expand.grid(
    delta_mode = c("empirical_context_delta", "local_context_delta", "oracle_context_delta", "oracle_per_edge_delta"),
    delta_scale = c(-1, -0.5, 0.5, 1),
    lambda_anchor = c(1, 10, 100),
    lambda_smooth = c(0, 0.1, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  list(base = base, grid = grid)
}

run_potential_prior_probe <- function(bundle, components, grf, task_info, dirs, abcd,
                                      force = FALSE, quick = "auto") {
  result_path <- file.path(dirs$results, "potential_prior_mean_probe.rds")
  if (!isTRUE(force) && file.exists(result_path)) {
    cached <- readRDS(result_path)
    if (!is.null(cached$metrics) && nrow(cached$metrics)) return(cached$metrics)
  }
  pc <- build_potential_configs(abcd)
  base <- pc$base
  grid <- pc$grid
  if (quick %in% c("auto", "true")) {
    base <- head(base, min(5L, nrow(base)))
    grid <- grid[grid$lambda_anchor %in% c(10, 100) & grid$lambda_smooth %in% c(0.1) &
                   grid$delta_scale %in% c(-1, -0.5, 0.5, 1), , drop = FALSE]
  }
  rows <- list()
  comp_rows <- list()
  idx <- 0L
  for (b in seq_len(nrow(base))) {
    for (g in seq_len(nrow(grid))) {
      cfg <- base[b, , drop = FALSE]
      gg <- grid[g, , drop = FALSE]
      cfg$experiment <- "E4"
      cfg$prior_mean_mode <- paste0("potential_", gg$delta_mode)
      cfg$prior_mean_scale <- gg$delta_scale
      cfg$candidate_id <- paste(
        config_id("E4", cfg$graph_edge_weight, cfg$lambda_l, cfg$lambda_e, cfg$sigma_obs,
                  cfg$anchor_var_mode, cfg$prior_mean_mode, cfg$prior_mean_scale,
                  cfg$anchor_count_reference_mode),
        paste0("la", gg$lambda_anchor),
        paste0("ls", gg$lambda_smooth),
        sep = "__"
      )
      message("[E4] ", b, "/", nrow(base), " grid ", g, "/", nrow(grid), " ", cfg$candidate_id)
      pm <- fit_potential_prior_mean(
        bundle$local, bundle$global_graph, gg$delta_mode, grf, task_info$lambda,
        delta_scale = gg$delta_scale,
        lambda_anchor = gg$lambda_anchor,
        lambda_smooth = gg$lambda_smooth,
        components = components
      )
      fit <- fit_cached(
        bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs,
        prior_mean = pm$mean,
        prior_mean_status = pm$status,
        force = force
      )
      idx <- idx + 1L
      rows[[idx]] <- transform(
        fit$metrics,
        potential_delta_mode = gg$delta_mode,
        potential_delta_scale = gg$delta_scale,
        potential_lambda_anchor = gg$lambda_anchor,
        potential_lambda_smooth = gg$lambda_smooth
      )
      cc <- pm$components
      cc$candidate_id <- cfg$candidate_id
      cc$graph_edge_weight <- cfg$graph_edge_weight
      cc$lambda_l <- cfg$lambda_l
      cc$lambda_e <- cfg$lambda_e
      cc$sigma_obs <- cfg$sigma_obs
      comp_rows[[idx]] <- cc
    }
  }
  out <- bind_rows_fill(rows)
  comp <- bind_rows_fill(comp_rows)
  write_tsv_safe(out, file.path(dirs$tables, "potential_prior_mean_probe.tsv"))
  write_tsv_safe(comp, file.path(dirs$tables, "potential_prior_mean_components.tsv"))
  far <- out[out$support_scope == "farfield" & out$metric_scale == "native", , drop = FALSE]
  top <- rbind(
    transform(head(far[order(far$shape_score, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "shape_score"),
    transform(head(far[order(-far$spearman, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "spearman"),
    transform(head(far[order(-far$pearson, far$centered_rmse), , drop = FALSE], 30L), ranking_mode = "pearson")
  )
  write_tsv_safe(top, file.path(dirs$tables, "potential_prior_mean_top.tsv"))
  saveRDS(list(metrics = out, components = comp, top = top, base = base, grid = grid),
          file.path(dirs$results, "potential_prior_mean_probe.rds"))
  out
}

run_oracle_controls <- function(bundle, components, grf, task_info, dirs, abcd, force = FALSE) {
  base <- base_config_table(abcd)
  forced <- data.frame(
    experiment = "E2",
    candidate_id = c("oracle_baseline", "oracle_normalized", "oracle_unit"),
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
  base <- dedupe_configs(rbind(forced, head(base, 5L)))
  truth_local <- clone_local_with_truth_direct_anchors(bundle$local, bundle$global_graph, grf, task_info$lambda)
  controls <- data.frame(
    oracle_control = c("O1_local_anchors_oracle_per_edge_delta",
                       "O2_truth_anchors_estimated_edge_delta",
                       "O3_truth_anchors_oracle_context_delta",
                       "O4_truth_anchors_zero_mean_graph",
                       "O5_truth_anchors_oracle_per_edge_delta"),
    anchor_source = c("local", "truth", "truth", "truth", "truth"),
    delta_mode = c("oracle_per_edge_delta", "empirical_context_delta", "oracle_context_delta", "zero", "oracle_per_edge_delta"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  idx <- 0L
  for (b in seq_len(nrow(base))) {
    for (o in seq_len(nrow(controls))) {
      cfg <- base[b, , drop = FALSE]
      ctrl <- controls[o, , drop = FALSE]
      local_use <- if (ctrl$anchor_source == "truth") truth_local else bundle$local
      cfg$experiment <- "E2"
      cfg$prior_mean_mode <- ctrl$delta_mode
      cfg$prior_mean_scale <- 1
      cfg$candidate_id <- paste0("E2__", ctrl$oracle_control, "__", cfg$graph_edge_weight,
                                 "__ll", cfg$lambda_l, "__le", cfg$lambda_e, "__so", cfg$sigma_obs)
      message("[E2] base ", b, "/", nrow(base), " control ", ctrl$oracle_control)
      if (identical(ctrl$delta_mode, "zero")) {
        prior <- build_prior_mean_zero(bundle$global_graph)
      } else {
        prior <- fit_potential_prior_mean(
          local_use, bundle$global_graph, ctrl$delta_mode, grf, task_info$lambda,
          delta_scale = 1, lambda_anchor = 10, lambda_smooth = 0.1, components = components
        )
      }
      fit <- fit_cached(local_use, bundle$global_graph, components, cfg, grf, task_info, dirs,
                        prior_mean = prior$mean, prior_mean_status = prior$status, force = force)
      idx <- idx + 1L
      rows[[idx]] <- transform(
        fit$metrics,
        oracle_control = ctrl$oracle_control,
        anchor_source = ctrl$anchor_source,
        delta_mode = ctrl$delta_mode
      )
    }
  }
  out <- bind_rows_fill(rows)
  by_scope <- out
  write_tsv_safe(out[out$support_scope %in% c("direct", "farfield", "all"), , drop = FALSE],
                 file.path(dirs$tables, "oracle_controls.tsv"))
  write_tsv_safe(by_scope, file.path(dirs$tables, "oracle_controls_by_scope.tsv"))
  saveRDS(list(metrics = out, controls = controls, base = base),
          file.path(dirs$results, "oracle_controls.rds"))
  out
}

local_grid <- function(quick) {
  if (quick %in% c("auto", "true")) {
    return(data.frame(
      local_shell_depth = c(1, 1, 1, 1, 0, 0, 1, 1),
      observation_model = c("dirichlet_multinomial", "dirichlet_multinomial", "dirichlet_multinomial",
                            "multinomial", "dirichlet_multinomial", "multinomial", "dirichlet_multinomial", "multinomial"),
      dm_concentration = c(50, 50, 50, NA, 50, NA, 200, NA),
      eval_max = c(500, 500, 500, 500, 500, 500, 2000, 2000),
      iter_max = c(500, 500, 500, 500, 500, 500, 2000, 2000),
      restart_id = c(1, 2, 3, 1, 1, 1, 1, 1),
      stringsAsFactors = FALSE
    ))
  }
  grid <- rbind(
    expand.grid(
      local_shell_depth = c(0, 1),
      observation_model = "multinomial",
      dm_concentration = NA_real_,
      eval_max = c(500, 2000, 5000),
      iter_max = c(500, 2000, 5000),
      restart_id = c(1, 2, 3),
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
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  grid
}

local_fit_cache_id <- function(row) {
  paste0("local_shell", row$local_shell_depth,
         "_", row$observation_model,
         "_dm", ifelse(is.na(row$dm_concentration), "NA", row$dm_concentration),
         "_eval", row$eval_max,
         "_iter", row$iter_max,
         "_restart", row$restart_id)
}

run_one_local_fit_probe <- function(bundle, row, dirs, force = FALSE) {
  cache_path <- file.path(dirs$cache, paste0(local_fit_cache_id(row), ".rds"))
  if (!isTRUE(force) && file.exists(cache_path)) return(readRDS(cache_path))
  graph <- if (as.integer(row$local_shell_depth) == as.integer(bundle$local$graph$shell_depth)) {
    bundle$local$graph
  } else {
    max_cn <- max(bundle$global_graph$karyotypes, na.rm = TRUE)
    alfak2::build_karyotype_graph(
      bundle$data,
      transition_kernel = "exact",
      shell_depth = as.integer(row$local_shell_depth),
      min_cn = 0L,
      max_cn = as.integer(max_cn),
      max_nodes = 150000L
    )
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
      retry_control = NULL
    )
  }, error = function(e) {
    list(error = conditionMessage(e), graph = graph, summary = data.frame(), diagnostics = list())
  })
  out <- list(fit = fit, config = row, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")))
  saveRDS(out, cache_path)
  out
}

local_edge_alignment <- function(fit_obj, grf, task_info, config_id_value) {
  if (!nrow(fit_obj$fit$summary)) return(data.frame())
  edges <- edge_context_table(fit_obj$fit$graph, fit_obj$fit, grf, task_info$lambda, "local_probe")
  out <- compute_edge_delta_metrics(edges, "local_delta", character())
  if (!nrow(out)) {
    return(data.frame(
      config_id = config_id_value,
      estimator = "local_delta",
      n_edges = 0L,
      n_scored_edges = 0L,
      truth_delta_mean = NA_real_,
      truth_delta_sd = NA_real_,
      estimated_delta_mean = NA_real_,
      estimated_delta_sd = NA_real_,
      delta_bias = NA_real_,
      delta_mae = NA_real_,
      delta_rmse = NA_real_,
      delta_pearson = NA_real_,
      delta_spearman = NA_real_,
      delta_sign_agreement = NA_real_,
      delta_sign_agreement_nonzero = NA_real_,
      fraction_truth_positive = NA_real_,
      fraction_estimated_positive = NA_real_,
      median_abs_delta_error = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  out$config_id <- config_id_value
  out
}

run_local_stability_grid <- function(bundle, grf, task_info, dirs, quick = "auto", force = FALSE) {
  grid <- local_grid(quick)
  rows <- list()
  metrics_rows <- list()
  edge_rows <- list()
  fits <- list()
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    cid <- local_fit_cache_id(row)
    message("[E3] ", i, "/", nrow(grid), " ", cid)
    res <- run_one_local_fit_probe(bundle, row, dirs, force = force)
    fits[[cid]] <- res$fit
    fit <- res$fit
    diag <- fit$diagnostics %||% list()
    rows[[i]] <- data.frame(
      config_id = cid,
      local_shell_depth = row$local_shell_depth,
      observation_model = row$observation_model,
      dm_concentration = row$dm_concentration,
      eval_max = row$eval_max,
      iter_max = row$iter_max,
      restart_id = row$restart_id,
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
    if (nrow(fit$summary %||% data.frame())) {
      cfg <- make_baseline_config("E3")
      cfg$candidate_id <- cid
      m <- score_summary_abcd(fit$summary, fit$graph, grf, task_info$lambda, task_info, cfg, "local_fit")
      m$config_id <- cid
      metrics_rows[[i]] <- m[m$support_scope %in% c("direct", "local_borrowed", "weakly_supported", "farfield", "all"), , drop = FALSE]
      edge_rows[[i]] <- local_edge_alignment(res, grf, task_info, cid)
    }
  }
  grid_tbl <- bind_rows_fill(rows)
  metric_tbl <- bind_rows_fill(metrics_rows)
  if (nrow(metric_tbl)) {
    keep <- metric_tbl$metric_scale == "native" & metric_tbl$support_scope %in% c("direct", "local_borrowed", "weakly_supported")
    wide <- metric_tbl[keep, c("config_id", "support_scope", "centered_rmse", "pearson", "spearman"), drop = FALSE]
    for (scope in unique(wide$support_scope)) {
      ss <- wide[wide$support_scope == scope, , drop = FALSE]
      names(ss)[names(ss) == "centered_rmse"] <- paste0(scope, "_native_centered_rmse")
      names(ss)[names(ss) == "pearson"] <- paste0(scope, "_native_pearson")
      names(ss)[names(ss) == "spearman"] <- paste0(scope, "_native_spearman")
      ss$support_scope <- NULL
      grid_tbl <- merge(grid_tbl, ss, by = "config_id", all.x = TRUE)
    }
  }
  edge_tbl <- bind_rows_fill(edge_rows)
  pair_rows <- list()
  pair_idx <- 0L
  fit_keys <- names(fits)
  if (length(fit_keys) >= 2L) {
    meta <- grid_tbl[, c("config_id", "local_shell_depth", "observation_model", "dm_concentration", "eval_max", "iter_max", "restart_id"), drop = FALSE]
    meta$key <- paste(meta$local_shell_depth, meta$observation_model, meta$dm_concentration, meta$eval_max, meta$iter_max, sep = "|")
    for (key in unique(meta$key)) {
      ids <- meta$config_id[meta$key == key]
      if (length(ids) < 2L) next
      cmb <- utils::combn(ids, 2L)
      for (j in seq_len(ncol(cmb))) {
        a <- fits[[cmb[1L, j]]]
        b <- fits[[cmb[2L, j]]]
        if (!nrow(a$summary %||% data.frame()) || !nrow(b$summary %||% data.frame())) next
        common <- intersect(as.character(a$summary$karyotype), as.character(b$summary$karyotype))
        fa <- setNames(a$summary$fitness_mean, a$summary$karyotype)[common]
        fb <- setNames(b$summary$fitness_mean, b$summary$karyotype)[common]
        tier <- setNames(as.character(a$summary$support_tier), a$summary$karyotype)[common]
        pair_idx <- pair_idx + 1L
        pair_rows[[pair_idx]] <- data.frame(
          key = key,
          config_id_a = cmb[1L, j],
          config_id_b = cmb[2L, j],
          n_common = length(common),
          fitness_mean_correlation = safe_cor2(fa, fb, "pearson"),
          direct_nodes_correlation = safe_cor2(fa[tier == "directly_informed"], fb[tier == "directly_informed"], "pearson"),
          borrowed_nodes_correlation = safe_cor2(fa[tier != "directly_informed"], fb[tier != "directly_informed"], "pearson"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pairwise <- bind_rows_fill(pair_rows)
  conv <- aggregate(
    grid_tbl[, c("gradient_norm", "elapsed_sec")],
    by = list(local_shell_depth = grid_tbl$local_shell_depth, observation_model = grid_tbl$observation_model,
              eval_max = grid_tbl$eval_max),
    FUN = function(x) median(as.numeric(x), na.rm = TRUE)
  )
  conv$n_configs <- as.integer(table(interaction(grid_tbl$local_shell_depth, grid_tbl$observation_model, grid_tbl$eval_max))[interaction(conv$local_shell_depth, conv$observation_model, conv$eval_max)])
  write_tsv_safe(grid_tbl, file.path(dirs$tables, "local_stability_grid.tsv"))
  write_tsv_safe(pairwise, file.path(dirs$tables, "local_stability_pairwise.tsv"))
  write_tsv_safe(edge_tbl, file.path(dirs$tables, "local_edge_delta_truth_alignment.tsv"))
  write_tsv_safe(conv, file.path(dirs$tables, "local_convergence_summary.tsv"))
  saveRDS(list(grid = grid_tbl, metrics = metric_tbl, edge_alignment = edge_tbl,
               pairwise = pairwise, convergence_summary = conv, quick_mode = quick),
          file.path(dirs$results, "local_stability_probe.rds"))
  list(grid = grid_tbl, metrics = metric_tbl, edge_alignment = edge_tbl, pairwise = pairwise)
}

classify_shape_config <- function(far, baseline_rmse = NA_real_) {
  sd_ok <- is.finite(far$estimate_sd_ratio) & far$estimate_sd_ratio >= 0.02
  valid <- sd_ok & is.finite(far$pearson) & is.finite(far$spearman) & far$pearson > 0 & far$spearman > 0
  wrong <- sd_ok & (far$pearson < 0 | far$spearman < 0)
  collapsed <- !sd_ok
  numeric_only <- !valid & !wrong & !collapsed
  if (is.finite(baseline_rmse)) {
    numeric_only <- numeric_only | (!valid & is.finite(far$centered_rmse) & far$centered_rmse < baseline_rmse)
  }
  out <- rep("numeric_only", nrow(far))
  out[wrong] <- "noncollapsed_wrong_direction"
  out[collapsed] <- "collapsed_shrinkage"
  out[valid] <- "valid_shape"
  out[numeric_only & !wrong & !collapsed & !valid] <- "numeric_only"
  out
}

dedupe_diag_configs <- function(x) {
  if (!nrow(x)) return(x)
  key_cols <- intersect(c(
    "candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
    "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode",
    "potential_delta_mode", "potential_delta_scale", "potential_lambda_anchor", "potential_lambda_smooth"
  ), names(x))
  key <- interaction(x[key_cols], drop = TRUE)
  x[!duplicated(key), , drop = FALSE]
}

configs_from_potential_metrics <- function(metrics, experiment = "CVdiag") {
  cols <- c(
    "candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
    "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode",
    "potential_delta_mode", "potential_delta_scale", "potential_lambda_anchor", "potential_lambda_smooth"
  )
  out <- unique(metrics[, intersect(cols, names(metrics)), drop = FALSE])
  out$experiment <- experiment
  out$solver <- "matrix_mean"
  for (nm in setdiff(cols, names(out))) out[[nm]] <- NA
  out[, c("experiment", cols, "solver"), drop = FALSE]
}

build_prior_for_config_diag <- function(cfg, bundle, components, grf, lambda) {
  mode <- as.character(cfg$prior_mean_mode[[1L]])
  if (startsWith(mode, "potential_")) {
    delta_mode <- as.character(cfg$potential_delta_mode[[1L]])
    if (!nzchar(delta_mode) || is.na(delta_mode)) delta_mode <- sub("^potential_", "", mode)
    delta_scale <- suppressWarnings(as.numeric(cfg$potential_delta_scale[[1L]]))
    if (!is.finite(delta_scale)) delta_scale <- suppressWarnings(as.numeric(cfg$prior_mean_scale[[1L]]))
    lambda_anchor <- suppressWarnings(as.numeric(cfg$potential_lambda_anchor[[1L]]))
    lambda_smooth <- suppressWarnings(as.numeric(cfg$potential_lambda_smooth[[1L]]))
    if (!is.finite(lambda_anchor)) lambda_anchor <- 10
    if (!is.finite(lambda_smooth)) lambda_smooth <- 0.1
    return(fit_potential_prior_mean(
      bundle$local, bundle$global_graph, delta_mode, grf, lambda,
      delta_scale = delta_scale,
      lambda_anchor = lambda_anchor,
      lambda_smooth = lambda_smooth,
      components = components
    ))
  }
  build_prior_mean(mode, bundle$local, bundle$global_graph, cfg$prior_mean_scale[[1L]])
}

run_direct_anchor_holdout_cv_diag <- function(bundle, components, configs, task_info, grf,
                                             cv_splits = 8L, seed = 89137L) {
  local <- bundle$local
  graph <- bundle$global_graph
  direct <- local$summary[
    as.character(local$summary$support_tier) == "directly_informed" &
      is.finite(local$summary$fitness_mean),
    ,
    drop = FALSE
  ]
  labels <- intersect(as.character(direct$karyotype), as.character(graph$labels))
  if (length(labels) < 4L) cv_splits <- min(cv_splits, length(labels))
  set.seed(seed)
  split_list <- vector("list", cv_splits)
  for (s in seq_len(cv_splits)) {
    n_hold <- if (length(labels) <= 6L) 1L else max(2L, min(length(labels) - 3L, round(length(labels) * 0.25)))
    split_list[[s]] <- sample(labels, n_hold)
  }
  truth_anchor <- setNames(as.numeric(local$summary$fitness_mean), as.character(local$summary$karyotype))
  rows <- list()
  idx <- 0L
  potential_cols <- c("potential_delta_mode", "potential_delta_scale", "potential_lambda_anchor", "potential_lambda_smooth")
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    pm <- build_prior_for_config_diag(cfg, bundle, components, grf, task_info$lambda)
    for (s in seq_along(split_list)) {
      hold <- split_list[[s]]
      fit <- fit_global_with_config(local, graph, components, cfg, task_info$minobs,
                                    prior_mean = pm$mean, anchor_exclude = hold)
      pred <- fit$summary$fitness_mean[match(hold, as.character(fit$summary$karyotype))]
      truth <- as.numeric(truth_anchor[hold])
      m <- cv_metric_row(pred, truth)
      base <- data.frame(
        split_id = s,
        n_holdout = length(hold),
        holdout_labels = paste(hold, collapse = ","),
        candidate_id = cfg$candidate_id,
        graph_edge_weight = cfg$graph_edge_weight,
        lambda_l = cfg$lambda_l,
        lambda_e = cfg$lambda_e,
        sigma_obs = cfg$sigma_obs,
        anchor_var_mode = cfg$anchor_var_mode,
        prior_mean_mode = cfg$prior_mean_mode,
        prior_mean_scale = cfg$prior_mean_scale,
        anchor_count_reference_mode = cfg$anchor_count_reference_mode,
        prior_mean_status = pm$status,
        stringsAsFactors = FALSE
      )
      for (nm in potential_cols) base[[nm]] <- if (nm %in% names(cfg)) cfg[[nm]] else NA
      idx <- idx + 1L
      rows[[idx]] <- cbind(base, m)
    }
  }
  bind_rows_fill(rows)
}

aggregate_cv_candidates_diag <- function(cv_splits) {
  key_cols <- intersect(c(
    "candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
    "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode",
    "potential_delta_mode", "potential_delta_scale", "potential_lambda_anchor", "potential_lambda_smooth"
  ), names(cv_splits))
  parts <- split(cv_splits, interaction(cv_splits[key_cols], drop = TRUE))
  rows <- lapply(parts, function(x) {
    data.frame(
      x[1L, key_cols, drop = FALSE],
      n_splits = nrow(x),
      cv_mse = median(x$cv_mse, na.rm = TRUE),
      cv_centered_rmse = median(x$cv_centered_rmse, na.rm = TRUE),
      cv_spearman = median(x$cv_spearman, na.rm = TRUE),
      cv_pairwise_rank_loss = median(x$cv_pairwise_rank_loss, na.rm = TRUE),
      cv_estimate_sd_ratio = median(x$cv_estimate_sd_ratio, na.rm = TRUE),
      cv_shape_score = median(x$cv_shape_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

run_cv_diagnostic <- function(bundle, components, task_info, dirs, abcd, potential_metrics, grf,
                              cv_splits = 8L, quick = "auto") {
  d_cfg <- if (nrow(abcd$d_candidates)) configs_from_metrics(abcd$d_candidates, "CVdiag") else data.frame()
  pot_far <- potential_metrics[potential_metrics$support_scope == "farfield" &
                                 potential_metrics$metric_scale == "native" &
                                 !grepl("oracle_", potential_metrics$prior_mean_mode), , drop = FALSE]
  pot_cfg <- if (nrow(pot_far)) configs_from_potential_metrics(head(pot_far[order(pot_far$shape_score), , drop = FALSE], 12L), "CVdiag") else data.frame()
  cfg <- dedupe_diag_configs(bind_rows_fill(list(d_cfg, pot_cfg, make_baseline_config("CVdiag"))))
  if (quick %in% c("auto", "true")) {
    cfg <- head(cfg, min(nrow(cfg), 24L))
    cv_splits <- min(cv_splits, 8L)
  }
  cv_splits_tbl <- run_direct_anchor_holdout_cv_diag(bundle, components, cfg, task_info, grf, cv_splits = cv_splits, seed = 89137L)
  cv_candidates <- aggregate_cv_candidates_diag(cv_splits_tbl)
  selected <- select_by_cv_objectives(cv_candidates)
  write_tsv_safe(selected, file.path(dirs$tables, "cv_diagnostic_selected_configs.tsv"))
  all_far <- bind_rows_fill(list(
    abcd$all[abcd$all$support_scope == "farfield" & abcd$all$metric_scale == "native", , drop = FALSE],
    potential_metrics[potential_metrics$support_scope == "farfield" & potential_metrics$metric_scale == "native", , drop = FALSE]
  ))
  baseline <- all_far[all_far$graph_edge_weight == "mutation" & all_far$lambda_l == 0.2 &
                        all_far$lambda_e == 1 & all_far$sigma_obs == 0.05 &
                        all_far$prior_mean_mode == "zero", , drop = FALSE]
  baseline_rmse <- if (nrow(baseline)) baseline$centered_rmse[[1L]] else NA_real_
  class_tbl <- unique(all_far[, intersect(c(
    "candidate_id", "experiment", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
    "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode",
    "centered_rmse", "pearson", "spearman", "estimate_sd_ratio", "shape_score",
    "amplitude_collapse"
  ), names(all_far)), drop = FALSE])
  class_tbl$shape_class <- classify_shape_config(class_tbl, baseline_rmse)
  class_tbl$selected_by_cv_objective <- vapply(seq_len(nrow(class_tbl)), function(i) {
    paste(selected$objective[selected$candidate_id == class_tbl$candidate_id[[i]]], collapse = ",")
  }, character(1))
  class_tbl$selected_by_cv_objective[class_tbl$selected_by_cv_objective == ""] <- NA_character_
  write_tsv_safe(class_tbl, file.path(dirs$tables, "config_shape_classification.tsv"))
  saveRDS(list(splits = cv_splits_tbl, candidates = cv_candidates, selected = selected,
               classification = class_tbl),
          file.path(dirs$results, "cv_diagnostic.rds"))
  list(selected = selected, classification = class_tbl)
}

make_diagnostic_summary <- function(edge_summary, neg, oracle, local, potential, cv) {
  far_neg <- neg[neg$support_scope == "farfield" & neg$metric_scale == "native", , drop = FALSE]
  far_oracle <- oracle[oracle$support_scope == "farfield" & oracle$metric_scale == "native", , drop = FALSE]
  far_pot <- potential[potential$support_scope == "farfield" & potential$metric_scale == "native", , drop = FALSE]
  rows <- list(
    data.frame(
      diagnostic = "E1_edge_delta",
      key_result = "global local_delta sign agreement",
      value = edge_summary$delta_sign_agreement[edge_summary$graph_layer == "global" & edge_summary$estimator == "local_delta"][1] %||% NA_real_,
      interpretation = "Values near or below 0.5 indicate unreliable local edge direction.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "E1b_negative_scale",
      key_result = "best farfield spearman",
      value = if (nrow(far_neg)) max(far_neg$spearman, na.rm = TRUE) else NA_real_,
      interpretation = "Compares negative and positive prior-mean scales.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "E2_oracle",
      key_result = "best oracle farfield spearman",
      value = if (nrow(far_oracle)) max(far_oracle$spearman, na.rm = TRUE) else NA_real_,
      interpretation = "Oracle controls separate architecture from estimated delta quality.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "E3_local_stability",
      key_result = "median gradient norm",
      value = if (nrow(local$grid)) stats::median(local$grid$gradient_norm, na.rm = TRUE) else NA_real_,
      interpretation = "Large gradients keep local covariance untrusted.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "E4_potential_prior",
      key_result = "best non-oracle farfield spearman",
      value = if (nrow(far_pot)) max(far_pot$spearman[!grepl("oracle_", far_pot$prior_mean_mode)], na.rm = TRUE) else NA_real_,
      interpretation = "Estimated potential priors must improve rank without collapse.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      diagnostic = "CV_followup",
      key_result = "n valid_shape configs",
      value = if (!is.null(cv$classification) && nrow(cv$classification)) sum(cv$classification$shape_class == "valid_shape", na.rm = TRUE) else NA_real_,
      interpretation = "Direct-anchor CV can select only among available valid candidates.",
      stringsAsFactors = FALSE
    )
  )
  bind_rows_fill(rows)
}

make_recommended_next_steps <- function(edge_summary, oracle, potential, local) {
  global_local <- edge_summary[edge_summary$graph_layer == "global" & edge_summary$estimator == "local_delta", , drop = FALSE]
  far_oracle <- oracle[oracle$support_scope == "farfield" & oracle$metric_scale == "native", , drop = FALSE]
  far_pot <- potential[potential$support_scope == "farfield" & potential$metric_scale == "native", , drop = FALSE]
  oracle_best <- if (nrow(far_oracle)) far_oracle[order(-far_oracle$spearman, far_oracle$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  estimated_best <- far_pot[!grepl("oracle_", far_pot$prior_mean_mode), , drop = FALSE]
  estimated_best <- if (nrow(estimated_best)) estimated_best[order(-estimated_best$spearman, estimated_best$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  data.frame(
    priority = 1:5,
    recommendation = c(
      "Do not implement C++ edge-gradient pseudo-observation yet.",
      "Fix or replace the edge-delta estimator before adding a production gradient prior.",
      "Prioritize local convergence/covariance diagnostics and robust anchor variance handling.",
      "Keep normalized as the benchmark/probe default, but keep amplitude-collapse failure labels in calibration.",
      "Extend CV from direct-anchor holdout to neighbor, shell-stratified, or edge-delta holdout."
    ),
    evidence = c(
      paste0("global local_delta sign agreement=", fmt_metric(global_local$delta_sign_agreement[1] %||% NA_real_)),
      paste0("best estimated potential spearman=", fmt_metric(estimated_best$spearman[1] %||% NA_real_),
             "; best oracle spearman=", fmt_metric(oracle_best$spearman[1] %||% NA_real_)),
      paste0("median local gradient=", fmt_metric(stats::median(local$grid$gradient_norm, na.rm = TRUE))),
      "normalized reduces numeric RMSE but often collapses amplitude without an explicit penalty.",
      "direct-anchor CV mostly distinguishes collapse from non-collapse and may not detect farfield direction errors."
    ),
    stringsAsFactors = FALSE
  )
}

write_diag_report <- function(dirs, ctx, task_info, abcd, edge, neg, oracle, local,
                              potential, cv, summary_tbl, next_steps, quick) {
  far_abcd <- abcd$all[abcd$all$support_scope == "farfield" & abcd$all$metric_scale == "native", , drop = FALSE]
  baseline <- far_abcd[far_abcd$candidate_id == "baseline_mutation_ll0p2_le1_so0p05", , drop = FALSE]
  if (!nrow(baseline)) baseline <- far_abcd[far_abcd$graph_edge_weight == "mutation" & far_abcd$lambda_l == 0.2 & far_abcd$lambda_e == 1, , drop = FALSE][1L, , drop = FALSE]
  e1_global <- edge$summary[edge$summary$graph_layer == "global", , drop = FALSE]
  local_delta <- e1_global[e1_global$estimator == "local_delta", , drop = FALSE]
  lctx <- e1_global[e1_global$estimator == "local_context_delta", , drop = FALSE]
  ectx <- e1_global[e1_global$estimator == "empirical_context_delta", , drop = FALSE]
  worst_context <- edge$by_context[edge$by_context$graph_layer == "global" &
                                     edge$by_context$estimator %in% c("local_context_delta", "empirical_context_delta"), , drop = FALSE]
  worst_context <- worst_context[order(worst_context$delta_sign_agreement_nonzero, -worst_context$n_scored_edges), , drop = FALSE]
  far_neg <- neg[neg$support_scope == "farfield" & neg$metric_scale == "native", , drop = FALSE]
  best_neg <- if (nrow(far_neg)) far_neg[order(-far_neg$spearman, far_neg$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_negative_scale <- far_neg[far_neg$prior_mean_scale < 0, , drop = FALSE]
  best_negative_scale <- if (nrow(best_negative_scale)) best_negative_scale[order(-best_negative_scale$spearman, best_negative_scale$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_positive_scale <- far_neg[far_neg$prior_mean_scale > 0, , drop = FALSE]
  best_positive_scale <- if (nrow(best_positive_scale)) best_positive_scale[order(-best_positive_scale$spearman, best_positive_scale$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  scale0_context <- far_neg[far_neg$prior_mean_mode != "zero" & far_neg$prior_mean_scale == 0, , drop = FALSE]
  scale0_context <- if (nrow(scale0_context)) scale0_context[order(scale0_context$shape_score, scale0_context$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  zero_mode <- far_neg[far_neg$prior_mean_mode == "zero", , drop = FALSE]
  zero_mode <- if (nrow(zero_mode)) zero_mode[order(zero_mode$shape_score, zero_mode$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  far_oracle <- oracle[oracle$support_scope == "farfield" & oracle$metric_scale == "native", , drop = FALSE]
  oracle_best_by_control <- if (nrow(far_oracle)) {
    do.call(rbind, lapply(split(far_oracle, far_oracle$oracle_control), function(x) {
      x[order(-x$spearman, x$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
    }))
  } else data.frame()
  far_pot <- potential[potential$support_scope == "farfield" & potential$metric_scale == "native", , drop = FALSE]
  best_pot <- if (nrow(far_pot)) far_pot[order(-far_pot$spearman, far_pot$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_pot_est <- far_pot[!grepl("oracle_", far_pot$prior_mean_mode), , drop = FALSE]
  best_pot_est <- if (nrow(best_pot_est)) best_pot_est[order(-best_pot_est$spearman, best_pot_est$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  cv_selected <- cv$selected
  cv_classes <- cv$classification
  valid_oracle <- if (nrow(cv_classes)) sum(cv_classes$shape_class == "valid_shape" & grepl("oracle_", cv_classes$prior_mean_mode), na.rm = TRUE) else 0L
  valid_nonoracle <- if (nrow(cv_classes)) sum(cv_classes$shape_class == "valid_shape" & !grepl("oracle_", cv_classes$prior_mean_mode), na.rm = TRUE) else 0L
  lines <- c(
    "# Farfield Shape Diagnostics Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", ctx$source_probe_dir %||% ctx$shared_input_dir, "`"),
    paste0("- abcd-dir: `", abcd$root, "`"),
    paste0("- simulation_id: ", task_info$simulation_id),
    paste0("- minobs: ", task_info$minobs),
    paste0("- input_policy: ", task_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    paste0("- quick mode: ", quick, " (auto limits the expensive E3 local TMB grid and E4 candidate fanout)."),
    "",
    "## Existing ABCD summary",
    paste0("- Baseline farfield/native: centered_rmse=", fmt_metric(baseline$centered_rmse[1]),
           ", pearson=", fmt_metric(baseline$pearson[1]),
           ", spearman=", fmt_metric(baseline$spearman[1]),
           ", estimate_sd_ratio=", fmt_metric(baseline$estimate_sd_ratio[1]), "."),
    "- ABCD showed normalized/unit can improve numeric RMSE and sometimes Spearman, but many such configs collapse amplitude. Constant anchor variance helps amplitude but not rank. Edge-delta prior means recover amplitude but tend to worsen direction.",
    "",
    "## E1 edge-delta direction diagnostics",
    paste0("- global local_delta vs truth_delta: sign_agreement=", fmt_metric(local_delta$delta_sign_agreement[1]),
           ", pearson=", fmt_metric(local_delta$delta_pearson[1]),
           ", spearman=", fmt_metric(local_delta$delta_spearman[1]), "."),
    paste0("- global local_context_delta vs truth_delta: sign_agreement=", fmt_metric(lctx$delta_sign_agreement[1]),
           ", pearson=", fmt_metric(lctx$delta_pearson[1]),
           ", spearman=", fmt_metric(lctx$delta_spearman[1]), "."),
    paste0("- global empirical_context_delta vs truth_delta: sign_agreement=", fmt_metric(ectx$delta_sign_agreement[1]),
           ", pearson=", fmt_metric(ectx$delta_pearson[1]),
           ", spearman=", fmt_metric(ectx$delta_spearman[1]), "."),
    paste0("- worst contexts by nonzero sign agreement: ",
           paste(utils::head(paste0(worst_context$context_label, "=", fmt_metric(worst_context$delta_sign_agreement_nonzero)), 5L), collapse = "; "), "."),
    "- Interpretation: edge-gradient pseudo-observation is unsafe unless these direction metrics are clearly better than chance and positive in correlation.",
    "",
    "## E1b negative-scale prior mean sanity check",
    paste0("- best farfield/native Spearman config: ", best_neg$candidate_id[1],
           " with prior_mean_mode=", best_neg$prior_mean_mode[1],
           ", scale=", best_neg$prior_mean_scale[1],
           ", pearson=", fmt_metric(best_neg$pearson[1]),
           ", spearman=", fmt_metric(best_neg$spearman[1]),
           ", sd_ratio=", fmt_metric(best_neg$estimate_sd_ratio[1]), "."),
    paste0("- best negative-scale config: ", best_negative_scale$candidate_id[1],
           " with pearson=", fmt_metric(best_negative_scale$pearson[1]),
           ", spearman=", fmt_metric(best_negative_scale$spearman[1]),
           ", sd_ratio=", fmt_metric(best_negative_scale$estimate_sd_ratio[1]), "."),
    paste0("- best positive-scale config: ", best_positive_scale$candidate_id[1],
           " with pearson=", fmt_metric(best_positive_scale$pearson[1]),
           ", spearman=", fmt_metric(best_positive_scale$spearman[1]),
           ", sd_ratio=", fmt_metric(best_positive_scale$estimate_sd_ratio[1]), "."),
    paste0("- scale=0 context-prior config is not equivalent to zero prior in this wrapper: best scale=0 context shape_score=",
           fmt_metric(scale0_context$shape_score[1]), " versus zero-mode shape_score=", fmt_metric(zero_mode$shape_score[1]),
           ". This happens because context-prior builders still seed and propagate direct-anchor means when delta_scale=0."),
    "- Negative scale does not rescue farfield rank here; the safest interpretation is unreliable/low-amplitude estimated delta rather than a simple global sign flip.",
    "",
    "## E2 oracle controls",
    paste0("- best O1/O5-style oracle result by control: ",
           paste(utils::head(paste0(oracle_best_by_control$oracle_control, ": sp=", fmt_metric(oracle_best_by_control$spearman),
                                    ", sd=", fmt_metric(oracle_best_by_control$estimate_sd_ratio)), 5L), collapse = "; "), "."),
    "- O1 tests architecture with local anchors and oracle per-edge deltas. O2/O3 isolate truth anchors with estimated or oracle context deltas. O4 tests truth anchors with zero-mean smoothing. O5 is the strongest oracle.",
    "- Oracle results are diagnostic only and were not used for production configuration selection.",
    "",
    "## E3 local stability / convergence",
    paste0("- local grid rows run: ", nrow(local$grid), ". Median gradient_norm=", fmt_metric(stats::median(local$grid$gradient_norm, na.rm = TRUE)), "."),
    paste0("- covariance statuses observed: ", paste(unique(local$grid$covariance_status), collapse = ", "), "."),
    paste0("- pairwise restart rows: ", nrow(local$pairwise), ". Mean fitness correlation median=", fmt_metric(stats::median(local$pairwise$fitness_mean_correlation, na.rm = TRUE)), "."),
    "- Auto mode is a sampled convergence probe, not the full 90-fit grid. It is enough to check whether the original nonconvergence is trivially fixed by small eval/model changes.",
    "",
    "## E4 potential-fit prior mean",
    paste0("- best overall potential prior: ", best_pot$candidate_id[1],
           " with mode=", best_pot$prior_mean_mode[1],
           ", scale=", best_pot$prior_mean_scale[1],
           ", pearson=", fmt_metric(best_pot$pearson[1]),
           ", spearman=", fmt_metric(best_pot$spearman[1]),
           ", sd_ratio=", fmt_metric(best_pot$estimate_sd_ratio[1]), "."),
    paste0("- best non-oracle estimated potential prior: ", best_pot_est$candidate_id[1],
           " with mode=", best_pot_est$prior_mean_mode[1],
           ", scale=", best_pot_est$prior_mean_scale[1],
           ", pearson=", fmt_metric(best_pot_est$pearson[1]),
           ", spearman=", fmt_metric(best_pot_est$spearman[1]),
           ", sd_ratio=", fmt_metric(best_pot_est$estimate_sd_ratio[1]), "."),
    "- If oracle potential helps but estimated potential does not, the next work is the delta estimator rather than C++ mechanics. If oracle also fails, graph support/shell depth is the bottleneck.",
    "",
    "## CV follow-up diagnostic",
    paste0("- selected objectives: ", paste(cv_selected$objective, cv_selected$candidate_id, sep = "=", collapse = "; "), "."),
    paste0("- shape classes: ", paste(names(table(cv_classes$shape_class)), as.integer(table(cv_classes$shape_class)), sep = "=", collapse = "; "), "."),
    paste0("- valid_shape configs are oracle-only in this run: valid_nonoracle=", valid_nonoracle, ", valid_oracle=", valid_oracle, "."),
    "- Direct-anchor CV remains limited: amplitude penalties avoid collapse, but direct anchors may not expose farfield wrong-direction errors. Neighbor-holdout, shell-stratified CV, or edge-delta CV should be tested next.",
    "",
    "## Default recommendation check",
    "- Keep benchmark/probe/calibration default graph_edge_weight as `normalized`, because mutation edge scale is not portable across graph support, but normalized must be paired with amplitude-collapse diagnostics.",
    "- Keep `unit` as a synthetic stress-test, not the default.",
    "- Do not default `anchor_count_reference=minobs` for full input; keep it as a candidate.",
    "- Calibration ranking should include amplitude-collapse as a failure state and may return no acceptable shape configuration.",
    "",
    "## Final conclusion",
    "- Current evidence does not support implementing C++ edge-gradient pseudo-observation yet. The blocker is not the lack of a C++ mechanism; it is unreliable delta direction and local fit/covariance diagnostics.",
    "- Highest priority is to fix local convergence/covariance and improve delta estimation/orientation checks. Only after estimated edge deltas are positively aligned with GRF truth in this diagnostic should a production edge-gradient prior be added.",
    "",
    "## Required outputs",
    "- `tables/all_diagnostics_long.tsv`",
    "- `tables/diagnostic_summary.tsv`",
    "- `tables/recommended_next_steps.tsv`",
    "- `tables/config_shape_classification.tsv`",
    "- `results/farfield_shape_diagnostics_all_results.rds`"
  )
  writeLines(lines, file.path(dirs$root, "farfield_shape_diagnostics_report.md"))
}

summarize_diagnostics <- function(dirs, ctx, task_info, abcd, quick) {
  edge <- readRDS(file.path(dirs$results, "edge_delta_diagnostics.rds"))
  neg <- read_tsv_safe(file.path(dirs$tables, "negative_scale_prior_mean_probe.tsv"))
  oracle <- read_tsv_safe(file.path(dirs$tables, "oracle_controls_by_scope.tsv"))
  local <- readRDS(file.path(dirs$results, "local_stability_probe.rds"))
  potential <- read_tsv_safe(file.path(dirs$tables, "potential_prior_mean_probe.tsv"))
  cv <- readRDS(file.path(dirs$results, "cv_diagnostic.rds"))
  all_long <- bind_rows_fill(list(
    transform(edge$summary, diagnostic = "E1_edge_delta"),
    transform(neg, diagnostic = "E1b_negative_scale"),
    transform(oracle, diagnostic = "E2_oracle"),
    transform(local$grid, diagnostic = "E3_local_stability"),
    transform(potential, diagnostic = "E4_potential_prior"),
    transform(cv$classification, diagnostic = "CV_shape_classification")
  ))
  write_tsv_safe(all_long, file.path(dirs$tables, "all_diagnostics_long.tsv"))
  summary_tbl <- make_diagnostic_summary(edge$summary, neg, oracle, local, potential, cv)
  next_steps <- make_recommended_next_steps(edge$summary, oracle, potential, local)
  write_tsv_safe(summary_tbl, file.path(dirs$tables, "diagnostic_summary.tsv"))
  write_tsv_safe(next_steps, file.path(dirs$tables, "recommended_next_steps.tsv"))
  saveRDS(list(
    context = ctx,
    task_info = task_info,
    abcd_summary = abcd,
    edge_delta = edge,
    negative_scale = neg,
    oracle = oracle,
    local_stability = local,
    potential = potential,
    cv = cv,
    all_diagnostics_long = all_long,
    diagnostic_summary = summary_tbl,
    recommended_next_steps = next_steps,
    quick_mode = quick
  ), file.path(dirs$results, "farfield_shape_diagnostics_all_results.rds"))
  write_diag_report(dirs, ctx, task_info, abcd, edge, neg, oracle, local, potential, cv,
                    summary_tbl, next_steps, quick)
  invisible(list(all = all_long, summary = summary_tbl, next_steps = next_steps))
}

main_diagnostics <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  mode <- match.arg(mode, c("prepare", "edge-delta", "oracle", "local-stability",
                            "potential-prior", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_shape_diagnostics"))
  simulation_id <- arg_integer(args, "simulation_id", 1L)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  quick <- as_bool_or_auto(arg_value(args, "quick", "auto"))
  cv_splits <- arg_integer(args, "cv_splits", 8L)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_diag_dirs(output_dir)
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
  saveRDS(list(context = ctx, task_info = task_info, local_diagnostics = local_diagnostics_table(bundle),
               abcd_dir = normalizePath(abcd_dir, winslash = "/", mustWork = TRUE)),
          file.path(dirs$results, "prepare_context.rds"))
  write_tsv_safe(local_diagnostics_table(bundle), file.path(dirs$tables, "local_fit_diagnostics.tsv"))
  if (identical(mode, "prepare")) return(invisible(dirs$root))

  if (mode %in% c("all", "edge-delta")) {
    edge <- run_edge_delta_diagnostics(bundle, grf, task_info, dirs)
    run_negative_scale_probe(bundle, components, grf, task_info, dirs, abcd, force = force)
    if (identical(mode, "edge-delta")) return(invisible(edge))
  }
  if (mode %in% c("all", "oracle")) {
    run_oracle_controls(bundle, components, grf, task_info, dirs, abcd, force = force)
    if (identical(mode, "oracle")) return(invisible(dirs$root))
  }
  if (mode %in% c("all", "local-stability")) {
    run_local_stability_grid(bundle, grf, task_info, dirs, quick = quick, force = force)
    if (identical(mode, "local-stability")) return(invisible(dirs$root))
  }
  potential <- data.frame()
  if (mode %in% c("all", "potential-prior")) {
    potential <- run_potential_prior_probe(bundle, components, grf, task_info, dirs, abcd,
                                           force = force, quick = quick)
    run_cv_diagnostic(bundle, components, task_info, dirs, abcd, potential, grf,
                      cv_splits = cv_splits, quick = quick)
    if (identical(mode, "potential-prior")) return(invisible(dirs$root))
  }
  if (mode %in% c("all", "summarize")) {
    summarize_diagnostics(dirs, ctx, task_info, abcd, quick)
  }
  message("Wrote farfield shape diagnostics under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) {
  main_diagnostics()
}
