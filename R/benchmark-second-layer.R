second_layer_rmse <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((pred[ok] - truth[ok])^2))
}

second_layer_mae <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  mean(abs(pred[ok] - truth[ok]))
}

second_layer_centered_rmse <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  pred <- pred[ok]
  truth <- truth[ok]
  sqrt(mean(((pred - mean(pred)) - (truth - mean(truth)))^2))
}

second_layer_affine_rmse <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L) return(NA_real_)
  pred <- pred[ok]
  truth <- truth[ok]
  if (!is.finite(stats::var(pred)) || stats::var(pred) <= 0) return(second_layer_rmse(pred, truth))
  fit <- stats::lm(truth ~ pred)
  sqrt(mean(stats::residuals(fit)^2))
}

second_layer_calibration_coef <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L || !is.finite(stats::var(pred[ok])) || stats::var(pred[ok]) <= 0) {
    return(c(intercept = NA_real_, slope = NA_real_))
  }
  cf <- stats::coef(stats::lm(truth[ok] ~ pred[ok]))
  c(intercept = unname(cf[[1]]), slope = unname(cf[[2]]))
}

second_layer_weighted_rmse <- function(pred, truth, weight) {
  ok <- is.finite(pred) & is.finite(truth) & is.finite(weight) & weight > 0
  if (!any(ok)) return(NA_real_)
  sqrt(stats::weighted.mean((pred[ok] - truth[ok])^2, weight[ok]))
}

second_layer_weighted_mae <- function(pred, truth, weight) {
  ok <- is.finite(pred) & is.finite(truth) & is.finite(weight) & weight > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(abs(pred[ok] - truth[ok]), weight[ok])
}

second_layer_cor <- function(pred, truth, method) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L) return(NA_real_)
  if (stats::sd(pred[ok]) <= 0 || stats::sd(truth[ok]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(pred[ok], truth[ok], method = method))
}

second_layer_sign_accuracy <- function(pred_gradient, truth_gradient, beneficial_only = FALSE, deleterious_only = FALSE) {
  ok <- is.finite(pred_gradient) & is.finite(truth_gradient)
  if (isTRUE(beneficial_only)) ok <- ok & truth_gradient > 0
  if (isTRUE(deleterious_only)) ok <- ok & truth_gradient < 0
  if (!any(ok)) return(NA_real_)
  mean(sign(pred_gradient[ok]) == sign(truth_gradient[ok]))
}

second_layer_top_k_overlap <- function(pred, truth, k = NULL) {
  ok <- is.finite(pred) & is.finite(truth)
  n <- sum(ok)
  if (!n) return(c(overlap_count = NA_real_, overlap_fraction = NA_real_, k = NA_real_))
  if (is.null(k)) k <- min(10L, max(1L, floor(0.1 * n)))
  k <- min(as.integer(k), n)
  idx <- which(ok)
  pred_top <- idx[order(pred[idx], decreasing = TRUE)][seq_len(k)]
  truth_top <- idx[order(truth[idx], decreasing = TRUE)][seq_len(k)]
  overlap <- length(intersect(pred_top, truth_top))
  c(overlap_count = overlap, overlap_fraction = overlap / k, k = k)
}

second_layer_shell_rows <- function(eval_nodes, shell) {
  d <- as.integer(eval_nodes$support_distance)
  switch(
    shell,
    all_nearfield = which(d %in% c(1L, 2L)),
    d0 = which(d == 0L),
    d1 = which(d == 1L),
    d2 = which(d == 2L),
    integer(0)
  )
}

second_layer_metric_values <- function(eval_nodes,
                                       eval_edges,
                                       shell = c("all_nearfield", "d0", "d1", "d2"),
                                       runtime_seconds = NA_real_,
                                       failure_status = "ok") {
  shell <- match.arg(shell)
  rows <- second_layer_shell_rows(eval_nodes, shell)
  df <- eval_nodes[rows, , drop = FALSE]
  ok <- is.finite(df$pred) & is.finite(df$truth)
  err <- df$pred - df$truth
  cal <- second_layer_calibration_coef(df$pred, df$truth)
  truth_sd <- stats::sd(df$truth[is.finite(df$truth)])
  weights <- if ("eval_weight" %in% names(df)) as.numeric(df$eval_weight) else rep(NA_real_, nrow(df))
  metrics <- c(
    rmse = second_layer_rmse(df$pred, df$truth),
    mae = second_layer_mae(df$pred, df$truth),
    bias = if (any(ok)) mean(err[ok]) else NA_real_,
    bias_abs = if (any(ok)) abs(mean(err[ok])) else NA_real_,
    median_absolute_error = if (any(ok)) stats::median(abs(err[ok])) else NA_real_,
    q90_absolute_error = if (any(ok)) stats::quantile(abs(err[ok]), 0.9, names = FALSE) else NA_real_,
    relative_rmse = if (is.finite(truth_sd) && truth_sd > 0) second_layer_rmse(df$pred, df$truth) / truth_sd else NA_real_,
    uncalibrated_r2 = if (any(ok) && sum((df$truth[ok] - mean(df$truth[ok]))^2) > 0) {
      1 - sum((df$pred[ok] - df$truth[ok])^2) / sum((df$truth[ok] - mean(df$truth[ok]))^2)
    } else NA_real_,
    calibration_intercept = cal[["intercept"]],
    calibration_slope = cal[["slope"]],
    centered_rmse = second_layer_centered_rmse(df$pred, df$truth),
    affine_rmse = second_layer_affine_rmse(df$pred, df$truth),
    count_weighted_rmse = second_layer_weighted_rmse(df$pred, df$truth, weights),
    count_weighted_mae = second_layer_weighted_mae(df$pred, df$truth, weights),
    pearson = second_layer_cor(df$pred, df$truth, "pearson"),
    spearman = second_layer_cor(df$pred, df$truth, "spearman")
  )
  if (!is.null(eval_edges) && nrow(eval_edges)) {
    edge_rows <- switch(
      shell,
      all_nearfield = which(eval_edges$child_distance %in% c(1L, 2L)),
      d0 = integer(0),
      d1 = which(eval_edges$child_distance == 1L),
      d2 = which(eval_edges$child_distance == 2L)
    )
    ed <- eval_edges[edge_rows, , drop = FALSE]
    edge_metrics <- c(
      edge_gradient_rmse = second_layer_rmse(ed$pred_gradient, ed$truth_gradient),
      edge_gradient_spearman = second_layer_cor(ed$pred_gradient, ed$truth_gradient, "spearman"),
      sign_accuracy = second_layer_sign_accuracy(ed$pred_gradient, ed$truth_gradient),
      beneficial_sign_accuracy = second_layer_sign_accuracy(ed$pred_gradient, ed$truth_gradient, beneficial_only = TRUE),
      deleterious_sign_accuracy = second_layer_sign_accuracy(ed$pred_gradient, ed$truth_gradient, deleterious_only = TRUE)
    )
  } else {
    edge_metrics <- c(
      edge_gradient_rmse = NA_real_,
      edge_gradient_spearman = NA_real_,
      sign_accuracy = NA_real_,
      beneficial_sign_accuracy = NA_real_,
      deleterious_sign_accuracy = NA_real_
    )
  }
  top <- if (identical(shell, "all_nearfield")) {
    second_layer_top_k_overlap(df$pred, df$truth)
  } else {
    c(overlap_count = NA_real_, overlap_fraction = NA_real_, k = NA_real_)
  }
  sd_ok <- is.finite(df$pred_sd) & df$pred_sd > 0 & is.finite(df$pred) & is.finite(df$truth)
  uncertainty <- if (any(sd_ok)) {
    c(
      interval_coverage_95 = mean(df$truth[sd_ok] >= df$pred[sd_ok] - 1.959963984540054 * df$pred_sd[sd_ok] &
                                   df$truth[sd_ok] <= df$pred[sd_ok] + 1.959963984540054 * df$pred_sd[sd_ok]),
      mean_pred_sd = mean(df$pred_sd[sd_ok]),
      median_pred_sd = stats::median(df$pred_sd[sd_ok]),
      standardized_rmse = sqrt(mean(((df$pred[sd_ok] - df$truth[sd_ok]) / df$pred_sd[sd_ok])^2))
    )
  } else {
    c(interval_coverage_95 = NA_real_, mean_pred_sd = NA_real_, median_pred_sd = NA_real_, standardized_rmse = NA_real_)
  }
  uncertainty <- c(
    uncertainty,
    interval_coverage_95_closeness = if (is.finite(uncertainty[["interval_coverage_95"]])) abs(uncertainty[["interval_coverage_95"]] - 0.95) else NA_real_
  )
  c(
    metrics,
    edge_metrics,
    top_k_overlap_count = top[["overlap_count"]],
    top_k_overlap_fraction = top[["overlap_fraction"]],
    top_k = top[["k"]],
    coverage = if (nrow(df)) sum(ok) / nrow(df) else NA_real_,
    uncertainty,
    runtime_seconds = as.numeric(runtime_seconds),
    failure_rate = ifelse(identical(failure_status, "ok"), 0, 1),
    failure_status_numeric = ifelse(identical(failure_status, "ok"), 0, ifelse(identical(failure_status, "partial"), 0.5, 1)),
    n_eval_nodes = nrow(df),
    n_valid_predictions = sum(ok),
    n_missing_predictions = nrow(df) - sum(ok)
  )
}

second_layer_metric_table <- function(eval_nodes,
                                      eval_edges,
                                      shells = c("all_nearfield", "d0", "d1", "d2"),
                                      runtime_seconds = NA_real_,
                                      failure_status = "ok") {
  scales <- c("raw", "anchor_calibrated")
  rows <- list()
  ri <- 0L
  for (scale in scales) {
    nodes <- eval_nodes
    edges <- eval_edges
    if (identical(scale, "anchor_calibrated")) {
      nodes$pred <- if ("pred_anchor_calibrated" %in% names(nodes)) nodes$pred_anchor_calibrated else nodes$pred
      if (nrow(edges)) {
        from_idx <- match(edges$parent_karyotype, nodes$karyotype)
        to_idx <- match(edges$child_karyotype, nodes$karyotype)
        edges$pred_gradient <- nodes$pred[to_idx] - nodes$pred[from_idx]
      }
    }
    for (shell in shells) {
      vals <- second_layer_metric_values(
        eval_nodes = nodes,
        eval_edges = edges,
        shell = shell,
        runtime_seconds = runtime_seconds,
        failure_status = failure_status
      )
      n_eval <- vals[["n_eval_nodes"]]
      n_valid <- vals[["n_valid_predictions"]]
      n_missing <- vals[["n_missing_predictions"]]
      vals <- vals[setdiff(names(vals), c("n_eval_nodes", "n_valid_predictions", "n_missing_predictions"))]
      ri <- ri + 1L
      rows[[ri]] <- data.frame(
        shell = shell,
        prediction_scale = scale,
        metric = names(vals),
        value = as.numeric(vals),
        n_eval_nodes = as.integer(n_eval),
        n_valid_predictions = as.integer(n_valid),
        n_missing_predictions = as.integer(n_missing),
        calibration_status = if (identical(scale, "anchor_calibrated")) "identity_or_internal" else "raw",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

second_layer_canonical_eval_graph <- function(counts,
                                              landscape,
                                              dt = 1,
                                              beta = 0.00005,
                                              transition_kernel = c("exact", "linear"),
                                              min_cn = NULL,
                                              max_cn = NULL,
                                              max_nodes = 5000) {
  transition_kernel <- match.arg(transition_kernel)
  if (is.null(min_cn)) min_cn <- landscape$min_cn
  if (is.null(max_cn)) max_cn <- landscape$max_cn
  data <- prepare_alfak2_data(counts, dt = dt, beta = beta)
  graph <- build_karyotype_graph(
    data,
    beta = beta,
    transition_kernel = transition_kernel,
    shell_depth = 2,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  truth <- predict_landscape_fitness(landscape, graph$karyotypes)
  observed_weight <- rowSums(data$counts)
  observed_weight <- observed_weight[match(as.character(graph$labels), names(observed_weight))]
  nodes <- data.frame(
    node_id = seq_along(graph$labels),
    karyotype = as.character(graph$labels),
    support_distance = as.integer(graph$support_distance),
    support_tier = as.character(graph$support_tier),
    truth = as.numeric(truth),
    eval_weight = as.numeric(observed_weight),
    stringsAsFactors = FALSE
  )
  edge_keep <- as.integer(graph$support_distance[graph$edge_to]) == as.integer(graph$support_distance[graph$edge_from]) + 1L &
    as.integer(graph$support_distance[graph$edge_to]) <= 2L
  edges <- data.frame(
    from = as.integer(graph$edge_from[edge_keep]),
    to = as.integer(graph$edge_to[edge_keep]),
    parent_karyotype = as.character(graph$labels[graph$edge_from[edge_keep]]),
    child_karyotype = as.character(graph$labels[graph$edge_to[edge_keep]]),
    parent_distance = as.integer(graph$support_distance[graph$edge_from[edge_keep]]),
    child_distance = as.integer(graph$support_distance[graph$edge_to[edge_keep]]),
    stringsAsFactors = FALSE
  )
  list(graph = graph, nodes = nodes, edges = edges)
}

second_layer_attach_predictions <- function(eval_graph, predictions) {
  nodes <- eval_graph$nodes
  pred_idx <- match(nodes$karyotype, as.character(predictions$karyotype))
  nodes$pred <- NA_real_
  nodes$pred_sd <- NA_real_
  nodes$pred_anchor_calibrated <- NA_real_
  nodes$prediction_status <- "missing"
  if (!"eval_weight" %in% names(nodes)) nodes$eval_weight <- NA_real_
  ok <- !is.na(pred_idx)
  if (any(ok)) {
    nodes$pred[ok] <- as.numeric(predictions$fitness_mean[pred_idx[ok]])
    nodes$pred_anchor_calibrated[ok] <- nodes$pred[ok]
    if ("fitness_sd" %in% names(predictions)) {
      nodes$pred_sd[ok] <- as.numeric(predictions$fitness_sd[pred_idx[ok]])
    }
    if ("prediction_status" %in% names(predictions)) {
      nodes$prediction_status[ok] <- as.character(predictions$prediction_status[pred_idx[ok]])
    }
  }
  edges <- eval_graph$edges
  if (nrow(edges)) {
    from_idx <- match(edges$parent_karyotype, nodes$karyotype)
    to_idx <- match(edges$child_karyotype, nodes$karyotype)
    edges$pred_gradient <- nodes$pred[to_idx] - nodes$pred[from_idx]
    edges$truth_gradient <- nodes$truth[to_idx] - nodes$truth[from_idx]
  } else {
    edges$pred_gradient <- numeric(0)
    edges$truth_gradient <- numeric(0)
  }
  list(nodes = nodes, edges = edges)
}

second_layer_alfak2_methods <- function() {
  alfak2_extrapolation_methods()
}

second_layer_alfakR_slots <- function() {
  out <- data.frame(
    NN_prior_slot = c(
      "none",
      "empirical",
      "empirical_censored",
      "empirical_censored_weighted_slot4",
      "empirical_censored_weighted_slot5"
    ),
    NN_prior = c(
      "none",
      "empirical",
      "empirical_censored",
      "empirical_censored_weighted",
      "empirical_censored_weighted"
    ),
    stringsAsFactors = FALSE
  )
  out$NN_prior_value <- out$NN_prior
  out
}

second_layer_build_run_index <- function(mode = c("full", "quick")) {
  mode <- match.arg(mode)
  lambdas <- if (identical(mode, "quick")) 0.6 else c(0.2, 0.6, 0.8)
  reps <- if (identical(mode, "quick")) 1L else 1:5
  landscape_grid <- expand.grid(
    lambda_index = seq_along(lambdas),
    landscape_rep = reps,
    stringsAsFactors = FALSE
  )
  landscape_grid$grf_lambda <- lambdas[landscape_grid$lambda_index]
  landscape_grid$landscape_id <- sprintf("lambda%s_rep%s", gsub("\\.", "p", as.character(landscape_grid$grf_lambda)), landscape_grid$landscape_rep)
  alfak2 <- merge(
    landscape_grid,
    expand.grid(
      package = "alfak2",
      input_mode = c("full", "minobs_matched", "soft_minobs"),
      extrapolation_method = second_layer_alfak2_methods(),
      stringsAsFactors = FALSE
    ),
    all = TRUE
  )
  alfak2$minobs <- NA_integer_
  alfak2$NN_prior_slot <- NA_character_
  alfak2$NN_prior <- NA_character_
  alfak2$NN_prior_value <- NA_character_
  slots <- second_layer_alfakR_slots()
  alfakR_modes <- expand.grid(
    package = "alfakR",
    minobs = if (identical(mode, "quick")) 5L else c(5L, 10L, 20L),
    slot_row = seq_len(nrow(slots)),
    stringsAsFactors = FALSE
  )
  alfakR_modes$NN_prior_slot <- slots$NN_prior_slot[alfakR_modes$slot_row]
  alfakR_modes$NN_prior <- slots$NN_prior[alfakR_modes$slot_row]
  alfakR_modes$NN_prior_value <- slots$NN_prior_value[alfakR_modes$slot_row]
  alfakR_modes$slot_row <- NULL
  alfakR <- merge(landscape_grid, alfakR_modes, all = TRUE)
  alfakR$input_mode <- NA_character_
  alfakR$extrapolation_method <- NA_character_
  out <- rbind(
    alfak2[, c("package", "grf_lambda", "lambda_index", "landscape_id", "landscape_rep", "input_mode", "extrapolation_method", "minobs", "NN_prior_slot", "NN_prior", "NN_prior_value")],
    alfakR[, c("package", "grf_lambda", "lambda_index", "landscape_id", "landscape_rep", "input_mode", "extrapolation_method", "minobs", "NN_prior_slot", "NN_prior", "NN_prior_value")]
  )
  out$run_id <- sprintf("run_%04d", seq_len(nrow(out)))
  out$landscape_seed <- 100000L + 1000L * as.integer(out$lambda_index) + as.integer(out$landscape_rep)
  out$count_seed <- 200000L + 1000L * as.integer(out$lambda_index) + as.integer(out$landscape_rep)
  out[, c("run_id", setdiff(names(out), "run_id"))]
}
