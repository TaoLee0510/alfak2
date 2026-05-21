alfak2_extrapolation_methods <- function() {
  c(
    "graph_gaussian_baseline",
    "edge_effect_empirical_bayes",
    "edge_effect_interaction_path_ensemble",
    "kronecker_or_graph_trend_filtering",
    "local_NNGP_or_GPnn",
    "delta_tree_ensemble",
    "tabpfn_nearfield_feature_model",
    "truncated_nearfield_gmrf",
    "local_polynomial_stencil"
  )
}

match_extrapolation_method <- function(method) {
  aliases <- c(
    graph_gaussian = "graph_gaussian_baseline",
    nearfield_local_edge_effect = "edge_effect_empirical_bayes",
    edge_interaction_path_ensemble = "edge_effect_interaction_path_ensemble",
    graph_trend_filtering = "kronecker_or_graph_trend_filtering"
  )
  if (length(method) > 1L) method <- method[[1L]]
  method <- match.arg(method, c(alfak2_extrapolation_methods(), names(aliases)))
  if (method %in% names(aliases)) unname(aliases[[method]]) else method
}

validate_prediction_distance <- function(max_prediction_distance) {
  validate_scalar(as.numeric(max_prediction_distance), "max_prediction_distance", lower = 0)
  as.integer(max_prediction_distance)
}

extrapolation_empty_summary <- function(graph, method, max_prediction_distance) {
  tier <- as.character(graph$support_tier)
  out_of_scope <- as.integer(graph$support_distance) > max_prediction_distance
  tier[out_of_scope] <- "out_of_scope"
  data.frame(
    node_id = seq_along(graph$labels),
    karyotype = as.character(graph$labels),
    support_tier = tier,
    support_distance = as.integer(graph$support_distance),
    fitness_mean = NA_real_,
    fitness_sd = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    extrapolation_method = method,
    prediction_status = ifelse(out_of_scope, "out_of_scope", "missing"),
    stringsAsFactors = FALSE
  )
}

extrapolation_add_intervals <- function(summary) {
  summary$conf_low <- as.numeric(summary$fitness_mean) - 1.959963984540054 * as.numeric(summary$fitness_sd)
  summary$conf_high <- as.numeric(summary$fitness_mean) + 1.959963984540054 * as.numeric(summary$fitness_sd)
  bad <- !is.finite(summary$fitness_mean) | !is.finite(summary$fitness_sd)
  summary$conf_low[bad] <- NA_real_
  summary$conf_high[bad] <- NA_real_
  summary
}

support_tier_anchor_weight <- function(tier) {
  tier <- as.character(tier)
  out <- rep(0.5, length(tier))
  out[tier == "directly_informed"] <- 1
  out[tier == "local_borrowed"] <- 0.6
  out[tier == "weakly_supported"] <- 0.3
  out
}

prepare_extrapolation_anchors <- function(local_fit,
                                          graph,
                                          anchor_support_tiers = "all",
                                          anchor_covariance_inflation = c(
                                            TMB_sdreport = 1,
                                            untrusted_gradient = 4,
                                            untrusted_nonconverged = 9,
                                            untrusted_sdreport_missing = 4,
                                            untrusted_sdreport_nonfinite = 4,
                                            unknown = 4
                                          ),
                                          anchor_count_reference = NULL,
                                          anchor_count_power = 1,
                                          anchor_min_effective_count = 0,
                                          anchor_exclude = character()) {
  if (!inherits(local_fit, "alfak2_local_fit")) {
    stop("`local_fit` must be an alfak2_local_fit object.", call. = FALSE)
  }
  if (!inherits(graph, "alfak2_graph")) stop("`graph` must be an alfak2_graph object.", call. = FALSE)
  anchor_tiers <- resolve_anchor_tiers(anchor_support_tiers)
  if (!is.null(anchor_min_effective_count)) {
    validate_scalar(as.numeric(anchor_min_effective_count), "anchor_min_effective_count", lower = 0)
    anchor_min_effective_count <- as.numeric(anchor_min_effective_count)
  }
  anchor_match <- match(local_fit$summary$karyotype, as.character(graph$labels))
  tier_ok <- rep(TRUE, nrow(local_fit$summary))
  if (!is.null(anchor_tiers)) {
    tier_ok <- as.character(local_fit$summary$support_tier) %in% anchor_tiers
  }
  if (!is.null(anchor_exclude) && length(anchor_exclude)) {
    tier_ok <- tier_ok & !(as.character(local_fit$summary$karyotype) %in% as.character(anchor_exclude))
  }
  anchor_count_all <- if ("effective_count_total" %in% names(local_fit$summary)) {
    as.numeric(local_fit$summary$effective_count_total)
  } else if ("count_total" %in% names(local_fit$summary)) {
    as.numeric(local_fit$summary$count_total)
  } else {
    rep(NA_real_, nrow(local_fit$summary))
  }
  count_ok <- rep(TRUE, nrow(local_fit$summary))
  if (!is.null(anchor_min_effective_count)) {
    if (!any(is.finite(anchor_count_all))) {
      alfak2_abort(
        "`anchor_min_effective_count` requires finite local effective counts.",
        diagnostics = list(stage = "extrapolation_anchor_filter")
      )
    }
    count_ok <- is.finite(anchor_count_all) & anchor_count_all > anchor_min_effective_count
  }
  keep <- which(!is.na(anchor_match) & is.finite(local_fit$summary$fitness_mean) & tier_ok & count_ok)
  if (!length(keep)) {
    alfak2_abort(
      "No local posterior anchors are present in the extrapolation graph.",
      diagnostics = list(
        stage = "extrapolation_anchor_filter",
        n_local = nrow(local_fit$summary),
        n_in_graph = sum(!is.na(anchor_match)),
        anchor_min_effective_count = anchor_min_effective_count
      )
    )
  }
  anchor_var_base <- as.numeric(local_fit$summary$fitness_sd[keep])^2
  if (any(!is.finite(anchor_var_base) | anchor_var_base <= 0)) {
    alfak2_abort(
      "Extrapolation anchors contain non-finite local fitness variances.",
      diagnostics = list(stage = "extrapolation_anchor_variance", bad_indices = keep[!is.finite(anchor_var_base) | anchor_var_base <= 0])
    )
  }
  if (!"covariance_status" %in% names(local_fit$summary)) {
    alfak2_abort(
      "Local fit summary is missing covariance status.",
      diagnostics = list(stage = "extrapolation_anchor_covariance_status")
    )
  }
  covariance_status <- as.character(local_fit$summary$covariance_status[keep])
  covariance_mult <- anchor_covariance_multiplier(covariance_status, anchor_covariance_inflation)
  anchor_count_for_weight <- anchor_count_all[keep]
  count_mult <- if (any(is.finite(anchor_count_for_weight))) {
    count_anchor_multiplier(anchor_count_for_weight, anchor_count_reference, anchor_count_power)
  } else {
    rep(1, length(keep))
  }
  anchor_var_multiplier <- covariance_mult * count_mult
  anchor_var <- anchor_var_base * anchor_var_multiplier
  data.frame(
    node_id = as.integer(anchor_match[keep]),
    karyotype = as.character(local_fit$summary$karyotype[keep]),
    mean = as.numeric(local_fit$summary$fitness_mean[keep]),
    variance_base = anchor_var_base,
    variance = anchor_var,
    variance_multiplier = anchor_var_multiplier,
    covariance_status = covariance_status,
    support_tier = as.character(local_fit$summary$support_tier[keep]),
    support_distance = as.integer(local_fit$summary$support_distance[keep]),
    count_total = if ("count_total" %in% names(local_fit$summary)) local_fit$summary$count_total[keep] else NA_integer_,
    effective_count_total = if ("effective_count_total" %in% names(local_fit$summary)) local_fit$summary$effective_count_total[keep] else NA_real_,
    anchor_count_for_weight = anchor_count_for_weight,
    stringsAsFactors = FALSE
  )
}

graph_edge_table <- function(graph, graph_edge_weight = c("mutation", "unit", "normalized")) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  data.frame(
    edge_id = seq_along(graph$edge_from),
    from = as.integer(graph$edge_from),
    to = as.integer(graph$edge_to),
    chr = as.integer(graph$edge_chr),
    direction = as.integer(graph$edge_direction),
    weight = as.numeric(graph_edge_weights(graph$edge_weight, graph_edge_weight)),
    stringsAsFactors = FALSE
  )
}

local_summary_map <- function(local_fit, graph) {
  idx <- match(as.character(graph$labels), as.character(local_fit$summary$karyotype))
  idx
}

fill_direct_anchor_predictions <- function(summary, anchors, status = "anchor") {
  idx <- anchors$node_id
  summary$fitness_mean[idx] <- anchors$mean
  summary$fitness_sd[idx] <- sqrt(pmax(anchors$variance, .Machine$double.eps))
  summary$prediction_status[idx] <- status
  extrapolation_add_intervals(summary)
}

fit_edge_effects <- function(graph,
                             anchors,
                             graph_edge_weight = c("mutation", "unit", "normalized"),
                             ridge_delta = 4,
                             ridge_gamma = 25,
                             max_anchor_distance = 1) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  validate_scalar(as.numeric(ridge_delta), "ridge_delta", lower = 0)
  validate_scalar(as.numeric(ridge_gamma), "ridge_gamma", lower = 0)
  edge <- graph_edge_table(graph, graph_edge_weight)
  anchor_idx <- anchors$node_id
  anchor_mean <- stats::setNames(anchors$mean, anchor_idx)
  anchor_var <- stats::setNames(anchors$variance, anchor_idx)
  anchor_ok <- anchors$support_distance <= max_anchor_distance
  anchor_nodes <- anchors$node_id[anchor_ok]
  edge$type <- paste(edge$chr, edge$direction, sep = ":")
  train <- edge[edge$from %in% anchor_nodes & edge$to %in% anchor_nodes, , drop = FALSE]
  warnings <- character()
  if (nrow(train)) {
    train$diff <- as.numeric(anchor_mean[as.character(train$to)] - anchor_mean[as.character(train$from)])
    train$var <- as.numeric(anchor_var[as.character(train$to)] + anchor_var[as.character(train$from)])
    train$w <- pmax(train$weight, 0) / pmax(train$var, .Machine$double.eps)
  }
  type_levels <- unique(edge$type)
  delta <- stats::setNames(rep(0, length(type_levels)), type_levels)
  delta_var <- stats::setNames(rep(1, length(type_levels)), type_levels)
  delta_n <- stats::setNames(integer(length(type_levels)), type_levels)
  if (nrow(train)) {
    for (tp in type_levels) {
      rows <- train[train$type == tp & is.finite(train$diff) & is.finite(train$w) & train$w > 0, , drop = FALSE]
      if (!nrow(rows)) next
      sw <- sum(rows$w)
      delta[[tp]] <- sum(rows$w * rows$diff) / (sw + ridge_delta)
      resid <- rows$diff - delta[[tp]]
      resid_var <- if (nrow(rows) >= 2L) stats::var(resid) else stats::weighted.mean(rows$var, rows$w)
      if (!is.finite(resid_var) || resid_var <= 0) resid_var <- 1 / (sw + ridge_delta)
      delta_var[[tp]] <- resid_var / max(1, nrow(rows)) + 1 / (sw + ridge_delta + .Machine$double.eps)
      delta_n[[tp]] <- nrow(rows)
    }
  } else {
    warnings <- c(warnings, "No anchor-to-anchor edges were available; edge deltas shrank to zero.")
  }

  gamma <- stats::setNames(numeric(0), character(0))
  gamma_var <- stats::setNames(numeric(0), character(0))
  gamma_n <- stats::setNames(integer(0), character(0))
  if (length(anchor_nodes) >= 3L) {
    e1 <- edge[edge$from %in% anchor_nodes & edge$to %in% anchor_nodes, , drop = FALSE]
    if (nrow(e1)) {
      rows <- list()
      ri <- 0L
      for (i in seq_len(nrow(e1))) {
        e2 <- e1[e1$from == e1$to[[i]] & e1$to %in% anchor_nodes, , drop = FALSE]
        if (!nrow(e2)) next
        for (j in seq_len(nrow(e2))) {
          if (identical(e1$from[[i]], e2$to[[j]])) next
          key <- paste(e1$type[[i]], e2$type[[j]], sep = "|")
          residual <- as.numeric(anchor_mean[as.character(e2$to[[j]])] - anchor_mean[as.character(e1$from[[i]])]) -
            delta[[e1$type[[i]]]] - delta[[e2$type[[j]]]]
          v <- as.numeric(anchor_var[as.character(e1$from[[i]])] + anchor_var[as.character(e1$to[[i]])] + anchor_var[as.character(e2$to[[j]])])
          ri <- ri + 1L
          rows[[ri]] <- data.frame(key = key, residual = residual, variance = v, stringsAsFactors = FALSE)
        }
      }
      if (length(rows)) {
        gtrain <- do.call(rbind, rows)
        for (key in unique(gtrain$key)) {
          rows <- gtrain[gtrain$key == key & is.finite(gtrain$residual), , drop = FALSE]
          w <- 1 / pmax(rows$variance, .Machine$double.eps)
          sw <- sum(w)
          gamma[[key]] <- sum(w * rows$residual) / (sw + ridge_gamma)
          rv <- if (nrow(rows) >= 2L) stats::var(rows$residual - gamma[[key]]) else stats::weighted.mean(rows$variance, w)
          if (!is.finite(rv) || rv <= 0) rv <- 1 / (sw + ridge_gamma + .Machine$double.eps)
          gamma_var[[key]] <- rv / max(1, nrow(rows)) + 1 / (sw + ridge_gamma + .Machine$double.eps)
          gamma_n[[key]] <- nrow(rows)
        }
      }
    }
  }
  list(
    delta = delta,
    delta_var = delta_var,
    delta_n = delta_n,
    gamma = gamma,
    gamma_var = gamma_var,
    gamma_n = gamma_n,
    training_edges = train,
    warnings = warnings
  )
}

combine_path_predictions <- function(mean, variance, path_weight, eps = 1e-8) {
  ok <- is.finite(mean) & is.finite(variance) & variance > 0 & is.finite(path_weight) & path_weight > 0
  if (!any(ok)) return(c(mean = NA_real_, variance = NA_real_, n_paths = 0))
  mean <- mean[ok]
  variance <- variance[ok]
  path_weight <- path_weight[ok]
  w <- path_weight / pmax(variance, eps)
  w <- w / sum(w)
  mu <- sum(w * mean)
  within <- sum(w * variance)
  between <- sum(w * (mean - mu)^2)
  c(mean = mu, variance = within + between, n_paths = length(mean))
}

path_ensemble_predict <- function(graph,
                                  anchors,
                                  edge_effects,
                                  graph_edge_weight = c("mutation", "unit", "normalized"),
                                  max_prediction_distance = 2,
                                  eps = 1e-8) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  edge <- graph_edge_table(graph, graph_edge_weight)
  edge$type <- paste(edge$chr, edge$direction, sep = ":")
  support_distance <- as.integer(graph$support_distance)
  anchor_mean <- stats::setNames(anchors$mean, anchors$node_id)
  anchor_var <- stats::setNames(anchors$variance, anchors$node_id)
  anchor_nodes <- anchors$node_id[anchors$support_distance == 0L]
  out_mean <- rep(NA_real_, length(graph$labels))
  out_var <- rep(NA_real_, length(graph$labels))
  n_paths <- integer(length(graph$labels))
  out_mean[anchors$node_id] <- anchors$mean
  out_var[anchors$node_id] <- anchors$variance
  n_paths[anchors$node_id] <- 1L
  for (target in which(support_distance > 0L & support_distance <= max_prediction_distance)) {
    if (support_distance[[target]] == 1L) {
      e1 <- edge[edge$to == target & edge$from %in% anchor_nodes, , drop = FALSE]
      if (!nrow(e1)) next
      pm <- pv <- pw <- numeric(nrow(e1))
      for (i in seq_len(nrow(e1))) {
        tp <- e1$type[[i]]
        from <- as.character(e1$from[[i]])
        pm[[i]] <- anchor_mean[[from]] + edge_effects$delta[[tp]]
        pv[[i]] <- anchor_var[[from]] + edge_effects$delta_var[[tp]]
        pw[[i]] <- pmax(e1$weight[[i]], 0)
      }
      cmb <- combine_path_predictions(pm, pv, pw, eps = eps)
    } else {
      last <- edge[edge$to == target & support_distance[edge$from] == 1L, , drop = FALSE]
      rows <- list()
      ri <- 0L
      if (nrow(last)) {
        for (i in seq_len(nrow(last))) {
          first <- edge[edge$to == last$from[[i]] & edge$from %in% anchor_nodes, , drop = FALSE]
          if (!nrow(first)) next
          for (j in seq_len(nrow(first))) {
            ri <- ri + 1L
            rows[[ri]] <- list(first = first[j, , drop = FALSE], last = last[i, , drop = FALSE])
          }
        }
      }
      if (!length(rows)) next
      pm <- pv <- pw <- numeric(length(rows))
      for (i in seq_along(rows)) {
        e0 <- rows[[i]]$first
        e1 <- rows[[i]]$last
        tp0 <- e0$type[[1]]
        tp1 <- e1$type[[1]]
        gkey <- paste(tp0, tp1, sep = "|")
        gamma <- if (gkey %in% names(edge_effects$gamma)) edge_effects$gamma[[gkey]] else 0
        gvar <- if (gkey %in% names(edge_effects$gamma_var)) edge_effects$gamma_var[[gkey]] else 1 / 25
        from <- as.character(e0$from[[1]])
        pm[[i]] <- anchor_mean[[from]] + edge_effects$delta[[tp0]] + edge_effects$delta[[tp1]] + gamma
        pv[[i]] <- anchor_var[[from]] + edge_effects$delta_var[[tp0]] + edge_effects$delta_var[[tp1]] + gvar
        pw[[i]] <- pmax(e0$weight[[1]], 0) * pmax(e1$weight[[1]], 0)
      }
      cmb <- combine_path_predictions(pm, pv, pw, eps = eps)
    }
    out_mean[[target]] <- cmb[["mean"]]
    out_var[[target]] <- cmb[["variance"]]
    n_paths[[target]] <- as.integer(cmb[["n_paths"]])
  }
  list(mean = out_mean, variance = out_var, n_paths = n_paths)
}

copy_local_summary_predictions <- function(summary, local_fit, graph, max_prediction_distance) {
  local_idx <- local_summary_map(local_fit, graph)
  ok <- which(!is.na(local_idx) &
                as.integer(graph$support_distance) <= max_prediction_distance &
                is.finite(local_fit$summary$fitness_mean[local_idx]))
  if (length(ok)) {
    summary$fitness_mean[ok] <- as.numeric(local_fit$summary$fitness_mean[local_idx[ok]])
    summary$fitness_sd[ok] <- as.numeric(local_fit$summary$fitness_sd[local_idx[ok]])
    summary$prediction_status[ok] <- ifelse(summary$support_distance[ok] == 0L, "anchor", "local_posterior")
  }
  extrapolation_add_intervals(summary)
}

build_extrapolation_diagnostics <- function(method,
                                            max_prediction_distance,
                                            anchors,
                                            summary,
                                            hyperparameters,
                                            tuning_grid,
                                            convergence_status,
                                            runtime_seconds,
                                            dependency_status = "ok",
                                            warnings = character(),
                                            extra = list()) {
  valid <- is.finite(summary$fitness_mean) & as.character(summary$prediction_status) != "out_of_scope"
  out_of_scope <- as.character(summary$prediction_status) == "out_of_scope"
  c(
    list(
      extrapolation_method = method,
      max_prediction_distance = max_prediction_distance,
      n_anchors = nrow(anchors),
      n_predicted = sum(valid),
      n_out_of_scope = sum(out_of_scope),
      hyperparameters = hyperparameters,
      tuning_grid = tuning_grid,
      convergence_status = convergence_status,
      runtime_seconds = as.numeric(runtime_seconds),
      warnings = as.character(warnings),
      dependency_status = as.character(dependency_status)
    ),
    extra
  )
}

finish_extrapolation_fit <- function(graph,
                                     summary,
                                     anchors,
                                     hyperparameters,
                                     tuning_grid,
                                     diagnostics) {
  new_alfak2_global_fit(list(
    graph = graph,
    summary = summary,
    anchors = anchors,
    hyperparameters = hyperparameters,
    tuning_grid = tuning_grid,
    diagnostics = diagnostics
  ))
}

nearfield_local_edge_effect_fit <- function(local_fit,
                                            graph,
                                            max_prediction_distance,
                                            anchors,
                                            graph_edge_weight = c("mutation", "unit", "normalized"),
                                            ridge_delta = 4,
                                            ...) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  started <- proc.time()[["elapsed"]]
  method <- "edge_effect_empirical_bayes"
  warnings <- character()
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  summary <- copy_local_summary_predictions(summary, local_fit, graph, max_prediction_distance)
  missing <- which(as.integer(graph$support_distance) <= max_prediction_distance &
                     !is.finite(summary$fitness_mean))
  effects <- fit_edge_effects(
    graph = graph,
    anchors = anchors,
    graph_edge_weight = graph_edge_weight,
    ridge_delta = ridge_delta,
    ridge_gamma = 1e6,
    max_anchor_distance = 1
  )
  warnings <- c(warnings, effects$warnings)
  if (length(missing)) {
    paths <- path_ensemble_predict(
      graph = graph,
      anchors = anchors,
      edge_effects = effects,
      graph_edge_weight = graph_edge_weight,
      max_prediction_distance = max_prediction_distance
    )
    ok <- missing[is.finite(paths$mean[missing])]
    summary$fitness_mean[ok] <- paths$mean[ok]
    summary$fitness_sd[ok] <- sqrt(pmax(paths$variance[ok], .Machine$double.eps))
    summary$prediction_status[ok] <- "local_edge_effect"
  }
  summary <- extrapolation_add_intervals(summary)
  hp <- list(graph_edge_weight = graph_edge_weight, ridge_delta = ridge_delta)
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp,
    tuning_grid = data.frame(parameter = "ridge_delta", value = ridge_delta),
    convergence_status = "ok",
    runtime_seconds = proc.time()[["elapsed"]] - started,
    warnings = warnings,
    extra = list(
      delta_bias = if (nrow(effects$training_edges)) mean(effects$training_edges$diff, na.rm = TRUE) else NA_real_,
      anchor_cv_rmse = NA_real_,
      anchor_cv_mae = NA_real_,
      anchor_cv_calibration_slope = NA_real_,
      anchor_cv_calibration_intercept = NA_real_,
      n_delta_contexts = length(effects$delta),
      delta_shrinkage = ridge_delta,
      interaction_status = "not_used"
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, diagnostics$tuning_grid, diagnostics)
}

edge_interaction_path_ensemble_fit <- function(local_fit,
                                               graph,
                                               max_prediction_distance,
                                               anchors,
                                               graph_edge_weight = c("mutation", "unit", "normalized"),
                                               ridge_delta = 4,
                                               ridge_gamma = 25,
                                               ...) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  started <- proc.time()[["elapsed"]]
  method <- "edge_effect_interaction_path_ensemble"
  effects <- fit_edge_effects(
    graph = graph,
    anchors = anchors,
    graph_edge_weight = graph_edge_weight,
    ridge_delta = ridge_delta,
    ridge_gamma = ridge_gamma,
    max_anchor_distance = 1
  )
  pred <- path_ensemble_predict(
    graph = graph,
    anchors = anchors,
    edge_effects = effects,
    graph_edge_weight = graph_edge_weight,
    max_prediction_distance = max_prediction_distance
  )
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  idx <- which(as.integer(graph$support_distance) <= max_prediction_distance & is.finite(pred$mean))
  summary$fitness_mean[idx] <- pred$mean[idx]
  summary$fitness_sd[idx] <- sqrt(pmax(pred$variance[idx], .Machine$double.eps))
  summary$prediction_status[idx] <- ifelse(summary$support_distance[idx] == 0L, "anchor", "path_ensemble")
  summary <- extrapolation_add_intervals(summary)
  hp <- list(
    graph_edge_weight = graph_edge_weight,
    ridge_delta = ridge_delta,
    ridge_gamma = ridge_gamma
  )
  grid <- rbind(
    data.frame(parameter = "ridge_delta", value = ridge_delta, stringsAsFactors = FALSE),
    data.frame(parameter = "ridge_gamma", value = ridge_gamma, stringsAsFactors = FALSE)
  )
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp, grid, "ok",
    proc.time()[["elapsed"]] - started,
    warnings = effects$warnings,
    extra = list(
      n_paths = sum(pred$n_paths, na.rm = TRUE),
      mean_paths_per_target = mean(pred$n_paths[as.integer(graph$support_distance) %in% c(1L, 2L)], na.rm = TRUE),
      max_paths_per_target = max(pred$n_paths[as.integer(graph$support_distance) %in% c(1L, 2L)], na.rm = TRUE),
      delta_shrinkage = ridge_delta,
      gamma_shrinkage = ridge_gamma,
      n_delta_contexts = length(effects$delta),
      n_gamma_contexts = length(effects$gamma),
      anchor_cv_rmse = NA_real_,
      anchor_cv_mae = NA_real_,
      anchor_cv_edge_gradient_rmse = NA_real_,
      anchor_cv_calibration_slope = NA_real_,
      anchor_cv_calibration_intercept = NA_real_,
      edge_delta_n = effects$delta_n,
      edge_interaction_n = effects$gamma_n,
      path_count_summary = summary(pred$n_paths)
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, grid, diagnostics)
}

anchor_lattice_distance <- function(target, anchor_mat) {
  rowSums(abs(sweep(anchor_mat, 2L, target, FUN = "-")))
}

fit_local_polynomial_target <- function(target,
                                        target_id,
                                        graph,
                                        anchors,
                                        ridge = 1,
                                        bandwidth = 2,
                                        min_train = 3,
                                        max_pairwise_terms = 20) {
  x0 <- as.numeric(graph$karyotypes[target_id, ])
  amat <- graph$karyotypes[anchors$node_id, , drop = FALSE]
  d <- anchor_lattice_distance(x0, amat)
  train <- which(is.finite(d) & d <= max(2, bandwidth * 2))
  if (length(train) < min_train) train <- seq_len(nrow(anchors))
  if (!length(train)) return(c(mean = NA_real_, variance = NA_real_, n_train = 0))
  dx <- sweep(amat[train, , drop = FALSE], 2L, x0, FUN = "-")
  p <- ncol(dx)
  X <- cbind("(Intercept)" = 1, dx, 0.5 * dx^2)
  colnames(X) <- c("(Intercept)", paste0("lin", seq_len(p)), paste0("quad", seq_len(p)))
  varying <- which(colSums(abs(dx)) > 0)
  pair_terms <- NULL
  if (length(train) >= ncol(X) + 5L && length(varying) >= 2L) {
    pairs <- utils::combn(varying, 2L)
    if (ncol(pairs) > max_pairwise_terms) pairs <- pairs[, seq_len(max_pairwise_terms), drop = FALSE]
    pair_terms <- vapply(seq_len(ncol(pairs)), function(i) dx[, pairs[1, i]] * dx[, pairs[2, i]], numeric(nrow(dx)))
    if (!is.matrix(pair_terms)) pair_terms <- matrix(pair_terms, ncol = 1L)
    colnames(pair_terms) <- paste0("pair", seq_len(ncol(pair_terms)))
    X <- cbind(X, pair_terms)
  }
  y <- anchors$mean[train]
  var <- pmax(anchors$variance[train], .Machine$double.eps)
  tier_w <- support_tier_anchor_weight(anchors$support_tier[train])
  w <- exp(-d[train] / max(bandwidth, .Machine$double.eps)) * tier_w / var
  w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0], 1)
  sw <- sqrt(w)
  Xw <- X * sw
  yw <- y * sw
  penalty <- diag(ridge, ncol(Xw))
  penalty[1, 1] <- 0
  beta <- tryCatch(solve(crossprod(Xw) + penalty, crossprod(Xw, yw)), error = function(e) NULL)
  if (is.null(beta)) return(c(mean = stats::weighted.mean(y, w), variance = stats::weighted.mean(var, w), n_train = length(train)))
  residual <- as.numeric(y - X %*% beta)
  rv <- if (length(residual) >= 2L) stats::weighted.mean(residual^2, w) else stats::weighted.mean(var, w)
  if (!is.finite(rv) || rv <= 0) rv <- stats::weighted.mean(var, w)
  leverage_var <- tryCatch(solve(crossprod(Xw) + penalty)[1, 1], error = function(e) NA_real_)
  pred_var <- rv * if (is.finite(leverage_var) && leverage_var > 0) leverage_var else 1 / max(1, sum(w))
  c(mean = as.numeric(beta[1]), variance = pmax(pred_var, .Machine$double.eps), n_train = length(train))
}

local_polynomial_stencil_fit <- function(local_fit,
                                         graph,
                                         max_prediction_distance,
                                         anchors,
                                         ridge = 1,
                                         bandwidth = 2,
                                         min_train = 3,
                                         max_pairwise_terms = 20,
                                         ...) {
  started <- proc.time()[["elapsed"]]
  method <- "local_polynomial_stencil"
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  summary <- fill_direct_anchor_predictions(summary, anchors, status = "anchor")
  scoped <- which(as.integer(graph$support_distance) > 0L &
                    as.integer(graph$support_distance) <= max_prediction_distance)
  n_train <- integer(length(graph$labels))
  for (target in scoped) {
    fit <- fit_local_polynomial_target(
      target = graph$karyotypes[target, ],
      target_id = target,
      graph = graph,
      anchors = anchors,
      ridge = ridge,
      bandwidth = bandwidth,
      min_train = min_train,
      max_pairwise_terms = max_pairwise_terms
    )
    if (is.finite(fit[["mean"]])) {
      summary$fitness_mean[target] <- fit[["mean"]]
      summary$fitness_sd[target] <- sqrt(pmax(fit[["variance"]], .Machine$double.eps))
      summary$prediction_status[target] <- "local_polynomial"
      n_train[target] <- as.integer(fit[["n_train"]])
    }
  }
  summary <- extrapolation_add_intervals(summary)
  hp <- list(ridge = ridge, bandwidth = bandwidth, min_train = min_train, max_pairwise_terms = max_pairwise_terms)
  grid <- data.frame(parameter = names(hp), value = unlist(hp, use.names = FALSE), stringsAsFactors = FALSE)
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp, grid, "ok",
    proc.time()[["elapsed"]] - started,
    extra = list(
      mean_training_neighbors = mean(n_train[n_train > 0], na.rm = TRUE),
      model_order = if (max_pairwise_terms > 0) "quadratic_with_limited_pairwise" else "quadratic_diagonal",
      ridge_lambda = ridge,
      pairwise_terms_used = max_pairwise_terms > 0,
      anchor_cv_rmse = NA_real_,
      anchor_cv_mae = NA_real_,
      fallback_count = 0L,
      n_train_summary = summary(n_train[n_train > 0])
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, grid, diagnostics)
}

local_gp_kernel <- function(x, y, lengthscale = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  out <- matrix(0, nrow(x), nrow(y))
  for (i in seq_len(nrow(x))) {
    out[i, ] <- exp(-rowSums(abs(sweep(y, 2L, x[i, ], FUN = "-"))) / max(lengthscale, .Machine$double.eps))
  }
  out
}

local_gp_predict_one <- function(target_x,
                                 anchors,
                                 graph,
                                 k = 30,
                                 lengthscale = 1,
                                 nugget = 0.05,
                                 jitter = 1e-8) {
  amat <- graph$karyotypes[anchors$node_id, , drop = FALSE]
  d <- anchor_lattice_distance(as.numeric(target_x), amat)
  keep <- order(d, anchors$variance)[seq_len(min(k, nrow(anchors)))]
  x <- amat[keep, , drop = FALSE]
  y <- anchors$mean[keep]
  avar <- pmax(anchors$variance[keep], .Machine$double.eps)
  mu0 <- stats::weighted.mean(y, 1 / avar)
  if (!is.finite(mu0)) mu0 <- mean(y)
  K <- local_gp_kernel(x, x, lengthscale = lengthscale)
  diag(K) <- diag(K) + nugget + avar + jitter
  kt <- as.numeric(local_gp_kernel(matrix(target_x, nrow = 1L), x, lengthscale = lengthscale))
  chol_status <- "ok"
  pred <- tryCatch({
    R <- chol(K)
    alpha <- backsolve(R, forwardsolve(t(R), y - mu0))
    v <- forwardsolve(t(R), kt)
    mean <- mu0 + sum(kt * alpha)
    variance <- pmax(1 + nugget - sum(v^2), .Machine$double.eps)
    c(mean = mean, variance = variance, n_neighbors = length(keep), cholesky_failed = 0)
  }, error = function(e) {
    w <- exp(-d[keep] / max(lengthscale, .Machine$double.eps)) / avar
    if (!any(is.finite(w) & w > 0)) w[] <- 1
    m <- stats::weighted.mean(y, w)
    rv <- stats::weighted.mean((y - m)^2 + avar, w)
    c(mean = m, variance = pmax(rv, .Machine$double.eps), n_neighbors = length(keep), cholesky_failed = 1)
  })
  pred
}

local_gp_anchor_cv <- function(anchors, graph, lengthscale = 1, nugget = 0.05, k = 30) {
  if (nrow(anchors) < 3L) {
    return(list(rmse = NA_real_, mae = NA_real_, nll = NA_real_, slope = NA_real_, intercept = NA_real_))
  }
  pred <- rep(NA_real_, nrow(anchors))
  psd <- rep(NA_real_, nrow(anchors))
  for (i in seq_len(nrow(anchors))) {
    tr <- anchors[-i, , drop = FALSE]
    p <- local_gp_predict_one(
      graph$karyotypes[anchors$node_id[[i]], ],
      tr,
      graph,
      k = min(k, nrow(tr)),
      lengthscale = lengthscale,
      nugget = nugget
    )
    pred[[i]] <- p[["mean"]]
    psd[[i]] <- sqrt(p[["variance"]])
  }
  ok <- is.finite(pred) & is.finite(anchors$mean)
  cal <- if (sum(ok) >= 2L && stats::sd(pred[ok]) > 0) stats::coef(stats::lm(anchors$mean[ok] ~ pred[ok])) else c(`(Intercept)` = NA_real_, `pred[ok]` = NA_real_)
  nll <- if (any(ok & is.finite(psd) & psd > 0)) {
    -mean(stats::dnorm(anchors$mean[ok], pred[ok], psd[ok], log = TRUE), na.rm = TRUE)
  } else NA_real_
  list(
    rmse = sqrt(mean((pred[ok] - anchors$mean[ok])^2)),
    mae = mean(abs(pred[ok] - anchors$mean[ok])),
    nll = nll,
    slope = unname(cal[[2]]),
    intercept = unname(cal[[1]])
  )
}

local_NNGP_or_GPnn_fit <- function(local_fit,
                                   graph,
                                   max_prediction_distance,
                                   anchors,
                                   k = 30,
                                   lengthscale_grid = c(0.75, 1.5, 3),
                                   nugget_grid = c(0.02, 0.05),
                                   ...) {
  started <- proc.time()[["elapsed"]]
  method <- "local_NNGP_or_GPnn"
  k <- min(as.integer(k), nrow(anchors))
  grid <- expand.grid(lengthscale = as.numeric(lengthscale_grid), nugget = as.numeric(nugget_grid))
  cv <- lapply(seq_len(nrow(grid)), function(i) {
    local_gp_anchor_cv(anchors, graph, lengthscale = grid$lengthscale[[i]], nugget = grid$nugget[[i]], k = k)
  })
  cv_rmse <- vapply(cv, `[[`, numeric(1), "rmse")
  best <- if (any(is.finite(cv_rmse))) which.min(cv_rmse) else 1L
  lengthscale <- grid$lengthscale[[best]]
  nugget <- grid$nugget[[best]]
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  scoped <- which(as.integer(graph$support_distance) <= max_prediction_distance)
  neighbor_count <- integer(length(graph$labels))
  cholesky_failures <- 0L
  for (target in scoped) {
    if (target %in% anchors$node_id) {
      ai <- match(target, anchors$node_id)
      summary$fitness_mean[target] <- anchors$mean[[ai]]
      summary$fitness_sd[target] <- sqrt(pmax(anchors$variance[[ai]], .Machine$double.eps))
      summary$prediction_status[target] <- "anchor"
      neighbor_count[target] <- 1L
      next
    }
    p <- local_gp_predict_one(graph$karyotypes[target, ], anchors, graph, k = k, lengthscale = lengthscale, nugget = nugget)
    summary$fitness_mean[target] <- p[["mean"]]
    summary$fitness_sd[target] <- sqrt(pmax(p[["variance"]], .Machine$double.eps))
    summary$prediction_status[target] <- ifelse(p[["cholesky_failed"]] > 0, "ridge_local_fallback", "local_gp")
    neighbor_count[target] <- as.integer(p[["n_neighbors"]])
    cholesky_failures <- cholesky_failures + as.integer(p[["cholesky_failed"]])
  }
  summary <- extrapolation_add_intervals(summary)
  hp <- list(kernel = "exp_l1_copy_number", lengthscale = lengthscale, nugget = nugget, k = k)
  tuning_grid <- cbind(grid, anchor_cv_rmse = cv_rmse)
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp, tuning_grid, "ok",
    proc.time()[["elapsed"]] - started,
    extra = list(
      mean_neighbor_count = mean(neighbor_count[neighbor_count > 0], na.rm = TRUE),
      min_neighbor_count = min(neighbor_count[neighbor_count > 0], na.rm = TRUE),
      max_neighbor_count = max(neighbor_count[neighbor_count > 0], na.rm = TRUE),
      kernel = "exp_l1_copy_number",
      lengthscale = lengthscale,
      nugget = nugget,
      anchor_cv_rmse = cv[[best]]$rmse,
      anchor_cv_mae = cv[[best]]$mae,
      anchor_cv_nll = cv[[best]]$nll,
      calibration_slope = cv[[best]]$slope,
      calibration_intercept = cv[[best]]$intercept,
      n_cholesky_failures = cholesky_failures,
      fallback_count = cholesky_failures
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, tuning_grid, diagnostics)
}

delta_feature_vector <- function(parent_node, target_node, graph, parent_anchor, path_length = NULL, edge_weight = 1) {
  px <- as.numeric(graph$karyotypes[parent_node, ])
  tx <- as.numeric(graph$karyotypes[target_node, ])
  delta <- tx - px
  if (is.null(path_length)) path_length <- sum(abs(delta))
  c(
    parent_fitness_mean = parent_anchor$mean,
    parent_fitness_sd = sqrt(pmax(parent_anchor$variance, .Machine$double.eps)),
    parent_support_weight = support_tier_anchor_weight(parent_anchor$support_tier),
    parent_count_total = if (is.finite(parent_anchor$count_total)) parent_anchor$count_total else 0,
    parent_effective_count = if (is.finite(parent_anchor$effective_count_total)) parent_anchor$effective_count_total else 0,
    support_distance = as.numeric(graph$support_distance[target_node]),
    path_length = path_length,
    edge_weight = edge_weight,
    local_anchor_density = NA_real_,
    candidate_parent_count = NA_real_,
    setNames(px, paste0("parent_cn_", seq_along(px))),
    setNames(tx, paste0("target_cn_", seq_along(tx))),
    setNames(delta, paste0("delta_cn_", seq_along(delta))),
    setNames(abs(delta), paste0("abs_delta_cn_", seq_along(delta)))
  )
}

delta_training_table <- function(graph, anchors, max_path_length = 2) {
  rows <- list()
  y <- numeric(0)
  ri <- 0L
  for (i in seq_len(nrow(anchors))) {
    for (j in seq_len(nrow(anchors))) {
      if (i == j) next
      d <- sum(abs(graph$karyotypes[anchors$node_id[[j]], ] - graph$karyotypes[anchors$node_id[[i]], ]))
      if (!is.finite(d) || d < 1 || d > max_path_length) next
      ri <- ri + 1L
      rows[[ri]] <- delta_feature_vector(anchors$node_id[[i]], anchors$node_id[[j]], graph, anchors[i, , drop = FALSE], path_length = d)
      y[[ri]] <- anchors$mean[[j]] - anchors$mean[[i]]
    }
  }
  if (!length(rows)) return(list(x = matrix(numeric(0), nrow = 0L), y = numeric(0)))
  x <- do.call(rbind, rows)
  x[, !is.finite(colMeans(x, na.rm = TRUE))] <- 0
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- 0
  list(x = x, y = y)
}

ridge_delta_fit <- function(x, y, lambda = 1) {
  if (!nrow(x) || length(y) < 2L) {
    return(list(beta = NULL, center = NULL, scale = NULL, intercept = if (length(y)) mean(y) else 0, residual_sd = 1))
  }
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  xs <- sweep(sweep(x, 2L, center), 2L, scale, "/")
  X <- cbind(1, xs)
  pen <- diag(lambda, ncol(X))
  pen[1, 1] <- 0
  beta <- tryCatch(as.numeric(solve(crossprod(X) + pen, crossprod(X, y))), error = function(e) NULL)
  if (is.null(beta)) {
    return(list(beta = NULL, center = center, scale = scale, intercept = mean(y), residual_sd = stats::sd(y)))
  }
  pred <- as.numeric(X %*% beta)
  list(beta = beta, center = center, scale = scale, intercept = beta[[1]], residual_sd = stats::sd(y - pred))
}

ridge_delta_predict <- function(model, x) {
  if (is.null(model$beta)) return(rep(model$intercept, nrow(x)))
  xs <- sweep(sweep(x, 2L, model$center), 2L, model$scale, "/")
  as.numeric(cbind(1, xs) %*% model$beta)
}

delta_tabular_fit <- function(x, y, backend, ridge_lambda = 1) {
  if (identical(backend, "xgboost") && requireNamespace("xgboost", quietly = TRUE) && nrow(x) >= 2L && length(y) >= 2L) {
    fit <- tryCatch({
      dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
      params <- list(
        objective = "reg:squarederror",
        max_depth = 2,
        eta = 0.1,
        subsample = 1,
        colsample_bytree = 1,
        min_child_weight = 1,
        nthread = 1,
        verbosity = 0
      )
      model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = 40,
        verbose = 0
      )
      pred <- as.numeric(stats::predict(model, x))
      list(
        model_type = "xgboost",
        model = model,
        residual_sd = stats::sd(y - pred),
        best_params = c(params, list(nrounds = 40))
      )
    }, error = function(e) NULL)
    if (!is.null(fit)) return(fit)
  }
  ridge <- ridge_delta_fit(x, y, lambda = ridge_lambda)
  list(
    model_type = "ridge",
    model = ridge,
    residual_sd = ridge$residual_sd,
    best_params = list(ridge_lambda = ridge_lambda)
  )
}

delta_tabular_predict <- function(model, x) {
  if (identical(model$model_type, "xgboost")) {
    return(as.numeric(stats::predict(model$model, x)))
  }
  ridge_delta_predict(model$model, x)
}

delta_model_predict_targets <- function(graph, anchors, model, max_prediction_distance, method_status) {
  out_mean <- rep(NA_real_, length(graph$labels))
  out_var <- rep(NA_real_, length(graph$labels))
  status <- rep("missing", length(graph$labels))
  out_mean[anchors$node_id] <- anchors$mean
  out_var[anchors$node_id] <- anchors$variance
  status[anchors$node_id] <- "anchor"
  for (target in which(as.integer(graph$support_distance) > 0L & as.integer(graph$support_distance) <= max_prediction_distance)) {
    amat <- graph$karyotypes[anchors$node_id, , drop = FALSE]
    d <- anchor_lattice_distance(graph$karyotypes[target, ], amat)
    parents <- which(d >= 1 & d <= max_prediction_distance)
    if (!length(parents)) next
    x <- do.call(rbind, lapply(parents, function(i) delta_feature_vector(anchors$node_id[[i]], target, graph, anchors[i, , drop = FALSE], path_length = d[[i]])))
    for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- 0
    delta <- delta_tabular_predict(model, x)
    pm <- anchors$mean[parents] + delta
    pv <- anchors$variance[parents] + pmax(model$residual_sd, 0.01)^2
    pw <- 1 / pmax(pv, .Machine$double.eps) / pmax(d[parents], 1)
    cmb <- combine_path_predictions(pm, pv, pw)
    out_mean[target] <- cmb[["mean"]]
    out_var[target] <- cmb[["variance"]]
    status[target] <- method_status
  }
  list(mean = out_mean, variance = out_var, status = status)
}

tree_backend_status <- function() {
  candidates <- c("xgboost", "lightgbm", "catboost", "xbart", "bart", "grf", "ranger", "randomForest", "glmnet")
  available <- vapply(candidates, requireNamespace, logical(1), quietly = TRUE)
  backend <- if (available[["xgboost"]]) "xgboost" else if (available[["lightgbm"]]) "lightgbm" else if (available[["catboost"]]) "catboost" else if (available[["xbart"]]) "xbart" else if (available[["bart"]]) "bart" else if (available[["grf"]]) "grf" else if (available[["ranger"]]) "ranger" else if (available[["randomForest"]]) "randomForest" else "ridge_fallback"
  list(backend = backend, available = available)
}

delta_tree_ensemble_fit <- function(local_fit,
                                    graph,
                                    max_prediction_distance,
                                    anchors,
                                    graph_edge_weight = c("mutation", "unit", "normalized"),
                                    ridge_lambda = 1,
                                    bootstrap_replicates = 0L,
                                    force_ridge = FALSE,
                                    method_name = "delta_tree_ensemble",
                                    dependency_status = NULL,
                                    ...) {
  started <- proc.time()[["elapsed"]]
  graph_edge_weight <- match.arg(graph_edge_weight)
  method <- method_name
  train <- delta_training_table(graph, anchors, max_path_length = 2)
  backend_status <- tree_backend_status()
  backend <- if (isTRUE(force_ridge)) "ridge_fallback" else backend_status$backend
  model <- delta_tabular_fit(train$x, train$y, backend = backend, ridge_lambda = ridge_lambda)
  method_status <- if (identical(model$model_type, "xgboost")) "delta_tree_xgboost" else "delta_tree_ridge"
  pred <- delta_model_predict_targets(graph, anchors, model, max_prediction_distance, method_status = method_status)
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  idx <- which(as.integer(graph$support_distance) <= max_prediction_distance & is.finite(pred$mean))
  summary$fitness_mean[idx] <- pred$mean[idx]
  summary$fitness_sd[idx] <- sqrt(pmax(pred$variance[idx], .Machine$double.eps))
  summary$prediction_status[idx] <- pred$status[idx]
  summary <- extrapolation_add_intervals(summary)
  hp <- list(backend = backend, model_type = model$model_type, ridge_lambda = ridge_lambda, graph_edge_weight = graph_edge_weight)
  grid <- data.frame(parameter = names(hp), value = as.character(unlist(hp, use.names = FALSE)), stringsAsFactors = FALSE)
  cv_pred <- if (nrow(train$x)) delta_tabular_predict(model, train$x) else numeric(0)
  abs_pred <- cv_pred
  cv_rmse <- if (length(abs_pred)) sqrt(mean((abs_pred - train$y)^2)) else NA_real_
  cv_mae <- if (length(abs_pred)) mean(abs(abs_pred - train$y)) else NA_real_
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp, grid, "ok",
    proc.time()[["elapsed"]] - started,
    dependency_status = if (is.null(dependency_status)) {
      if (identical(model$model_type, "xgboost")) {
        "xgboost"
      } else if (identical(backend, "ridge_fallback")) {
        "tree_backend_unavailable_used_ridge_fallback"
      } else {
        paste0(backend, "_available_used_ridge_fallback")
      }
    } else {
      dependency_status
    },
    extra = list(
      backend = backend,
      backend_available = backend_status$available,
      n_training_edges = length(train$y),
      n_features = if (length(dim(train$x))) ncol(train$x) else 0L,
      cv_score = cv_rmse,
      anchor_cv_absolute_f_rmse = cv_rmse,
      anchor_cv_absolute_f_mae = cv_mae,
      best_params = model$best_params,
      bootstrap_replicates = as.integer(bootstrap_replicates),
      fallback_status = if (identical(model$model_type, "xgboost")) "none" else if (identical(backend, "ridge_fallback")) "ridge_fallback" else "ridge_fallback_after_backend_failure"
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, grid, diagnostics)
}

tabpfn_nearfield_feature_model_fit <- function(local_fit,
                                               graph,
                                               max_prediction_distance,
                                               anchors,
                                               graph_edge_weight = c("mutation", "unit", "normalized"),
                                               ...) {
  python_available <- requireNamespace("reticulate", quietly = TRUE)
  tabpfn_available <- FALSE
  tabpfn_version <- NA_character_
  if (python_available) {
    tabpfn_available <- tryCatch(reticulate::py_module_available("tabpfn"), error = function(e) FALSE)
    if (isTRUE(tabpfn_available)) {
      tabpfn_version <- tryCatch(as.character(reticulate::py_get_attr(reticulate::import("tabpfn"), "__version__")), error = function(e) NA_character_)
    }
  }
  dependency_status <- if (isTRUE(tabpfn_available)) "tabpfn_available_ridge_adapter" else "tabpfn_unavailable_used_tree_fallback"
  fit <- delta_tree_ensemble_fit(
    local_fit = local_fit,
    graph = graph,
    max_prediction_distance = max_prediction_distance,
    anchors = anchors,
    graph_edge_weight = graph_edge_weight,
    force_ridge = FALSE,
    method_name = "tabpfn_nearfield_feature_model",
    dependency_status = dependency_status,
    ...
  )
  fit$diagnostics$python_available <- python_available
  fit$diagnostics$tabpfn_available <- tabpfn_available
  fit$diagnostics$tabpfn_version <- tabpfn_version
  fit$diagnostics$fallback_used <- !isTRUE(tabpfn_available)
  if (!isTRUE(tabpfn_available)) {
    fit$diagnostics$fallback_status <- "tree_fallback_used"
  }
  fit
}

induce_graph_nodes <- function(graph, keep) {
  keep <- sort(unique(as.integer(keep)))
  map <- integer(length(graph$labels))
  map[keep] <- seq_along(keep)
  edge_keep <- graph$edge_from %in% keep & graph$edge_to %in% keep
  out <- list(
    labels = graph$labels[keep],
    karyotypes = graph$karyotypes[keep, , drop = FALSE],
    support_distance = graph$support_distance[keep],
    support_tier = graph$support_tier[keep],
    observed_index = map[graph$observed_index[graph$observed_index %in% keep]],
    edge_from = map[graph$edge_from[edge_keep]],
    edge_to = map[graph$edge_to[edge_keep]],
    edge_chr = graph$edge_chr[edge_keep],
    edge_direction = graph$edge_direction[edge_keep],
    edge_weight = graph$edge_weight[edge_keep],
    transition_from0 = integer(0),
    transition_to0 = integer(0),
    transition_weight = numeric(0),
    parent_from0 = integer(0),
    parent_to0 = integer(0),
    parent_weight = numeric(0),
    parent_context0 = integer(0),
    context_label = character(0),
    context_group0 = integer(0),
    n_chr = graph$n_chr,
    beta = graph$beta,
    transition_kernel = graph$transition_kernel,
    shell_depth = min(max(graph$support_distance[keep]), graph$shell_depth),
    min_cn = graph$min_cn,
    max_cn = graph$max_cn,
    original_node_id = keep
  )
  new_alfak2_graph(out)
}

dense_graph_penalty_matrix <- function(n, rows) {
  if (!length(rows)) return(matrix(0, nrow = 0L, ncol = n))
  out <- matrix(0, nrow = length(rows), ncol = n)
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    out[i, row$idx] <- row$coef
  }
  out
}

trend_filter_solve <- function(n,
                               anchor_sub_idx,
                               anchor_y,
                               anchor_var,
                               edge,
                               lambda1,
                               lambda2,
                               admm_rho,
                               max_iter,
                               tol) {
  A <- matrix(0, n, n)
  b <- numeric(n)
  aw <- 1 / pmax(anchor_var, .Machine$double.eps)
  for (i in seq_along(anchor_sub_idx)) {
    idx <- anchor_sub_idx[[i]]
    A[idx, idx] <- A[idx, idx] + aw[[i]]
    b[idx] <- b[idx] + aw[[i]] * anchor_y[[i]]
  }
  d1_rows <- lapply(seq_len(nrow(edge)), function(i) list(idx = c(edge$from[[i]], edge$to[[i]]), coef = c(-1, 1)))
  D1 <- dense_graph_penalty_matrix(n, d1_rows)
  d2_rows <- list()
  if (nrow(edge)) {
    ri <- 0L
    for (i in seq_len(nrow(edge))) {
      next_edges <- edge[edge$from == edge$to[[i]], , drop = FALSE]
      if (!nrow(next_edges)) next
      for (j in seq_len(nrow(next_edges))) {
        if (next_edges$to[[j]] == edge$from[[i]]) next
        ri <- ri + 1L
        d2_rows[[ri]] <- list(idx = c(edge$from[[i]], edge$to[[i]], next_edges$to[[j]]), coef = c(1, -2, 1))
      }
    }
  }
  D2 <- dense_graph_penalty_matrix(n, d2_rows)
  f <- rep(stats::weighted.mean(anchor_y, aw), n)
  if (!is.finite(f[[1]])) f[] <- 0
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    z1 <- if (nrow(D1)) as.numeric(D1 %*% f) else numeric(0)
    z2 <- if (nrow(D2)) as.numeric(D2 %*% f) else numeric(0)
    w1 <- if (length(z1)) lambda1 / pmax(abs(z1), tol) else numeric(0)
    w2 <- if (length(z2)) lambda2 / pmax(abs(z2), tol) else numeric(0)
    P <- A + diag(admm_rho, n)
    if (nrow(D1)) P <- P + crossprod(D1, D1 * w1)
    if (nrow(D2)) P <- P + crossprod(D2, D2 * w2)
    f_new <- tryCatch(as.numeric(solve(P, b)), error = function(e) rep(NA_real_, n))
    if (!all(is.finite(f_new))) break
    if (max(abs(f_new - f), na.rm = TRUE) <= tol) {
      f <- f_new
      converged <- TRUE
      break
    }
    f <- f_new
  }
  residual <- anchor_y - f[anchor_sub_idx]
  sigma2 <- if (length(residual) >= 2L) stats::weighted.mean(residual^2, aw) else stats::weighted.mean(anchor_var, aw)
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- stats::weighted.mean(anchor_var, aw)
  list(mean = f, sd = rep(sqrt(pmax(sigma2, .Machine$double.eps)), n), converged = converged, iterations = iter)
}

graph_trend_filtering_fit <- function(local_fit,
                                      graph,
                                      max_prediction_distance,
                                      anchors,
                                      graph_edge_weight = c("mutation", "unit", "normalized"),
                                      lambda1 = 0.25,
                                      lambda2 = 0.1,
                                      admm_rho = 1e-6,
                                      max_iter = 80,
                                      tol = 1e-4,
                                      max_dense_nodes = 1500,
                                      ...) {
  graph_edge_weight <- match.arg(graph_edge_weight)
  started <- proc.time()[["elapsed"]]
  method <- "kronecker_or_graph_trend_filtering"
  scoped <- which(as.integer(graph$support_distance) <= max_prediction_distance)
  subgraph <- induce_graph_nodes(graph, scoped)
  anchor_sub_idx <- match(anchors$node_id, scoped)
  keep_anchor <- !is.na(anchor_sub_idx)
  warnings <- character()
  if (length(scoped) > max_dense_nodes) {
    warnings <- c(warnings, sprintf("Induced nearfield graph has %d nodes; using path-ensemble deterministic fallback.", length(scoped)))
    fallback <- edge_interaction_path_ensemble_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight
    )
    fallback$summary$extrapolation_method <- method
    fallback$diagnostics$extrapolation_method <- method
    fallback$diagnostics$convergence_status <- "fallback_path_ensemble"
    fallback$diagnostics$warnings <- c(fallback$diagnostics$warnings, warnings)
    return(new_alfak2_global_fit(fallback))
  }
  edge <- graph_edge_table(subgraph, graph_edge_weight)
  sol <- trend_filter_solve(
    n = length(scoped),
    anchor_sub_idx = anchor_sub_idx[keep_anchor],
    anchor_y = anchors$mean[keep_anchor],
    anchor_var = anchors$variance[keep_anchor],
    edge = edge,
    lambda1 = lambda1,
    lambda2 = lambda2,
    admm_rho = admm_rho,
    max_iter = max_iter,
    tol = tol
  )
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  summary$fitness_mean[scoped] <- sol$mean
  summary$fitness_sd[scoped] <- sol$sd
  summary$prediction_status[scoped] <- ifelse(summary$support_distance[scoped] == 0L, "anchor_smoothed", "trend_filtered")
  summary <- extrapolation_add_intervals(summary)
  hp <- list(
    graph_edge_weight = graph_edge_weight,
    lambda1 = lambda1,
    lambda2 = lambda2,
    admm_rho = admm_rho,
    max_iter = max_iter,
    tol = tol
  )
  grid <- data.frame(parameter = names(hp), value = as.character(unlist(hp, use.names = FALSE)), stringsAsFactors = FALSE)
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, anchors, summary, hp, grid,
    convergence_status = if (isTRUE(sol$converged)) "converged" else "max_iter",
    runtime_seconds = proc.time()[["elapsed"]] - started,
    warnings = warnings,
    extra = list(
      trend_filtering_backend = "graph",
      lambda1 = lambda1,
      lambda2 = lambda2,
      admm_iterations = sol$iterations,
      primal_residual = NA_real_,
      dual_residual = NA_real_,
      anchor_cv_rmse = NA_real_,
      anchor_cv_edge_gradient_rmse = NA_real_,
      iterations = sol$iterations,
      induced_nodes = length(scoped),
      induced_edges = nrow(edge)
    )
  )
  finish_extrapolation_fit(graph, summary, anchors, hp, grid, diagnostics)
}

truncated_nearfield_gmrf_fit <- function(local_fit,
                                         graph,
                                         max_prediction_distance,
                                         anchors,
                                         lambda_l_grid = c(0.2, 1, 5),
                                         lambda_e_grid = c(0.05, 0.25, 1),
                                         sigma_obs_grid = c(0.02, 0.05, 0.1),
                                         eps = 1e-5,
                                         graph_edge_weight = c("mutation", "unit", "normalized"),
                                         anchor_support_tiers = "all",
                                         anchor_covariance_inflation = c(
                                           TMB_sdreport = 1,
                                           untrusted_gradient = 4,
                                           untrusted_nonconverged = 9,
                                           untrusted_sdreport_missing = 4,
                                           untrusted_sdreport_nonfinite = 4,
                                           unknown = 4
                                         ),
                                         anchor_count_reference = NULL,
                                         anchor_count_power = 1,
                                         anchor_min_effective_count = 0,
                                         anchor_exclude = character(),
                                         compute_sd = TRUE,
                                         ...) {
  started <- proc.time()[["elapsed"]]
  method <- "truncated_nearfield_gmrf"
  scoped <- which(as.integer(graph$support_distance) <= max_prediction_distance)
  subgraph <- induce_graph_nodes(graph, scoped)
  subfit <- fit_graph_posterior(
    local_fit = local_fit,
    graph = subgraph,
    lambda_l_grid = lambda_l_grid,
    lambda_e_grid = lambda_e_grid,
    sigma_obs_grid = sigma_obs_grid,
    eps = eps,
    graph_edge_weight = graph_edge_weight,
    anchor_support_tiers = anchor_support_tiers,
    anchor_covariance_inflation = anchor_covariance_inflation,
    anchor_count_reference = anchor_count_reference,
    anchor_count_power = anchor_count_power,
    anchor_min_effective_count = anchor_min_effective_count,
    anchor_exclude = anchor_exclude,
    compute_sd = compute_sd
  )
  summary <- extrapolation_empty_summary(graph, method, max_prediction_distance)
  summary$fitness_mean[scoped] <- subfit$summary$fitness_mean
  summary$fitness_sd[scoped] <- subfit$summary$fitness_sd
  summary$conf_low[scoped] <- subfit$summary$conf_low
  summary$conf_high[scoped] <- subfit$summary$conf_high
  summary$prediction_status[scoped] <- ifelse(summary$support_distance[scoped] == 0L, "anchor_gmrf", "truncated_gmrf")
  hp <- subfit$hyperparameters
  hp$truncated_scope <- paste0("support_distance<=", max_prediction_distance)
  diagnostics <- build_extrapolation_diagnostics(
    method, max_prediction_distance, subfit$anchors, summary, hp, subfit$tuning_grid,
    convergence_status = subfit$diagnostics$factorization_status,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    warnings = character(),
    extra = c(
      subfit$diagnostics,
      list(
        n_induced_nodes = length(scoped),
        n_induced_edges = length(subgraph$edge_from),
        lambda_l = hp$lambda_l,
        lambda_e = hp$lambda_e,
        sigma_obs = hp$sigma_obs,
        cv_score = hp$cv_score,
        anchor_cv_rmse = hp$cv_score,
        anchor_cv_mae = NA_real_,
        factorization_status = subfit$diagnostics$factorization_status
      )
    )
  )
  finish_extrapolation_fit(graph, summary, subfit$anchors, hp, subfit$tuning_grid, diagnostics)
}

wrap_graph_gaussian_fit <- function(local_fit,
                                    graph,
                                    max_prediction_distance,
                                    ...) {
  started <- proc.time()[["elapsed"]]
  out <- fit_graph_posterior(local_fit = local_fit, graph = graph, ...)
  out$summary$extrapolation_method <- "graph_gaussian_baseline"
  if (!"prediction_status" %in% names(out$summary)) {
    out$summary$prediction_status <- ifelse(is.finite(out$summary$fitness_mean), "graph_gaussian", "missing")
  }
  out$diagnostics$extrapolation_method <- "graph_gaussian_baseline"
  out$diagnostics$max_prediction_distance <- max_prediction_distance
  out$diagnostics$n_anchors <- nrow(out$anchors)
  out$diagnostics$n_predicted <- sum(is.finite(out$summary$fitness_mean))
  out$diagnostics$n_out_of_scope <- sum(as.character(out$summary$support_tier) %in% c("out_of_scope"))
  out$diagnostics$hyperparameters <- out$hyperparameters
  out$diagnostics$tuning_grid <- out$tuning_grid
  out$diagnostics$convergence_status <- out$diagnostics$factorization_status
  out$diagnostics$runtime_seconds <- as.numeric(proc.time()[["elapsed"]] - started)
  out$diagnostics$warnings <- character()
  out$diagnostics$dependency_status <- "ok"
  new_alfak2_global_fit(out)
}

#' Fit a second-layer extrapolation model
#'
#' @param local_fit An `alfak2_local_fit` object.
#' @param graph Optional global graph. Defaults to `local_fit$graph`.
#' @param method Extrapolation method.
#' @param max_prediction_distance Maximum near-field support distance reported as
#'   a valid prediction for near-field methods.
#' @param ... Method-specific controls.
#'
#' @return An `alfak2_global_fit` compatible list with `graph`, `summary`,
#'   `anchors`, `hyperparameters`, `tuning_grid`, and `diagnostics`.
#' @export
fit_extrapolation_layer <- function(local_fit,
                                    graph = NULL,
                                    method = c(
                                      "graph_gaussian_baseline",
                                      "edge_effect_empirical_bayes",
                                      "edge_effect_interaction_path_ensemble",
                                      "kronecker_or_graph_trend_filtering",
                                      "local_NNGP_or_GPnn",
                                      "delta_tree_ensemble",
                                      "tabpfn_nearfield_feature_model",
                                      "truncated_nearfield_gmrf",
                                      "local_polynomial_stencil"
                                    ),
                                    max_prediction_distance = 2,
                                    lambda_l_grid = c(0.2, 1, 5),
                                    lambda_e_grid = c(0.05, 0.25, 1),
                                    sigma_obs_grid = c(0.02, 0.05, 0.1),
                                    eps = 1e-5,
                                    graph_edge_weight = c("mutation", "unit", "normalized"),
                                    anchor_support_tiers = "all",
                                    anchor_covariance_inflation = c(
                                      TMB_sdreport = 1,
                                      untrusted_gradient = 4,
                                      untrusted_nonconverged = 9,
                                      untrusted_sdreport_missing = 4,
                                      untrusted_sdreport_nonfinite = 4,
                                      unknown = 4
                                    ),
                                    anchor_count_reference = NULL,
                                    anchor_count_power = 1,
                                    anchor_min_effective_count = 0,
                                    anchor_exclude = character(),
                                    compute_sd = TRUE,
                                    ...) {
  if (!inherits(local_fit, "alfak2_local_fit")) {
    stop("`local_fit` must be an alfak2_local_fit object.", call. = FALSE)
  }
  if (is.null(graph)) graph <- local_fit$graph
  if (!inherits(graph, "alfak2_graph")) stop("`graph` must be an alfak2_graph object.", call. = FALSE)
  method <- match_extrapolation_method(method)
  max_prediction_distance <- validate_prediction_distance(max_prediction_distance)
  graph_edge_weight <- match.arg(graph_edge_weight)
  if (identical(method, "graph_gaussian_baseline")) {
    return(wrap_graph_gaussian_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      lambda_l_grid = lambda_l_grid,
      lambda_e_grid = lambda_e_grid,
      sigma_obs_grid = sigma_obs_grid,
      eps = eps,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = anchor_support_tiers,
      anchor_covariance_inflation = anchor_covariance_inflation,
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power,
      anchor_min_effective_count = anchor_min_effective_count,
      anchor_exclude = anchor_exclude,
      compute_sd = compute_sd,
      ...
    ))
  }
  anchors <- prepare_extrapolation_anchors(
    local_fit = local_fit,
    graph = graph,
    anchor_support_tiers = anchor_support_tiers,
    anchor_covariance_inflation = anchor_covariance_inflation,
    anchor_count_reference = anchor_count_reference,
    anchor_count_power = anchor_count_power,
    anchor_min_effective_count = anchor_min_effective_count,
    anchor_exclude = anchor_exclude
  )
  switch(
    method,
    edge_effect_empirical_bayes = nearfield_local_edge_effect_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight,
      ...
    ),
    edge_effect_interaction_path_ensemble = edge_interaction_path_ensemble_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight,
      ...
    ),
    local_polynomial_stencil = local_polynomial_stencil_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      ...
    ),
    kronecker_or_graph_trend_filtering = graph_trend_filtering_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight,
      ...
    ),
    local_NNGP_or_GPnn = local_NNGP_or_GPnn_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      ...
    ),
    delta_tree_ensemble = delta_tree_ensemble_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight,
      ...
    ),
    tabpfn_nearfield_feature_model = tabpfn_nearfield_feature_model_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      graph_edge_weight = graph_edge_weight,
      ...
    ),
    truncated_nearfield_gmrf = truncated_nearfield_gmrf_fit(
      local_fit = local_fit,
      graph = graph,
      max_prediction_distance = max_prediction_distance,
      anchors = anchors,
      lambda_l_grid = lambda_l_grid,
      lambda_e_grid = lambda_e_grid,
      sigma_obs_grid = sigma_obs_grid,
      eps = eps,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = anchor_support_tiers,
      anchor_covariance_inflation = anchor_covariance_inflation,
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power,
      anchor_min_effective_count = anchor_min_effective_count,
      anchor_exclude = anchor_exclude,
      compute_sd = compute_sd,
      ...
    )
  )
}
