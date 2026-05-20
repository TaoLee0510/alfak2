alfak2_first_nonnull <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

alfak2_finite_sd <- function(x, fallback = 0.05) {
  x <- as.numeric(x)
  out <- if (sum(is.finite(x)) >= 2L) stats::sd(x, na.rm = TRUE) else NA_real_
  if (!is.finite(out) || out <= 0) fallback else out
}

alfak2_bootstrap_counts <- function(counts) {
  counts <- validate_count_matrix(counts)
  out <- counts
  for (j in seq_len(ncol(counts))) {
    total <- sum(counts[, j], na.rm = TRUE)
    if (!is.finite(total) || total <= 0) {
      out[, j] <- 0
      next
    }
    prob <- as.numeric(counts[, j]) / total
    prob[!is.finite(prob) | prob < 0] <- 0
    if (sum(prob) <= 0) {
      out[, j] <- 0
      next
    }
    out[, j] <- as.integer(stats::rmultinom(1L, size = as.integer(round(total)), prob = prob)[, 1L])
  }
  dimnames(out) <- dimnames(counts)
  out
}

alfak2_graph_edges_1based <- function(graph) {
  from <- as.integer(graph$edge_from)
  to <- as.integer(graph$edge_to)
  if (!length(from) || !length(to)) {
    return(data.frame(from = integer(), to = integer(), weight = numeric()))
  }
  if (min(c(from, to), na.rm = TRUE) == 0L) {
    from <- from + 1L
    to <- to + 1L
  }
  data.frame(from = from, to = to, weight = as.numeric(graph$edge_weight))
}

alfak2_graph_adjacency <- function(graph) {
  n <- length(graph$labels)
  edges <- alfak2_graph_edges_1based(graph)
  adj <- vector("list", n)
  if (!nrow(edges)) return(adj)
  for (i in seq_len(nrow(edges))) {
    a <- edges$from[i]
    b <- edges$to[i]
    if (a >= 1L && a <= n && b >= 1L && b <= n) {
      adj[[a]] <- c(adj[[a]], b)
      adj[[b]] <- c(adj[[b]], a)
    }
  }
  lapply(adj, unique)
}

alfak2_bfs_distances <- function(adj, source, max_depth = Inf) {
  n <- length(adj)
  dist <- rep(Inf, n)
  if (!length(source) || is.na(source) || source < 1L || source > n) return(dist)
  source <- as.integer(source)
  dist[source] <- 0
  queue <- integer(n)
  queue[1L] <- source
  head <- 1L
  tail <- 1L
  while (head <= tail) {
    u <- queue[head]
    head <- head + 1L
    if (dist[u] >= max_depth) next
    nb <- adj[[u]]
    if (!length(nb)) next
    new <- nb[is.infinite(dist[nb])]
    if (length(new)) {
      dist[new] <- dist[u] + 1L
      n_new <- length(new)
      queue[(tail + 1L):(tail + n_new)] <- new
      tail <- tail + n_new
    }
  }
  dist
}

alfak2_distance_matrix_from_sources <- function(graph, sources, max_depth = Inf) {
  adj <- alfak2_graph_adjacency(graph)
  sources <- as.integer(sources)
  out <- matrix(Inf, nrow = length(graph$labels), ncol = length(sources))
  for (j in seq_along(sources)) out[, j] <- alfak2_bfs_distances(adj, sources[j], max_depth = max_depth)
  colnames(out) <- as.character(graph$labels)[sources]
  rownames(out) <- as.character(graph$labels)
  out
}

alfak2_parent_edges_1based <- function(graph) {
  from <- as.integer(unlist(graph$parent_from0))
  to <- as.integer(unlist(graph$parent_to0))
  weight <- as.numeric(unlist(graph$parent_weight))
  if (!length(from) || !length(to)) {
    return(data.frame(from = integer(), to = integer(), weight = numeric()))
  }
  data.frame(from = from + 1L, to = to + 1L, weight = weight)
}

fit_alfak2_direct_local_experimental <- function(counts,
                                                 dt = 1,
                                                 beta = 0.00005,
                                                 transition_kernel = c("exact", "linear"),
                                                 min_cn = 0,
                                                 max_cn = 5,
                                                 max_nodes = 5000,
                                                 input_depth = c("raw", "effective"),
                                                 effective_depth = NULL,
                                                 effective_depth_mode = c("min", "cap", "fixed"),
                                                 effective_depth_rounding = c("hash", "largest_remainder", "stochastic"),
                                                 effective_depth_seed = NULL,
                                                 observation_model = NULL,
                                                 dm_concentration = NULL,
                                                 observation_weight_mode = c("likelihood", "fractional_count"),
                                                 control = list(eval.max = 500, iter.max = 500),
                                                 retry_control = list(eval.max = 2000, iter.max = 2000),
                                                 ...) {
  modes <- validate_effective_depth_mode(input_depth, effective_depth_mode)
  input_depth <- modes$input_depth
  effective_depth_mode <- modes$effective_depth_mode
  transition_kernel <- match.arg(transition_kernel)
  effective_depth_rounding <- match.arg(effective_depth_rounding)
  observation_weight_mode <- match_observation_weight_mode(observation_weight_mode)
  obs_controls <- resolve_fit_observation_controls(input_depth, observation_model, dm_concentration)
  data <- prepare_counts_for_input_depth(
    counts,
    dt = dt,
    beta = beta,
    input_depth = input_depth,
    effective_depth = effective_depth,
    effective_depth_mode = effective_depth_mode,
    effective_depth_rounding = effective_depth_rounding,
    effective_depth_seed = effective_depth_seed
  )
  graph <- build_karyotype_graph(
    data,
    transition_kernel = transition_kernel,
    shell_depth = 0L,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  fit_local_posterior(
    data,
    graph,
    observation_model = obs_controls$observation_model,
    dm_concentration = obs_controls$dm_concentration,
    observation_weight_mode = observation_weight_mode,
    control = control,
    retry_control = retry_control,
    ...
  )
}

extract_alfak2_direct_state_experimental <- function(fit,
                                                     fitness_boot = NULL,
                                                     bootstrap_diagnostics = NULL) {
  if (!inherits(fit, "alfak2_fit")) {
    stop("`fit` must be an alfak2_fit object.", call. = FALSE)
  }
  local_summary <- fit$local$summary
  graph <- fit$global$graph
  graph_labels <- as.character(graph$labels)
  direct <- as.character(local_summary$support_tier) == "directly_informed" &
    is.finite(local_summary$fitness_mean) &
    is.finite(local_summary$fitness_sd)
  if ("effective_count_total" %in% names(local_summary)) {
    direct <- direct & is.finite(local_summary$effective_count_total) &
      local_summary$effective_count_total > 0
  } else if ("count_total" %in% names(local_summary)) {
    direct <- direct & is.finite(local_summary$count_total) & local_summary$count_total > 0
  }
  local_direct <- local_summary[direct, , drop = FALSE]
  if (!nrow(local_direct)) {
    stop("No finite directly informed local nodes were available.", call. = FALSE)
  }
  direct_labels <- as.character(local_direct$karyotype)
  direct_index <- match(direct_labels, graph_labels)
  keep <- !is.na(direct_index)
  local_direct <- local_direct[keep, , drop = FALSE]
  direct_labels <- direct_labels[keep]
  direct_index <- direct_index[keep]
  node_table <- data.frame(
    node_id = seq_along(graph_labels),
    karyotype = graph_labels,
    support_tier = as.character(graph$support_tier),
    support_distance = as.integer(graph$support_distance),
    is_direct = seq_along(graph_labels) %in% direct_index,
    stringsAsFactors = FALSE
  )
  fitness_hat <- stats::setNames(as.numeric(local_direct$fitness_mean), direct_labels)
  fitness_sd <- stats::setNames(as.numeric(local_direct$fitness_sd), direct_labels)
  x0_hat <- if ("pi0" %in% names(local_direct)) {
    stats::setNames(as.numeric(local_direct$pi0), direct_labels)
  } else {
    stats::setNames(rep(NA_real_, length(direct_labels)), direct_labels)
  }
  x0_hat[!is.finite(x0_hat) | x0_hat < 0] <- 0
  if (sum(x0_hat) > 0) x0_hat <- x0_hat / sum(x0_hat)
  if (!is.null(fitness_boot)) {
    fitness_boot <- as.matrix(fitness_boot)
    missing <- setdiff(direct_labels, colnames(fitness_boot))
    if (length(missing)) {
      add <- matrix(NA_real_, nrow = nrow(fitness_boot), ncol = length(missing),
                    dimnames = list(rownames(fitness_boot), missing))
      fitness_boot <- cbind(fitness_boot, add)
    }
    fitness_boot <- fitness_boot[, direct_labels, drop = FALSE]
  }
  diagnostics <- list(
    local_convergence = fit$local$diagnostics$convergence,
    local_gradient_norm = fit$local$diagnostics$gradient_norm,
    local_covariance_status = fit$local$diagnostics$covariance_status,
    local_objective = fit$local$diagnostics$objective,
    n_direct = length(direct_labels),
    bootstrap = bootstrap_diagnostics
  )
  out <- list(
    node_table = node_table,
    direct_node_ids = direct_index,
    direct_labels = direct_labels,
    fitness_hat = fitness_hat,
    fitness_sd = fitness_sd,
    fitness_boot = fitness_boot,
    x0_hat = x0_hat,
    counts = fit$data$counts,
    ntot = colSums(fit$data$counts),
    timepoints = seq(0, by = fit$data$dt, length.out = ncol(fit$data$counts)),
    dt = fit$data$dt,
    beta = fit$data$beta,
    local_summary = local_summary,
    graph = graph,
    graph_edges = alfak2_graph_edges_1based(graph),
    parent_edges = alfak2_parent_edges_1based(graph),
    diagnostics = diagnostics
  )
  class(out) <- "alfak2_direct_state"
  out
}

alfak2_direct_delta_prior <- function(direct_state, direct_f, nn_prior = "empirical") {
  direct_labels <- direct_state$direct_labels
  direct_idx <- direct_state$direct_node_ids
  f_by_idx <- rep(NA_real_, length(direct_state$graph$labels))
  f_by_idx[direct_idx] <- as.numeric(direct_f[direct_labels])
  pe <- direct_state$parent_edges
  if (!nrow(pe)) {
    return(list(mu = 0, sigma = alfak2_finite_sd(direct_f), n_delta = 0L, mode_used = "none"))
  }
  ok <- pe$from %in% direct_idx & pe$to %in% direct_idx
  delta <- f_by_idx[pe$to[ok]] - f_by_idx[pe$from[ok]]
  delta <- delta[is.finite(delta)]
  if (!length(delta)) {
    return(list(mu = 0, sigma = alfak2_finite_sd(direct_f), n_delta = 0L, mode_used = "none"))
  }
  sigma_floor <- max(alfak2_finite_sd(delta, fallback = 0.05), 0.01)
  if (identical(nn_prior, "none")) {
    return(list(mu = 0, sigma = sigma_floor, n_delta = length(delta), mode_used = "none"))
  }
  if (identical(nn_prior, "empirical")) {
    return(list(mu = mean(delta), sigma = sigma_floor, n_delta = length(delta), mode_used = "empirical"))
  }
  if (identical(nn_prior, "empirical_censored")) {
    delta2 <- c(delta, rep(0, max(1L, ceiling(length(delta) / 2))))
    return(list(mu = mean(delta2), sigma = max(stats::sd(delta2), 0.01),
                n_delta = length(delta), mode_used = "empirical_censored"))
  }
  if (identical(nn_prior, "empirical_censored_weighted")) {
    counts <- direct_state$local_summary$count_total[
      match(names(direct_f), as.character(direct_state$local_summary$karyotype))
    ]
    count_by_idx <- rep(NA_real_, length(direct_state$graph$labels))
    count_by_idx[direct_idx] <- as.numeric(counts)
    w <- sqrt(pmax(1, count_by_idx[pe$from[ok]] + count_by_idx[pe$to[ok]]))
    w <- w[is.finite(delta)]
    if (!length(w) || sum(w) <= 0) w <- rep(1, length(delta))
    mu <- stats::weighted.mean(delta, w)
    return(list(mu = mu, sigma = sigma_floor, n_delta = length(delta),
                mode_used = "empirical_censored_weighted"))
  }
  if (identical(nn_prior, "empirical_two_step")) {
    return(list(mu = mean(delta), sigma = max(sigma_floor, 0.02),
                n_delta = length(delta), mode_used = "empirical_two_step"))
  }
  stop("Unsupported `nn_prior`: ", nn_prior, call. = FALSE)
}

alfak2_predict_native_nn_once <- function(direct_state,
                                          direct_f,
                                          nn_prior = "empirical_censored",
                                          nn_shell_depth = 1L,
                                          prior = NULL) {
  graph <- direct_state$graph
  labels <- as.character(graph$labels)
  direct_idx <- direct_state$direct_node_ids
  direct_labels <- direct_state$direct_labels
  direct_f <- stats::setNames(as.numeric(direct_f[direct_labels]), direct_labels)
  prior <- alfak2_first_nonnull(prior, alfak2_direct_delta_prior(direct_state, direct_f, nn_prior))
  dmat <- alfak2_distance_matrix_from_sources(graph, direct_idx, max_depth = nn_shell_depth)
  min_dist <- apply(dmat, 1L, min, na.rm = TRUE)
  target_idx <- which(!seq_along(labels) %in% direct_idx & is.finite(min_dist) & min_dist <= nn_shell_depth)
  if (!length(target_idx)) {
    return(list(summary = data.frame(), fitness = stats::setNames(numeric(), character()),
                diagnostics = data.frame()))
  }
  pe <- direct_state$parent_edges
  rows <- vector("list", length(target_idx))
  diag <- vector("list", length(target_idx))
  f_by_idx <- rep(NA_real_, length(labels))
  f_by_idx[direct_idx] <- as.numeric(direct_f[direct_labels])
  names(f_by_idx) <- labels
  for (ii in seq_along(target_idx)) {
    target <- target_idx[ii]
    dist <- dmat[target, ]
    src <- which(is.finite(dist) & dist <= min_dist[target])
    src_idx <- direct_idx[src]
    src_labels <- labels[src_idx]
    w <- exp(-as.numeric(dist[src]))
    if (!length(w) || sum(w) <= 0) w <- rep(1, length(src_idx))
    direction <- rep(0, length(src_idx))
    for (j in seq_along(src_idx)) {
      if (nrow(pe)) {
        if (any(pe$from == src_idx[j] & pe$to == target)) direction[j] <- 1
        if (any(pe$from == target & pe$to == src_idx[j])) direction[j] <- -1
      }
    }
    direction[direction == 0] <- 1
    step_mult <- if (identical(nn_prior, "empirical_two_step")) pmax(1, dist[src]) else 1
    pred_each <- as.numeric(direct_f[src_labels]) + direction * prior$mu * step_mult
    pred <- stats::weighted.mean(pred_each, w)
    rows[[ii]] <- data.frame(
      node_id = target,
      karyotype = labels[target],
      support_scope = "nn",
      support_distance = as.integer(min_dist[target]),
      fitness_mean = pred,
      fitness_sd = prior$sigma / sqrt(max(1, length(src_idx))),
      nn_prior_used = prior$mode_used,
      n_parent_anchors = length(src_idx),
      stringsAsFactors = FALSE
    )
    diag[[ii]] <- data.frame(
      karyotype = labels[target],
      nn_prior_mode_requested = nn_prior,
      nn_prior_mode_used = prior$mode_used,
      prior_mu_hat = prior$mu,
      prior_sigma_hat = prior$sigma,
      n_delta = prior$n_delta,
      n_parent_anchors = length(src_idx),
      nearest_distance = min_dist[target],
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, rows)
  fitness <- stats::setNames(summary$fitness_mean, summary$karyotype)
  list(summary = summary, fitness = fitness, diagnostics = do.call(rbind, diag), prior = prior)
}

fit_alfak2_native_nn_experimental <- function(direct_state,
                                              nn_prior = c("empirical_censored", "empirical_censored_weighted", "empirical_two_step", "none", "empirical"),
                                              nn_shell_depth = 1L) {
  nn_prior <- match.arg(nn_prior)
  point <- alfak2_predict_native_nn_once(
    direct_state,
    direct_state$fitness_hat,
    nn_prior = nn_prior,
    nn_shell_depth = nn_shell_depth
  )
  boot <- NULL
  if (!is.null(direct_state$fitness_boot) && nrow(direct_state$fitness_boot) > 0L && length(point$fitness)) {
    boot <- matrix(NA_real_, nrow = nrow(direct_state$fitness_boot), ncol = length(point$fitness),
                   dimnames = list(rownames(direct_state$fitness_boot), names(point$fitness)))
    for (b in seq_len(nrow(direct_state$fitness_boot))) {
      pred <- alfak2_predict_native_nn_once(
        direct_state,
        direct_state$fitness_boot[b, , drop = TRUE],
        nn_prior = nn_prior,
        nn_shell_depth = nn_shell_depth,
        prior = point$prior
      )
      common <- intersect(colnames(boot), names(pred$fitness))
      boot[b, common] <- pred$fitness[common]
    }
    if (nrow(boot) >= 2L) {
      boot_sd <- apply(boot, 2L, stats::sd, na.rm = TRUE)
      point$summary$fitness_sd <- pmax(point$summary$fitness_sd, boot_sd[point$summary$karyotype], na.rm = TRUE)
    }
  }
  out <- list(
    summary = point$summary,
    fitness = point$fitness,
    fitness_boot = boot,
    diagnostics = point$diagnostics,
    prior = point$prior,
    nn_prior = nn_prior,
    nn_shell_depth = nn_shell_depth
  )
  class(out) <- "alfak2_native_nn_state"
  out
}

fit_alfak2_native_kriging_experimental <- function(direct_state,
                                                   nn_state = NULL,
                                                   range = NULL,
                                                   nugget = 1e-4,
                                                   max_anchors = 300L) {
  graph <- direct_state$graph
  labels <- as.character(graph$labels)
  direct_idx <- direct_state$direct_node_ids
  obs_labels <- direct_state$direct_labels
  obs_idx <- direct_idx
  obs_y <- as.numeric(direct_state$fitness_hat[obs_labels])
  obs_sd <- as.numeric(direct_state$fitness_sd[obs_labels])
  if (!is.null(nn_state) && nrow(nn_state$summary)) {
    obs_labels <- c(obs_labels, as.character(nn_state$summary$karyotype))
    obs_idx <- c(obs_idx, as.integer(nn_state$summary$node_id))
    obs_y <- c(obs_y, as.numeric(nn_state$summary$fitness_mean))
    obs_sd <- c(obs_sd, as.numeric(nn_state$summary$fitness_sd))
  }
  ok <- is.finite(obs_idx) & is.finite(obs_y)
  obs_labels <- obs_labels[ok]
  obs_idx <- obs_idx[ok]
  obs_y <- obs_y[ok]
  obs_sd <- obs_sd[ok]
  unique_anchor <- !duplicated(obs_idx)
  obs_labels <- obs_labels[unique_anchor]
  obs_idx <- obs_idx[unique_anchor]
  obs_y <- obs_y[unique_anchor]
  obs_sd <- obs_sd[unique_anchor]
  if (length(obs_y) < 2L) {
    stop("Native Kriging requires at least two finite direct/NN anchors.", call. = FALSE)
  }
  if (length(obs_y) > max_anchors) {
    direct_keep <- seq_along(obs_y)[obs_idx %in% direct_idx]
    other <- setdiff(seq_along(obs_y), direct_keep)
    other <- other[order(obs_sd[other], na.last = TRUE)]
    keep <- c(direct_keep, utils::head(other, max_anchors - length(direct_keep)))
    keep <- keep[seq_len(min(length(keep), max_anchors))]
    obs_labels <- obs_labels[keep]
    obs_idx <- obs_idx[keep]
    obs_y <- obs_y[keep]
    obs_sd <- obs_sd[keep]
  }
  sigma2 <- stats::var(obs_y)
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- NA_real_
  range_arg <- if (is.null(range)) NA_real_ else as.numeric(range)[1L]
  anchor_values <- matrix(obs_y, nrow = 1L,
                          dimnames = list("point", obs_labels))
  if (!is.null(direct_state$fitness_boot) && nrow(direct_state$fitness_boot) > 0L) {
    direct_boot <- direct_state$fitness_boot
    nn_boot <- if (!is.null(nn_state)) nn_state$fitness_boot else NULL
    boot_values <- matrix(rep(obs_y, each = nrow(direct_boot)),
                          nrow = nrow(direct_boot),
                          ncol = length(obs_y),
                          dimnames = list(rownames(direct_boot), obs_labels))
    for (b in seq_len(nrow(direct_boot))) {
      direct_common <- intersect(obs_labels, colnames(direct_boot))
      if (length(direct_common)) {
        boot_values[b, match(direct_common, obs_labels)] <- direct_boot[b, direct_common]
      }
      if (!is.null(nn_boot) && b <= nrow(nn_boot)) {
        nn_common <- intersect(obs_labels, colnames(nn_boot))
        if (length(nn_common)) {
          boot_values[b, match(nn_common, obs_labels)] <- nn_boot[b, nn_common]
        }
      }
    }
    anchor_values <- rbind(anchor_values, boot_values)
  }
  edges <- alfak2_graph_edges_1based(graph)
  cpp <- alfak2_native_kriging_cpp(
    n_nodes = length(labels),
    edge_from = as.integer(edges$from),
    edge_to = as.integer(edges$to),
    anchor_index = as.integer(obs_idx),
    anchor_values = anchor_values,
    anchor_sd = as.numeric(obs_sd),
    range = range_arg,
    sigma2 = sigma2,
    nugget = nugget,
    compute_variance = TRUE
  )
  predictions <- as.matrix(cpp$predictions)
  rownames(predictions) <- rownames(anchor_values)
  colnames(predictions) <- labels
  samples <- if (nrow(predictions) > 1L) predictions[-1L, , drop = FALSE] else NULL
  summary <- data.frame(
    node_id = seq_along(labels),
    karyotype = labels,
    fitness_mean = as.numeric(cpp$mean),
    fitness_sd = as.numeric(cpp$sd),
    support_scope = "kriging",
    stringsAsFactors = FALSE
  )
  out <- list(
    summary = summary,
    posterior_samples = samples,
    anchors = data.frame(karyotype = obs_labels, node_id = obs_idx, fitness = obs_y, fitness_sd = obs_sd,
                         stringsAsFactors = FALSE),
    diagnostics = list(
      engine = cpp$engine,
      range = cpp$range,
      sigma2 = cpp$sigma2,
      nugget = cpp$nugget,
      auto_range = cpp$auto_range,
      n_anchors = length(obs_y),
      max_anchors = max_anchors,
      solve_status = cpp$solve_status
    )
  )
  class(out) <- "alfak2_native_kriging_state"
  out
}

alfak2_native_summary <- function(direct_state, nn_state, kriging_state) {
  out <- kriging_state$summary
  out$support_scope <- "kriging"
  direct_idx <- direct_state$direct_node_ids
  out$support_scope[direct_idx] <- "direct"
  out$fitness_mean[direct_idx] <- as.numeric(direct_state$fitness_hat[direct_state$direct_labels])
  out$fitness_sd[direct_idx] <- as.numeric(direct_state$fitness_sd[direct_state$direct_labels])
  if (!is.null(nn_state) && nrow(nn_state$summary)) {
    nn_idx <- as.integer(nn_state$summary$node_id)
    out$support_scope[nn_idx] <- "nn"
    out$fitness_mean[nn_idx] <- as.numeric(nn_state$summary$fitness_mean)
    out$fitness_sd[nn_idx] <- as.numeric(nn_state$summary$fitness_sd)
  }
  out$fq <- out$support_scope == "direct"
  out$nn <- out$support_scope == "nn"
  out
}

fit_alfak2_nn_kriging_experimental <- function(counts,
                                               dt = 1,
                                               beta = 0.00005,
                                               transition_kernel = c("exact", "linear"),
                                               local_shell_depth = 0,
                                               global_extra_shell = 1,
                                               min_cn = 0,
                                               max_cn = 5,
                                               max_nodes = 5000,
                                               lambda_l_grid = c(0.2, 1, 5),
                                               lambda_e_grid = c(0.05, 0.25, 1),
                                               sigma_obs_grid = c(0.02, 0.05, 0.1),
                                               graph_edge_weight = c("mutation", "unit", "normalized"),
                                               input_depth = c("raw", "effective"),
                                               effective_depth = NULL,
                                               effective_depth_mode = c("min", "cap", "fixed"),
                                               effective_depth_rounding = c("hash", "largest_remainder", "stochastic"),
                                               effective_depth_seed = NULL,
                                               observation_model = NULL,
                                               dm_concentration = NULL,
                                               observation_weight_mode = c("likelihood", "fractional_count"),
                                               nn_prior = c("empirical_censored", "empirical_censored_weighted", "empirical_two_step", "none", "empirical"),
                                               nn_shell_depth = 1L,
                                               kriging_range = NULL,
                                               kriging_nugget = 1e-4,
                                               kriging_max_anchors = 300L,
                                               nboot = 20L,
                                               seed = NULL,
                                               control = list(eval.max = 500, iter.max = 500),
                                               retry_control = list(eval.max = 2000, iter.max = 2000),
                                               ...) {
  transition_kernel <- match.arg(transition_kernel)
  graph_edge_weight <- match.arg(graph_edge_weight)
  input_depth <- match.arg(input_depth)
  effective_depth_mode <- match.arg(effective_depth_mode)
  effective_depth_rounding <- match.arg(effective_depth_rounding)
  observation_weight_mode <- match_observation_weight_mode(observation_weight_mode)
  nn_prior <- match.arg(nn_prior)
  nboot <- as.integer(nboot)
  if (!is.finite(nboot) || nboot < 0L) stop("`nboot` must be a non-negative integer.", call. = FALSE)
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  local_shell_depth <- as.integer(local_shell_depth)
  if (local_shell_depth != 0L) {
    warning("Experimental native NN/Kriging is intended for direct-informed `local_shell_depth = 0`.")
  }
  graphgp_fit <- fit_alfak2(
    counts,
    dt = dt,
    beta = beta,
    transition_kernel = transition_kernel,
    local_shell_depth = local_shell_depth,
    global_extra_shell = global_extra_shell,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes,
    lambda_l_grid = lambda_l_grid,
    lambda_e_grid = lambda_e_grid,
    sigma_obs_grid = sigma_obs_grid,
    graph_edge_weight = graph_edge_weight,
    anchor_support_tiers = "directly_informed",
    input_depth = input_depth,
    effective_depth = effective_depth,
    effective_depth_mode = effective_depth_mode,
    effective_depth_rounding = effective_depth_rounding,
    effective_depth_seed = effective_depth_seed,
    observation_model = observation_model,
    dm_concentration = dm_concentration,
    observation_weight_mode = observation_weight_mode,
    control = control,
    retry_control = retry_control,
    ...
  )
  direct_labels <- graphgp_fit$local$summary$karyotype[
    as.character(graphgp_fit$local$summary$support_tier) == "directly_informed"
  ]
  direct_labels <- as.character(direct_labels)
  fitness_boot <- NULL
  boot_diag <- data.frame()
  if (nboot > 0L && length(direct_labels)) {
    fitness_boot <- matrix(NA_real_, nrow = nboot, ncol = length(direct_labels),
                           dimnames = list(paste0("boot_", seq_len(nboot)), direct_labels))
    boot_rows <- vector("list", nboot)
    for (b in seq_len(nboot)) {
      bc <- alfak2_bootstrap_counts(counts)
      lf <- tryCatch(
        fit_alfak2_direct_local_experimental(
          bc,
          dt = dt,
          beta = beta,
          transition_kernel = transition_kernel,
          min_cn = min_cn,
          max_cn = max_cn,
          max_nodes = max_nodes,
          input_depth = input_depth,
          effective_depth = effective_depth,
          effective_depth_mode = effective_depth_mode,
          effective_depth_rounding = effective_depth_rounding,
          effective_depth_seed = effective_depth_seed,
          observation_model = observation_model,
          dm_concentration = dm_concentration,
          observation_weight_mode = observation_weight_mode,
          control = control,
          retry_control = retry_control,
          ...
        ),
        error = function(e) e
      )
      if (inherits(lf, "error")) {
        boot_rows[[b]] <- data.frame(
          replicate = b,
          status = "error",
          error_message = conditionMessage(lf),
          local_convergence = NA_integer_,
          local_gradient_norm = NA_real_,
          local_covariance_status = NA_character_,
          stringsAsFactors = FALSE
        )
        next
      }
      s <- lf$summary
      m <- match(direct_labels, as.character(s$karyotype))
      ok <- !is.na(m) & is.finite(s$fitness_mean[m])
      fitness_boot[b, ok] <- as.numeric(s$fitness_mean[m[ok]])
      boot_rows[[b]] <- data.frame(
        replicate = b,
        status = "ok",
        error_message = NA_character_,
        local_convergence = lf$diagnostics$convergence,
        local_gradient_norm = lf$diagnostics$gradient_norm,
        local_covariance_status = lf$diagnostics$covariance_status,
        stringsAsFactors = FALSE
      )
    }
    boot_diag <- do.call(rbind, boot_rows)
  }
  direct_state <- extract_alfak2_direct_state_experimental(
    graphgp_fit,
    fitness_boot = fitness_boot,
    bootstrap_diagnostics = boot_diag
  )
  nn_state <- fit_alfak2_native_nn_experimental(
    direct_state,
    nn_prior = nn_prior,
    nn_shell_depth = nn_shell_depth
  )
  kriging_state <- fit_alfak2_native_kriging_experimental(
    direct_state,
    nn_state,
    range = kriging_range,
    nugget = kriging_nugget,
    max_anchors = kriging_max_anchors
  )
  out <- list(
    data = graphgp_fit$data,
    graphgp = graphgp_fit$global,
    local = graphgp_fit$local,
    direct_state = direct_state,
    nn = nn_state,
    kriging = kriging_state,
    summary = alfak2_native_summary(direct_state, nn_state, kriging_state),
    posterior_samples = kriging_state$posterior_samples,
    diagnostics = list(
      backend = "native_nn_kriging_experimental",
      nn_prior = nn_prior,
      nn_shell_depth = nn_shell_depth,
      local = graphgp_fit$local$diagnostics,
      graphgp = graphgp_fit$global$diagnostics,
      direct = direct_state$diagnostics,
      nn = nn_state$diagnostics,
      kriging = kriging_state$diagnostics
    )
  )
  class(out) <- "alfak2_native_nn_kriging_fit"
  out
}
