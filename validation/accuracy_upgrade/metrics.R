# These helpers never fit an affine map or choose a reference using truth.
finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
rmse <- function(x) sqrt(finite_mean(x^2))

heldout_count_metrics <- function(fit, dataset, model, fit_status, covariance_status) {
  frequency <- fit$posterior_predictive
  labels <- fit$global$summary$karyotype
  if (is.null(frequency)) {
    frequency <- fit$local$posterior_predictive
    labels <- fit$local$summary$karyotype
  }
  if (!is.null(frequency$karyotype)) labels <- frequency$karyotype
  if (length(frequency$pi0) != length(labels) || length(frequency$pi1) != length(labels)) return(NULL)
  idx <- match(dataset$target_labels, labels)
  p0 <- p1 <- numeric(length(dataset$target_labels))
  ok <- !is.na(idx)
  p0[ok] <- frequency$pi0[idx[ok]]; p1[ok] <- frequency$pi1[idx[ok]]
  if (any(!is.finite(c(p0, p1))) || any(c(p0, p1) < 0) || sum(p0) <= 0 || sum(p1) <= 0) return(NULL)
  set.seed(dataset$sampling_seed + 50000000L)
  heldout0 <- as.numeric(stats::rmultinom(1, dataset$config$n0, dataset$truth_pi0))
  heldout1 <- as.numeric(stats::rmultinom(1, dataset$config$n1, dataset$truth_pi1))
  floor <- 1e-12
  log_score <- function(y, p) {
    pp <- pmax(p, floor); pp <- pp / sum(pp)
    -sum(y * log(pp)) / sum(y)
  }
  vals <- c(heldout_t0_nll_per_cell_floor_1e12 = log_score(heldout0, p0),
            heldout_t1_nll_per_cell_floor_1e12 = log_score(heldout1, p1),
            heldout_t0_positive_count_on_zero_prediction = sum(heldout0[p0 == 0]),
            heldout_t1_positive_count_on_zero_prediction = sum(heldout1[p1 == 0]),
            expected_t0_cross_entropy_floor_1e12 = log_score(dataset$truth_pi0, p0),
            expected_t1_cross_entropy_floor_1e12 = log_score(dataset$truth_pi1, p1),
            pi0_rmse_full_target = rmse(p0 - dataset$truth_pi0),
            pi1_rmse_full_target = rmse(p1 - dataset$truth_pi1),
            predicted_pi0_target_mass = sum(p0), predicted_pi1_target_mass = sum(p1),
            truth_pi0_mass_with_prediction_support = sum(dataset$truth_pi0[p0 > 0]),
            truth_pi1_mass_with_prediction_support = sum(dataset$truth_pi1[p1 > 0]))
  data.frame(dataset_id = dataset$dataset_id, family = dataset$family, landscape_id = dataset$landscape_id,
             split = dataset$split, model = model, distance_group = "heldout_counts",
             fit_status = fit_status, covariance_status = covariance_status,
             metric = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)
}

# Recover the frozen baseline's complete Gaussian covariance without changing
# its predictions. Both its marginal variances and posterior means must match.
baseline_graph_covariance <- function(global, eps = 1e-5) {
  g <- global$graph
  h <- global$hyperparameters
  if (!identical(h$graph_edge_weight, "unit")) stop("Baseline covariance recovery requires unit edges")
  n <- length(g$labels)
  q <- diag(eps, n)
  add <- function(idx, coeff, scale) {
    q[idx, idx] <<- q[idx, idx, drop = FALSE] + scale * tcrossprod(coeff)
  }
  for (e in seq_along(g$edge_from)) {
    idx <- c(g$edge_from[e], g$edge_to[e])
    if (idx[1] != idx[2]) add(idx, c(1, -1), h$lambda_l)
  }
  label <- function(x) paste(x, collapse = ".")
  k <- g$karyotypes
  for (i in seq_len(n)) {
    x <- k[i, ]
    for (c in seq_len(ncol(k))) {
      lo <- hi <- x; lo[c] <- lo[c] - 1L; hi[c] <- hi[c] + 1L
      idx <- match(c(label(lo), label(x), label(hi)), g$labels)
      if (!anyNA(idx)) add(idx, c(1, -2, 1), h$lambda_e)
    }
    for (c1 in seq_len(ncol(k) - 1L)) for (c2 in (c1 + 1L):ncol(k)) {
      xc <- xd <- xcd <- x
      xc[c1] <- xc[c1] + 1L; xd[c2] <- xd[c2] + 1L
      xcd[c1] <- xcd[c1] + 1L; xcd[c2] <- xcd[c2] + 1L
      idx <- match(c(label(x), label(xc), label(xd), label(xcd)), g$labels)
      if (!anyNA(idx)) add(idx, c(1, -1, -1, 1), 0.5 * h$lambda_e)
    }
  }
  anchor <- global$anchors
  idx <- match(anchor$karyotype, g$labels)
  precision <- 1 / pmax(1e-10, anchor$variance + h$sigma_obs^2)
  rhs <- numeric(n)
  for (a in seq_along(idx)) {
    q[idx[a], idx[a]] <- q[idx[a], idx[a]] + precision[a]
    rhs[idx[a]] <- rhs[idx[a]] + precision[a] * anchor$mean[a]
  }
  covariance <- solve(q)
  recovered_mean <- as.numeric(covariance %*% rhs)
  mean_error <- max(abs(recovered_mean - global$summary$fitness_mean))
  variance_error <- max(abs(diag(covariance) - global$summary$fitness_sd^2))
  if (!is.finite(mean_error) || mean_error > 1e-7 || variance_error > 1e-7) {
    stop(sprintf("Baseline covariance recovery failed: mean=%g variance=%g", mean_error, variance_error))
  }
  dimnames(covariance) <- list(as.character(g$labels), as.character(g$labels))
  attr(covariance, "verification") <- c(mean_error = mean_error, variance_error = variance_error)
  covariance
}

extract_target_predictions <- function(fit, dataset, model, recover_baseline = FALSE) {
  s <- fit$global$summary
  idx <- match(dataset$target_labels, as.character(s$karyotype))
  pred <- s$fitness_mean[idx]
  sd <- s$fitness_sd[idx]
  w <- dataset$gauge_weight
  ref_ok <- all(is.finite(pred[w > 0]))
  relative <- if (ref_ok) pred - sum(w[w > 0] * pred[w > 0]) else rep(NA_real_, length(pred))
  covariance <- fit$global$fitness_covariance
  if (is.null(covariance)) covariance <- fit$global$covariance
  covariance_source <- if (!is.null(covariance)) "returned_full_covariance" else "unavailable"
  if (is.null(covariance) && recover_baseline) {
    covariance <- baseline_graph_covariance(fit$global)
    covariance_source <- "verified_reconstructed_baseline_precision"
  }
  relative_sd <- rep(NA_real_, length(pred))
  if (!is.null(covariance)) {
    covariance <- as.matrix(covariance)
    cov_labels <- rownames(covariance)
    if (is.null(cov_labels)) cov_labels <- fit$global$fitness_covariance_labels
    if (is.null(cov_labels)) cov_labels <- as.character(s$karyotype)
    ii <- match(dataset$target_labels, cov_labels)
    if (!anyNA(ii)) {
      covariance <- covariance[ii, ii, drop = FALSE]
      a <- diag(length(w)) - matrix(w, nrow = length(w), ncol = length(w), byrow = TRUE)
      relative_sd <- sqrt(pmax(0, diag(a %*% covariance %*% t(a))))
    } else covariance_source <- "covariance_missing_target_nodes"
  }
  nodes <- data.frame(dataset_id = dataset$dataset_id, family = dataset$family,
             landscape_id = dataset$landscape_id, split = dataset$split,
             model = model, karyotype = dataset$target_labels, distance = dataset$target_distance,
             truth_raw = dataset$truth, truth_relative = dataset$truth_relative,
             pred_raw = pred, pred_relative = relative, sd_raw = sd, sd_relative = relative_sd,
             covariance_source = covariance_source, stringsAsFactors = FALSE)
  interval_source <- if (all(is.finite(relative_sd))) "Gaussian_contrast_full_covariance" else "unvalidated_missing_covariance"
  for (level in dataset$config$interval_levels) {
    suffix <- sprintf("%02d", round(100 * level))
    q <- stats::qnorm((1 + level) / 2)
    nodes[[paste0("relative_low_", suffix)]] <- relative - q * relative_sd
    nodes[[paste0("relative_high_", suffix)]] <- relative + q * relative_sd
  }
  latent_mean <- fit$local$latent_fitness_mean
  latent_cov <- fit$local$latent_fitness_covariance
  if (identical(fit$local$diagnostics$forward_model, "division_process") &&
      !is.null(latent_mean) && !is.null(latent_cov)) {
    latent_idx <- match(dataset$target_labels, fit$local$summary$karyotype)
    if (!anyNA(latent_idx)) {
      latent_cov <- as.matrix(latent_cov)[latent_idx, latent_idx, drop = FALSE]
      latent_mean <- as.numeric(latent_mean)[latent_idx]
      eig <- eigen((latent_cov + t(latent_cov)) / 2, symmetric = TRUE)
      if (min(eig$values) > -1e-7 && all(is.finite(eig$values))) {
        set.seed(dataset$sampling_seed + 60000000L)
        z <- matrix(stats::rnorm(nrow(latent_cov) * 8000L), nrow = nrow(latent_cov))
        transform <- sweep(eig$vectors, 2L, sqrt(pmax(0, eig$values)), "*")
        draws <- exp(sweep(transform %*% z, 1L, latent_mean, "+"))
        reference <- as.numeric(crossprod(w, draws))
        draws <- sweep(draws, 2L, reference, "-")
        for (level in dataset$config$interval_levels) {
          suffix <- sprintf("%02d", round(100 * level))
          quantiles <- t(apply(draws, 1L, stats::quantile, probs = c((1 - level) / 2, (1 + level) / 2), names = FALSE))
          nodes[[paste0("relative_low_", suffix)]] <- quantiles[, 1]
          nodes[[paste0("relative_high_", suffix)]] <- quantiles[, 2]
        }
        interval_source <- "8000_joint_lognormal_Laplace_contrast_draws"
      }
    }
  }
  nodes$interval_source <- interval_source
  nodes
}

score_dataset <- function(nodes, dataset, runtime_seconds, fit_status, covariance_status) {
  rows <- list()
  distances <- sort(unique(dataset$target_distance))
  groups <- c(setNames(lapply(distances, function(d) which(nodes$distance == d)), paste0("d", distances)),
              list(all_node_micro = seq_len(nrow(nodes))))
  for (group in names(groups)) {
    z <- nodes[groups[[group]], , drop = FALSE]
    ok <- is.finite(z$pred_relative)
    err <- z$pred_relative - z$truth_relative
    sd_ok <- ok & is.finite(z$sd_relative) & z$sd_relative >= 0
    vals <- c(n_target = nrow(z), prediction_coverage = mean(ok),
              relative_rmse = rmse(err), relative_mae = finite_mean(abs(err)),
              relative_bias = finite_mean(err), raw_rmse = rmse(z$pred_raw - z$truth_raw),
              valid_interval_fraction = mean(sd_ok), runtime_seconds = runtime_seconds)
    for (level in dataset$config$interval_levels) {
      alpha <- 1 - level
      q <- stats::qnorm(1 - alpha / 2)
      suffix <- sprintf("%02d", round(100 * level))
      lo <- if (paste0("relative_low_", suffix) %in% names(z)) z[[paste0("relative_low_", suffix)]] else z$pred_relative - q * z$sd_relative
      hi <- if (paste0("relative_high_", suffix) %in% names(z)) z[[paste0("relative_high_", suffix)]] else z$pred_relative + q * z$sd_relative
      score <- hi - lo + 2 / alpha * pmax(lo - z$truth_relative, 0) +
        2 / alpha * pmax(z$truth_relative - hi, 0)
      vals[paste0("coverage_", suffix)] <- if (any(sd_ok)) mean(z$truth_relative[sd_ok] >= lo[sd_ok] & z$truth_relative[sd_ok] <= hi[sd_ok]) else NA_real_
      vals[paste0("interval_score_", suffix)] <- if (any(sd_ok)) mean(score[sd_ok]) else NA_real_
      vals[paste0("interval_width_", suffix)] <- if (any(sd_ok)) mean((hi - lo)[sd_ok]) else NA_real_
    }
    edges <- dataset$target_edges
    edge_distance <- pmax(nodes$distance[edges[, 1]], nodes$distance[edges[, 2]])
    ekeep <- if (group == "all_node_micro") rep(TRUE, nrow(edges)) else edge_distance == as.integer(sub("d", "", group))
    e <- edges[ekeep, , drop = FALSE]
    pg <- nodes$pred_relative[e[, 2]] - nodes$pred_relative[e[, 1]]
    tg <- nodes$truth_relative[e[, 2]] - nodes$truth_relative[e[, 1]]
    vals["edge_gradient_rmse"] <- rmse(pg - tg)
    vals["edge_sign_accuracy"] <- finite_mean(as.numeric(sign(pg) == sign(tg)))
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_id = dataset$dataset_id, family = dataset$family, landscape_id = dataset$landscape_id,
      split = dataset$split, model = nodes$model[1], distance_group = group,
      fit_status = fit_status, covariance_status = covariance_status,
      metric = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  exact <- out[grepl("^d[0-9]+$", out$distance_group), , drop = FALSE]
  macro <- aggregate(value ~ dataset_id + family + landscape_id + split + model + fit_status + covariance_status + metric,
                     exact, finite_mean, na.action = na.pass)
  macro$distance_group <- "all_distance_macro"
  rbind(out, macro[, names(out)])
}

paired_landscape_summary <- function(metrics, baseline_model = "baseline", bootstrap_replicates = 1000L, seed = 911204L) {
  set.seed(seed)
  keys <- c("dataset_id", "family", "landscape_id", "split", "distance_group", "metric")
  baseline <- metrics[metrics$model == baseline_model, c(keys, "value")]
  names(baseline)[names(baseline) == "value"] <- "baseline"
  candidate <- metrics[metrics$model != baseline_model, c(keys, "model", "value")]
  paired <- merge(candidate, baseline, by = keys, all.x = TRUE)
  paired$delta <- paired$value - paired$baseline
  groups <- split(paired, interaction(paired$model, paired$distance_group, paired$metric, drop = TRUE))
  rows <- lapply(groups, function(x) {
    valid <- is.finite(x$value) & is.finite(x$baseline)
    x <- x[valid, , drop = FALSE]
    if (!nrow(x)) return(NULL)
    land <- aggregate(cbind(value, baseline, delta) ~ landscape_id, x, mean)
    boot <- replicate(bootstrap_replicates, mean(sample(land$delta, nrow(land), replace = TRUE)))
    data.frame(model = x$model[1], distance_group = x$distance_group[1], metric = x$metric[1],
               n_landscapes = nrow(land), n_paired_datasets = nrow(x), baseline_mean = mean(land$baseline),
               candidate_mean = mean(land$value), paired_delta = mean(land$delta),
               delta_low = unname(stats::quantile(boot, 0.025)), delta_high = unname(stats::quantile(boot, 0.975)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}
