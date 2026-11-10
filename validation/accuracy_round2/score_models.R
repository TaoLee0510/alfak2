round2_extract_nodes <- function(fit, dataset, model, legacy_scoring) {
  nodes <- legacy_scoring$extract_target_predictions(fit, dataset, model,
    recover_baseline = identical(model, "original_baseline"))
  components <- fit$global$posterior_components
  if (length(components)) {
    draws <- tryCatch(draw_alfak2_posterior_components(components, dataset$target_labels,
      n_draws = 8000L, seed = dataset$sampling_seed + 60000000L,
      default_labels = as.character(fit$global$summary$karyotype)), error = function(e) e)
    if (inherits(draws, "error")) {
      for (level in dataset$config$interval_levels) {
        suffix <- sprintf("%02d", round(100 * level))
        nodes[[paste0("relative_low_", suffix)]] <- NA_real_
        nodes[[paste0("relative_high_", suffix)]] <- NA_real_
      }
      nodes$interval_source <- paste("unvalidated_mixture_distribution", conditionMessage(draws))
    } else {
      draws <- sweep(draws, 2L, as.numeric(crossprod(dataset$gauge_weight, draws)), "-")
      for (level in dataset$config$interval_levels) {
        suffix <- sprintf("%02d", round(100 * level))
        q <- t(apply(draws, 1L, stats::quantile, probs = c((1 - level) / 2, (1 + level) / 2), names = FALSE))
        nodes[[paste0("relative_low_", suffix)]] <- q[, 1]
        nodes[[paste0("relative_high_", suffix)]] <- q[, 2]
      }
      nodes$interval_source <- "8000_actual_component_mixture_reference_contrasts"
    }
  }
  nodes
}

round2_score_nodes <- function(nodes, dataset) {
  e <- nodes$pred_relative - nodes$truth_relative
  centered <- e - round2_mean(e)
  edges <- dataset$target_edges
  d <- dataset$target_distance
  edge_d <- pmax(d[edges[, 1]], d[edges[, 2]])
  low <- nodes$relative_low_95; high <- nodes$relative_high_95
  if (is.null(low)) low <- high <- rep(NA_real_, nrow(nodes))
  safe_cor <- function(a, b, method) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5L || stats::sd(a[ok]) < 1e-12 || stats::sd(b[ok]) < 1e-12) return(NA_real_)
    suppressWarnings(stats::cor(a[ok], b[ok], method = method))
  }
  rows <- list()
  for (distance in c(-1L, sort(unique(d)))) {
    ii <- if (distance == -1L) seq_len(nrow(nodes)) else which(d == distance)
    ee <- if (distance == -1L) seq_len(nrow(edges)) else which(edge_d == distance)
    gradient_error <- e[edges[ee, 2]] - e[edges[ee, 1]]
    valid <- is.finite(e[ii]); iv <- valid & is.finite(low[ii]) & is.finite(high[ii])
    truth_sd <- if (length(ii) > 1L) stats::sd(nodes$truth_relative[ii]) else NA_real_
    pred_sd <- if (length(ii) > 1L) stats::sd(nodes$pred_relative[ii]) else NA_real_
    amplitude <- if (is.finite(truth_sd) && truth_sd > 0) pred_sd / truth_sd else NA_real_
    score <- high[ii] - low[ii] + 40 * pmax(low[ii] - nodes$truth_relative[ii], 0) +
      40 * pmax(nodes$truth_relative[ii] - high[ii], 0)
    rows[[length(rows) + 1L]] <- data.frame(dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
      family = dataset$family, split = dataset$split, model = nodes$model[1], distance = distance,
      n_target = length(ii), n_edges = length(ee), prediction_fraction = mean(valid), interval_fraction = mean(iv),
      rmse = sqrt(round2_mean(e[ii]^2)), mse = round2_mean(e[ii]^2), mae = round2_mean(abs(e[ii])),
      spatially_centered_rmse = sqrt(round2_mean(centered[ii]^2)), spatially_centered_mse = round2_mean(centered[ii]^2),
      global_offset_squared = round2_mean(e)^2,
      edge_gradient_rmse = sqrt(round2_mean(gradient_error^2)),
      spearman = safe_cor(nodes$pred_relative[ii], nodes$truth_relative[ii], "spearman"),
      pearson = safe_cor(nodes$pred_relative[ii], nodes$truth_relative[ii], "pearson"),
      amplitude_ratio = amplitude, variance_ratio = amplitude^2, amplitude_fidelity_error = abs(amplitude - 1),
      coverage95 = if (any(iv)) mean(nodes$truth_relative[ii][iv] >= low[ii][iv] & nodes$truth_relative[ii][iv] <= high[ii][iv]) else NA_real_,
      interval_score95 = if (any(iv)) mean(score[iv]) else NA_real_,
      width95 = if (any(iv)) mean((high[ii] - low[ii])[iv]) else NA_real_)
  }
  do.call(rbind, rows)
}
