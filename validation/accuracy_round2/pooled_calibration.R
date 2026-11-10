# Landscape-balanced empirical calibration. Native point estimates and posterior
# SDs are retained. A separate calibration scale can flatten unreliable far-field
# SD ordering; a shared CDF regularizes distance-specific quantiles.
round2_calibration_weights <- function(x, balance_distance = FALSE) {
  w <- numeric(nrow(x))
  landscapes <- unique(x$landscape_id)
  for (landscape in landscapes) {
    datasets <- unique(x$dataset_id[x$landscape_id == landscape])
    for (dataset in datasets) {
      rows <- which(x$landscape_id == landscape & x$dataset_id == dataset)
      groups <- if (balance_distance) split(rows, x$distance[rows]) else list(rows)
      for (ii in groups) w[ii] <- 1 / length(landscapes) / length(datasets) / length(groups) / length(ii)
    }
  }
  w
}

round2_weighted_quantile <- function(x, w, p) {
  order <- order(x); x <- x[order]; w <- w[order] / sum(w)
  x[which(cumsum(w) >= p - 1e-14)[1L]]
}

round2_calibration_scale <- function(sd, median_sd, distance, gamma_far) {
  gamma <- ifelse(distance <= 1, 1, gamma_far)
  out <- rep(NA_real_, length(sd))
  good <- is.finite(sd) & sd >= 0 & is.finite(median_sd) & median_sd > 0
  positive <- good & sd > 0
  out[positive] <- median_sd[positive] * exp(gamma[positive] * log(sd[positive] / median_sd[positive]))
  # Preserve a deterministic reference contrast; do not manufacture its variance.
  out[good & sd == 0] <- 0
  out
}

round2_fit_calibration_internal <- function(data, model_key, scenario, gamma_far, kappa,
                                            levels, min_landscapes, role) {
  required <- c("landscape_id", "dataset_id", "split", "distance", "pred_relative", "sd_relative", "truth_relative")
  if (!all(required %in% names(data)) || !nrow(data)) stop("Missing calibration input columns or rows")
  if (identical(role, "calibration") && any(is.na(data$split) | data$split != "calibration"))
    stop("Final calibration requires only independently reserved calibration data")
  if (identical(role, "development_fold") && any(is.na(data$split) | data$split != "development"))
    stop("Development selection requires explicit development data")
  if (!is.numeric(gamma_far) || length(gamma_far) != 1 || !is.finite(gamma_far) || gamma_far < 0 || gamma_far > 1)
    stop("gamma_far must lie in [0,1]")
  if (length(kappa) != 1 || is.na(kappa) || kappa < 0) stop("kappa must be nonnegative")
  if (any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) stop("Invalid interval levels")
  valid <- is.finite(data$truth_relative) & is.finite(data$pred_relative) &
    is.finite(data$sd_relative) & data$sd_relative >= 0 & is.finite(data$distance)
  valid <- valid & !(data$sd_relative == 0 & abs(data$truth_relative - data$pred_relative) > 1e-12)
  x <- data[valid, , drop = FALSE]
  if (length(unique(x$landscape_id)) < min_landscapes) stop("Insufficient independent calibration landscapes")
  distances <- sort(unique(x$distance))
  medians <- vapply(distances, function(d) {
    z <- x[x$distance == d & x$sd_relative > 0, , drop = FALSE]
    if (!nrow(z)) return(1)
    round2_weighted_quantile(z$sd_relative, round2_calibration_weights(z), .5)
  }, numeric(1))
  names(medians) <- as.character(distances)
  x$calibration_scale <- round2_calibration_scale(x$sd_relative,
    medians[match(as.character(x$distance), names(medians))], x$distance, gamma_far)
  residual <- ifelse(x$calibration_scale > 0,
    (x$truth_relative - x$pred_relative) / x$calibration_scale, 0)
  pool_weight <- round2_calibration_weights(x, balance_distance = TRUE)
  table <- list()
  for (d in distances) {
    rows <- which(x$distance == d)
    n_land <- length(unique(x$landscape_id[rows]))
    local_weight <- if (is.infinite(kappa)) 0 else n_land / (n_land + kappa)
    weights <- c(local_weight * round2_calibration_weights(x[rows, , drop = FALSE]),
                 (1 - local_weight) * pool_weight)
    values <- c(residual[rows], residual)
    supported <- n_land >= min_landscapes
    for (level in levels) table[[length(table) + 1L]] <- data.frame(
      distance = d, level = level,
      q_lower = if (supported) round2_weighted_quantile(values, weights, (1 - level) / 2) else NA_real_,
      q_upper = if (supported) round2_weighted_quantile(values, weights, (1 + level) / 2) else NA_real_,
      n_landscapes = n_land, n_datasets = length(unique(x$dataset_id[rows])),
      local_cdf_weight = local_weight,
      status = if (supported) "empirical_partial_pooling" else "unsupported_distance")
  }
  structure(list(schema_version = 1L, role = role, model_key = model_key,
    scenario = scenario, gamma_far = gamma_far, kappa = kappa, median_sd = medians,
    table = do.call(rbind, table), levels = sort(levels),
    calibration_landscapes = sort(unique(x$landscape_id)),
    n_invalid_rows = sum(!valid), point_estimates_changed = FALSE,
    posterior_standard_errors_changed = FALSE,
    interpretation = "empirical error intervals; no distribution-free finite-sample guarantee"),
    class = "alfak2_round2_calibration")
}

fit_round2_pooled_calibration <- function(data, model_key, scenario, gamma_far = 1,
                                         kappa = 0, levels = c(.5, .8, .95), min_landscapes = 2L) {
  round2_fit_calibration_internal(data, model_key, scenario, gamma_far, kappa,
    levels, min_landscapes, role = "calibration")
}

apply_round2_pooled_calibration <- function(data, calibration, model_key, scenario) {
  if (!inherits(calibration, "alfak2_round2_calibration")) stop("Invalid calibration object")
  if (!identical(model_key, calibration$model_key)) stop("Model fingerprint mismatch")
  if (!identical(scenario, calibration$scenario)) stop("Scenario fingerprint mismatch")
  if (!all(c("pred_relative", "sd_relative", "distance") %in% names(data))) stop("Invalid predictions")
  out <- data
  out$calibration_scale <- round2_calibration_scale(out$sd_relative,
    calibration$median_sd[match(as.character(out$distance), names(calibration$median_sd))],
    out$distance, calibration$gamma_far)
  out$calibration_status <- "unsupported_distance_or_scale"
  for (level in calibration$levels) {
    suffix <- sprintf("%02d", round(level * 100))
    q <- calibration$table[calibration$table$level == level, , drop = FALSE]
    index <- match(out$distance, q$distance)
    lower <- out$pred_relative + out$calibration_scale * q$q_lower[index]
    upper <- out$pred_relative + out$calibration_scale * q$q_upper[index]
    out[[paste0("calibrated_low_", suffix)]] <- lower
    out[[paste0("calibrated_high_", suffix)]] <- upper
    good <- is.finite(lower) & is.finite(upper)
    out$calibration_status[good] <- "empirical_partial_pooling"
  }
  out
}

select_round2_calibration <- function(development, model_key, scenario,
                                      gammas = c(0, .5, 1), kappas = c(0, 6, 18, Inf),
                                      levels = c(.5, .8, .95), required_distances = NULL) {
  required <- c("landscape_id", "dataset_id", "split", "distance", "pred_relative", "sd_relative", "truth_relative")
  if (!all(required %in% names(development)) || !nrow(development)) stop("Missing development columns or rows")
  if (anyNA(development$split) || any(development$split != "development")) stop("Selection is restricted to development data")
  if (anyNA(development$landscape_id) || anyNA(development$dataset_id) ||
      any(!is.finite(development$distance)) || any(!is.finite(development$pred_relative)) ||
      any(!is.finite(development$truth_relative)) || any(!is.finite(development$sd_relative)) ||
      any(development$sd_relative < 0) ||
      any(development$sd_relative == 0 & abs(development$truth_relative - development$pred_relative) > 1e-12))
    stop("Incomplete or invalid development predictions cannot be excluded from calibration selection")
  distance_scope <- if (is.null(required_distances)) "observed_development_distances_only" else "explicit_predeclared_distances"
  if (is.null(required_distances)) required_distances <- sort(unique(development$distance))
  if (!length(required_distances) || any(!is.finite(required_distances)) || anyDuplicated(required_distances))
    stop("required_distances must be a nonempty unique finite vector")
  if (!setequal(sort(unique(development$distance)), required_distances))
    stop("Development distance strata do not match required_distances")
  counts <- vapply(required_distances, function(d)
    length(unique(development$landscape_id[development$distance == d])), integer(1))
  if (any(counts < 3L)) stop("Unsupported distance: landscape-LOO requires at least three independent landscapes per distance")
  landscapes <- unique(development$landscape_id)
  if (length(landscapes) < 4) stop("At least four independent development landscapes are required")
  grid <- expand.grid(gamma_far = gammas, kappa = kappas)
  scores <- list()
  for (g in seq_len(nrow(grid))) for (held in landscapes) {
    train <- development[development$landscape_id != held, , drop = FALSE]
    test <- development[development$landscape_id == held, , drop = FALSE]
    cal <- round2_fit_calibration_internal(train, model_key, scenario,
      grid$gamma_far[g], grid$kappa[g], levels, 2L, "development_fold")
    pred <- apply_round2_pooled_calibration(test, cal, model_key, scenario)
    for (level in levels) {
      suffix <- sprintf("%02d", round(level * 100))
      low <- pred[[paste0("calibrated_low_", suffix)]]
      high <- pred[[paste0("calibrated_high_", suffix)]]
      alpha <- 1 - level; truth <- pred$truth_relative
      loss <- high - low + 2 / alpha * pmax(low - truth, 0) + 2 / alpha * pmax(truth - high, 0)
      if (any(!is.finite(loss)))
        stop("LOO calibration has unsupported held-out rows; missing intervals cannot be silently dropped")
      z <- data.frame(dataset_id = pred$dataset_id, distance = pred$distance, loss = loss)
      case <- stats::aggregate(loss ~ dataset_id + distance, z, mean)
      distance <- stats::aggregate(loss ~ distance, case, mean)
      scores[[length(scores) + 1L]] <- data.frame(grid_id = g, landscape_id = held,
        level = level, distance = distance$distance, loss = distance$loss)
    }
  }
  scores <- do.call(rbind, scores)
  risk <- stats::aggregate(loss ~ grid_id + distance + level, scores, mean)
  expected <- expand.grid(distance = required_distances, level = levels)
  expected_keys <- paste(expected$distance, expected$level)
  for (g in seq_len(nrow(grid))) {
    x <- risk[risk$grid_id == g, ]
    if (nrow(x) != nrow(expected) || !setequal(paste(x$distance, x$level), expected_keys) ||
        any(!is.finite(x$loss))) stop("Incomplete distance-by-level calibration risk grid")
  }
  baseline <- which(grid$gamma_far == 1 & grid$kappa == 0)
  if (length(baseline) != 1L) stop("Selection grid must include the original gamma=1,kappa=0 reference")
  reference <- risk[risk$grid_id == baseline, ]
  summary <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
    x <- risk[risk$grid_id == g, ]; index <- match(paste(x$distance, x$level), paste(reference$distance, reference$level))
    data.frame(grid_id = g, gamma_far = grid$gamma_far[g], kappa = grid$kappa[g],
      macro_interval_score = mean(x$loss),
      all_distance_nonworse = all(x$loss <= reference$loss[index] + 1e-12))
  }))
  eligible <- which(summary$all_distance_nonworse)
  selected <- eligible[which.min(summary$macro_interval_score[eligible])]
  list(selected = summary[selected, , drop = FALSE], grid_summary = summary,
       distance_risk = risk, heldout_scores = scores,
       required_distances = required_distances, distance_scope = distance_scope,
       landscapes_per_distance = stats::setNames(counts, required_distances),
       rule = "landscape-LOO development only; no distance mean interval score worse than original calibration")
}
