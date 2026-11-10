# Empirical uncertainty calibration. Fitness point estimates and SDs are never
# modified. This is not a finite-sample conformal guarantee.
alfak2_calibration_scenario <- function(dataset) {
  config <- dataset$config
  list(beta = as.double(config$beta), dt = as.double(config$dt),
       depth_t0 = as.integer(config$n0), depth_t1 = as.integer(config$n1),
       n_chr = as.integer(config$n_chr), min_cn = as.integer(config$min_cn),
       max_cn = as.integer(config$max_cn),
       generating_protocol = if (is.null(config$generating_protocol)) "v1_tau_leap_division" else as.character(config$generating_protocol),
       simulation_time_step = if (identical(config$simulator_backend, "gillespie")) NA_real_ else as.double(config$time_step),
       reference_definition = as.character(config$gauge),
       reference_depth = if (is.null(config$depth_sensitivity)) NA_integer_ else as.integer(config$depth_sensitivity$reference_depth))
}

alfak2_calibration_model_key <- function(provenance) {
  required <- c("model", "implementation_hash", "arguments_hash")
  if (!all(required %in% names(provenance))) stop("Model provenance is missing required fingerprints")
  provenance[required]
}

alfak2_weighted_quantile <- function(x, weight, probability) {
  order <- order(x)
  x <- x[order]; weight <- weight[order] / sum(weight)
  x[which(cumsum(weight) + 1e-14 >= probability)[1L]]
}

fit_alfak2_interval_calibration <- function(calibration_table, model_key, scenario,
                                           levels = c(0.5, 0.8, 0.95), min_landscapes = 2L) {
  required <- c("landscape_id", "dataset_id", "split", "distance", "pred_relative", "sd_relative", "truth_relative")
  if (!all(required %in% names(calibration_table))) stop("Calibration table is missing required columns")
  if (!nrow(calibration_table) || anyNA(calibration_table$split) || any(calibration_table$split != "calibration"))
    stop("Only an independently reserved calibration split may fit uncertainty calibration")
  if (any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) stop("Levels must lie strictly between 0 and 1")
  if (length(unique(calibration_table$landscape_id)) < min_landscapes) stop("Insufficient independent calibration landscapes")
  table <- list()
  for (distance in sort(unique(calibration_table$distance))) {
    x <- calibration_table[calibration_table$distance == distance, , drop = FALSE]
    error <- x$truth_relative - x$pred_relative
    valid <- is.finite(error) & is.finite(x$sd_relative) & x$sd_relative >= 0
    # An exact reference contrast can have zero posterior variance and zero error.
    contradiction <- valid & x$sd_relative == 0 & abs(error) > 1e-12
    valid <- valid & !contradiction
    residual <- rep(NA_real_, nrow(x))
    residual[valid & x$sd_relative > 0] <- error[valid & x$sd_relative > 0] / x$sd_relative[valid & x$sd_relative > 0]
    residual[valid & x$sd_relative == 0] <- 0
    valid_landscapes <- unique(x$landscape_id[valid])
    weight <- rep(0, nrow(x))
    for (landscape in valid_landscapes) {
      datasets <- unique(x$dataset_id[valid & x$landscape_id == landscape])
      for (dataset in datasets) {
        ii <- which(valid & x$landscape_id == landscape & x$dataset_id == dataset)
        weight[ii] <- 1 / length(valid_landscapes) / length(datasets) / length(ii)
      }
    }
    supported <- length(valid_landscapes) >= min_landscapes
    for (level in levels) {
      lower <- upper <- NA_real_
      if (supported) {
        lower <- alfak2_weighted_quantile(residual[valid], weight[valid], (1 - level) / 2)
        upper <- alfak2_weighted_quantile(residual[valid], weight[valid], (1 + level) / 2)
      }
      table[[length(table) + 1L]] <- data.frame(distance = as.character(distance), level = level,
        q_lower = lower, q_upper = upper, n_landscapes = length(valid_landscapes),
        n_datasets = length(unique(x$dataset_id[valid])), n_nodes = nrow(x),
        n_valid = sum(valid), n_invalid = sum(!valid), n_zero_scale_contradictions = sum(contradiction),
        status = if (supported) "empirically_calibrated" else "unvalidated_insufficient_landscapes")
    }
  }
  structure(list(schema_version = 1L, model_key = model_key, scenario = scenario,
                 levels = levels, table = do.call(rbind, table),
                 calibration_landscapes = sort(unique(calibration_table$landscape_id)),
                 calibration_datasets = sort(unique(calibration_table$dataset_id)),
                 weighting = "equal landscape; equal dataset within landscape; equal node within dataset/distance",
                 interpretation = "empirical signed standardized-error intervals; no finite-sample coverage guarantee",
                 point_estimates_changed = FALSE, standard_errors_changed = FALSE),
            class = "alfak2_interval_calibration")
}

apply_alfak2_interval_calibration <- function(predictions, calibration, model_key, scenario) {
  if (!inherits(calibration, "alfak2_interval_calibration")) stop("Invalid calibration object")
  if (!identical(model_key, calibration$model_key)) stop("Calibration model fingerprint mismatch")
  if (!identical(scenario, calibration$scenario)) stop("Calibration scenario fingerprint mismatch (depth/beta/dt/chromosomes/CN/reference)")
  if (!all(c("distance", "pred_relative", "sd_relative") %in% names(predictions))) stop("Predictions are missing required columns")
  out <- predictions
  out$calibration_status <- "unvalidated_distance"
  valid_prediction <- is.finite(out$pred_relative) & is.finite(out$sd_relative) & out$sd_relative >= 0
  for (level in calibration$levels) {
    suffix <- sprintf("%02d", round(100 * level))
    rows <- calibration$table[calibration$table$level == level, , drop = FALSE]
    ii <- match(as.character(out$distance), rows$distance)
    supported <- !is.na(ii) & is.finite(rows$q_lower[ii]) & is.finite(rows$q_upper[ii])
    valid <- supported & valid_prediction
    lower <- upper <- rep(NA_real_, nrow(out))
    lower[valid] <- out$pred_relative[valid] + out$sd_relative[valid] * rows$q_lower[ii[valid]]
    upper[valid] <- out$pred_relative[valid] + out$sd_relative[valid] * rows$q_upper[ii[valid]]
    out[[paste0("calibrated_low_", suffix)]] <- lower
    out[[paste0("calibrated_high_", suffix)]] <- upper
    out$calibration_status[valid] <- "empirically_calibrated"
    out$calibration_status[supported & !valid_prediction] <- "invalid_prediction_or_scale"
  }
  out
}
