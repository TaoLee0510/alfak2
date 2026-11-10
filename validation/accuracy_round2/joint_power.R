# Prospective planning from actual out-of-landscape paired loss vectors.
# The simulated event covers the declared core risk/shape conjunction only,
# using one-sided alpha=.05 per component. It does not certify the separate
# rank/amplitude/availability gates or the full promotion decision.
round2_joint_power <- function(oof_pairs, n_total = c(36L, 45L, 60L), variance_multiplier = 1.5,
                                repetitions = 5000L, alpha = 0.05, seed = 610041L,
                                effect_mode = c("observed_oof", "common_5percent_target"),
                                effect_shrink = 0.5, include_shape = TRUE,
                                required_distances = NULL) {
  effect_mode <- match.arg(effect_mode)
  required <- c("candidate_id", "family", "landscape_id", "distance", "metric", "reference_loss", "candidate_loss", "oof_verified")
  if (!all(required %in% names(oof_pairs)) || !nrow(oof_pairs) ||
      anyNA(oof_pairs$oof_verified) || any(!oof_pairs$oof_verified)) stop("Verified OOF paired losses are required")
  if (length(unique(oof_pairs$candidate_id)) != 1L) stop("Freeze one winner/selection procedure before joint power planning")
  families <- sort(unique(oof_pairs$family)); stopifnot(length(families) == 3L)
  benefits <- c("rmse", "interval_score95", if (include_shape) c("spatially_centered_rmse", "edge_gradient_rmse"))
  safety <- benefits
  if (is.null(required_distances) || !length(required_distances) ||
      any(!is.finite(required_distances)) || anyDuplicated(required_distances))
    stop("Supply explicit predeclared required_distances; observed input alone cannot define the promotion gates")
  if (length(variance_multiplier) != 1L || !is.finite(variance_multiplier) || variance_multiplier <= 0 ||
      length(effect_shrink) != 1L || !is.finite(effect_shrink) || effect_shrink < 0 || effect_shrink > 1)
    stop("Invalid planning variance multiplier or effect shrinkage")
  oof_pairs <- oof_pairs[oof_pairs$metric %in% benefits, ]
  if (any(!is.finite(oof_pairs$reference_loss)) || any(!is.finite(oof_pairs$candidate_loss))) stop("Incomplete OOF strata cannot be silently dropped or imputed")
  required_keys <- expand.grid(metric = benefits, distance = required_distances, stringsAsFactors = FALSE)
  keys <- sort(paste(required_keys$metric, required_keys$distance, sep = "/"))
  if (!setequal(unique(paste(oof_pairs$metric, oof_pairs$distance, sep = "/")), keys))
    stop("Missing or unexpected predeclared metric-by-distance promotion gates")
  row_key <- paste(oof_pairs$family, oof_pairs$landscape_id, oof_pairs$metric, oof_pairs$distance, sep = "/")
  if (anyNA(oof_pairs$family) || anyNA(oof_pairs$landscape_id) || anyDuplicated(row_key))
    stop("OOF gate rows must have unique nonmissing family/landscape/metric/distance keys")
  key_risk <- sub("/.*$", "", keys)
  data <- vector("list", 3L)
  for (fi in seq_along(families)) {
    x <- oof_pairs[oof_pairs$family == families[fi], ]
    lands <- sort(unique(x$landscape_id))
    if (length(lands) < 4L) stop("At least four independent OOF landscapes per family required for planning; no guarantee is implied")
    reference <- candidate <- matrix(NA_real_, length(lands), length(keys))
    indices <- cbind(match(x$landscape_id, lands), match(paste(x$metric, x$distance, sep = "/"), keys))
    reference[indices] <- x$reference_loss; candidate[indices] <- x$candidate_loss
    if (anyNA(reference) || anyNA(candidate)) stop("Missing OOF gate vectors; planning remains unvalidated")
    means_r <- colMeans(reference); means_c <- colMeans(candidate)
    target_c <- if (effect_mode == "common_5percent_target") 0.95 * means_r else means_r + effect_shrink * (means_c - means_r)
    joint <- cbind(reference, candidate)
    centered <- sweep(joint, 2L, colMeans(joint), "-") * sqrt(variance_multiplier)
    data[[fi]] <- list(centered = centered, mean = c(means_r, target_c),
      superiority_variance = apply(candidate - reference, 2L, stats::var) * variance_multiplier,
      safety_variance = apply(candidate - 1.05 * reference, 2L, stats::var) * variance_multiplier)
  }
  set.seed(seed); results <- list(); z <- stats::qnorm(1 - alpha); k <- length(keys)
  for (n in n_total) {
    if (n %% 3L || n > 60L) stop("Test budget must be family-balanced and at most 60 landscapes per domain")
    per_family <- n / 3L
    global_se <- sqrt(Reduce(`+`, lapply(data, `[[`, "superiority_variance")) / per_family) / 3
    safety_se <- lapply(data, function(x) sqrt(x$safety_variance / per_family))
    gate_count <- k + length(families) * sum(key_risk %in% safety)
    pass <- matrix(FALSE, repetitions, gate_count)
    for (b in seq_len(repetitions)) {
      means <- lapply(data, function(x) {
        chosen <- sample.int(nrow(x$centered), per_family, replace = TRUE)
        colMeans(x$centered[chosen, , drop = FALSE]) + x$mean
      })
      differences <- lapply(means, function(x) x[(k + 1L):(2L * k)] - x[seq_len(k)])
      global <- Reduce(`+`, differences) / 3
      flags <- global + z * global_se < 0
      for (fi in seq_along(families)) {
        delta <- means[[fi]][(k + 1L):(2L * k)] - 1.05 * means[[fi]][seq_len(k)]
        flags <- c(flags, (delta + z * safety_se[[fi]] < 0)[key_risk %in% safety])
      }
      pass[b, ] <- flags
    }
    joint_pass <- rowMeans(pass) == 1
    probability <- mean(joint_pass)
    results[[length(results) + 1L]] <- data.frame(n_test_landscapes = n, n_per_family = per_family,
      n_gates = gate_count, core_test_rule = "IUT_one_sided_0p05", effect_mode = effect_mode,
      full_promotion_validated = FALSE,
      observed_effect_shrink = if (effect_mode == "observed_oof") effect_shrink else NA_real_,
      variance_multiplier = variance_multiplier, bootstrap_planning_repetitions = repetitions,
      required_distances = paste(sort(required_distances), collapse = ","),
      approximate_core_gate_pass_probability = probability, monte_carlo_se_only = sqrt(probability * (1 - probability) / repetitions),
      weakest_marginal_pass_probability = min(colMeans(pass)),
      limitations = paste("OOF empirical-vector/normal-bound core-gate planning only; excludes rank/amplitude/availability gates;",
        "requires nested calibration/selection OOF; not a power or coverage guarantee"))
  }
  do.call(rbind, results)
}
