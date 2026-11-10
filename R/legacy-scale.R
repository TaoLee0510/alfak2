validate_legacy_scale_controls <- function(n0, nb, dt, beta, correct_efflux) {
  validate_scalar(n0, "n0", lower = .Machine$double.eps)
  validate_scalar(nb, "nb", lower = .Machine$double.eps)
  validate_scalar(dt, "dt", lower = .Machine$double.eps)
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  if (!is.logical(correct_efflux) || length(correct_efflux) != 1L || is.na(correct_efflux)) {
    stop("`correct_efflux` must be TRUE or FALSE.", call. = FALSE)
  }
}

normalize_legacy_weights <- function(weights, fitness, label) {
  if (length(weights) != length(fitness)) {
    stop(sprintf("`%s` weights must match the number of summary rows.", label), call. = FALSE)
  }
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0 | !is.finite(fitness)] <- 0
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop(sprintf("`%s` weights do not contain positive mass on finite fitness values.", label), call. = FALSE)
  }
  weights / total
}

resolve_legacy_scale_weights <- function(local_fit,
                                         target_summary,
                                         legacy_weight = c("pi0", "directly_informed", "uniform")) {
  legacy_weight <- match.arg(legacy_weight)
  labels <- as.character(target_summary$karyotype)
  weights <- numeric(length(labels))

  if (legacy_weight == "pi0") {
    pi0 <- local_fit$posterior_predictive$pi0
    if (is.null(pi0) || length(pi0) != nrow(local_fit$summary)) {
      stop("Local fit does not contain a usable `pi0` posterior predictive vector.", call. = FALSE)
    }
    pi0 <- as.numeric(pi0)
    names(pi0) <- as.character(local_fit$summary$karyotype)
    idx <- match(labels, names(pi0))
    ok <- !is.na(idx)
    weights[ok] <- pi0[idx[ok]]
  } else if (legacy_weight == "directly_informed") {
    weights[as.character(target_summary$support_tier) == "directly_informed"] <- 1
  } else {
    weights[] <- 1
  }

  normalize_legacy_weights(weights, target_summary$fitness_mean, legacy_weight)
}

add_alfakR_scale_to_summary <- function(summary,
                                        karyotypes,
                                        weights,
                                        dt,
                                        beta,
                                        n0,
                                        nb,
                                        correct_efflux = FALSE,
                                        fitness_covariance = NULL) {
  validate_legacy_scale_controls(n0, nb, dt, beta, correct_efflux)
  if (!is.matrix(karyotypes) || nrow(karyotypes) != nrow(summary)) {
    stop("`karyotypes` must be a matrix with one row per summary row.", call. = FALSE)
  }
  if (!all(c("fitness_mean", "fitness_sd") %in% names(summary))) {
    stop("`summary` must contain `fitness_mean` and `fitness_sd`.", call. = FALSE)
  }

  fitness <- as.numeric(summary$fitness_mean)
  weights <- normalize_legacy_weights(weights, fitness, "legacy_scale")
  g0 <- log(nb / n0) / dt
  viability <- 2 * (1 - beta)^rowSums(karyotypes) - 1
  if (any(!is.finite(viability))) {
    stop("Computed efflux viability contains non-finite values.", call. = FALSE)
  }

  if (isTRUE(correct_efflux)) {
    if (any(viability <= 0)) {
      affected <- as.character(summary$karyotype[viability <= 0])
      stop(
        sprintf(
          "`correct_efflux = TRUE` requires positive viability for every graph node; affected: %s",
          paste(utils::head(affected, 5), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  # Match finite-time population growth, not the instantaneous arithmetic mean.
  exponent <- dt * fitness
  exponent[!is.finite(fitness)] <- -Inf
  shift <- max(exponent[weights > 0])
  growth_weights <- weights * exp(exponent - shift)
  log_growth <- (log(sum(growth_weights)) + shift) / dt
  growth_weights <- growth_weights / sum(growth_weights)
  offset <- log_growth - g0
  scale <- if (isTRUE(correct_efflux)) 1 / viability else rep(1, length(fitness))
  legacy_mean <- (fitness - offset) * scale
  if (is.null(fitness_covariance)) {
    marginal_variance <- as.numeric(summary$fitness_sd)^2
    marginal_variance[!is.finite(marginal_variance)] <- 0
    covariance <- diag(marginal_variance, nrow(summary))
  } else {
    covariance <- as.matrix(fitness_covariance)
    if (!identical(dim(covariance), c(nrow(summary), nrow(summary))) || any(!is.finite(covariance)))
      stop("Legacy calibration covariance must be finite and aligned with summary.", call. = FALSE)
  }
  jacobian <- diag(nrow(summary)) - matrix(growth_weights, nrow(summary), nrow(summary), byrow = TRUE)
  jacobian <- jacobian * scale
  calibrated_covariance <- jacobian %*% covariance %*% t(jacobian)
  legacy_sd <- sqrt(pmax(0, diag(calibrated_covariance)))
  legacy_sd[!is.finite(summary$fitness_sd) | !is.finite(fitness)] <- NA_real_

  summary$fitness_mean_alfakR_scale <- legacy_mean
  summary$fitness_sd_alfakR_scale <- legacy_sd
  summary$conf_low_alfakR_scale <- legacy_mean - 1.959963984540054 * legacy_sd
  summary$conf_high_alfakR_scale <- legacy_mean + 1.959963984540054 * legacy_sd
  if ("interval_valid" %in% names(summary)) {
    invalid <- !summary$interval_valid
    summary$conf_low_alfakR_scale[invalid] <- NA_real_
    summary$conf_high_alfakR_scale[invalid] <- NA_real_
  }
  summary$efflux_viability <- viability
  attr(summary, "fitness_covariance_alfakR_scale") <- calibrated_covariance
  summary
}

add_alfakR_scale_to_fit <- function(fit,
                                    n0,
                                    nb,
                                    correct_efflux = FALSE,
                                    legacy_weight = c("pi0", "directly_informed", "uniform")) {
  legacy_weight <- match.arg(legacy_weight)
  validate_legacy_scale_controls(n0, nb, fit$data$dt, fit$data$beta, correct_efflux)
  if (identical(fit$local$diagnostics$forward_model, "division_process"))
    stop("Division-process growth is calibrated during fitting: supply `population_sizes = c(N0, Nb)` instead of `alfakR_scale`.", call. = FALSE)

  local_weights <- resolve_legacy_scale_weights(fit$local, fit$local$summary, legacy_weight)
  global_weights <- resolve_legacy_scale_weights(fit$local, fit$global$summary, legacy_weight)

  fit$local$summary <- add_alfakR_scale_to_summary(
    fit$local$summary,
    fit$local$graph$karyotypes,
    weights = local_weights,
    dt = fit$data$dt,
    beta = fit$data$beta,
    n0 = n0,
    nb = nb,
    correct_efflux = correct_efflux,
    fitness_covariance = fit$local$fitness_covariance
  )
  fit$global$summary <- add_alfakR_scale_to_summary(
    fit$global$summary,
    fit$global$graph$karyotypes,
    weights = global_weights,
    dt = fit$data$dt,
    beta = fit$data$beta,
    n0 = n0,
    nb = nb,
    correct_efflux = correct_efflux,
    fitness_covariance = fit$global$fitness_covariance
  )
  fit$diagnostics$alfakR_scale <- list(
    enabled = TRUE,
    n0 = n0,
    nb = nb,
    g0 = log(nb / n0) / fit$data$dt,
    correct_efflux = correct_efflux,
    legacy_weight = legacy_weight,
    note = "finite-time population growth calibration; covariance conditional on pi0 weights and supplied population sizes"
  )
  fit
}
