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
                                        correct_efflux = FALSE) {
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
    offset <- (sum(weights * fitness / viability) - g0) / sum(weights / viability)
    legacy_mean <- (fitness - offset) / viability
    legacy_sd <- as.numeric(summary$fitness_sd) / viability
  } else {
    center <- sum(weights * fitness)
    legacy_mean <- fitness + g0 - center
    legacy_sd <- as.numeric(summary$fitness_sd)
  }

  summary$fitness_mean_alfakR_scale <- legacy_mean
  summary$fitness_sd_alfakR_scale <- legacy_sd
  summary$conf_low_alfakR_scale <- legacy_mean - 1.959963984540054 * legacy_sd
  summary$conf_high_alfakR_scale <- legacy_mean + 1.959963984540054 * legacy_sd
  summary$efflux_viability <- viability
  summary
}

add_alfakR_scale_to_fit <- function(fit,
                                    n0,
                                    nb,
                                    correct_efflux = FALSE,
                                    legacy_weight = c("pi0", "directly_informed", "uniform")) {
  legacy_weight <- match.arg(legacy_weight)
  validate_legacy_scale_controls(n0, nb, fit$data$dt, fit$data$beta, correct_efflux)

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
    correct_efflux = correct_efflux
  )
  fit$global$summary <- add_alfakR_scale_to_summary(
    fit$global$summary,
    fit$global$graph$karyotypes,
    weights = global_weights,
    dt = fit$data$dt,
    beta = fit$data$beta,
    n0 = n0,
    nb = nb,
    correct_efflux = correct_efflux
  )
  fit$diagnostics$alfakR_scale <- list(
    enabled = TRUE,
    n0 = n0,
    nb = nb,
    g0 = log(nb / n0) / fit$data$dt,
    correct_efflux = correct_efflux,
    legacy_weight = legacy_weight,
    note = "alfak2 posterior fitness expressed on the alfakR absolute-growth scale"
  )
  fit
}
