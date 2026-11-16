# The sparse joint model when a rate it assumes is not known.
#
# The baseline birth rate cannot be learned from two time points and has to be
# supplied; neither a single dataset's nor a cohort's marginal likelihood
# locates it. What it changes is known: the fitness difference between two
# well-sampled clones comes from their frequency change and does not depend on
# it, while the deficit of a rare or unobserved karyotype comes from its level
# under mutation-selection balance, which fixes the deficit per division. Per
# unit of time that deficit is proportional to the assumed rate. Rankings
# barely move. Fitting over a range of rates shows all of this for a given
# dataset.
#
# The missegregation rate acts on the same deficits through the influx it
# assumes. Its marginal likelihood points somewhere, but with few cells per
# time point, short intervals or strong drift it points too high
# (validation/sparse_joint/beta_choice), so it is not chosen from the data
# either. A range of missegregation rates is fitted the same way.

#' Fit the sparse joint model over a range of birth or missegregation rates
#'
#' Use this when the baseline birth rate or the missegregation rate is unknown
#' and only a plausible range can be given. The model is fitted at each of a
#' few rates across the range, with the other rate held at its given value.
#' [stats::predict()] then reports relative fitness at one of them, an
#' interval that covers every rate in the range, and relative fitness per
#' division.
#'
#' The two directions of error are not alike. Assuming too low a rate shrinks
#' the predicted deficits of rare and unobserved karyotypes towards zero; on
#' simulated landscapes a birth rate 4 to 10 times too low still predicts
#' better than no differences. Assuming too high a rate inflates them, and a
#' birth rate 4 times too high predicts worse than no differences. Rankings
#' change little either way. The point prediction is therefore taken at the
#' lowest rate of the range by default, which is never too high when the range
#' holds the true rate. The same holds for the missegregation rate: on
#' simulated passaged cultures with true rates of 1, 4 and 16 times 5e-5, a
#' rate 4 times too high raised the median error 1.9 to 3.7 times, and one 4
#' times too low 1.2 to 1.6 times.
#'
#' The missegregation rate is not chosen by marginal likelihood: with
#' few cells per time point, short intervals or strong drift between passages,
#' the likelihood prefers rates several times too high, and no rule tried
#' (gates on interval and cells, evidence thresholds, a prior) separated those
#' cases from real differences. A range shows how much the conclusions depend
#' on it.
#'
#' @param ... Arguments of the fitting function other than the rate that is
#'   varied: for [fit_sparse_joint()], `counts`, `dt`, `beta` or
#'   `birth_rate`, and so on.
#' @param birth_rate_range Lowest and highest plausible baseline birth rate
#'   per unit of `dt`.
#' @param beta_range Lowest and highest plausible missegregation rate (per
#'   chromosome per division), in place of `birth_rate_range`. Give one of
#'   the two.
#' @param n_rates Number of rates fitted, spaced evenly on the log scale. An
#'   odd number puts one at the centre of the range.
#' @param point Rate of the point prediction: `"lowest"` of the range, or its
#'   geometric `"centre"`, which is more accurate when the true rate is near
#'   the centre and much less accurate when it is near the lowest.
#' @param fit_function [fit_sparse_joint()] or [fit_sparse_joint_cohort()].
#'
#' @return An `alfak2_sparse_joint_rate_fit` object: `rates`, the `fits` at
#'   each of them, `parameter` (`"birth_rate"` or `"beta"`, the rate that
#'   varies), `point`, the index of the fit used for point predictions, and
#'   `central`, the index of the fit at the centre.
#' @export
fit_sparse_joint_rates <- function(..., birth_rate_range = NULL, beta_range = NULL, n_rates = 3L,
                                   point = c("lowest", "centre"), fit_function = fit_sparse_joint) {
  if (is.null(birth_rate_range) == is.null(beta_range)) {
    stop("Give one of `birth_rate_range` and `beta_range`.", call. = FALSE)
  }
  parameter <- if (is.null(beta_range)) "birth_rate" else "beta"
  range_name <- paste0(parameter, "_range")
  range <- if (parameter == "beta") beta_range else birth_rate_range
  if (!is.numeric(range) || length(range) != 2L || any(!is.finite(range)) ||
      range[[1L]] <= 0 || range[[2L]] <= range[[1L]]) {
    stop(sprintf("`%s` must give a positive lowest and a larger highest rate.", range_name), call. = FALSE)
  }
  validate_scalar(n_rates, "n_rates", lower = 2)
  point <- match.arg(point)
  if (parameter %in% names(list(...))) {
    stop(sprintf("Give `%s`, not `%s`.", range_name, parameter), call. = FALSE)
  }
  rates <- exp(seq(log(range[[1L]]), log(range[[2L]]), length.out = as.integer(n_rates)))
  fits <- lapply(rates, function(rate) {
    if (parameter == "beta") fit_function(..., beta = rate) else fit_function(..., birth_rate = rate)
  })
  central <- (length(rates) + 1L) %/% 2L
  structure(list(rates = rates, fits = fits, parameter = parameter,
                 point = if (point == "lowest") 1L else central, central = central),
            class = "alfak2_sparse_joint_rate_fit")
}

#' Predict fitness over a range of birth or missegregation rates
#'
#' @param object An `alfak2_sparse_joint_rate_fit`.
#' @param karyotypes Character labels or an integer matrix of karyotypes.
#' @param ... Passed to the prediction method of the fits, for instance
#'   `sample` for a cohort.
#'
#' @return Data frame with one row per karyotype. `relative_fitness` is
#'   fitness per unit of `dt`, relative to the major karyotypes, at the rate
#'   chosen by `point` in [fit_sparse_joint_rates()]; `relative_low` and
#'   `relative_high` enclose its 95% intervals at every fitted rate, and
#'   `relative_sd` is the standard deviation whose symmetric 95% interval
#'   about `relative_fitness` contains them. `per_division` is relative
#'   fitness divided by the birth rate of the fit, with `per_division_low`
#'   and `per_division_high` its range over the rates. The attribute
#'   `rank_agreement` is the lowest Spearman correlation between the
#'   predictions at any two rates.
#' @export
predict.alfak2_sparse_joint_rate_fit <- function(object, karyotypes, ...) {
  predictions <- lapply(object$fits, function(fit) stats::predict(fit, karyotypes, ...))
  z <- stats::qnorm(0.975)
  mean_at <- vapply(predictions, function(p) p$relative_fitness, numeric(nrow(predictions[[1L]])))
  sd_at <- vapply(predictions, function(p) p$relative_sd, numeric(nrow(predictions[[1L]])))
  mean_at <- matrix(mean_at, ncol = length(object$rates))
  sd_at <- matrix(sd_at, ncol = length(object$rates))
  birth_rates <- if (identical(object$parameter, "beta")) {
    vapply(object$fits, function(fit) fit$diagnostics$birth_rate, numeric(1))
  } else object$rates
  per_division <- sweep(mean_at, 2L, birth_rates, "/")
  low <- apply(mean_at - z * sd_at, 1L, min)
  high <- apply(mean_at + z * sd_at, 1L, max)
  at <- object$point
  out <- data.frame(
    karyotype = predictions[[1L]]$karyotype,
    relative_fitness = mean_at[, at],
    relative_low = low, relative_high = high,
    # The point prediction is not in the middle of the interval: the
    # standard deviation is the one whose symmetric interval contains it.
    relative_sd = pmax(high - mean_at[, at], mean_at[, at] - low) / z,
    per_division = per_division[, at],
    per_division_low = apply(per_division, 1L, min), per_division_high = apply(per_division, 1L, max),
    stringsAsFactors = FALSE
  )
  agreement <- if (nrow(out) > 2L) {
    correlation <- suppressWarnings(stats::cor(mean_at, method = "spearman"))
    min(correlation[upper.tri(correlation)], na.rm = TRUE)
  } else NA_real_
  attr(out, "rank_agreement") <- agreement
  out
}

#' @export
print.alfak2_sparse_joint_rate_fit <- function(x, ...) {
  cat("<alfak2_sparse_joint_rate_fit>\n")
  label <- if (identical(x$parameter, "beta")) "missegregation rates:" else "baseline birth rates:"
  cat("  ", label, " ", paste(signif(x$rates, 3), collapse = ", "),
      "  point predictions at ", signif(x$rates[[x$point]], 3), "\n", sep = "")
  invisible(x)
}
