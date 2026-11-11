# The sparse joint model when the baseline birth rate is not known.
#
# The rate cannot be learned from two samples and has to be supplied; neither
# a single sample's nor a cohort's marginal likelihood locates it. What it
# changes is known: the fitness difference between two well-sampled clones
# comes from their frequency change and does not depend on it, while the
# deficit of a rare or unobserved karyotype comes from its level under
# mutation-selection balance, which fixes the deficit per division. Per unit
# of time that deficit is proportional to the assumed rate. Rankings barely
# move. Fitting over a range of rates shows all of this for a given dataset.

#' Fit the sparse joint model over a range of baseline birth rates
#'
#' Use this when the baseline birth rate is unknown and only a plausible range
#' can be given. The model is fitted at each of a few rates across the range.
#' [stats::predict()] then reports relative fitness at one of them, an
#' interval that covers every rate in the range, and relative fitness per
#' division.
#'
#' The two directions of error are not alike. Assuming too low a rate shrinks
#' the predicted deficits of rare and unobserved karyotypes towards zero; on
#' simulated landscapes a rate 4 to 10 times too low still predicts better
#' than no differences. Assuming too high a rate inflates them, and a rate 4
#' times too high predicts worse than no differences. Rankings change little
#' either way. The point prediction is therefore taken at the lowest rate of
#' the range by default, which is never too high when the range holds the
#' true rate.
#'
#' @param ... Arguments of the fitting function other than `birth_rate`:
#'   for [fit_sparse_joint()], `counts`, `dt`, `beta` and so on.
#' @param birth_rate_range Lowest and highest plausible baseline birth rate
#'   per unit of `dt`.
#' @param n_rates Number of rates fitted, spaced evenly on the log scale. An
#'   odd number puts one at the centre of the range.
#' @param point Rate of the point prediction: `"lowest"` of the range, or its
#'   geometric `"centre"`, which is more accurate when the true rate is near
#'   the centre and much less accurate when it is near the lowest.
#' @param fit_function [fit_sparse_joint()] or [fit_sparse_joint_cohort()].
#'
#' @return An `alfak2_sparse_joint_rate_fit` object: `rates`, the `fits` at
#'   each of them, `point`, the index of the fit used for point predictions,
#'   and `central`, the index of the fit at the centre.
#' @export
fit_sparse_joint_rates <- function(..., birth_rate_range, n_rates = 3L, point = c("lowest", "centre"),
                                   fit_function = fit_sparse_joint) {
  if (!is.numeric(birth_rate_range) || length(birth_rate_range) != 2L || any(!is.finite(birth_rate_range)) ||
      birth_rate_range[[1L]] <= 0 || birth_rate_range[[2L]] <= birth_rate_range[[1L]]) {
    stop("`birth_rate_range` must give a positive lowest and a larger highest rate.", call. = FALSE)
  }
  validate_scalar(n_rates, "n_rates", lower = 2)
  point <- match.arg(point)
  if ("birth_rate" %in% names(list(...))) {
    stop("Give `birth_rate_range`, not `birth_rate`.", call. = FALSE)
  }
  rates <- exp(seq(log(birth_rate_range[[1L]]), log(birth_rate_range[[2L]]), length.out = as.integer(n_rates)))
  fits <- lapply(rates, function(rate) fit_function(..., birth_rate = rate))
  central <- (length(rates) + 1L) %/% 2L
  structure(list(rates = rates, fits = fits, point = if (point == "lowest") 1L else central, central = central),
            class = "alfak2_sparse_joint_rate_fit")
}

#' Predict fitness over a range of baseline birth rates
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
#'   fitness divided by the assumed rate, with `per_division_low` and
#'   `per_division_high` its range over the rates. The attribute
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
  per_division <- sweep(mean_at, 2L, object$rates, "/")
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
  cat("  baseline birth rates:", paste(signif(x$rates, 3), collapse = ", "),
      " point predictions at", signif(x$rates[[x$point]], 3), "\n")
  invisible(x)
}
