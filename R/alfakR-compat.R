# Outputs in the format of alfakR, so that analyses written for alfakR run on
# alfak2 fits.
#
# alfakR::alfak() (branch dev) returns its cross-validation R2R invisibly and
# writes four files: bootstrap_res.Rds (fitness of the frequent karyotypes and
# their one-step neighbours, one row per bootstrap replicate),
# landscape.Rds (mean, median and sd of fitness over those karyotypes and
# their neighbours' neighbours), landscape_posterior_samples.Rds (one column
# per bootstrap replicate of that landscape) and xval.Rds. Here the replicates
# are draws from the posterior of the sparse joint model, the landscape is
# its prediction, and the cross-validation holds out the counts of each
# frequent karyotype and its neighbours in turn.

# One-step neighbours of `ids` in alfakR's order (gen_all_neighbours()): for
# each id, chromosome by chromosome, one copy fewer then one more; first
# occurrence kept, the ids themselves and nullisomic karyotypes removed.
alfakR_neighbours <- function(ids) {
  if (!length(ids)) return(character())
  k <- parse_karyotypes(ids)
  n_chr <- ncol(k)
  steps <- matrix(0L, 2L * n_chr, n_chr)
  steps[cbind(seq(1L, 2L * n_chr, 2L), seq_len(n_chr))] <- -1L
  steps[cbind(seq(2L, 2L * n_chr, 2L), seq_len(n_chr))] <- 1L
  moved <- k[rep(seq_len(nrow(k)), each = 2L * n_chr), , drop = FALSE] +
    steps[rep(seq_len(2L * n_chr), nrow(k)), , drop = FALSE]
  moved <- moved[rowSums(moved < 1L) == 0L, , drop = FALSE]
  labels <- unique(format_karyotypes(moved))
  labels[!labels %in% ids]
}

# R-squared of alfakR's cross-validation: both series centred first.
alfakR_r2r <- function(observation, prediction) {
  observation <- observation - mean(observation)
  prediction <- prediction - mean(prediction)
  1 - sum((prediction - observation)^2) / sum(observation^2)
}

# Draws of relative fitness at `labels` about its posterior mean, one column
# per draw. Domain karyotypes take their log birth rates from the joint
# posterior, karyotypes predicted from the trend take its coefficients from it
# and add the residual, and observed karyotypes outside the domain vary
# independently about the observed level. The reference, the count-weighted
# mean of the majors, is drawn with them. Each karyotype then receives
# independent noise for the part of its calibrated standard deviation that the
# model's own does not explain, so that differences between karyotypes keep
# the model's correlations.
sparse_joint_relative_draws <- function(fit, solver, labels, prediction, k) {
  latent <- solver$draw(k)
  g_index <- which(solver$names == "g")
  coef_index <- which(solver$names == "coef")
  inside <- match(labels, fit$summary$karyotype)
  log_birth <- matrix(0, length(labels), k)
  rows <- which(!is.na(inside))
  log_birth[rows, ] <- latent[g_index[inside[rows]], , drop = FALSE]
  outside <- which(is.na(inside))
  if (length(outside)) {
    trend <- fit$trend
    left_out <- labels[outside] %in% fit$dropped & !is.null(trend$observed_level)
    from_trend <- outside[!left_out]
    if (length(from_trend)) {
      design <- sparse_joint_design(parse_karyotypes(labels[from_trend]), trend$centre, trend$major,
                                    trend$shell_power)$design
      log_birth[from_trend, ] <- design %*% latent[coef_index, , drop = FALSE] +
        trend$scales[["residual"]] * matrix(stats::rnorm(length(from_trend) * k), length(from_trend), k)
    }
    level <- outside[left_out]
    if (length(level)) {
      log_birth[level, ] <- trend$observed_level[["sd"]] * matrix(stats::rnorm(length(level) * k), length(level), k)
    }
  }
  fitness <- prediction$fitness_mean * log_birth
  major <- seq_len(fit$diagnostics$n_major)
  weight <- (fit$summary$count_t0 + fit$summary$count_t1)[major]
  weight <- weight / sum(weight)
  reference <- colSums(weight * fit$summary$fitness_mean[major] * latent[g_index[major], , drop = FALSE])
  relative <- sweep(fitness, 2L, reference, "-")
  extra <- sqrt(pmax(prediction$relative_sd^2 - prediction$relative_sd_model^2, 0))
  relative <- relative + extra * matrix(stats::rnorm(length(labels) * k), length(labels), k)
  list(relative = relative, latent = latent)
}

# `nboot` columns about `mean` from `deviation` (one column per pair): each
# deviation and its mirror image, interleaved, and the mean itself when
# `nboot` is odd. Each row is then scaled to standard deviation `sd`. Row
# means and medians equal `mean` and row standard deviations equal `sd`
# exactly, while correlations between rows are those of the draws.
alfakR_symmetric_columns <- function(mean, sd, deviation, nboot) {
  pairs <- nboot %/% 2L
  columns <- matrix(0, length(mean), nboot)
  if (pairs > 0L) {
    columns[, seq(1L, 2L * pairs, 2L)] <- deviation[, seq_len(pairs)]
    columns[, seq(2L, 2L * pairs, 2L)] <- -deviation[, seq_len(pairs)]
    spread <- sqrt(rowSums(columns^2) / (nboot - 1L))
    columns <- columns * ifelse(spread > 0, sd / spread, 0)
  }
  mean + columns
}

alfakR_row_summary <- function(x, f) {
  apply(x, 1L, function(v) {
    v <- v[is.finite(v)]
    if (!length(v) || (identical(f, stats::sd) && length(v) < 2L)) NA_real_ else f(v)
  })
}

# Cross-validation in alfakR's folds: each frequent karyotype with those of
# its one-step neighbours not already in an earlier fold. A fold's karyotypes
# inside the model domain are treated as never observed: their counts leave
# the likelihood, their initial abundance gets the prior of an unseen
# karyotype, and the trend's centre and major karyotypes, from which the
# distance covariate is measured, are recomputed without them. The posterior
# mode is then found again with the scales fixed, and the prediction for each
# held-out karyotype is compared with its fitness from all the data.
# Karyotypes outside the domain are not tested: their fitness never came from
# their own counts.
sparse_joint_alfakR_xval <- function(fit, fq, observation, level, divisor) {
  empty <- list(R2R = NA_real_, predictions = numeric(0), observations = numeric(0),
                xval_data = data.frame(k = character(0), observation = numeric(0), prediction = numeric(0),
                                       stringsAsFactors = FALSE),
                tmp = matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("test_f", "est_f"))))
  if (!length(fq)) return(empty)
  summary <- fit$summary
  n <- nrow(summary)
  n_major <- fit$diagnostics$n_major
  fold <- unlist(lapply(seq_along(fq), function(i) {
    members <- c(fq[[i]], alfakR_neighbours(fq[[i]]))
    stats::setNames(rep(i, length(members)), members)
  }))
  fold <- fold[!duplicated(names(fold))]
  problem <- sparse_joint_problem(fit$data, fit$settings)
  never_seen <- log(sum(problem$tmb_data$y0)) + fit$settings$initial_frequency[[1L]]
  karyotypes <- problem$domain$karyotypes
  self <- 2 * exp(problem$kernel$log_self) - 1
  weight_major <- (summary$count_t0 + summary$count_t1)[seq_len(n_major)]
  reference <- fit$trend$reference$fitness
  pairs <- lapply(unique(fold), function(i) {
    members <- names(fold)[fold == i]
    state <- match(members, summary$karyotype)
    tested <- members[!is.na(state)]
    state <- state[!is.na(state)]
    if (!length(state)) return(NULL)
    weight <- rep(1, n)
    weight[state] <- 0
    prior_mean <- problem$tmb_data$zeta_prior_mean
    prior_sd <- problem$tmb_data$zeta_prior_sd
    prior_mean[state] <- never_seen
    prior_sd[state] <- fit$settings$initial_frequency[[2L]]
    # A fold holding every major karyotype leaves nothing to predict from.
    major <- setdiff(seq_len(n_major), state)
    if (!length(major)) return(NULL)
    total <- problem$total * weight
    design <- sparse_joint_design(karyotypes, colSums(karyotypes * total) / sum(total),
                                  karyotypes[major, , drop = FALSE], fit$settings$shell_power)$design
    mode <- tryCatch(sparse_joint_rebuild(fit, weight, prior_mean, prior_sd, design = design, posterior = FALSE),
                     error = function(e) NULL)
    if (is.null(mode)) return(NULL)
    held_out <- self * exp(mode$g)
    # Anchor the refit on the majors outside the fold, whose relative fitness
    # it shares with the full fit up to a constant.
    w <- weight_major[major] / sum(weight_major[major])
    shift <- sum(w * (held_out[major] - summary$fitness_mean[major]))
    relative <- held_out[state] - shift - reference
    cbind(test_f = observation[tested], est_f = (relative + level) / divisor[tested])
  })
  tmp <- do.call(rbind, pairs)
  if (is.null(tmp)) return(empty)
  tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]
  observations <- tmp[, "test_f"]
  predictions <- tmp[, "est_f"]
  list(R2R = if (nrow(tmp) >= 2L) alfakR_r2r(observations, predictions) else NA_real_,
       predictions = predictions, observations = observations,
       xval_data = data.frame(k = rownames(tmp), observation = unname(observations),
                              prediction = unname(predictions), stringsAsFactors = FALSE),
       tmp = tmp)
}

#' Convert a sparse joint fit to the outputs of alfakR
#'
#' Produces, from a fit of [fit_sparse_joint()], the objects that
#' `alfakR::alfak()` (branch dev) saves, with the same names, columns and
#' layout, so that code written for alfakR results can read them.
#' [write_alfakR_outputs()] saves them under alfakR's file names, and
#' [alfak_compat()] does both in one call with alfakR's arguments.
#'
#' The karyotypes are chosen as alfakR chooses them. The frequent karyotypes
#' (`fq`) have at least `minobs` cells over both samples; `nn` are their
#' one-step neighbours without nullisomy; the landscape covers both and the
#' one-step neighbours of `nn`. `minobs` affects only this choice. The fit
#' itself uses every cell, with its own major karyotypes.
#'
#' Where alfakR has bootstrap replicates, these are draws from the posterior
#' of the model: jointly for all karyotypes, so that a difference between two
#' karyotypes has the uncertainty the model gives it, and widened to the
#' calibrated standard deviation of [predict.alfak2_sparse_joint_fit()]. They
#' come in mirrored pairs. As a result the landscape's `mean` and `median` are
#' the model's prediction and its `sd` the calibrated standard deviation,
#' while each is also the summary of `landscape_posterior_samples`, as in
#' alfakR. `initial_fitness` equals `final_fitness`: the model has no
#' separate first estimate.
#'
#' The frequencies of the frequent karyotypes follow alfakR's convention:
#' they are the frequencies at time 0 of alfakR's time axis, from which
#' alfakR projects to time `t` as `x0 * exp(f * t)`. The model estimates them
#' at the first sample, so each draw is moved back from `first_time`, the
#' time of the first sample on that axis, with the draw's `final_fitness`;
#' projecting the reported frequencies to `first_time` in alfakR's way gives
#' back the model's frequencies at the first sample. With `first_time = 0`
#' the two coincide.
#'
#' The scale is alfakR's. Fitness is the net growth rate per unit of time,
#' shifted so that the count-weighted mean of the major karyotypes is
#' `log(nb / n0) / dt`, the growth rate alfakR assigns to the population.
#' alfakR weights by the estimated initial frequencies of the frequent
#' karyotypes instead. With `n0 = NULL` and `nb = NULL` the model's own level
#' at its assumed birth rate is kept. `correct_efflux = TRUE` divides by
#' `2 P_ii - 1`, the share of divisions that return the parent's karyotype,
#' and so reports birth rates, as alfakR's efflux correction does.
#'
#' The cross-validation uses alfakR's folds: each frequent karyotype with its
#' one-step neighbours, each karyotype in the first fold that contains it. The
#' counts of a fold are held out and treated as never observed, the posterior
#' mode is found again with the fitted scales, and the predictions for the
#' fold are compared with the fitness estimated from all the data. Held-out
#' karyotypes also leave the trend's reference points: the distance
#' covariate is measured to the remaining major karyotypes. Only karyotypes
#' inside the model domain are tested, and a fold holding every major
#' karyotype is skipped. `R2R` is alfakR's statistic on these pairs.
#'
#' @param fit An `alfak2_sparse_joint_fit`, or an
#'   `alfak2_sparse_joint_rate_fit`, whose point fit is then used.
#' @param minobs Minimum number of cells, over both samples, of a frequent
#'   karyotype.
#' @param nboot Number of posterior draws, alfakR's bootstrap replicates.
#' @param n0,nb Population sizes that set the scale as in alfakR, or both
#'   `NULL` to keep the model's level.
#' @param correct_efflux Logical; report birth rates instead of net growth
#'   rates.
#' @param xval Logical; run the cross-validation. With `FALSE`, `xval` holds
#'   no pairs and `R2R` is `NA`.
#' @param landscape_data_output Logical; also return `landscape_data`, the
#'   counterpart of alfakR's Kriging models: a list with elements `mean` and
#'   `median`, each an object whose [predict()] method returns the fitness of
#'   any karyotypes on the same scale as `landscape`.
#' @param first_time Time of the first sample on alfakR's time axis (the
#'   first passage time, or column name times `dt`); sets the time 0 of the
#'   reported frequencies.
#'
#' @return An `alfak2_alfakR_output` list with `bootstrap_res`, `landscape`,
#'   `landscape_posterior_samples`, `xval` (`R2R`, `predictions`,
#'   `observations`, `xval_data`, and `tmp`, the paired matrix of earlier
#'   alfakR versions), `landscape_data` (`NULL` unless requested) and `fit`.
#' @export
as_alfakR <- function(fit, minobs = 20, nboot = 45, n0 = 1e5, nb = 1e7, correct_efflux = FALSE,
                      xval = TRUE, landscape_data_output = FALSE, first_time = 0) {
  if (inherits(fit, "alfak2_sparse_joint_rate_fit")) fit <- fit$fits[[fit$point]]
  if (!inherits(fit, "alfak2_sparse_joint_fit")) {
    stop("`fit` must come from fit_sparse_joint() or fit_sparse_joint_rates().", call. = FALSE)
  }
  validate_scalar(minobs, "minobs", lower = 0)
  validate_scalar(nboot, "nboot", lower = 1)
  validate_scalar(first_time, "first_time")
  if (nboot != round(nboot)) stop("`nboot` must be a whole number.", call. = FALSE)
  nboot <- as.integer(nboot)
  if (xor(is.null(n0), is.null(nb))) stop("Give both `n0` and `nb`, or neither.", call. = FALSE)
  if (!is.null(n0)) {
    validate_scalar(n0, "n0", lower = .Machine$double.eps)
    validate_scalar(nb, "nb", lower = .Machine$double.eps)
  }
  for (flag in c("correct_efflux", "xval", "landscape_data_output")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", flag), call. = FALSE)
    }
  }
  counts <- fit$data$counts
  fq <- rownames(counts)[rowSums(counts) >= minobs]
  if (!length(fq)) {
    stop(sprintf("no frequent karyotypes detected for minobs = %s; frequent karyotypes require rowSums(x) >= minobs",
                 format(minobs, trim = TRUE)), call. = FALSE)
  }
  nn <- alfakR_neighbours(fq)
  k <- unique(c(fq, nn, alfakR_neighbours(nn)))
  prediction <- stats::predict(fit, k)
  level <- if (is.null(n0)) fit$trend$reference$fitness else log(nb / n0) / fit$data$dt
  divisor <- if (correct_efflux) 2 * sparse_joint_self_probability(parse_karyotypes(k), fit$data$beta) - 1 else
    rep(1, length(k))
  names(divisor) <- k

  solver <- sparse_joint_rebuild(fit)$solver
  draws <- sparse_joint_relative_draws(fit, solver, k, prediction, max(nboot %/% 2L, 1L))
  samples <- alfakR_symmetric_columns(prediction$relative_fitness + level, prediction$relative_sd,
                                      draws$relative, nboot) / divisor
  landscape <- data.frame(k = k,
                          mean = alfakR_row_summary(samples, mean),
                          median = alfakR_row_summary(samples, stats::median),
                          sd = alfakR_row_summary(samples, stats::sd),
                          fq = k %in% fq, nn = k %in% nn, stringsAsFactors = FALSE)

  final_fitness <- t(samples[match(fq, k), , drop = FALSE])
  colnames(final_fitness) <- fq

  # Frequencies of the frequent karyotypes among themselves: at the first
  # sample, drawn with the fitness for those in the domain, the observed count
  # for the others; then moved back to time 0 with each draw's fitness.
  zeta_index <- which(solver$names == "zeta")
  inside <- match(fq, fit$summary$karyotype)
  log_count <- matrix(log(pmax(counts[fq, 1L], .Machine$double.xmin)), length(fq), nboot)
  rows <- which(!is.na(inside))
  if (length(rows)) {
    mode <- fit$posterior$zeta[inside[rows]]
    deviation <- draws$latent[zeta_index[inside[rows]], , drop = FALSE]
    pairs <- nboot %/% 2L
    shift <- matrix(0, length(rows), nboot)
    if (pairs > 0L) {
      shift[, seq(1L, 2L * pairs, 2L)] <- deviation[, seq_len(pairs)]
      shift[, seq(2L, 2L * pairs, 2L)] <- -deviation[, seq_len(pairs)]
    }
    log_count[rows, ] <- mode + shift
  }
  log_count <- log_count - first_time * t(final_fitness)
  frequency <- t(exp(sweep(log_count, 2L, apply(log_count, 2L, max), "-")))
  frequency <- frequency / rowSums(frequency)
  colnames(frequency) <- fq

  nn_fitness <- if (length(nn)) t(samples[match(nn, k), , drop = FALSE]) else matrix(numeric(0), nboot, 0L)
  if (length(nn)) colnames(nn_fitness) <- nn
  bootstrap_res <- list(initial_fitness = final_fitness, final_fitness = final_fitness,
                        initial_frequencies = frequency, final_frequencies = frequency, nn_fitness = nn_fitness)
  attr(bootstrap_res, "bootstrap_retry_log") <- data.frame(
    target_replicate = integer(0), attempt = integer(0), error_message = character(0),
    warning_message = character(0), stringsAsFactors = FALSE
  )
  attr(bootstrap_res, "bootstrap_retry_summary") <- list(
    requested_nboot = nboot, successful_replicates = nboot, total_attempts = nboot,
    failed_attempts = 0L, max_attempts = nboot * 20L
  )

  observation <- stats::setNames(landscape$mean, k)
  cross_validation <- if (xval) {
    sparse_joint_alfakR_xval(fit, fq, observation, level, divisor)
  } else {
    sparse_joint_alfakR_xval(fit, character(0), observation, level, divisor)
  }

  landscape_data <- if (landscape_data_output) {
    model <- structure(list(fit = fit, level = level, correct_efflux = correct_efflux),
                       class = "alfak2_alfakR_landscape")
    list(mean = model, median = model)
  } else NULL

  structure(list(bootstrap_res = bootstrap_res, landscape = landscape, landscape_posterior_samples = unname(samples),
                 xval = cross_validation, landscape_data = landscape_data, fit = fit),
            class = "alfak2_alfakR_output")
}

#' Save alfakR outputs
#'
#' Writes the objects of [as_alfakR()] under the file names of
#' `alfakR::alfak()`: `bootstrap_res.Rds`, `landscape.Rds`,
#' `landscape_posterior_samples.Rds`, `xval.Rds` and, when present,
#' `landscape_data.Rds`. The alfak2 fit is saved beside them as
#' `alfak2_fit.Rds`.
#'
#' @param x An `alfak2_alfakR_output` from [as_alfakR()].
#' @param outdir Directory, created if needed.
#'
#' @return The paths written, invisibly.
#' @export
write_alfakR_outputs <- function(x, outdir) {
  if (!inherits(x, "alfak2_alfakR_output")) stop("`x` must come from as_alfakR().", call. = FALSE)
  if (!is.character(outdir) || length(outdir) != 1L || !nzchar(outdir)) {
    stop("`outdir` must be one directory path.", call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  objects <- list(bootstrap_res.Rds = x$bootstrap_res, landscape.Rds = x$landscape,
                  landscape_posterior_samples.Rds = x$landscape_posterior_samples, xval.Rds = x$xval)
  if (!is.null(x$landscape_data)) objects$landscape_data.Rds <- x$landscape_data
  objects$alfak2_fit.Rds <- x$fit
  paths <- file.path(outdir, names(objects))
  for (i in seq_along(objects)) saveRDS(objects[[i]], file = paths[[i]])
  invisible(paths)
}

#' Fit a landscape with the arguments and outputs of alfakR
#'
#' A replacement for `alfakR::alfak()` (branch dev) that fits the sparse
#' joint model, [fit_sparse_joint()], instead. It takes the same input,
#' writes the same files to `outdir` (see [as_alfakR()] for how each is
#' obtained) and returns the cross-validation R-squared invisibly, so that a
#' pipeline built on alfakR can switch by calling this function. The fit
#' itself is also saved, as `alfak2_fit.Rds`.
#'
#' alfak2 fits two samples. `yi$x` must have two columns; their times are
#' `passage_times` or, without `passage_times`,
#' `as.numeric(colnames(yi$x)) * yi$dt`, as in alfakR. The fit uses the time
#' between them; the time of the first sets the time 0 of the reported
#' frequencies (`first_time` of [as_alfakR()]). `pm` is the
#' missegregation rate, `beta` of [fit_sparse_joint()]. alfakR's arguments
#' for its own neighbour prior and Kriging (`nn_prior`, `nn_prior_sd`,
#' `nn_prior_sd_floor`, `krig_bootstrap_mode`) have no counterpart; giving
#' one is an error.
#'
#' The model needs a baseline birth rate per unit of time, `birth_rate`
#' (default 1), passed through `...` like any other argument of
#' [fit_sparse_joint()]. It is not the population growth `log(nb / n0)`,
#' which only sets the level of the reported fitness.
#'
#' @param yi List with `x`, a count matrix with karyotype row names and one
#'   column per sample, and `dt`, the time per unit of the column names.
#' @param outdir Directory for the output files.
#' @param passage_times Optional times of the two samples.
#' @param minobs,nboot,n0,nb,correct_efflux,landscape_data_output See
#'   [as_alfakR()].
#' @param pm Missegregation rate per chromosome copy and division.
#' @param allow_noninteger_counts Logical; round non-integer counts with a
#'   warning instead of stopping.
#' @param ... Further arguments of [fit_sparse_joint()], such as
#'   `birth_rate`, `min_cn` and `max_cn`.
#'
#' @return The cross-validation R-squared, invisibly.
#' @export
alfak_compat <- function(yi, outdir, passage_times = NULL, minobs = 20, nboot = 45, n0 = 1e5, nb = 1e7,
                  pm = 0.00005, allow_noninteger_counts = FALSE, correct_efflux = FALSE,
                  landscape_data_output = FALSE, ...) {
  extra <- names(list(...))
  unsupported <- intersect(extra, c("nn_prior", "nn_prior_sd", "nn_prior_sd_floor", "krig_bootstrap_mode"))
  if (length(unsupported)) {
    stop(sprintf("alfak2 has no counterpart of alfakR's %s.", paste0("`", unsupported, "`", collapse = ", ")),
         call. = FALSE)
  }
  if (!is.list(yi) || is.null(yi$x)) stop("`yi` must be a list containing a count matrix `x`.", call. = FALSE)
  x <- yi$x
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) stop("`yi$x` must be a numeric matrix.", call. = FALSE)
  if (ncol(x) != 2L) {
    stop(sprintf("alfak2 fits two samples; `yi$x` has %d columns.", ncol(x)), call. = FALSE)
  }
  if (any(abs(x - round(x)) > sqrt(.Machine$double.eps), na.rm = TRUE)) {
    if (!isTRUE(allow_noninteger_counts)) {
      stop("`yi$x` contains non-integer counts; set `allow_noninteger_counts = TRUE` to round them.", call. = FALSE)
    }
    warning("Non-integer counts in `yi$x` were rounded.", call. = FALSE)
    x <- round(x)
  }
  times <- if (is.null(passage_times)) {
    if (is.null(yi$dt) || !is.numeric(yi$dt) || length(yi$dt) != 1L || !is.finite(yi$dt)) {
      stop("`yi$dt` must be a single finite number when `passage_times` is NULL.", call. = FALSE)
    }
    axis <- suppressWarnings(as.numeric(colnames(x)))
    if (length(axis) != 2L || any(!is.finite(axis))) {
      stop("Without `passage_times`, `colnames(yi$x)` must be numeric.", call. = FALSE)
    }
    axis * yi$dt
  } else {
    if (!is.numeric(passage_times) || length(passage_times) != 2L || any(!is.finite(passage_times))) {
      stop("`passage_times` must give two finite times.", call. = FALSE)
    }
    passage_times
  }
  interval <- diff(as.numeric(times))
  if (interval <= 0) stop("The second sample must come after the first.", call. = FALSE)
  fit <- fit_sparse_joint(x, dt = interval, beta = pm, ...)
  out <- as_alfakR(fit, minobs = minobs, nboot = nboot, n0 = n0, nb = nb, correct_efflux = correct_efflux,
                   landscape_data_output = landscape_data_output, first_time = as.numeric(times)[1L])
  write_alfakR_outputs(out, outdir)
  invisible(out$xval$R2R)
}

#' Fitness from an alfakR-style landscape model
#'
#' @param object An element of `landscape_data` from [as_alfakR()].
#' @param x Karyotypes, as a matrix with one row per karyotype or as labels.
#' @param ... Unused.
#'
#' @return A one-column matrix of fitness on the scale of the landscape.
#' @export
predict.alfak2_alfakR_landscape <- function(object, x, ...) {
  labels <- if (is.character(x)) x else {
    karyotypes <- round(as.matrix(x))
    storage.mode(karyotypes) <- "integer"
    format_karyotypes(karyotypes)
  }
  value <- stats::predict(object$fit, labels)$relative_fitness + object$level
  if (isTRUE(object$correct_efflux)) {
    value <- value / (2 * sparse_joint_self_probability(parse_karyotypes(labels), object$fit$data$beta) - 1)
  }
  matrix(value, ncol = 1L)
}

#' @export
print.alfak2_alfakR_output <- function(x, ...) {
  cat("<alfak2_alfakR_output>\n")
  cat("  frequent karyotypes:", sum(x$landscape$fq), " neighbours:", sum(x$landscape$nn),
      " landscape rows:", nrow(x$landscape), "\n")
  cat("  posterior draws:", ncol(x$landscape_posterior_samples), " xval R2R:", format(x$xval$R2R), "\n")
  invisible(x)
}
