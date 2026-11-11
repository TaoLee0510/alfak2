# Sparse joint count model: domain, design, fit and prediction.
#
# The model is anchored on the well-sampled ("major") karyotypes. Its domain
# is the majors and every one-step neighbour of a major, observed or not, so
# that a neighbour seen in one cell is judged against all its unseen siblings.
# Fitness further away is predicted from the fitted trend.

one_step_neighbours <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- as.matrix(karyotypes)
  n_chr <- ncol(karyotypes)
  out <- vector("list", 2L * n_chr)
  for (chr in seq_len(n_chr)) {
    for (direction in c(-1L, 1L)) {
      moved <- karyotypes
      moved[, chr] <- moved[, chr] + direction
      keep <- moved[, chr] >= min_cn & moved[, chr] <= max_cn
      out[[2L * (chr - 1L) + (direction + 3L) %/% 2L]] <- moved[keep, , drop = FALSE]
    }
  }
  out <- do.call(rbind, out)
  storage.mode(out) <- "integer"
  out[!duplicated(format_karyotypes(out)), , drop = FALSE]
}

# Major karyotypes anchor the domain: at least `major_frequency` of either
# sample and at least `major_count` cells. The share adapts the rule to depth:
# one cell in 3000 may be a member of the mutational cloud around a clone.
# The count is needed as well. A major's neighbours enter the likelihood
# whether observed or not, and an unobserved neighbour is evidence only if a
# cell falling on it would have left it a neighbour. Were a single cell enough
# to make a karyotype major, all cells of a small sample would sit on majors
# by construction and every unobserved neighbour would look strongly unfit.
# Without any such karyotype the most abundant ones are used. If the domain
# would exceed `max_states`, the share is raised to the lowest value that fits.
sparse_joint_major <- function(counts, major_frequency, major_count = 3, max_states = Inf,
                               min_cn = 1L, max_cn = 8L) {
  share <- pmax(counts[, 1L] / sum(counts[, 1L]), counts[, 2L] / sum(counts[, 2L]))
  total <- rowSums(counts)
  labels <- rownames(counts)
  eligible <- share >= major_frequency & total >= major_count
  if (!any(eligible)) eligible <- total == max(total)
  tiers <- sort(unique(share[eligible]))
  domain_size <- function(threshold) {
    major <- labels[eligible & share >= threshold]
    length(union(major, format_karyotypes(one_step_neighbours(parse_karyotypes(major), min_cn, max_cn))))
  }
  tier <- 1L
  while (tier < length(tiers) && domain_size(tiers[tier]) > max_states) tier <- tier + 1L
  labels[eligible & share >= tiers[tier]]
}

# Lattice distance from each karyotype to the nearest major karyotype.
sparse_joint_major_distance <- function(karyotypes, major) {
  out <- rep(Inf, nrow(karyotypes))
  for (m in seq_len(nrow(major))) {
    out <- pmin(out, rowSums(abs(sweep(karyotypes, 2L, major[m, ], "-"))))
  }
  out
}

# Domain: the major karyotypes, listed first, and every one-step neighbour of
# one. Observed karyotypes outside it are left out with their cells.
sparse_joint_domain <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn) {
  observed <- rownames(counts)
  major_labels <- sparse_joint_major(counts, major_frequency, major_count, max_states, min_cn, max_cn)
  major <- parse_karyotypes(major_labels)
  labels <- union(major_labels, format_karyotypes(one_step_neighbours(major, min_cn, max_cn)))
  in_data <- match(labels, observed)
  y <- matrix(0, length(labels), 2L)
  y[!is.na(in_data), ] <- counts[in_data[!is.na(in_data)], , drop = FALSE]
  outside <- !observed %in% labels
  list(labels = labels, karyotypes = parse_karyotypes(labels), major = major,
       n_major = length(major_labels), y0 = y[, 1L], y1 = y[, 2L],
       dropped = observed[outside], dropped_cells = sum(counts[outside, ]))
}

# Trend design: intercept, chromosome dosage and curvature around the
# population's mean karyotype, and a shell term. The shell term is the
# distance to the nearest major karyotype raised to `shell_power`. It is
# measured to majors only: whether a neighbour was observed must not enter its
# own covariates, because its count is what the likelihood explains.
sparse_joint_design <- function(karyotypes, centre, major, shell_power) {
  dx <- sweep(karyotypes, 2L, centre, "-")
  design <- cbind(1, dx, dx^2, sparse_joint_major_distance(karyotypes, major)^shell_power)
  list(design = unname(design),
       group = c(0L, rep(1L, ncol(dx)), rep(2L, ncol(dx)), 4L))
}

# Probability that a daughter has exactly its parent's karyotype.
sparse_joint_self_probability <- function(karyotypes, beta) {
  stay <- vapply(0:max(karyotypes), function(cn) alfak2_pij_cpp(cn, cn, beta), numeric(1))
  exp(rowSums(matrix(log(stay)[karyotypes + 1L], nrow = nrow(karyotypes))))
}

#' Fit the sparse joint count model
#'
#' Fits fitness to two sparse karyotype samples with a single count likelihood
#' on a sparse domain: the well-sampled ("major") karyotypes and all their
#' one-step neighbours. Division dynamics use a closed form in which each
#' karyotype grows on its own and receives mutational influx from its
#' neighbours; an unobserved neighbour of an abundant karyotype is thereby
#' evidence that the neighbour is less fit. Fitness away from the domain is
#' predicted with [predict()] from a fitted trend in chromosome dosage,
#' curvature and distance to the nearest major karyotype.
#'
#' Observed karyotypes that are not next to a major karyotype are left out of
#' the likelihood, together with their cells. They would enter only because
#' they were observed, while their unobserved siblings would not, which
#' biases their fitness upwards. [predict()] places them at the level of the
#' observed neighbours inside the domain.
#'
#' The baseline birth rate is supplied rather than estimated. Two samples
#' identify it only weakly, and a free baseline drifts upwards and inflates
#' every fitness difference. Rerun with other values to see how much a
#' conclusion depends on it.
#'
#' @param counts Two-column count matrix with karyotype row names, or an
#'   `alfak2_data` object.
#' @param dt Time between the samples, if `counts` is a matrix.
#' @param beta Missegregation rate per chromosome copy and division, if
#'   `counts` is a matrix.
#' @param birth_rate Assumed baseline birth (division) rate per unit of `dt`.
#' @param birth_rate_log_sd Prior standard deviation of the log baseline birth
#'   rate around `log(birth_rate)`.
#' @param major_frequency Minimum share of either sample for a karyotype to
#'   anchor the domain as a major karyotype.
#' @param major_count Minimum number of cells, over both samples, for a major
#'   karyotype. A single cell must not be enough.
#' @param max_states Largest domain to fit. If the domain would be larger,
#'   the share required of a major karyotype is raised until it fits.
#' @param min_cn,max_cn Copy-number bounds of viable karyotypes. `max_cn`
#'   defaults to one above the largest observed copy number.
#' @param shell_power Exponent of the distance to the nearest major karyotype
#'   in the trend. Mean fitness falls slightly less than linearly with that
#'   distance.
#' @param shell_sd Prior standard deviation of the shell coefficient.
#' @param initial_frequency Prior mean and standard deviation of the log
#'   frequency at the first timepoint of a karyotype seen at neither
#'   timepoint. It is a frequency, not a count: a domain has thousands of
#'   unseen karyotypes, and a prior in counts would give them a large share of
#'   a small sample. A karyotype seen only at the second timepoint is centred
#'   on its frequency there, with the same standard deviation: absent other
#'   evidence its frequency is taken not to have changed.
#' @param control List passed to [stats::nlminb()] for the scale parameters.
#' @param inner_control TMB inner Newton controls.
#'
#' @return An `alfak2_sparse_joint_fit` object. `summary` has one row per
#'   domain karyotype; `fitness_mean` is the net growth rate per unit of `dt`.
#' @export
fit_sparse_joint <- function(counts,
                             dt = 1,
                             beta = 0.00005,
                             birth_rate = 1,
                             birth_rate_log_sd = 0.05,
                             major_frequency = 0.01,
                             major_count = 3,
                             max_states = 5000,
                             min_cn = 1,
                             max_cn = NULL,
                             shell_power = 0.85,
                             shell_sd = 0.5,
                             initial_frequency = c(mean = log(1e-5), sd = 3),
                             control = list(eval.max = 200, iter.max = 100),
                             inner_control = list(maxit = 200, tol = 1e-10)) {
  data <- if (inherits(counts, "alfak2_data")) counts else prepare_alfak2_data(counts, dt = dt, beta = beta)
  validate_scalar(birth_rate, "birth_rate", lower = .Machine$double.eps)
  validate_scalar(birth_rate_log_sd, "birth_rate_log_sd", lower = .Machine$double.eps)
  validate_scalar(major_frequency, "major_frequency", lower = 0, upper = 1)
  validate_scalar(major_count, "major_count", lower = 2)
  validate_scalar(max_states, "max_states", lower = 1)
  validate_scalar(min_cn, "min_cn", lower = 0)
  validate_scalar(shell_power, "shell_power", lower = .Machine$double.eps)
  validate_scalar(shell_sd, "shell_sd", lower = .Machine$double.eps)
  if (!is.numeric(initial_frequency) || length(initial_frequency) != 2L ||
      any(!is.finite(initial_frequency)) || initial_frequency[[2L]] <= 0) {
    stop("`initial_frequency` must give a finite mean and a positive standard deviation.", call. = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  counts <- data$counts
  observed_cn <- parse_karyotypes(rownames(counts))
  if (is.null(max_cn)) max_cn <- max(observed_cn) + 1L
  validate_scalar(max_cn, "max_cn", lower = min_cn)
  if (any(observed_cn < min_cn | observed_cn > max_cn)) {
    stop("Observed karyotypes must lie within the copy-number bounds.", call. = FALSE)
  }
  if (any(colSums(counts) <= 0)) stop("Both timepoints need at least one cell.", call. = FALSE)

  domain <- sparse_joint_domain(counts, major_frequency, major_count, max_states,
                                as.integer(min_cn), as.integer(max_cn))
  n <- length(domain$labels)
  if (any(c(sum(domain$y0), sum(domain$y1)) <= 0)) {
    stop("A timepoint has no cell inside the model domain.", call. = FALSE)
  }
  kernel <- sparse_division_kernel(domain$karyotypes, data$beta)
  total <- domain$y0 + domain$y1
  centre <- colSums(domain$karyotypes * total) / sum(total)
  design <- sparse_joint_design(domain$karyotypes, centre, domain$major, shell_power)
  seen0 <- domain$y0 > 0
  scale_names <- c("dosage", "curvature", "interaction", "residual")
  tmb_data <- list(
    model = "sparse_joint",
    y0 = domain$y0, y1 = domain$y1, dt = data$dt,
    edge_from = kernel$from - 1L, edge_to = kernel$to - 1L,
    edge_log_prob = kernel$log_prob, log_self = kernel$log_self,
    design = design$design, coef_group = design$group,
    reference = which.max(total) - 1L,
    # Log expected count at the first timepoint. Seen then: the count itself,
    # essentially free. Seen only later: the later frequency. Never seen: rare.
    zeta_prior_mean = ifelse(
      seen0, log(pmax(domain$y0, 1)),
      ifelse(domain$y1 > 0, log(pmax(domain$y1, 1) * sum(domain$y0) / sum(domain$y1)),
             log(sum(domain$y0)) + initial_frequency[[1L]])
    ),
    zeta_prior_sd = ifelse(seen0, 10, initial_frequency[[2L]]),
    intercept_mean = log(birth_rate), intercept_sd = birth_rate_log_sd, shell_sd = shell_sd,
    # Centred on a 1% effect per step; wide enough to reach 0.05% and 20%.
    log_scale_prior_mean = rep(log(0.01), 4L),
    log_scale_prior_sd = rep(1.5, 4L)
  )
  parameters <- list(
    zeta = tmb_data$zeta_prior_mean,
    g = rep(log(birth_rate), n),
    coef = c(log(birth_rate), rep(0, ncol(design$design) - 1L)),
    theta1 = log(sum(domain$y1) / sum(domain$y0)),
    log_scale = tmb_data$log_scale_prior_mean
  )
  # There are no interaction terms in the design, so their scale is not free.
  free_scale <- c(TRUE, TRUE, FALSE, TRUE)
  objective <- TMB::MakeADFun(
    tmb_data, parameters,
    map = list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA))),
    random = c("zeta", "g", "coef", "theta1"), DLL = "alfak2", silent = TRUE,
    inner.control = inner_control
  )
  # Scales stay in a finite box: outside it the inner problem can become
  # degenerate and return non-finite gradients.
  optimum <- tryCatch(
    suppressWarnings(stats::nlminb(objective$par, objective$fn, objective$gr,
                                   lower = log(1e-4), upper = log(0.5), control = control)),
    error = function(e) NULL
  )
  if (is.null(optimum)) {
    # The optimiser stepped to scales where the inner problem has no finite
    # gradient. The fit is still usable at the prior-centred scales.
    warning("Scale optimisation failed; the prior-centred scales are used.", call. = FALSE)
    optimum <- list(par = objective$par, objective = objective$fn(objective$par), convergence = 99L,
                    message = "scale optimisation failed; prior-centred scales used")
  }
  # Uncertainty is conditional on the fitted scales. The marginal likelihood
  # is often flat in a scale, and propagating that through a Hessian that is
  # not positive definite gives unusable standard deviations.
  sdrep <- TMB::sdreport(objective, par.fixed = optimum$par, ignore.parm.uncertainty = TRUE)
  mode <- objective$env$last.par.best
  random <- summary(sdrep, "random")
  pick <- function(name) random[rownames(random) == name, , drop = FALSE]
  g_hat <- pick("g")
  coef_hat <- pick("coef")
  report <- objective$report(mode)
  hessian <- objective$env$spHess(mode, random = TRUE)
  coef_index <- which(names(mode)[objective$env$random] == "coef")
  selector <- Matrix::sparseMatrix(i = coef_index, j = seq_along(coef_index), x = 1,
                                   dims = c(nrow(hessian), length(coef_index)))
  coef_covariance <- as.matrix(Matrix::solve(hessian, selector)[coef_index, , drop = FALSE])
  scale <- tmb_data$log_scale_prior_mean
  scale[free_scale] <- optimum$par
  scale <- stats::setNames(exp(scale), scale_names)

  net <- as.numeric(report$net)
  self <- 2 * exp(kernel$log_self) - 1
  summary <- data.frame(
    karyotype = domain$labels,
    count_t0 = domain$y0, count_t1 = domain$y1,
    major = seq_len(n) <= domain$n_major,
    shell = sparse_joint_major_distance(domain$karyotypes, domain$major),
    log_birth_mean = g_hat[, 1L], log_birth_sd = g_hat[, 2L],
    fitness_mean = net, fitness_sd = net * g_hat[, 2L],
    conf_low = self * exp(g_hat[, 1L] - stats::qnorm(0.975) * g_hat[, 2L]),
    conf_high = self * exp(g_hat[, 1L] + stats::qnorm(0.975) * g_hat[, 2L]),
    stringsAsFactors = FALSE
  )
  # Level of the observed karyotypes that are not major: used to predict the
  # observed karyotypes left out of the domain. An observed karyotype is much
  # closer in fitness to the majors than an unobserved one at the same distance.
  observed_minor <- !summary$major & (summary$count_t0 + summary$count_t1) > 0
  observed_level <- if (any(observed_minor)) {
    c(log_birth = mean(summary$log_birth_mean[observed_minor]),
      sd = sqrt(mean(summary$log_birth_sd[observed_minor]^2) +
                  if (sum(observed_minor) > 1L) stats::var(summary$log_birth_mean[observed_minor]) else 0))
  } else NULL
  structure(list(
    data = data,
    summary = summary,
    dropped = domain$dropped,
    trend = list(coefficients = coef_hat[, 1L], covariance = coef_covariance, centre = centre,
                 major = domain$major, shell_power = shell_power, scales = scale,
                 observed_level = observed_level),
    diagnostics = list(
      convergence = optimum$convergence, message = optimum$message, objective = optimum$objective,
      finite_uncertainty = all(is.finite(g_hat[, 2L])),
      n_states = n, n_major = domain$n_major, n_edges = length(kernel$from),
      dropped_karyotypes = length(domain$dropped), dropped_cells = domain$dropped_cells,
      birth_rate = birth_rate, fitted_birth_rate = exp(coef_hat[1L, 1L]),
      min_cn = as.integer(min_cn), max_cn = as.integer(max_cn),
      seconds = proc.time()[["elapsed"]] - started
    )
  ), class = "alfak2_sparse_joint_fit")
}

#' Predict fitness from a sparse joint fit
#'
#' Karyotypes inside the fitted domain return their posterior. An observed
#' karyotype that was left out of the domain is placed at the level of the
#' observed, non-major karyotypes inside it. Any other karyotype is predicted
#' from the trend: chromosome dosage and curvature around the population's
#' mean karyotype, and distance to the nearest major karyotype, with the
#' residual scale as additional uncertainty.
#'
#' @param object An `alfak2_sparse_joint_fit`.
#' @param karyotypes Character labels or an integer matrix of karyotypes.
#' @param ... Unused.
#'
#' @return Data frame with `karyotype`, `fitness_mean` (net growth rate per
#'   unit of `dt`), `fitness_sd`, `conf_low`, `conf_high`, `shell` and
#'   `in_domain`.
#' @export
predict.alfak2_sparse_joint_fit <- function(object, karyotypes, ...) {
  labels <- if (is.character(karyotypes)) karyotypes else format_karyotypes(as.matrix(karyotypes))
  inside <- match(labels, object$summary$karyotype)
  out <- data.frame(karyotype = labels, fitness_mean = NA_real_, fitness_sd = NA_real_,
                    conf_low = NA_real_, conf_high = NA_real_, shell = NA_real_,
                    in_domain = !is.na(inside), stringsAsFactors = FALSE)
  columns <- c("fitness_mean", "fitness_sd", "conf_low", "conf_high", "shell")
  out[!is.na(inside), columns] <- object$summary[inside[!is.na(inside)], columns]
  outside <- which(is.na(inside))
  if (length(outside)) {
    trend <- object$trend
    cn <- parse_karyotypes(labels[outside])
    design <- sparse_joint_design(cn, trend$centre, trend$major, trend$shell_power)$design
    log_birth <- as.numeric(design %*% trend$coefficients)
    log_birth_sd <- sqrt(rowSums((design %*% trend$covariance) * design) + trend$scales[["residual"]]^2)
    left_out <- labels[outside] %in% object$dropped
    if (any(left_out) && !is.null(trend$observed_level)) {
      log_birth[left_out] <- trend$observed_level[["log_birth"]]
      log_birth_sd[left_out] <- trend$observed_level[["sd"]]
    }
    self <- 2 * sparse_joint_self_probability(cn, object$data$beta) - 1
    out$fitness_mean[outside] <- self * exp(log_birth)
    out$fitness_sd[outside] <- out$fitness_mean[outside] * log_birth_sd
    out$conf_low[outside] <- self * exp(log_birth - stats::qnorm(0.975) * log_birth_sd)
    out$conf_high[outside] <- self * exp(log_birth + stats::qnorm(0.975) * log_birth_sd)
    out$shell[outside] <- sparse_joint_major_distance(cn, trend$major)
  }
  out
}

#' @export
print.alfak2_sparse_joint_fit <- function(x, ...) {
  d <- x$diagnostics
  cat("<alfak2_sparse_joint_fit>\n")
  cat("  states:", d$n_states, "(", d$n_major, "major )\n")
  cat("  observed karyotypes outside the domain:", d$dropped_karyotypes, "(", d$dropped_cells, "cells )\n")
  cat("  baseline birth rate:", d$birth_rate, " convergence:", d$convergence, "\n")
  invisible(x)
}
