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

# Prior on the log expected count at the first timepoint. Seen then: the
# count itself, essentially free. Seen only later: the later frequency. Never
# seen: rare, on the frequency scale given by `initial_frequency`.
sparse_joint_initial_prior <- function(y0, y1, initial_frequency) {
  seen_first <- y0 > 0
  list(
    mean = ifelse(seen_first, log(pmax(y0, 1)),
                  ifelse(y1 > 0, log(pmax(y1, 1) * sum(y0) / sum(y1)),
                         log(sum(y0)) + initial_frequency[[1L]])),
    sd = ifelse(seen_first, 10, initial_frequency[[2L]])
  )
}

# Scale parameters by Laplace marginal likelihood, kept inside a finite box:
# outside it the inner problem can become degenerate and return non-finite
# gradients. If the optimiser still steps somewhere without a finite
# gradient, the fit is usable at the starting, prior-centred scales.
sparse_joint_optimise_scales <- function(objective, control, lower = log(1e-4), upper = log(0.5)) {
  optimum <- tryCatch(
    suppressWarnings(stats::nlminb(objective$par, objective$fn, objective$gr,
                                   lower = lower, upper = upper, control = control)),
    error = function(e) NULL
  )
  if (is.null(optimum)) {
    warning("Scale optimisation failed; the prior-centred scales are used.", call. = FALSE)
    optimum <- list(par = objective$par, objective = as.numeric(objective$fn(objective$par)), convergence = 99L,
                    message = "scale optimisation failed; prior-centred scales used")
  }
  optimum
}

# Posterior covariance of the latent variables given the scales, as a solver.
#
# The covariance is the inverse of the expected information at the mode: the
# Hessian of the same objective with the counts replaced by their fitted
# expectations, and the residuals of the balance prior set to zero. The observed Hessian is not used. The intensities are not
# log-linear in the birth rates, so where the fit to the counts is poor the
# observed Hessian has directions of negative curvature, and its inverse then
# gives standard deviations far wider than the prior or none at all. The
# expected information is positive definite by construction.
#
# Returns the mode, a solver for right-hand sides, the variances of selected
# latent variables, draws of all latent variables about the mode, the names of
# the latent variables, and the fitted expected counts.
sparse_joint_posterior <- function(objective, tmb_data, parameters, map, random, par_fixed, inner_control) {
  objective$fn(par_fixed)
  mode <- objective$env$last.par
  latent <- objective$env$random
  report <- objective$report(mode)
  expected <- tmb_data
  expected$y0 <- exp(unname(mode[names(mode) == "zeta"]))
  expected$y1 <- exp(as.numeric(report$log_intensity1))
  # The balance prior on initial abundances is not log-linear either: its
  # curvature is taken at zero residual in the same way.
  expected$balance_offset <- as.numeric(report$balance_residual)
  information <- TMB::MakeADFun(expected, parameters, map = map, random = random, DLL = "alfak2",
                                silent = TRUE, inner.control = inner_control)$env$spHess(mode, random = TRUE)
  factor <- Matrix::Cholesky(information, LDL = FALSE)
  n <- nrow(information)
  solve <- function(rhs) as.matrix(Matrix::solve(factor, rhs))
  list(
    mode = unname(mode[latent]), names = names(mode)[latent], solve = solve,
    # Fitted expected counts, for comparison with the observed ones.
    expected_t0 = expected$y0, expected_t1 = expected$y1,
    # `k` draws of the latent variables from the Gaussian posterior, centred on
    # zero, as columns. With P A P' = L L', P' L^-T z has covariance A^-1.
    draw = function(k) {
      z <- matrix(stats::rnorm(n * k), n, k)
      as.matrix(Matrix::solve(factor, Matrix::solve(factor, z, system = "Lt"), system = "Pt"))
    },
    # Diagonal of the inverse at `index`, a few hundred columns at a time.
    variance = function(index) {
      out <- numeric(length(index))
      for (chunk in split(seq_along(index), ceiling(seq_along(index) / 400))) {
        rhs <- Matrix::sparseMatrix(i = index[chunk], j = seq_along(chunk), x = 1, dims = c(n, length(chunk)))
        out[chunk] <- solve(rhs)[cbind(index[chunk], seq_along(chunk))]
      }
      out
    }
  )
}

# Posterior covariance of latent variables at `coef_index`, by default all
# trend coefficients.
sparse_joint_coef_covariance <- function(solver, coef_index = which(solver$names == "coef")) {
  selector <- Matrix::sparseMatrix(i = coef_index, j = seq_along(coef_index), x = 1,
                                   dims = c(length(solver$names), length(coef_index)))
  solver$solve(selector)[coef_index, , drop = FALSE]
}

# Calibration of the reported uncertainty of relative fitness.
#
# The model's own posterior standard deviation is too small. It is conditional
# on the fitted scales, and it takes the model to be true: the dynamics are
# first order in the mutational flux, and the spread of fitness among
# unobserved karyotypes is set by the observed ones, which selection has made
# alike. On simulated populations its 95% intervals cover about three in four
# observed karyotypes and half of the unobserved ones. The reported standard
# deviation adds two terms of model error to it:
#
#   sqrt(model sd^2 + (deficit * relative fitness)^2 +
#        (floor * baseline birth rate * (1 + shell))^2).
#
# The first says that a karyotype predicted far from the reference is known
# in proportion to that distance. The second is an error per division that no
# amount of data removes, growing with `shell`, the distance to the nearest
# major karyotype. Both constants were fitted to the development landscapes
# of `benchmark/sparse_two_timepoint` (validation/sparse_joint/calibration_fit.R).
sparse_joint_default_calibration <- c(floor = 0.009, deficit = 0.6)

sparse_joint_validate_calibration <- function(calibration) {
  if (is.null(calibration)) return(c(floor = 0, deficit = 0))
  if (!is.numeric(calibration) || !setequal(names(calibration), c("floor", "deficit")) ||
      any(!is.finite(calibration)) || any(calibration < 0)) {
    stop("`calibration` must be NULL or non-negative numbers named floor and deficit.", call. = FALSE)
  }
  calibration
}

sparse_joint_calibrate <- function(model_sd, relative, shell, birth_rate, calibration) {
  sqrt(model_sd^2 + (calibration[["deficit"]] * relative)^2 +
         (calibration[["floor"]] * birth_rate * (1 + shell))^2)
}

# Fitness relative to a reference: the count-weighted mean fitness of the
# major karyotypes. Two samples say nothing about absolute fitness, and the
# posterior of every log birth rate carries the uncertainty of the baseline,
# which is common to all karyotypes and cancels in a difference. `latent` is
# the position of each row's log birth rate among the latent variables,
# `weight` its weight in the reference (zero outside it).
sparse_joint_relative <- function(solver, latent, fitness, log_birth_sd, weight) {
  weight <- weight / sum(weight)
  rhs <- numeric(length(solver$names))
  contribution <- weight * fitness
  rhs[unique(latent)] <- as.numeric(tapply(contribution, factor(latent, levels = unique(latent)), sum))
  covariance <- as.numeric(solver$solve(rhs))
  reference_variance <- sum(rhs * covariance)
  list(
    mean = fitness - sum(contribution),
    sd = sqrt(pmax(fitness^2 * log_birth_sd^2 + reference_variance - 2 * fitness * covariance[latent], 0)),
    # What `predict()` needs for karyotypes outside the domain.
    reference = list(fitness = sum(contribution), variance = reference_variance,
                     coef_covariance = covariance[solver$names == "coef"]),
    covariance = covariance
  )
}

# Relative fitness of karyotypes predicted from the trend or from the level
# of the observed neighbours, neither of which is a latent variable.
sparse_joint_relative_outside <- function(reference, fitness, log_birth_sd, design = NULL) {
  covariance <- if (is.null(design)) 0 else as.numeric(design %*% reference$coef_covariance)
  list(mean = fitness - reference$fitness,
       sd = sqrt(pmax(fitness^2 * log_birth_sd^2 + reference$variance - 2 * fitness * covariance, 0)))
}

# Level of the observed karyotypes that are not major: used to predict the
# observed karyotypes left out of the domain. An observed karyotype is much
# closer in fitness to the majors than an unobserved one at the same distance.
sparse_joint_observed_level <- function(summary) {
  observed_minor <- !summary$major & (summary$count_t0 + summary$count_t1) > 0
  if (!any(observed_minor)) return(NULL)
  c(log_birth = mean(summary$log_birth_mean[observed_minor]),
    sd = sqrt(mean(summary$log_birth_sd[observed_minor]^2) +
                if (sum(observed_minor) > 1L) stats::var(summary$log_birth_mean[observed_minor]) else 0))
}

# Fitness columns from a log birth rate and its standard deviation. `self` is
# 2 P_ii - 1, the factor between birth rate and net growth rate.
sparse_joint_fitness_columns <- function(self, log_birth, log_birth_sd) {
  z <- stats::qnorm(0.975)
  data.frame(
    log_birth_mean = log_birth, log_birth_sd = log_birth_sd,
    fitness_mean = self * exp(log_birth), fitness_sd = self * exp(log_birth) * log_birth_sd,
    conf_low = self * exp(log_birth - z * log_birth_sd),
    conf_high = self * exp(log_birth + z * log_birth_sd)
  )
}

# Probability that a daughter has exactly its parent's karyotype.
sparse_joint_self_probability <- function(karyotypes, beta) {
  stay <- vapply(0:max(karyotypes), function(cn) alfak2_pij_cpp(cn, cn, beta), numeric(1))
  exp(rowSums(matrix(log(stay)[karyotypes + 1L], nrow = nrow(karyotypes))))
}

# Candidate histories of the balance prior: the interval, and 30, 100, 300,
# 1000 and 3000 divisions at the baseline birth rate where longer than the
# interval. Freely chosen by marginal likelihood, the history went to 30-100
# divisions even where the interval was 180-1080 days and predicted better;
# no candidate is therefore shorter than the interval
# (validation/sparse_joint/short_interval).
sparse_joint_history_candidates <- function(dt, birth_rate) {
  divisions <- c(30, 100, 300, 1000, 3000) / birth_rate
  c(dt, divisions[divisions > dt])
}

# Everything a fit needs that follows from the data and the settings alone:
# the domain, the division kernel, the trend design and the TMB inputs. Used by
# the fit and to rebuild a fit's posterior later from its stored settings.
sparse_joint_problem <- function(data, settings) {
  domain <- sparse_joint_domain(data$counts, settings$major_frequency, settings$major_count, settings$max_states,
                                as.integer(settings$min_cn), as.integer(settings$max_cn))
  n <- length(domain$labels)
  kernel <- sparse_division_kernel(domain$karyotypes, data$beta)
  total <- domain$y0 + domain$y1
  centre <- colSums(domain$karyotypes * total) / sum(total)
  design <- sparse_joint_design(domain$karyotypes, centre, domain$major, settings$shell_power)
  initial <- sparse_joint_initial_prior(domain$y0, domain$y1, settings$initial_frequency)
  tmb_data <- list(
    model = "sparse_joint",
    y0 = domain$y0, y1 = domain$y1, dt = data$dt,
    edge_from = kernel$from - 1L, edge_to = kernel$to - 1L,
    edge_log_prob = kernel$log_prob, log_self = kernel$log_self,
    design = design$design, coef_group = design$group,
    reference = which.max(total) - 1L,
    zeta_prior_mean = initial$mean, zeta_prior_sd = initial$sd,
    intercept_mean = log(settings$birth_rate), intercept_sd = settings$birth_rate_log_sd,
    shell_sd = settings$shell_sd,
    # Centred on a 1% effect per step; wide enough to reach 0.05% and 20%.
    log_scale_prior_mean = rep(log(0.01), 4L),
    log_scale_prior_sd = rep(1.5, 4L),
    n_major = domain$n_major, balance = as.integer(settings$initial_state == "balance"),
    history = settings$history, balance_sd = settings$balance_sd, balance_offset = rep(0, n),
    weight = rep(1, n)
  )
  parameters <- list(
    zeta = tmb_data$zeta_prior_mean,
    g = rep(log(settings$birth_rate), n),
    coef = c(log(settings$birth_rate), rep(0, ncol(design$design) - 1L)),
    theta1 = log(sum(domain$y1) / sum(domain$y0)),
    log_scale = tmb_data$log_scale_prior_mean
  )
  # There are no interaction terms in the design, so their scale is not free.
  free_scale <- c(TRUE, TRUE, FALSE, TRUE)
  list(domain = domain, kernel = kernel, total = total, centre = centre, design = design,
       tmb_data = tmb_data, parameters = parameters, free_scale = free_scale,
       map = list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA))),
       random = c("zeta", "g", "coef", "theta1"))
}

# The posterior of a fitted model, rebuilt from its stored settings and mode
# with the scales fixed at their fitted values. `weight`,
# `zeta_prior_mean`/`zeta_prior_sd` and `design` replace the defaults, for
# instance to hold out the counts of some states. Returns the problem, the solver of
# `sparse_joint_posterior()` and the TMB objective; with `posterior = FALSE`
# only the mode, as a list of the latent parameters.
sparse_joint_rebuild <- function(fit, weight = NULL, zeta_prior_mean = NULL, zeta_prior_sd = NULL,
                                 design = NULL, posterior = TRUE) {
  if (is.null(fit$settings) || is.null(fit$posterior)) {
    stop("This fit predates stored settings; refit it with the current alfak2.", call. = FALSE)
  }
  problem <- sparse_joint_problem(fit$data, fit$settings)
  tmb_data <- problem$tmb_data
  if (!is.null(weight)) tmb_data$weight <- as.numeric(weight)
  if (!is.null(zeta_prior_mean)) tmb_data$zeta_prior_mean <- as.numeric(zeta_prior_mean)
  if (!is.null(zeta_prior_sd)) tmb_data$zeta_prior_sd <- as.numeric(zeta_prior_sd)
  if (!is.null(design)) tmb_data$design <- design
  parameters <- fit$posterior[c("zeta", "g", "coef", "theta1", "log_scale")]
  map <- list(log_scale = factor(rep(NA, length(parameters$log_scale))))
  inner_control <- fit$settings$inner_control
  objective <- TMB::MakeADFun(tmb_data, parameters, map = map, random = problem$random, DLL = "alfak2",
                              silent = TRUE, inner.control = inner_control)
  if (!posterior) {
    objective$fn(numeric(0))
    mode <- objective$env$last.par
    return(lapply(stats::setNames(nm = problem$random), function(name) unname(mode[names(mode) == name])))
  }
  solver <- sparse_joint_posterior(objective, tmb_data, parameters, map, problem$random, numeric(0), inner_control)
  list(problem = problem, solver = solver, objective = objective, tmb_data = tmb_data)
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
#' `relative_sd` is not the model's own posterior standard deviation, which
#' is kept as `relative_sd_model`. That one takes the model to be exact, and
#' on simulated populations its 95% intervals cover about three in four
#' observed karyotypes and half of the unobserved ones. The reported value is
#' `sqrt(model sd^2 + (deficit * relative fitness)^2 +
#' (floor * birth_rate * (1 + shell))^2)`, where `shell` is the distance to
#' the nearest major karyotype. The two constants were fitted to simulated
#' landscapes, on which the resulting 95% intervals cover 91-97%. They carry
#' over to real samples only as far as those landscapes resemble them.
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
#' @param initial_state How abundances at the first timepoint are treated.
#'   `"balance"`, the default, starts every karyotype that is not major at
#'   mutation-selection balance with the major karyotypes next to it, so that
#'   its count in the first sample also informs its fitness. `"free"`
#'   estimates one for every karyotype and so uses the first sample only as a
#'   starting point. Balance takes fitness before the first sample to have
#'   resembled fitness between the samples. On simulated landscapes it was
#'   still more accurate than `"free"` when the landscape changed at the first
#'   sample, as under a new treatment.
#' @param history Time, in units of `dt`, over which the balance built up
#'   before the first sample. By default it is chosen among the interval
#'   between the samples and 30, 100, 300, 1000 and 3000 divisions at
#'   `birth_rate` that are longer than the interval: the model is fitted at
#'   each and the fit with the highest marginal likelihood is kept, so a fit
#'   takes up to six times as long; `diagnostics$history_profile` lists the
#'   candidates. A single value fixes the history; several values are the
#'   candidates. Up to alfak2 8327282 the history was the interval itself
#'   (`history = dt`). With samples a few weeks apart that is far shorter than
#'   the population has been evolving: the balance prior then holds the
#'   neighbours of the major karyotypes far below their observed abundance, and
#'   their fitness is overestimated (validation/sparse_joint/short_interval).
#' @param balance_sd Spread, on the log scale, of initial abundances around
#'   the balance.
#' @param calibration Constants `floor` and `deficit` that turn the model's
#'   posterior standard deviation of relative fitness into the reported
#'   `relative_sd`; see Details. `NULL` reports the model's own.
#' @param control List passed to [stats::nlminb()] for the scale parameters.
#' @param inner_control TMB inner Newton controls.
#'
#' @return An `alfak2_sparse_joint_fit` object. `summary` has one row per
#'   domain karyotype, with its counts and the counts the fit expects;
#'   `fitness_mean` is the net growth rate per unit of `dt`,
#'   and `relative_fitness` is that rate minus the count-weighted mean over
#'   the major karyotypes, with `relative_sd` its uncertainty.
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
                             initial_state = c("balance", "free"),
                             history = NULL,
                             balance_sd = 1,
                             calibration = sparse_joint_default_calibration,
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
  initial_state <- match.arg(initial_state)
  calibration <- sparse_joint_validate_calibration(calibration)
  if (is.null(history)) {
    history <- if (initial_state == "balance") sparse_joint_history_candidates(data$dt, birth_rate) else data$dt
  }
  if (!is.numeric(history) || !length(history)) stop("`history` must be one or more positive times.", call. = FALSE)
  for (h in history) validate_scalar(h, "history", lower = .Machine$double.eps)
  # Without the balance prior the history does not enter the model.
  if (initial_state == "free") history <- history[[1L]]
  validate_scalar(balance_sd, "balance_sd", lower = .Machine$double.eps)
  if (length(history) > 1L) {
    # Candidate histories differ only in the means of the balance prior, so
    # their marginal likelihoods compare directly: fit at each, keep the best.
    started <- proc.time()[["elapsed"]]
    arguments <- mget(names(formals(sys.function())), envir = environment())
    arguments$counts <- data
    fits <- lapply(history, function(h) {
      arguments$history <- h
      do.call(fit_sparse_joint, arguments)
    })
    objective <- vapply(fits, function(f) f$diagnostics$objective, numeric(1))
    fit <- fits[[which.min(objective)]]
    fit$diagnostics$history_profile <- data.frame(
      history = history, objective = objective,
      convergence = vapply(fits, function(f) f$diagnostics$convergence, numeric(1))
    )
    fit$diagnostics$seconds <- proc.time()[["elapsed"]] - started
    return(fit)
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

  settings <- list(birth_rate = birth_rate, birth_rate_log_sd = birth_rate_log_sd,
                   major_frequency = major_frequency, major_count = major_count, max_states = max_states,
                   min_cn = as.integer(min_cn), max_cn = as.integer(max_cn), shell_power = shell_power,
                   shell_sd = shell_sd, initial_frequency = initial_frequency, initial_state = initial_state,
                   history = history, balance_sd = balance_sd, calibration = calibration,
                   inner_control = inner_control)
  problem <- sparse_joint_problem(data, settings)
  domain <- problem$domain
  n <- length(domain$labels)
  if (any(c(sum(domain$y0), sum(domain$y1)) <= 0)) {
    stop("A timepoint has no cell inside the model domain.", call. = FALSE)
  }
  kernel <- problem$kernel
  total <- problem$total
  centre <- problem$centre
  design <- problem$design
  scale_names <- c("dosage", "curvature", "interaction", "residual")
  tmb_data <- problem$tmb_data
  parameters <- problem$parameters
  free_scale <- problem$free_scale
  map <- problem$map
  random <- problem$random
  objective <- TMB::MakeADFun(tmb_data, parameters, map = map, random = random, DLL = "alfak2",
                              silent = TRUE, inner.control = inner_control)
  optimum <- sparse_joint_optimise_scales(objective, control)
  # Uncertainty is conditional on the fitted scales: the marginal likelihood
  # is often flat in a scale.
  solver <- sparse_joint_posterior(objective, tmb_data, parameters, map, random, optimum$par, inner_control)
  g_index <- which(solver$names == "g")
  g_hat <- cbind(solver$mode[g_index], sqrt(solver$variance(g_index)))
  coef_hat <- cbind(solver$mode[solver$names == "coef"])
  coef_covariance <- sparse_joint_coef_covariance(solver)
  log_scale <- tmb_data$log_scale_prior_mean
  log_scale[free_scale] <- optimum$par
  scale <- stats::setNames(exp(log_scale), scale_names)

  self <- 2 * exp(kernel$log_self) - 1
  summary <- cbind(
    data.frame(karyotype = domain$labels, count_t0 = domain$y0, count_t1 = domain$y1,
               expected_t0 = solver$expected_t0, expected_t1 = solver$expected_t1,
               major = seq_len(n) <= domain$n_major,
               shell = sparse_joint_major_distance(domain$karyotypes, domain$major),
               stringsAsFactors = FALSE),
    sparse_joint_fitness_columns(self, g_hat[, 1L], g_hat[, 2L])
  )
  # Relative to the count-weighted mean fitness of the majors.
  relative <- sparse_joint_relative(solver, g_index, summary$fitness_mean,
                                    summary$log_birth_sd, ifelse(summary$major, total, 0))
  summary$relative_fitness <- relative$mean
  summary$relative_sd <- sparse_joint_calibrate(relative$sd, relative$mean, summary$shell, birth_rate, calibration)
  summary$relative_sd_model <- relative$sd
  observed_level <- sparse_joint_observed_level(summary)
  structure(list(
    data = data,
    summary = summary,
    dropped = domain$dropped,
    # What rebuilding the posterior needs: see sparse_joint_rebuild().
    settings = settings,
    posterior = list(zeta = solver$mode[solver$names == "zeta"], g = solver$mode[solver$names == "g"],
                     coef = solver$mode[solver$names == "coef"], theta1 = solver$mode[solver$names == "theta1"],
                     log_scale = log_scale),
    trend = list(coefficients = coef_hat[, 1L], covariance = coef_covariance, centre = centre,
                 major = domain$major, shell_power = shell_power, scales = scale,
                 observed_level = observed_level, reference = relative$reference,
                 calibration = calibration, birth_rate = birth_rate),
    diagnostics = list(
      convergence = optimum$convergence, message = optimum$message, objective = optimum$objective,
      finite_uncertainty = all(is.finite(g_hat[, 2L])) && all(is.finite(relative$sd)),
      n_states = n, n_major = domain$n_major, n_edges = length(kernel$from),
      dropped_karyotypes = length(domain$dropped), dropped_cells = domain$dropped_cells,
      birth_rate = birth_rate, fitted_birth_rate = exp(coef_hat[1L, 1L]),
      initial_state = initial_state, history = history,
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
#'   unit of `dt`), `fitness_sd`, `conf_low`, `conf_high`, `relative_fitness`,
#'   `relative_sd`, `relative_sd_model`, `shell` and `in_domain`.
#'   `relative_fitness` is fitness minus the count-weighted mean fitness of
#'   the major karyotypes, which is what two samples identify. `relative_sd`
#'   is its calibrated uncertainty and `relative_sd_model` the model's own;
#'   see [fit_sparse_joint()]. Neither includes the uncertainty of the
#'   baseline birth rate, which is common to all karyotypes.
#' @export
predict.alfak2_sparse_joint_fit <- function(object, karyotypes, ...) {
  labels <- if (is.character(karyotypes)) karyotypes else format_karyotypes(as.matrix(karyotypes))
  inside <- match(labels, object$summary$karyotype)
  out <- data.frame(karyotype = labels, fitness_mean = NA_real_, fitness_sd = NA_real_,
                    conf_low = NA_real_, conf_high = NA_real_, relative_fitness = NA_real_,
                    relative_sd = NA_real_, relative_sd_model = NA_real_, shell = NA_real_,
                    in_domain = !is.na(inside), stringsAsFactors = FALSE)
  columns <- c("fitness_mean", "fitness_sd", "conf_low", "conf_high", "relative_fitness", "relative_sd",
               "relative_sd_model", "shell")
  out[!is.na(inside), columns] <- object$summary[inside[!is.na(inside)], columns]
  outside <- which(is.na(inside))
  if (length(outside)) {
    trend <- object$trend
    cn <- parse_karyotypes(labels[outside])
    design <- sparse_joint_design(cn, trend$centre, trend$major, trend$shell_power)$design
    log_birth <- as.numeric(design %*% trend$coefficients)
    log_birth_sd <- sqrt(rowSums((design %*% trend$covariance) * design) + trend$scales[["residual"]]^2)
    left_out <- labels[outside] %in% object$dropped & !is.null(trend$observed_level)
    if (any(left_out)) {
      log_birth[left_out] <- trend$observed_level[["log_birth"]]
      log_birth_sd[left_out] <- trend$observed_level[["sd"]]
    }
    self <- 2 * sparse_joint_self_probability(cn, object$data$beta) - 1
    fitness <- sparse_joint_fitness_columns(self, log_birth, log_birth_sd)
    out[outside, c("fitness_mean", "fitness_sd", "conf_low", "conf_high")] <-
      fitness[c("fitness_mean", "fitness_sd", "conf_low", "conf_high")]
    relative <- sparse_joint_relative_outside(trend$reference, fitness$fitness_mean, log_birth_sd, design)
    level <- sparse_joint_relative_outside(trend$reference, fitness$fitness_mean, log_birth_sd)
    out$relative_fitness[outside] <- relative$mean
    out$relative_sd_model[outside] <- ifelse(left_out, level$sd, relative$sd)
    out$shell[outside] <- sparse_joint_major_distance(cn, trend$major)
    out$relative_sd[outside] <- sparse_joint_calibrate(out$relative_sd_model[outside], relative$mean,
                                                       out$shell[outside], trend$birth_rate, trend$calibration)
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
