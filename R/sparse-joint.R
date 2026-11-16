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

# Major karyotypes anchor the domain: at least `major_frequency` of any
# sample and at least `major_count` cells over all samples. The share adapts the rule to depth:
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
  share <- do.call(pmax, lapply(seq_len(ncol(counts)), function(k) counts[, k] / sum(counts[, k])))
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
# one. Observed karyotypes outside it are left out with their cells. `y` has
# one column per timepoint; `y0` and `y1` are its first two.
sparse_joint_domain <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn,
                                major_columns = NULL) {
  observed <- rownames(counts)
  # `major_columns`: the samples the majors are chosen from (all by default).
  from <- if (is.null(major_columns)) counts else counts[, major_columns, drop = FALSE]
  major_labels <- sparse_joint_major(from, major_frequency, major_count, max_states, min_cn, max_cn)
  major <- parse_karyotypes(major_labels)
  labels <- union(major_labels, format_karyotypes(one_step_neighbours(major, min_cn, max_cn)))
  in_data <- match(labels, observed)
  y <- matrix(0, length(labels), ncol(counts))
  y[!is.na(in_data), ] <- counts[in_data[!is.na(in_data)], , drop = FALSE]
  outside <- !observed %in% labels
  list(labels = labels, karyotypes = parse_karyotypes(labels), major = major,
       n_major = length(major_labels), y = y, y0 = y[, 1L], y1 = y[, 2L],
       dropped = observed[outside], dropped_cells = sum(counts[outside, ]))
}

# How much of each sample the domain holds: per timepoint, its cells and those
# inside the domain, its karyotypes, those that meet the share and count of a
# major in that sample (`major_candidates`) and the majors seen in it, the
# karyotypes it shares with another timepoint (with two timepoints, those both
# share), and the median distance from its cells outside the domain to the
# nearest major.
sparse_joint_domain_coverage <- function(counts, domain, major_frequency, major_count) {
  labels <- rownames(counts)
  inside <- labels %in% domain$labels
  majors <- labels %in% domain$labels[seq_len(domain$n_major)]
  total <- rowSums(counts)
  distance <- rep(NA_real_, length(labels))
  if (any(!inside)) distance[!inside] <- sparse_joint_major_distance(parse_karyotypes(labels[!inside]), domain$major)
  seen <- counts > 0
  rows <- lapply(seq_len(ncol(counts)), function(t) {
    cells <- counts[, t]
    out <- !inside & cells > 0
    data.frame(
      timepoint = colnames(counts)[[t]], cells = sum(cells), cells_in_domain = sum(cells[inside]),
      fraction = sum(cells[inside]) / sum(cells), karyotypes = sum(cells > 0),
      major_candidates = sum(cells / sum(cells) >= major_frequency & total >= major_count),
      majors = sum(majors & cells > 0), shared = sum(seen[, t] & rowSums(seen[, -t, drop = FALSE]) > 0),
      median_distance_outside = if (any(out)) stats::median(rep(distance[out], cells[out])) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Error for a sample with a timepoint that has no cell inside its domain. The
# first sentence of the message is fixed, since callers match it; the rest
# says why. `coverage` is attached to the condition; with `by_sample` it has
# a `sample` column (cohort fits).
sparse_joint_coverage_error <- function(coverage, first, by_sample = FALSE) {
  plural <- function(n, word) sprintf("%d %s%s", as.integer(n), word, if (n == 1) "" else "s")
  empty <- coverage[coverage$cells_in_domain <= 0, , drop = FALSE]
  reasons <- vapply(seq_len(nrow(empty)), function(i) {
    row <- empty[i, ]
    same <- if (by_sample) coverage$sample == row$sample else rep(TRUE, nrow(coverage))
    other <- coverage[same & coverage$timepoint != row$timepoint, , drop = FALSE]
    where <- if (by_sample) sprintf("Sample %s, timepoint %s", row$sample, row$timepoint) else
      sprintf("Timepoint %s", row$timepoint)
    if (nrow(other) == 1L) {
      shared <- if (row$shared == 0) "the samples share no karyotype" else
        sprintf("the samples share %s", plural(row$shared, "karyotype"))
      majors <- other$majors[[1L]]
      source <- if (majors == 1) sprintf("the only major comes from %s", other$timepoint[[1L]]) else
        sprintf("all %d majors come from %s", as.integer(majors), other$timepoint[[1L]])
    } else {
      shared <- if (row$shared == 0) "it shares no karyotype with another timepoint" else
        sprintf("it shares %s with other timepoints", plural(row$shared, "karyotype"))
      source <- sprintf("the majors come from %s", paste(other$timepoint[other$majors > 0], collapse = ", "))
    }
    sprintf("%s: none of its %s lies within one step of a major karyotype; %s, %s, and the median distance to the nearest major is %s.",
            where, plural(row$cells, "cell"), source, shared, format(row$median_distance_outside))
  }, character(1))
  message <- paste(c(first, reasons, "The samples may be too far apart, or one of them may be noisy."), collapse = " ")
  stop(errorCondition(message, class = "alfak2_domain_coverage_error", coverage = coverage, call = NULL))
}

# Coverage of a sample's domain; an error if a timepoint in `required` (all by
# default) has no cell in it.
sparse_joint_check_coverage <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn,
                                        major_columns = NULL, required = NULL) {
  domain <- sparse_joint_domain(counts, major_frequency, major_count, max_states, min_cn, max_cn, major_columns)
  coverage <- sparse_joint_domain_coverage(counts, domain, major_frequency, major_count)
  if (is.null(required)) required <- seq_len(nrow(coverage))
  if (any(coverage$cells_in_domain[required] <= 0)) {
    sparse_joint_coverage_error(coverage, "A timepoint has no cell inside the model domain.")
  }
  coverage
}

# "t0 812 of 1000 (81%), t1 640 of 1000 (64%)", optionally prefixed by sample.
sparse_joint_format_coverage <- function(coverage) {
  label <- if (is.null(coverage$sample)) coverage$timepoint else paste(coverage$sample, coverage$timepoint)
  paste(sprintf("%s %d of %d (%.0f%%)", label, as.integer(coverage$cells_in_domain), as.integer(coverage$cells),
                100 * coverage$fraction), collapse = ", ")
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
# count itself, essentially free. Seen only later: its frequency at the
# earliest timepoint it was seen, scaled to the first sample's depth. Never
# seen: rare, on the frequency scale given by `initial_frequency`.
sparse_joint_initial_prior <- function(y0, y1, initial_frequency) {
  sparse_joint_initial_prior_matrix(cbind(y0, y1), initial_frequency)
}

sparse_joint_initial_prior_matrix <- function(y, initial_frequency, depth = colSums(y)) {
  y0 <- y[, 1L]
  seen_first <- y0 > 0
  mean <- rep(log(depth[[1L]]) + initial_frequency[[1L]], length(y0))
  # Latest timepoint first, so that the earliest one seen is written last.
  for (k in rev(seq_len(ncol(y))[-1L])) {
    yk <- y[, k]
    mean <- ifelse(yk > 0, log(pmax(yk, 1) * depth[[1L]] / depth[[k]]), mean)
  }
  list(
    mean = ifelse(seen_first, log(pmax(y0, 1)), mean),
    sd = ifelse(seen_first, 10, initial_frequency[[2L]])
  )
}

# Scale parameters: sparse_joint_optimise_scales() in sparse-joint-convergence.R.

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
  series_method <- tmb_data$model %in% c("sparse_joint_pairwise", "sparse_joint_state_space")
  if (series_method) {
    # The prior residuals of these templates are set to zero there as well.
    computed <- sparse_joint_series_method_expected(tmb_data, parameters, map, mode, report)
    expected <- computed$data
  } else if (identical(tmb_data$model, "sparse_joint_series")) {
    expected$y <- cbind(exp(unname(mode[names(mode) == "zeta"])), exp(report$log_intensity))
  } else {
    expected$y0 <- exp(unname(mode[names(mode) == "zeta"]))
    expected$y1 <- exp(as.numeric(report$log_intensity1))
  }
  # The balance prior on initial abundances is not log-linear either: its
  # curvature is taken at zero residual in the same way.
  if (!series_method) expected$balance_offset <- as.numeric(report$balance_residual)
  information <- TMB::MakeADFun(expected, parameters, map = map, random = random, DLL = "alfak2",
                                silent = TRUE, inner.control = inner_control)$env$spHess(mode, random = TRUE)
  factor <- Matrix::Cholesky(information, LDL = FALSE)
  n <- nrow(information)
  solve <- function(rhs) as.matrix(Matrix::solve(factor, rhs))
  expected_counts <- if (series_method) computed$counts else if (is.null(expected$y)) cbind(expected$y0, expected$y1) else expected$y
  list(
    mode = unname(mode[latent]), names = names(mode)[latent], solve = solve,
    # Fitted expected counts, for comparison with the observed ones: one
    # column per timepoint, and the first two separately.
    expected = expected_counts, expected_t0 = expected_counts[, 1L], expected_t1 = expected_counts[, 2L],
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
# major karyotype. Floor and deficit were fitted to the development landscapes
# of `benchmark/sparse_two_timepoint` (validation/sparse_joint/calibration_fit.R).
#
# Karyotypes seen in neither sample have a deficit of their own. How far below
# the majors they lie is not identified by zero counts, and the predicted
# drop is too shallow on moderately rugged landscapes, by up to 40% in young
# populations or at intervals of days to weeks; with the deficit of observed
# karyotypes their 95% intervals covered as few as 79%. No observable told
# those regimes apart from smoother or more rugged ones, so the deficit is the
# smallest that keeps every simulated regime at 90% or more
# (validation/sparse_joint/unobserved_bias/calibration_forms.R).
sparse_joint_default_calibration <- c(floor = 0.009, deficit = 0.6, deficit_unobserved = 1.1)

# Calibration constants in full. Without `deficit_unobserved`, unobserved
# karyotypes take the deficit of observed ones, as up to alfak2 33a045c.
sparse_joint_validate_calibration <- function(calibration) {
  if (is.null(calibration)) return(c(floor = 0, deficit = 0, deficit_unobserved = 0))
  two <- c("floor", "deficit")
  if (!is.numeric(calibration) || is.null(names(calibration)) ||
      !(setequal(names(calibration), two) || setequal(names(calibration), c(two, "deficit_unobserved"))) ||
      any(!is.finite(calibration)) || any(calibration < 0)) {
    stop("`calibration` must be NULL or non-negative numbers named floor, deficit and optionally ",
         "deficit_unobserved.", call. = FALSE)
  }
  unobserved <- if ("deficit_unobserved" %in% names(calibration)) calibration[["deficit_unobserved"]] else
    calibration[["deficit"]]
  c(floor = calibration[["floor"]], deficit = calibration[["deficit"]], deficit_unobserved = unobserved)
}

# `observed`: whether each karyotype was seen in a sample. Fits made before
# `deficit_unobserved` existed carry two constants and use `deficit` for all.
sparse_joint_calibrate <- function(model_sd, relative, shell, birth_rate, calibration, observed) {
  unobserved <- if ("deficit_unobserved" %in% names(calibration)) calibration[["deficit_unobserved"]] else
    calibration[["deficit"]]
  deficit <- ifelse(observed, calibration[["deficit"]], unobserved)
  sqrt(model_sd^2 + (deficit * relative)^2 + (calibration[["floor"]] * birth_rate * (1 + shell))^2)
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
  observed_minor <- !summary$major & sparse_joint_summary_total(summary) > 0
  if (!any(observed_minor)) return(NULL)
  c(log_birth = mean(summary$log_birth_mean[observed_minor]),
    sd = sqrt(mean(summary$log_birth_sd[observed_minor]^2) +
                if (sum(observed_minor) > 1L) stats::var(summary$log_birth_mean[observed_minor]) else 0))
}

# Cells of each karyotype of a fit's summary over all timepoints (columns
# count_t0, count_t1, ...).
sparse_joint_summary_total <- function(summary) {
  columns <- grep("^count_t[0-9]+$", names(summary), value = TRUE)
  if (length(columns) == 2L) summary[[columns[[1L]]]] + summary[[columns[[2L]]]] else
    rowSums(as.matrix(summary[columns]))
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

# Template of a fit: "pair" (src/sparse_joint_model.h) for two timepoints,
# "series" (src/sparse_joint_series_model.h) for more. The series template
# gives the same objective for two timepoints; `settings$template = "series"`
# forces it, to check exactly that.
sparse_joint_template <- function(data, settings) {
  if (ncol(data$counts) > 2L && !is.null(settings$series) && settings$series != "anchored") return(settings$series)
  if (ncol(data$counts) > 2L || identical(settings$template, "series")) "series" else "pair"
}

# The samples the majors are chosen from: all of them, or with more than two
# samples and `settings$intermediate` "majors" the first and the last. A clone
# that rises and falls between them, or appears in between only, would
# otherwise become a major, and the domain grow around it.
sparse_joint_major_columns <- function(data, settings) {
  n_time <- ncol(data$counts)
  if (n_time > 2L && identical(settings$intermediate, "majors")) c(1L, n_time) else NULL
}

# Samples between the first and the last with no cell on the states they
# count for: the majors with `intermediate` "majors", the domain with "all".
# Such a sample says nothing under the model. Its depth offset has a flat
# prior, so its likelihood is largest as the offset goes to minus infinity,
# where the sample contributes nothing, and the Laplace approximation over
# that offset is improper. Leaving it out is that limit. Leaving it out
# changes neither the majors nor the domain: it has no cell on them.
sparse_joint_uncounted_samples <- function(coverage, intermediate) {
  n_time <- nrow(coverage)
  if (n_time <= 2L) return(integer(0))
  middle <- seq(2L, n_time - 1L)
  counted <- if (identical(intermediate, "majors")) coverage$majors[middle] else coverage$cells_in_domain[middle]
  middle[counted <= 0]
}

# The data with only the samples in `kept` (in order, the first among them),
# as prepare_alfak2_data() would give them: karyotypes with no cell left are
# dropped.
sparse_joint_keep_samples <- function(data, kept) {
  times <- alfak2_data_times(data)[kept]
  times <- times - times[[1L]]
  counts <- data$counts[, kept, drop = FALSE]
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  new_alfak2_data(counts = counts, dt = times[[length(times)]], beta = data$beta,
                  metadata = data$metadata, times = if (length(kept) > 2L) times else NULL)
}

# Dispersion of the majors' counts in the samples between the first and the
# last, against the fit's expectations scaled to each sample's total over the
# majors: about 1 for sampling noise alone, more when the population drifts
# between samples (Pearson statistic per count).
sparse_joint_intermediate_dispersion <- function(fit) {
  summary <- fit$summary
  counts <- as.matrix(summary[grep("^count_t[0-9]+$", names(summary))])
  expected <- as.matrix(summary[grep("^expected_t[0-9]+$", names(summary))])
  n_time <- ncol(counts)
  if (n_time < 3L || !any(summary$major)) return(NA_real_)
  middle <- seq(2L, n_time - 1L)
  y <- counts[summary$major, middle, drop = FALSE]
  e <- expected[summary$major, middle, drop = FALSE]
  e <- sweep(e, 2L, colSums(y) / colSums(e), "*")
  mean((y - e)^2 / pmax(e, 1e-8))
}

# The domain's counts as the likelihood uses them, with `sample_weight` (one
# row per state, one column per sample after the first: 1 counts, 0 left out)
# as an attribute. With more than two samples and `settings$intermediate`
# "majors", the counts of the other states in the samples between the first
# and the last are left out (and taken as zero for the priors): a rare clone
# seen in one of them owes its count largely to the chance with which such a
# clone arises, drifts and is lost, which the model does not describe, and the
# same chance event would otherwise count once per sample.
sparse_joint_used_counts <- function(domain, settings) {
  used <- domain$y
  n_time <- ncol(used)
  sample_weight <- matrix(1, nrow(used), n_time - 1L)
  if (n_time > 2L && identical(settings$intermediate, "majors")) {
    minor <- seq_len(nrow(used)) > domain$n_major
    middle <- seq(2L, n_time - 1L)
    used[minor, middle] <- 0
    sample_weight[minor, middle - 1L] <- 0
    weight <- if (is.null(settings$intermediate_weight)) 1 else settings$intermediate_weight
    sample_weight[!minor, middle - 1L] <- weight
  }
  attr(used, "sample_weight") <- sample_weight
  used
}

# Everything a fit needs that follows from the data and the settings alone:
# the domain, the division kernel, the trend design and the TMB inputs. Used by
# the fit and to rebuild a fit's posterior later from its stored settings.
sparse_joint_problem <- function(data, settings) {
  template <- sparse_joint_template(data, settings)
  if (template == "pairwise") return(sparse_joint_pairwise_problem(data, settings))
  if (template == "state_space") return(sparse_joint_state_space_problem(data, settings))
  domain <- sparse_joint_domain(data$counts, settings$major_frequency, settings$major_count, settings$max_states,
                                as.integer(settings$min_cn), as.integer(settings$max_cn),
                                sparse_joint_major_columns(data, settings))
  n <- length(domain$labels)
  kernel <- sparse_division_kernel(domain$karyotypes, data$beta)
  # The counts the likelihood uses: with `intermediate = "majors"`, the
  # samples between the first and the last count for the majors only.
  used <- sparse_joint_used_counts(domain, settings)
  total <- if (ncol(domain$y) == 2L) domain$y0 + domain$y1 else rowSums(used)
  centre <- colSums(domain$karyotypes * total) / sum(total)
  design <- sparse_joint_design(domain$karyotypes, centre, domain$major, settings$shell_power)
  initial <- sparse_joint_initial_prior_matrix(used, settings$initial_frequency, depth = colSums(domain$y))
  counts_data <- if (template == "pair") {
    list(model = "sparse_joint", y0 = domain$y0, y1 = domain$y1, dt = data$dt)
  } else {
    list(model = "sparse_joint_series", y = domain$y, tau = alfak2_data_times(data)[-1L],
         sample_weight = attr(used, "sample_weight"))
  }
  tmb_data <- c(counts_data, list(
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
  ))
  # Log depth ratio of each later sample to the first.
  depth <- if (template == "pair") list(theta1 = log(sum(domain$y1) / sum(domain$y0))) else
    list(theta = log(colSums(domain$y[, -1L, drop = FALSE]) / sum(domain$y[, 1L])))
  parameters <- c(list(
    zeta = tmb_data$zeta_prior_mean,
    g = rep(log(settings$birth_rate), n),
    coef = c(log(settings$birth_rate), rep(0, ncol(design$design) - 1L))),
    depth,
    list(log_scale = tmb_data$log_scale_prior_mean)
  )
  # There are no interaction terms in the design, so their scale is not free.
  free_scale <- c(TRUE, TRUE, FALSE, TRUE)
  list(domain = domain, kernel = kernel, total = total, centre = centre, design = design, used = used,
       tmb_data = tmb_data, parameters = parameters, free_scale = free_scale, template = template,
       map = list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA))),
       random = c("zeta", "g", "coef", names(depth)))
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
  # A matrix keeps its shape (pairwise method: one column per pair).
  keep_shape <- function(new, old) if (is.matrix(old)) matrix(as.numeric(new), nrow(old)) else as.numeric(new)
  if (!is.null(zeta_prior_mean)) tmb_data$zeta_prior_mean <- keep_shape(zeta_prior_mean, tmb_data$zeta_prior_mean)
  if (!is.null(zeta_prior_sd)) tmb_data$zeta_prior_sd <- keep_shape(zeta_prior_sd, tmb_data$zeta_prior_sd)
  if (!is.null(design)) tmb_data$design <- design
  parameters <- fit$posterior[names(problem$parameters)]
  map <- problem$map
  map$log_scale <- factor(rep(NA, length(parameters$log_scale)))
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
#' Fits fitness to sparse karyotype samples of one population, two or more
#' timepoints, with a single count likelihood on a sparse domain: the
#' well-sampled ("major") karyotypes and all their one-step neighbours. Division dynamics use a closed form in which each
#' karyotype grows on its own and receives mutational influx from its
#' neighbours; an unobserved neighbour of an abundant karyotype is thereby
#' evidence that the neighbour is less fit. Fitness away from the domain is
#' predicted with [predict()] from a fitted trend in chromosome dosage,
#' curvature and distance to the nearest major karyotype.
#'
#' With more than two timepoints every later sample is compared with the
#' same closed-form dynamics started at the first sample, each with its own
#' depth; the dynamics are not restarted at intermediate samples. Fitness is
#' taken to be constant over the whole series, and the history of the balance
#' prior is at least the time from the first to the last sample. Two-timepoint
#' fits are unchanged. By default (`intermediate = "majors"`) the majors and
#' the domain come from the first and the last sample, and the samples in
#' between count for the majors only, weighted by the inverse of their
#' dispersion in a first fit (`diagnostics$intermediate_weight`,
#' `intermediate_dispersion`). Clones that arise, grow by chance and are lost
#' between samples, or are thinned by passaging, would otherwise count once
#' per sample and read as fit. On simulated passaged cultures the default was
#' as accurate as the first and the last sample alone; every sample for every
#' karyotype (`"all"`) was 4-6% more accurate with large bottlenecks (20000
#' cells) and up to 2.4 times less accurate with small ones (200)
#' (validation/sparse_joint/multi_timepoint).
#'
#' The first and the last sample must each have a cell in the domain. A
#' sample in between with no cell on the states it counts for (the majors by
#' default, the domain with `"all"`) is left out with a warning
#' (`diagnostics$dropped_samples`). Its depth is free, so such a sample says
#' nothing under the model, and the fit is the fit to the other samples.
#' Leaving it out is the limit of its own likelihood. The input as given is
#' kept in `input`, for [as_alfakR()].
#'
#' With `series` `"pairwise"` or `"state_space"`, each pair of consecutive
#' samples has its own domain, and both samples of a pair need a cell in it.
#' A sample in between without one is left out in the same way, and its
#' neighbours become a pair. The first and the last sample are kept. When one
#' of them has no cell in its pair's domain, its neighbour in between is left
#' out instead, so that it pairs with a sample further in. This repeats until
#' every pair is covered. At worst the first and the last sample remain, and
#' the fit is the two-sample fit.
#'
#' Observed karyotypes that are not next to a major karyotype are left out of
#' the likelihood, together with their cells. They would enter only because
#' they were observed, while their unobserved siblings would not, which
#' biases their fitness upwards. [predict()] places them at the level of the
#' observed neighbours inside the domain.
#'
#' `diagnostics$domain_coverage` reports, for each timepoint, how many of its
#' cells lie inside the domain. If a timepoint has none, typically because
#' every major karyotype comes from the other sample and the two samples share
#' no karyotype, there is nothing to fit and the function stops with an error
#' of class `alfak2_domain_coverage_error`. Its message begins "A timepoint
#' has no cell inside the model domain." and its `coverage` element holds the
#' same table. [check_sparse_joint_counts()] summarises a count matrix before
#' fitting.
#'
#' The scales of the trend and the residual are fitted by Laplace marginal
#' likelihood. `nlminb`'s own verdict on convergence is not taken as it is:
#' where the inner Hessian is nearly singular the Laplace approximation drops
#' for no reason in the data, and `nlminb` stops there with "false
#' convergence", or reports that at a true optimum. The fit checks the point
#' itself (healthy inner Hessian, the same objective from its own mode, no
#' better inner mode from the initial values, a small finite-difference
#' gradient). `diagnostics$convergence` is 0 when the point was accepted and 1
#' when not, with a warning of class `alfak2_convergence_warning`;
#' `diagnostics$optimisation` tells how (`status`, `stage`, `reasons`,
#' `restarts`, `evaluations`, `nlminb_convergence`, `nlminb_message`,
#' `min_eigen`, `projected_gradient`, `objective_difference`,
#' `joint_difference`, `flat_direction`, `flat_composition`). The fit itself
#' is the same as before this check (alfak2 ff13cf4); so is the choice of
#' history, by the lowest objective.
#'
#' Most fits that are not accepted have a karyotype with very few cells, whose
#' initial abundance or birth rate the data hardly constrain. A count in the
#' second sample, for one, has two explanations of nearly equal weight: the
#' karyotype was there and not sampled, or its neighbours made it during the
#' interval. The Laplace approximation is then nearly flat along that
#' quantity, or jumps from one explanation to the other as the scales change.
#' For a fit that is not accepted, `flat_direction` lists the latent variables
#' with the largest shares of the flattest direction of the inner Hessian
#' (`variable`, `karyotype`, its `cells` in each sample, `major`, `share`),
#' `flat_composition` its make-up by kind of latent variable (`variable`,
#' `major`, `share`), and the warning names both when that Hessian is nearly
#' singular. On simulated data, fits that searching again turned into
#' accepted ones predicted fitness as well as before
#' (validation/sparse_joint/convergence).
#'
#' `relative_sd` is not the model's own posterior standard deviation, which
#' is kept as `relative_sd_model`. That one takes the model to be exact, and
#' on simulated populations its 95% intervals cover about three in four
#' observed karyotypes and half of the unobserved ones. The reported value is
#' `sqrt(model sd^2 + (deficit * relative fitness)^2 +
#' (floor * birth_rate * (1 + shell))^2)`, where `shell` is the distance to
#' the nearest major karyotype. The deficit is 0.6 for karyotypes seen in a
#' sample and 1.1 for those seen in neither: zero counts do not tell how far
#' below the major karyotypes an unobserved karyotype lies, and on moderately
#' rugged simulated landscapes the predicted drop was up to 40% too shallow.
#' On simulated landscapes the resulting 95% intervals covered at least 90%
#' of the karyotypes in every regime tested, and close to all of them on
#' smooth or very rugged landscapes, where they are wider than needed. The
#' constants carry over to real samples only as far as those landscapes
#' resemble them.
#'
#' The baseline birth rate is supplied rather than estimated. Two samples
#' identify it only weakly, and a free baseline drifts upwards and inflates
#' every fitness difference. Rerun with other values to see how much a
#' conclusion depends on it.
#'
#' @param counts Count matrix with karyotype row names and one column per
#'   timepoint, at least two, in time order; or an `alfak2_data` object.
#' @param dt Time between the samples, if `counts` is a matrix: one interval
#'   for two timepoints, or one per pair of consecutive timepoints. See also
#'   `times`.
#' @param beta Missegregation rate per chromosome copy and division, if
#'   `counts` is a matrix.
#' @param birth_rate Assumed baseline birth (division) rate per unit of `dt`.
#' @param birth_rate_log_sd Prior standard deviation of the log baseline birth
#'   rate around `log(birth_rate)`.
#' @param major_frequency Minimum share of any sample for a karyotype to
#'   anchor the domain as a major karyotype.
#' @param major_count Minimum number of cells, over all samples, for a major
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
#'   frequency at the first timepoint of a karyotype seen at no timepoint. It
#'   is a frequency, not a count: a domain has thousands of unseen
#'   karyotypes, and a prior in counts would give them a large share of a
#'   small sample. A karyotype seen only later is centred on its frequency at
#'   the earliest timepoint it was seen, with the same standard deviation:
#'   absent other evidence its frequency is taken not to have changed.
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
#'   between the first and the last sample and 30, 100, 300, 1000 and 3000
#'   divisions at `birth_rate` that are longer than that interval: the model is fitted at
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
#' @param calibration Constants `floor`, `deficit` and `deficit_unobserved`
#'   that turn the model's posterior standard deviation of relative fitness
#'   into the reported `relative_sd`; see Details. Without
#'   `deficit_unobserved`, unobserved karyotypes take `deficit`, as up to
#'   alfak2 33a045c (`c(floor = 0.009, deficit = 0.6)`). `NULL` reports the
#'   model's own.
#' @param control List passed to [stats::nlminb()] for the first search of
#'   the scale parameters; see Details for what follows it.
#' @param inner_control TMB inner Newton controls.
#' @param times Sampling times, one per column of `counts`, instead of `dt`;
#'   only their differences matter. When given, `dt` is not used.
#' @param series With more than two samples, how they are tied together.
#'   `"anchored"`, the default: one trajectory from the first sample, with
#'   the majors and the domain from the first and the last sample (see
#'   `intermediate`). `"pairwise"`: every pair of consecutive samples is
#'   fitted as two samples are, with its own majors, domain, initial
#'   abundances and depth, and all pairs share one landscape. `"state_space"`:
#'   every sample has latent abundances, linked from one sample to the next
#'   by the same dynamics plus a drift whose scale is estimated
#'   (`diagnostics$drift`); the pairs' domains are joined. `"auto"`: one of
#'   the three, chosen from the data:
#'   - the anchored method when every sample in between has at least 30% of
#'     its cells on the first and last samples' majors;
#'   - otherwise (a clonal sweep) the pairwise method when the anchored
#'     intermediate weight is at least 0.5 (little drift), or when no sample
#'     in between remains to measure it;
#'   - else the state-space method.
#'
#'   The fit is that method's fit, and `diagnostics$series_choice` gives the
#'   choice and the two measures. Two samples give the same fit with all
#'   methods. The methods and the choice are compared in
#'   validation/sparse_joint/series_methods.
#' @param intermediate With more than two samples, how the samples between
#'   the first and the last are used. `"majors"`, the default: the majors and
#'   the domain are chosen from the first and the last sample, as for two
#'   samples, and the samples in between count for the majors only, with
#'   weight `intermediate_weight`; the other karyotypes are fitted to the first
#'   and the last sample. `"all"`: like the first and the last, for every
#'   state, with the majors chosen from every sample. See Details.
#' @param intermediate_weight With `intermediate = "majors"`, the weight of the
#'   samples in between, in (0, 1]. `NULL` estimates it: the inverse of the
#'   dispersion of the majors' counts in those samples in a first fit (at the
#'   shortest history), at most 1, so that a population that drifts between
#'   samples (small passaging bottlenecks) leans on the first and the last.
#'
#' @return An `alfak2_sparse_joint_fit` object. `summary` has one row per
#'   domain karyotype, with its counts (`count_t0`, `count_t1`, ..., one
#'   column per timepoint) and the counts the fit expects (`expected_t0`, ...);
#'   `fitness_mean` is the net growth rate per unit of `dt`,
#'   and `relative_fitness` is that rate minus the count-weighted mean over
#'   the major karyotypes, with `relative_sd` its uncertainty.
#'   `diagnostics$domain_coverage` has one row per timepoint: `cells`,
#'   `cells_in_domain` and their ratio `fraction`, `karyotypes`,
#'   `major_candidates` (karyotypes with at least `major_frequency` of that
#'   sample and `major_count` cells over both), `majors` (major karyotypes seen
#'   at that timepoint), `shared` (karyotypes seen at both timepoints) and
#'   `median_distance_outside` (median distance from its cells outside the
#'   domain to the nearest major karyotype). With more than two timepoints,
#'   `diagnostics$times` are the sample times fitted and
#'   `diagnostics$dropped_samples` the samples left out. It has one row per
#'   sample left out (`sample`, its column in the input; `time`; and the
#'   coverage columns), and no rows when none was. For the pairwise and
#'   state-space methods its columns are `sample`, `time`, `timepoint`,
#'   `cells` and `uncovered`: the sample without a cell in its pair's domain,
#'   which is the sample itself, or a first or last sample next to it.
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
                             inner_control = list(maxit = 200, tol = 1e-10),
                             times = NULL,
                             series = c("anchored", "pairwise", "state_space", "auto"),
                             intermediate = c("majors", "all"),
                             intermediate_weight = NULL) {
  data <- if (inherits(counts, "alfak2_data")) {
    counts
  } else if (is.null(times)) {
    prepare_alfak2_data(counts, dt = dt, beta = beta)
  } else {
    prepare_alfak2_data(counts, beta = beta, times = times)
  }
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
  intermediate <- match.arg(intermediate)
  series <- match.arg(series)
  calibration <- sparse_joint_validate_calibration(calibration)
  if (is.null(history)) {
    history <- if (initial_state == "balance") sparse_joint_history_candidates(data$dt, birth_rate) else data$dt
  }
  if (!is.numeric(history) || !length(history)) stop("`history` must be one or more positive times.", call. = FALSE)
  for (h in history) validate_scalar(h, "history", lower = .Machine$double.eps)
  # Without the balance prior the history does not enter the model.
  if (initial_state == "free") history <- history[[1L]]
  validate_scalar(balance_sd, "balance_sd", lower = .Machine$double.eps)
  counts <- data$counts
  observed_cn <- parse_karyotypes(rownames(counts))
  if (is.null(max_cn)) max_cn <- max(observed_cn) + 1L
  validate_scalar(max_cn, "max_cn", lower = min_cn)
  if (any(observed_cn < min_cn | observed_cn > max_cn)) {
    stop("Observed karyotypes must lie within the copy-number bounds.", call. = FALSE)
  }
  if (any(colSums(counts) <= 0)) {
    stop(if (ncol(counts) == 2L) "Both timepoints need at least one cell." else "Every timepoint needs at least one cell.",
         call. = FALSE)
  }
  # Checked once, before any candidate history is fitted: a first or last
  # sample with no cell next to a major karyotype leaves nothing to fit.
  n_input <- ncol(counts)
  # Automatic choice of the series method (sparse_joint_auto_series()). Two
  # samples give the two-sample fit with every method.
  if (series == "auto" && n_input <= 2L) series <- "anchored"
  if (series == "auto") {
    if (intermediate != "majors") {
      stop("`series = \"auto\"` chooses with `intermediate = \"majors\"`.", call. = FALSE)
    }
    share <- sparse_joint_ends_major_share(counts, major_frequency, major_count, max_states,
                                           as.integer(min_cn), as.integer(max_cn))
    weight <- NA_real_
    if (share < sparse_joint_auto_thresholds[["ends_major"]]) {
      weight <- if (!is.null(intermediate_weight)) intermediate_weight else {
        # The anchored fit's own weight: from a first fit at the shortest
        # history, on the samples the anchored method keeps.
        coverage <- sparse_joint_check_coverage(counts, major_frequency, major_count, max_states,
                                                as.integer(min_cn), as.integer(max_cn),
                                                sparse_joint_major_columns(data, list(intermediate = "majors")),
                                                required = c(1L, n_input))
        left_out <- sparse_joint_uncounted_samples(coverage, "majors")
        kept <- if (length(left_out)) sparse_joint_keep_samples(data, setdiff(seq_len(n_input), left_out)) else data
        if (ncol(kept$counts) <= 2L) NA_real_ else {
          arguments <- mget(names(formals(sys.function())), envir = environment())
          arguments$counts <- kept
          arguments$series <- "anchored"
          arguments$history <- history[[1L]]
          arguments$intermediate_weight <- 1
          first <- withCallingHandlers(do.call(fit_sparse_joint, arguments),
                                       alfak2_convergence_warning = function(w) invokeRestart("muffleWarning"))
          dispersion <- sparse_joint_intermediate_dispersion(first)
          if (is.finite(dispersion)) 1 / max(1, dispersion) else 1
        }
      }
    }
    chosen <- sparse_joint_auto_series(share, weight)
    arguments <- mget(names(formals(sys.function())), envir = environment())
    arguments$counts <- data
    arguments$series <- chosen
    fit <- do.call(fit_sparse_joint, arguments)
    fit$diagnostics$series_choice <- list(requested = "auto", series = chosen, ends_major = share, weight = weight,
                                          thresholds = sparse_joint_auto_thresholds)
    return(fit)
  }
  # The pairwise and state-space methods: every pair of consecutive samples
  # with its own domain, each sample with a cell in it. Samples are left out
  # until that holds (sparse_joint_pair_left_out()).
  other_method <- n_input > 2L && series != "anchored"
  if (other_method) {
    pairs <- sparse_joint_pair_left_out(counts, major_frequency, major_count, max_states,
                                        as.integer(min_cn), as.integer(max_cn))
    left_out <- pairs$left_out
    if (length(left_out)) {
      # `uncovered`: the sample without a cell in its pair's domain, the left
      # out one itself or a first or last sample next to it.
      detail <- data.frame(sample = left_out, time = alfak2_data_times(data)[left_out],
                           timepoint = colnames(counts)[left_out], cells = colSums(counts)[left_out],
                           uncovered = pairs$uncovered, stringsAsFactors = FALSE)
      reason <- paste("a pair of consecutive samples had a sample with no cell in its domain",
                      "(the sample itself, or a first or last sample next to it)")
    } else {
      coverage <- sparse_joint_check_pair_coverage(counts, major_frequency, major_count, max_states,
                                                   as.integer(min_cn), as.integer(max_cn))
    }
  } else {
    coverage <- sparse_joint_check_coverage(counts, major_frequency, major_count, max_states,
                                            as.integer(min_cn), as.integer(max_cn),
                                            sparse_joint_major_columns(data, list(intermediate = intermediate)),
                                            required = c(1L, n_input))
    # Anchored: a sample in between with no cell on the states it counts for
    # is left out (see sparse_joint_uncounted_samples()), and the rest fitted.
    left_out <- sparse_joint_uncounted_samples(coverage, intermediate)
    if (length(left_out)) {
      detail <- cbind(data.frame(sample = left_out, time = alfak2_data_times(data)[left_out]),
                      coverage[left_out, , drop = FALSE])
      reason <- sprintf("no cell on the %s", if (intermediate == "majors") "major karyotypes" else "model domain")
    }
  }
  if (length(left_out)) {
    dropped_samples <- detail
    rownames(dropped_samples) <- NULL
    warning(warningCondition(
      sprintf(paste0("Left out %d of %d samples (%s): %s. Such a sample carries no information ",
                     "under the model; see diagnostics$dropped_samples."),
              length(left_out), n_input, paste(colnames(counts)[left_out], collapse = ", "), reason),
      class = "alfak2_samples_left_out_warning", call = NULL))
    arguments <- mget(names(formals(sys.function())), envir = environment())
    arguments$counts <- sparse_joint_keep_samples(data, setdiff(seq_len(n_input), left_out))
    fit <- do.call(fit_sparse_joint, arguments)
    fit$diagnostics$dropped_samples <- dropped_samples
    # The input as given, for as_alfakR(): alfakR picks its karyotypes from
    # every sample.
    fit$input <- data
    return(fit)
  }
  # With intermediate samples for the majors only, their weight: the inverse of
  # the majors' dispersion in a first fit (at the shortest history), at most 1.
  dispersion <- NA_real_
  if (ncol(counts) > 2L && intermediate == "majors" && series == "anchored") {
    if (is.null(intermediate_weight)) {
      arguments <- mget(names(formals(sys.function())), envir = environment())
      arguments$counts <- data
      arguments$history <- history[[1L]]
      arguments$intermediate_weight <- 1
      first <- withCallingHandlers(do.call(fit_sparse_joint, arguments),
                                   alfak2_convergence_warning = function(w) invokeRestart("muffleWarning"))
      dispersion <- sparse_joint_intermediate_dispersion(first)
      intermediate_weight <- if (is.finite(dispersion)) 1 / max(1, dispersion) else 1
    }
    validate_scalar(intermediate_weight, "intermediate_weight", lower = .Machine$double.eps, upper = 1)
  }
  if (length(history) > 1L) {
    # Candidate histories differ only in the means of the balance prior, so
    # their marginal likelihoods compare directly: fit at each, keep the best.
    started <- proc.time()[["elapsed"]]
    arguments <- mget(names(formals(sys.function())), envir = environment())
    arguments$counts <- data
    # Only the chosen candidate's convergence is warned about.
    fits <- lapply(history, function(h) {
      arguments$history <- h
      withCallingHandlers(do.call(fit_sparse_joint, arguments),
                          alfak2_convergence_warning = function(w) invokeRestart("muffleWarning"))
    })
    objective <- vapply(fits, function(f) f$diagnostics$objective, numeric(1))
    convergence <- vapply(fits, function(f) f$diagnostics$convergence, numeric(1))
    fit <- fits[[sparse_joint_choose_candidate(objective, convergence)]]
    fit$diagnostics$history_profile <- data.frame(history = history, objective = objective,
                                                  convergence = convergence)
    if (is.finite(dispersion)) fit$diagnostics$intermediate_dispersion <- dispersion
    fit$diagnostics$seconds <- proc.time()[["elapsed"]] - started
    sparse_joint_convergence_warning(fit$diagnostics$optimisation)
    return(fit)
  }
  started <- proc.time()[["elapsed"]]

  settings <- list(birth_rate = birth_rate, birth_rate_log_sd = birth_rate_log_sd,
                   major_frequency = major_frequency, major_count = major_count, max_states = max_states,
                   min_cn = as.integer(min_cn), max_cn = as.integer(max_cn), shell_power = shell_power,
                   shell_sd = shell_sd, initial_frequency = initial_frequency, initial_state = initial_state,
                   history = history, balance_sd = balance_sd, calibration = calibration,
                   inner_control = inner_control)
  # Only with more than two samples, so that two-sample fits stay as they were.
  if (ncol(data$counts) > 2L && series == "anchored") {
    settings$intermediate <- intermediate
    if (intermediate == "majors") settings$intermediate_weight <- intermediate_weight
  } else if (ncol(data$counts) > 2L) {
    settings$series <- series
    # Internal options of the evaluation (validation/sparse_joint/series_methods).
    if (series == "pairwise") settings$pair_history <- getOption("alfak2.pair_history", "age")
    if (series == "state_space") settings$drift <- getOption("alfak2.state_space_drift")
  }
  # Internal: fit two timepoints with the template for more, to check that it
  # gives the same fit.
  if (identical(getOption("alfak2.sparse_joint_template"), "series")) settings$template <- "series"
  problem <- sparse_joint_problem(data, settings)
  domain <- problem$domain
  n <- length(domain$labels)
  kernel <- problem$kernel
  total <- problem$total
  centre <- problem$centre
  design <- problem$design
  scale_names <- if (is.null(problem$scale_names)) c("dosage", "curvature", "interaction", "residual") else
    problem$scale_names
  tmb_data <- problem$tmb_data
  parameters <- problem$parameters
  free_scale <- problem$free_scale
  map <- problem$map
  random <- problem$random
  objective <- TMB::MakeADFun(tmb_data, parameters, map = map, random = random, DLL = "alfak2",
                              silent = TRUE, inner.control = inner_control)
  latent_names <- names(objective$env$par)[objective$env$random]
  labels <- sparse_joint_latent_labels(
    latent_names,
    states = data.frame(karyotype = domain$labels, cells = apply(domain$y, 1L, paste, collapse = " -> "),
                        major = seq_len(n) <= domain$n_major, stringsAsFactors = FALSE),
    state_index = sparse_joint_latent_states(problem, latent_names)
  )
  optimum <- if (is.null(problem$scale_bounds)) {
    sparse_joint_optimise_scales(objective, control, labels = labels)
  } else {
    sparse_joint_optimise_scales(objective, control, lower = problem$scale_bounds$lower,
                                 upper = problem$scale_bounds$upper, labels = labels)
  }
  # Uncertainty is conditional on the fitted scales: the marginal likelihood
  # is often flat in a scale.
  solver <- sparse_joint_posterior(objective, tmb_data, parameters, map, random, optimum$par, inner_control)
  g_index <- which(solver$names == "g")
  g_hat <- cbind(solver$mode[g_index], sqrt(solver$variance(g_index)))
  coef_hat <- cbind(solver$mode[solver$names == "coef"])
  coef_covariance <- sparse_joint_coef_covariance(solver)
  # Scales not free keep their values (the prior mean, or a fixed drift).
  log_scale <- parameters$log_scale
  log_scale[free_scale] <- optimum$par
  scale <- stats::setNames(exp(log_scale), scale_names)

  self <- 2 * exp(kernel$log_self) - 1
  n_time <- ncol(domain$y)
  timepoint <- seq_len(n_time) - 1L
  summary <- cbind(
    data.frame(karyotype = domain$labels, stringsAsFactors = FALSE),
    stats::setNames(as.data.frame(domain$y), paste0("count_t", timepoint)),
    stats::setNames(as.data.frame(solver$expected), paste0("expected_t", timepoint)),
    data.frame(major = seq_len(n) <= domain$n_major,
               shell = sparse_joint_major_distance(domain$karyotypes, domain$major)),
    sparse_joint_fitness_columns(self, g_hat[, 1L], g_hat[, 2L])
  )
  # Relative to the count-weighted mean fitness of the majors.
  relative <- sparse_joint_relative(solver, g_index, summary$fitness_mean,
                                    summary$log_birth_sd, ifelse(summary$major, total, 0))
  summary$relative_fitness <- relative$mean
  summary$relative_sd <- sparse_joint_calibrate(relative$sd, relative$mean, summary$shell, birth_rate, calibration,
                                                observed = rowSums(problem$used) > 0)
  summary$relative_sd_model <- relative$sd
  observed_level <- sparse_joint_observed_level(summary)
  fit <- structure(list(
    data = data,
    summary = summary,
    dropped = domain$dropped,
    # What rebuilding the posterior needs: see sparse_joint_rebuild().
    settings = settings,
    posterior = c(sparse_joint_fill_latent(parameters, map, solver$mode, solver$names, random),
                  list(log_scale = log_scale)),
    trend = list(coefficients = coef_hat[, 1L], covariance = coef_covariance, centre = centre,
                 major = domain$major, shell_power = shell_power, scales = scale,
                 observed_level = observed_level, reference = relative$reference,
                 calibration = calibration, birth_rate = birth_rate),
    diagnostics = list(
      convergence = optimum$convergence, message = optimum$message, objective = optimum$objective,
      optimisation = optimum$check,
      finite_uncertainty = all(is.finite(g_hat[, 2L])) && all(is.finite(relative$sd)),
      n_states = n, n_major = domain$n_major, n_edges = length(kernel$from),
      dropped_karyotypes = length(domain$dropped), dropped_cells = domain$dropped_cells,
      domain_coverage = coverage,
      birth_rate = birth_rate, fitted_birth_rate = exp(coef_hat[1L, 1L]),
      initial_state = initial_state, history = history,
      min_cn = as.integer(min_cn), max_cn = as.integer(max_cn),
      seconds = proc.time()[["elapsed"]] - started
    )
  ), class = "alfak2_sparse_joint_fit")
  if (n_time > 2L) {
    if (!is.null(settings$series)) {
      fit$diagnostics$series <- settings$series
      if (settings$series == "pairwise") fit$diagnostics$pair_history <- settings$pair_history
      if (settings$series == "state_space") fit$diagnostics$drift <- scale[["drift"]]
    }
    fit$diagnostics$times <- alfak2_data_times(data)
    # None left out here; a caller that left samples out replaces it.
    fit$diagnostics$dropped_samples <- cbind(data.frame(sample = integer(0), time = numeric(0)), coverage[0L, ])
  }
  if (!is.null(settings$intermediate_weight)) fit$diagnostics$intermediate_weight <- settings$intermediate_weight
  if (is.finite(dispersion)) fit$diagnostics$intermediate_dispersion <- dispersion
  sparse_joint_convergence_warning(fit$diagnostics$optimisation)
  fit
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
    # Observed karyotypes left out of the domain count as observed; the rest
    # are predicted from the trend and were seen in neither sample.
    out$relative_sd[outside] <- sparse_joint_calibrate(out$relative_sd_model[outside], relative$mean,
                                                       out$shell[outside], trend$birth_rate, trend$calibration,
                                                       observed = left_out)
  }
  out
}

#' @export
print.alfak2_sparse_joint_fit <- function(x, ...) {
  d <- x$diagnostics
  cat("<alfak2_sparse_joint_fit>\n")
  cat("  states:", d$n_states, "(", d$n_major, "major )\n")
  cat("  observed karyotypes outside the domain:", d$dropped_karyotypes, "(", d$dropped_cells, "cells )\n")
  if (!is.null(d$domain_coverage)) cat("  cells in the domain:", sparse_joint_format_coverage(d$domain_coverage), "\n")
  cat("  baseline birth rate:", d$birth_rate, " convergence:", d$convergence, "\n")
  invisible(x)
}
