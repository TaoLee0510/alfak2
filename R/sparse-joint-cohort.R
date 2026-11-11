# Sparse joint count model for a cohort: several two-timepoint samples that
# share one landscape.
#
# A "node" is a karyotype of the shared landscape; a "state" is a node in one
# sample's domain. Each sample has its own domain, initial abundances and
# depth, as in `fit_sparse_joint()`. Samples are tied through the log birth
# rate of the nodes they have in common and through one trend.

#' Fit the sparse joint count model to a cohort
#'
#' Fits several two-timepoint samples at once. Each sample keeps its own
#' domain, as in [fit_sparse_joint()], and its own initial abundances; the
#' samples share the fitness of the karyotypes they have in common and one
#' trend. A karyotype that is abundant in one sample and an unobserved
#' neighbour in another is informed by both.
#'
#' With `sharing = "partial"` a karyotype's birth rate in a sample deviates
#' from the shared value by an amount whose scale is estimated, so that
#' samples which do not follow one landscape are not forced onto it. With
#' `sharing = "full"` it is the shared value.
#'
#' With `trend = "sample"` every sample has its own trend: dosage and
#' curvature about its own mean karyotype, and distance to its own major
#' karyotypes. Samples then share only the karyotypes they have in common and
#' the prior scales, and samples with no karyotype in common are fitted as if
#' alone. `trend = "shared"` fits one trend for the cohort, which helps when
#' the samples follow one landscape and on simulated cohorts from different
#' landscapes predicted much worse than fitting each sample alone.
#' `trend = "select"`, the default, fits both and keeps the one with the
#' higher marginal likelihood; on simulated cohorts it chose the shared trend
#' for samples of one landscape and separate trends otherwise. It takes twice
#' as long. The marginal likelihoods are in `diagnostics$trend_objective`.
#'
#' @param samples List of two-column count matrices with karyotype row names,
#'   or of `alfak2_data` objects. Names identify the samples.
#' @param dt Time between the two samples: one value, or one per sample.
#'   Used for count matrices only.
#' @param sharing `"partial"` or `"full"`; see above.
#' @param trend `"select"`, `"sample"` or `"shared"`; see above.
#' @param initial_state `"free"` or `"balance"`; see [fit_sparse_joint()].
#'   Unlike a single sample, a cohort starts free by default: with balance,
#'   the scale optimisation of some simulated cohorts failed and left no
#'   marginal likelihood to choose the trend by.
#' @inheritParams fit_sparse_joint
#'
#' @return An `alfak2_sparse_joint_cohort_fit` object. `landscape` has one row
#'   per karyotype in any sample's domain, with the shared fitness; `samples`
#'   has one summary per sample, as returned by [fit_sparse_joint()].
#' @export
fit_sparse_joint_cohort <- function(samples,
                                    dt = NULL,
                                    beta = 0.00005,
                                    birth_rate = 1,
                                    birth_rate_log_sd = 0.05,
                                    sharing = c("partial", "full"),
                                    trend = c("select", "sample", "shared"),
                                    major_frequency = 0.01,
                                    major_count = 3,
                                    max_states = 5000,
                                    min_cn = 1,
                                    max_cn = NULL,
                                    shell_power = 0.85,
                                    shell_sd = 0.5,
                                    initial_frequency = c(mean = log(1e-5), sd = 3),
                                    initial_state = c("free", "balance"),
                                    history = NULL,
                                    balance_sd = 1,
                                    calibration = sparse_joint_default_calibration,
                                    control = list(eval.max = 200, iter.max = 100),
                                    inner_control = list(maxit = 200, tol = 1e-10)) {
  if (!is.list(samples) || inherits(samples, "alfak2_data") || !length(samples)) {
    stop("`samples` must be a list of count matrices or `alfak2_data` objects.", call. = FALSE)
  }
  trend_mode <- match.arg(trend)
  if (trend_mode == "select") {
    # Fit with each trend and keep the one the data support better. Both
    # models have the same hyperparameters, so their Laplace marginal
    # likelihoods are comparable.
    arguments <- mget(names(formals(sys.function())), envir = environment())
    fits <- lapply(c(shared = "shared", sample = "sample"), function(mode) {
      tryCatch(do.call(fit_sparse_joint_cohort, utils::modifyList(arguments, list(trend = mode))),
               error = function(e) e)
    })
    failed <- vapply(fits, inherits, logical(1), what = "error")
    if (all(failed)) stop(conditionMessage(fits[[1L]]), call. = FALSE)
    objective <- vapply(fits, function(fit) if (inherits(fit, "error")) Inf else fit$diagnostics$objective, numeric(1))
    # A fit whose scale optimisation failed has no marginal likelihood to
    # compare; without one on either side there is nothing to choose by.
    usable <- is.finite(objective)
    if (!any(usable)) {
      stop("Neither trend has a finite marginal likelihood (the scale optimisation failed for both); ",
           "choose one with `trend = \"shared\"` or `trend = \"sample\"`.", call. = FALSE)
    }
    best <- fits[[which(usable)[which.min(objective[usable])]]]
    best$diagnostics$trend_objective <- objective
    best$diagnostics$seconds <- sum(vapply(fits[!failed], function(fit) fit$diagnostics$seconds, numeric(1)))
    return(best)
  }
  sharing <- match.arg(sharing)
  initial_state <- match.arg(initial_state)
  calibration <- sparse_joint_validate_calibration(calibration)
  n_sample <- length(samples)
  if (is.null(names(samples))) names(samples) <- paste0("sample", seq_len(n_sample))
  if (anyDuplicated(names(samples)) || any(!nzchar(names(samples)))) {
    stop("`samples` need distinct, non-empty names.", call. = FALSE)
  }
  is_data <- vapply(samples, inherits, logical(1), what = "alfak2_data")
  if (any(!is_data)) {
    if (is.null(dt) || !length(dt) %in% c(1L, n_sample)) {
      stop("`dt` must give one interval, or one per sample.", call. = FALSE)
    }
    dt <- rep_len(dt, n_sample)
  }
  data <- lapply(seq_len(n_sample), function(s) {
    if (is_data[[s]]) samples[[s]] else prepare_alfak2_data(samples[[s]], dt = dt[[s]], beta = beta)
  })
  beta <- data[[1L]]$beta
  if (any(vapply(data, function(d) !isTRUE(all.equal(d$beta, beta)), logical(1)))) {
    stop("All samples must have the same missegregation rate.", call. = FALSE)
  }
  validate_scalar(birth_rate, "birth_rate", lower = .Machine$double.eps)
  validate_scalar(birth_rate_log_sd, "birth_rate_log_sd", lower = .Machine$double.eps)
  validate_scalar(major_frequency, "major_frequency", lower = 0, upper = 1)
  validate_scalar(major_count, "major_count", lower = 2)
  validate_scalar(max_states, "max_states", lower = 1)
  validate_scalar(min_cn, "min_cn", lower = 0)
  validate_scalar(shell_power, "shell_power", lower = .Machine$double.eps)
  validate_scalar(shell_sd, "shell_sd", lower = .Machine$double.eps)
  validate_scalar(balance_sd, "balance_sd", lower = .Machine$double.eps)
  if (!is.numeric(initial_frequency) || length(initial_frequency) != 2L ||
      any(!is.finite(initial_frequency)) || initial_frequency[[2L]] <= 0) {
    stop("`initial_frequency` must give a finite mean and a positive standard deviation.", call. = FALSE)
  }
  intervals <- vapply(data, function(d) d$dt, numeric(1))
  history <- if (is.null(history)) intervals else rep_len(history, n_sample)
  if (any(!is.finite(history)) || any(history <= 0)) stop("`history` must be positive.", call. = FALSE)
  started <- proc.time()[["elapsed"]]

  observed_cn <- lapply(data, function(d) parse_karyotypes(rownames(d$counts)))
  if (is.null(max_cn)) max_cn <- max(unlist(observed_cn)) + 1L
  validate_scalar(max_cn, "max_cn", lower = min_cn)
  if (any(unlist(observed_cn) < min_cn | unlist(observed_cn) > max_cn)) {
    stop("Observed karyotypes must lie within the copy-number bounds.", call. = FALSE)
  }
  if (any(vapply(data, function(d) any(colSums(d$counts) <= 0), logical(1)))) {
    stop("Both timepoints of every sample need at least one cell.", call. = FALSE)
  }
  domains <- lapply(data, function(d) {
    sparse_joint_domain(d$counts, major_frequency, major_count, max_states, as.integer(min_cn), as.integer(max_cn))
  })
  if (any(vapply(domains, function(d) any(c(sum(d$y0), sum(d$y1)) <= 0), logical(1)))) {
    stop("A timepoint has no cell inside its sample's model domain.", call. = FALSE)
  }

  # States of all samples, concatenated; nodes are their distinct karyotypes.
  sizes <- vapply(domains, function(d) length(d$labels), integer(1))
  offset <- c(0L, cumsum(sizes))
  n_state <- sum(sizes)
  state_sample <- rep(seq_len(n_sample), sizes)
  state_label <- unlist(lapply(domains, `[[`, "labels"), use.names = FALSE)
  state_major <- unlist(lapply(domains, function(d) seq_along(d$labels) <= d$n_major), use.names = FALSE)
  y0 <- unlist(lapply(domains, `[[`, "y0"), use.names = FALSE)
  y1 <- unlist(lapply(domains, `[[`, "y1"), use.names = FALSE)
  node_labels <- unique(state_label)
  n_node <- length(node_labels)
  state_node <- match(state_label, node_labels)
  node_karyotypes <- parse_karyotypes(node_labels)
  kernels <- lapply(domains, function(d) sparse_division_kernel(d$karyotypes, beta))
  bind <- function(name, shift = FALSE) {
    unlist(lapply(seq_len(n_sample), function(s) kernels[[s]][[name]] + if (shift) offset[[s]] else 0L),
           use.names = FALSE)
  }
  log_self <- bind("log_self")

  # Trend designs. A sample's own trend is centred on its count-weighted mean
  # karyotype and measures distance to its own majors. The cohort's trend is
  # centred on the average of those means, weighting samples equally, and
  # measures distance to the nearest major of any sample.
  major_labels <- unique(state_label[state_major])
  major <- parse_karyotypes(major_labels)
  sample_centre <- lapply(domains, function(d) {
    total <- d$y0 + d$y1
    colSums(d$karyotypes * total) / sum(total)
  })
  centre <- Reduce(`+`, sample_centre) / n_sample
  shared_trend <- trend_mode == "shared"
  trend_centre <- if (shared_trend) rep(list(centre), n_sample) else sample_centre
  trend_major <- if (shared_trend) rep(list(major), n_sample) else lapply(domains, `[[`, "major")
  designs <- lapply(seq_len(n_sample), function(s) {
    sparse_joint_design(domains[[s]]$karyotypes, trend_centre[[s]], trend_major[[s]], shell_power)
  })
  design <- list(design = do.call(rbind, lapply(designs, `[[`, "design")), group = designs[[1L]]$group)
  n_coef <- ncol(design$design)
  initial <- lapply(domains, function(d) sparse_joint_initial_prior(d$y0, d$y1, initial_frequency))
  reference <- vapply(seq_len(n_sample), function(s) {
    offset[[s]] + which.max(domains[[s]]$y0 + domains[[s]]$y1)
  }, numeric(1))

  deviations <- sharing == "partial"
  scale_names <- c("dosage", "curvature", "interaction", "residual", "deviation")
  tmb_data <- list(
    model = "sparse_joint_cohort",
    y0 = y0, y1 = y1,
    state_sample = state_sample - 1L, state_node = state_node - 1L,
    state_major = as.integer(state_major), dt = intervals,
    edge_from = bind("from", TRUE) - 1L, edge_to = bind("to", TRUE) - 1L,
    edge_log_prob = bind("log_prob"), log_self = log_self,
    design = design$design, coef_group = design$group,
    reference = as.integer(reference) - 1L,
    zeta_prior_mean = unlist(lapply(initial, `[[`, "mean"), use.names = FALSE),
    zeta_prior_sd = unlist(lapply(initial, `[[`, "sd"), use.names = FALSE),
    intercept_mean = log(birth_rate), intercept_sd = birth_rate_log_sd, shell_sd = shell_sd,
    log_scale_prior_mean = rep(log(0.01), 5L),
    log_scale_prior_sd = rep(1.5, 5L),
    deviations = as.integer(deviations),
    balance = as.integer(initial_state == "balance"), history = history, balance_sd = balance_sd,
    balance_offset = rep(0, n_state), shared_trend = as.integer(shared_trend)
  )
  parameters <- list(
    zeta = tmb_data$zeta_prior_mean,
    m = rep(log(birth_rate), n_node),
    g = rep(log(birth_rate), n_state),
    coef = matrix(c(log(birth_rate), rep(0, n_coef - 1L)), n_coef, n_sample),
    theta1 = vapply(domains, function(d) log(sum(d$y1) / sum(d$y0)), numeric(1)),
    log_scale = tmb_data$log_scale_prior_mean
  )
  # No interaction terms in the design; no deviations under full sharing.
  free_scale <- c(TRUE, TRUE, FALSE, TRUE, deviations)
  map <- list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA)))
  if (!deviations) map$g <- factor(rep(NA, n_state))
  # One trend: every column of `coef` is the same parameters.
  if (shared_trend) map$coef <- factor(rep(seq_len(n_coef), n_sample))
  random <- c("zeta", "m", if (deviations) "g", "coef", "theta1")
  objective <- TMB::MakeADFun(tmb_data, parameters, map = map, random = random, DLL = "alfak2",
                              silent = TRUE, inner.control = inner_control)
  optimum <- sparse_joint_optimise_scales(objective, control)
  solver <- sparse_joint_posterior(objective, tmb_data, parameters, map, random, optimum$par, inner_control)
  node_latent <- which(solver$names == "m")
  state_latent <- if (deviations) which(solver$names == "g") else node_latent[state_node]
  m_hat <- cbind(solver$mode[node_latent], sqrt(solver$variance(node_latent)))
  # Positions of each sample's coefficients among the latent variables.
  coef_latent <- which(solver$names == "coef")
  coef_index <- lapply(seq_len(n_sample), function(s) {
    if (shared_trend) coef_latent else coef_latent[(s - 1L) * n_coef + seq_len(n_coef)]
  })
  state_hat <- if (deviations) {
    cbind(solver$mode[state_latent], sqrt(solver$variance(state_latent)))
  } else m_hat[state_node, , drop = FALSE]
  sample_trend <- lapply(seq_len(n_sample), function(s) {
    list(centre = trend_centre[[s]], major = trend_major[[s]], coefficients = solver$mode[coef_index[[s]]],
         covariance = sparse_joint_coef_covariance(solver, coef_index[[s]]))
  })
  names(sample_trend) <- names(samples)
  scale <- tmb_data$log_scale_prior_mean
  scale[free_scale] <- optimum$par
  scale <- stats::setNames(exp(scale), scale_names)
  if (!deviations) scale[["deviation"]] <- 0

  node_self <- 2 * sparse_joint_self_probability(node_karyotypes, beta) - 1
  in_domain_of <- tabulate(state_node, n_node)
  observed_in_domain <- tabulate(state_node[y0 + y1 > 0], n_node)
  landscape <- cbind(
    data.frame(karyotype = node_labels, n_samples = in_domain_of,
               n_major = tabulate(state_node[state_major], n_node),
               n_observed = observed_in_domain,
               cells = as.numeric(tapply(y0 + y1, factor(state_node, levels = seq_len(n_node)), sum)),
               shell = sparse_joint_major_distance(node_karyotypes, major),
               stringsAsFactors = FALSE),
    sparse_joint_fitness_columns(node_self, m_hat[, 1L], m_hat[, 2L])
  )
  # The cohort's reference weights every sample's majors by their share of
  # that sample's major cells, and every sample equally.
  major_share <- unlist(lapply(seq_len(n_sample), function(s) {
    rows <- offset[[s]] + seq_len(sizes[[s]])
    weight <- ifelse(state_major[rows], y0[rows] + y1[rows], 0)
    weight / sum(weight)
  }), use.names = FALSE)
  node_weight <- as.numeric(tapply(major_share, factor(state_node, levels = seq_len(n_node)), sum))
  relative <- sparse_joint_relative(solver, node_latent, landscape$fitness_mean, landscape$log_birth_sd, node_weight)
  # Covariance of each sample's trend coefficients with the cohort reference.
  relative$reference$sample_coef_covariance <- lapply(coef_index, function(index) relative$covariance[index])
  landscape$relative_fitness <- relative$mean
  landscape$relative_sd <- sparse_joint_calibrate(relative$sd, relative$mean, landscape$shell, birth_rate,
                                                  calibration)
  landscape$relative_sd_model <- relative$sd
  references <- vector("list", n_sample)
  summaries <- lapply(seq_len(n_sample), function(s) {
    rows <- offset[[s]] + seq_len(sizes[[s]])
    domain <- domains[[s]]
    summary <- cbind(
      data.frame(karyotype = domain$labels, count_t0 = domain$y0, count_t1 = domain$y1,
                 expected_t0 = solver$expected_t0[rows], expected_t1 = solver$expected_t1[rows],
                 major = state_major[rows],
                 shell = sparse_joint_major_distance(domain$karyotypes, domain$major),
                 stringsAsFactors = FALSE),
      sparse_joint_fitness_columns(2 * exp(log_self[rows]) - 1, state_hat[rows, 1L], state_hat[rows, 2L])
    )
    # Relative to the sample's own majors.
    own <- sparse_joint_relative(solver, state_latent[rows], summary$fitness_mean, summary$log_birth_sd,
                                 major_share[rows])
    summary$relative_fitness <- own$mean
    summary$relative_sd <- sparse_joint_calibrate(own$sd, own$mean, summary$shell, birth_rate, calibration)
    summary$relative_sd_model <- own$sd
    # Covariance of every node's shared value, and of the sample's trend
    # coefficients, with this sample's reference.
    own$reference$coef_covariance <- own$covariance[coef_index[[s]]]
    references[[s]] <<- c(own$reference, list(node_covariance = own$covariance[node_latent]))
    summary
  })
  names(summaries) <- names(samples)
  structure(list(
    data = stats::setNames(data, names(samples)),
    landscape = landscape,
    samples = summaries,
    dropped = stats::setNames(lapply(domains, `[[`, "dropped"), names(samples)),
    trend = list(mode = trend_mode, samples = sample_trend, centre = centre,
                 major = major, shell_power = shell_power, scales = scale,
                 observed_level = lapply(summaries, sparse_joint_observed_level),
                 reference = relative$reference, sample_reference = references,
                 calibration = calibration, birth_rate = birth_rate),
    diagnostics = list(
      convergence = optimum$convergence, message = optimum$message, objective = optimum$objective,
      finite_uncertainty = all(is.finite(state_hat[, 2L])) && all(is.finite(m_hat[, 2L])) &&
        all(is.finite(landscape$relative_sd)),
      n_samples = n_sample, n_nodes = n_node, n_states = n_state, n_major = length(major_labels),
      shared_nodes = sum(in_domain_of > 1L),
      dropped_karyotypes = sum(lengths(lapply(domains, `[[`, "dropped"))),
      dropped_cells = sum(vapply(domains, `[[`, numeric(1), "dropped_cells")),
      sharing = sharing, trend = trend_mode, birth_rate = birth_rate,
      fitted_birth_rate = exp(mean(vapply(sample_trend, function(t) t$coefficients[[1L]], numeric(1)))),
      initial_state = initial_state, beta = beta,
      min_cn = as.integer(min_cn), max_cn = as.integer(max_cn),
      seconds = proc.time()[["elapsed"]] - started
    )
  ), class = "alfak2_sparse_joint_cohort_fit")
}

#' Predict fitness from a cohort fit
#'
#' Without `sample`, returns the shared landscape: the posterior of a
#' karyotype in any sample's domain, and the trend elsewhere. With `sample`,
#' returns fitness in that sample: its own posterior inside its domain; the
#' shared value, widened by the between-sample deviation, for a karyotype that
#' only other samples have in their domain; and otherwise as
#' [predict.alfak2_sparse_joint_fit()] does for a single sample.
#'
#' @param object An `alfak2_sparse_joint_cohort_fit`.
#' @param karyotypes Character labels or an integer matrix of karyotypes.
#' @param sample Name or index of a sample, or `NULL` for the shared landscape.
#' @param ... Unused.
#'
#' @return Data frame with `karyotype`, `fitness_mean`, `fitness_sd`,
#'   `conf_low`, `conf_high`, `relative_fitness`, `relative_sd`,
#'   `relative_sd_model`, `shell` and `source`, which is `"sample"`,
#'   `"landscape"`, `"observed_level"` or `"trend"`. Relative fitness is
#'   measured from the count-weighted mean of the sample's major karyotypes,
#'   or of all samples' when `sample` is `NULL`; `shell` is the distance to the
#'   nearest major karyotype of the sample, or of any sample.
#' @export
predict.alfak2_sparse_joint_cohort_fit <- function(object, karyotypes, sample = NULL, ...) {
  labels <- if (is.character(karyotypes)) karyotypes else format_karyotypes(as.matrix(karyotypes))
  trend <- object$trend
  shared <- is.null(sample)
  if (!shared) {
    if (is.character(sample)) sample <- match(sample, names(object$samples))
    if (length(sample) != 1L || is.na(sample) || sample < 1L || sample > length(object$samples)) {
      stop("`sample` must name one sample of the cohort.", call. = FALSE)
    }
  }
  reference <- if (shared) trend$reference else trend$sample_reference[[sample]]
  cn <- parse_karyotypes(labels)
  self <- 2 * sparse_joint_self_probability(cn, object$diagnostics$beta) - 1
  # In a sample, a karyotype deviates from the shared value.
  deviation <- if (shared) 0 else trend$scales[["deviation"]]

  # The trend everywhere, then overwritten by better-informed sources. In a
  # sample it is that sample's trend. For the shared landscape it is the
  # average of the samples' trends; the variance taken is the average of
  # theirs, which is at least the variance of the average.
  trend_part <- function(s, reference_covariance) {
    part <- trend$samples[[s]]
    design <- sparse_joint_design(cn, part$centre, part$major, trend$shell_power)$design
    list(mean = as.numeric(design %*% part$coefficients),
         variance = rowSums((design %*% part$covariance) * design),
         covariance = as.numeric(design %*% reference_covariance))
  }
  parts <- if (shared) {
    lapply(seq_along(trend$samples), function(s) trend_part(s, reference$sample_coef_covariance[[s]]))
  } else list(trend_part(sample, reference$coef_covariance))
  average <- function(name) Reduce(`+`, lapply(parts, `[[`, name)) / length(parts)
  log_birth <- average("mean")
  log_birth_sd <- sqrt(average("variance") + trend$scales[["residual"]]^2 + deviation^2)
  covariance <- average("covariance")
  source <- rep("trend", length(labels))

  node <- match(labels, object$landscape$karyotype)
  known <- !is.na(node)
  log_birth[known] <- object$landscape$log_birth_mean[node[known]]
  log_birth_sd[known] <- sqrt(object$landscape$log_birth_sd[node[known]]^2 + deviation^2)
  # Covariance of the node's log birth rate with the reference: stored for a
  # sample; for the shared landscape it is implied by the stored relative sd.
  if (!shared) covariance[known] <- reference$node_covariance[node[known]]
  source[known] <- "landscape"

  if (!shared) {
    # An observed karyotype left out of this sample's domain sits at the level
    # of the sample's observed neighbours, unless another sample observed it
    # inside its own domain and so informs it directly.
    level <- trend$observed_level[[sample]]
    left_out <- labels %in% object$dropped[[sample]] & !is.null(level) &
      !(known & object$landscape$n_observed[ifelse(known, node, 1L)] > 0)
    if (any(left_out)) {
      log_birth[left_out] <- level[["log_birth"]]
      log_birth_sd[left_out] <- level[["sd"]]
      covariance[left_out] <- 0
      source[left_out] <- "observed_level"
    }
  }
  fitness <- sparse_joint_fitness_columns(self, log_birth, log_birth_sd)
  out <- cbind(
    data.frame(karyotype = labels, stringsAsFactors = FALSE),
    fitness[c("fitness_mean", "fitness_sd", "conf_low", "conf_high")],
    data.frame(
      relative_fitness = fitness$fitness_mean - reference$fitness,
      relative_sd = NA_real_,
      relative_sd_model = sqrt(pmax(fitness$fitness_mean^2 * log_birth_sd^2 + reference$variance -
                                      2 * fitness$fitness_mean * covariance, 0)),
      # Distance to the nearest major of the sample, or of any sample.
      shell = sparse_joint_major_distance(cn, if (shared) trend$major else {
        own <- object$samples[[sample]]
        parse_karyotypes(own$karyotype[own$major])
      }),
      source = source, stringsAsFactors = FALSE
    )
  )
  out$relative_sd <- sparse_joint_calibrate(out$relative_sd_model, out$relative_fitness, out$shell,
                                            trend$birth_rate, trend$calibration)
  columns <- c("fitness_mean", "fitness_sd", "conf_low", "conf_high", "relative_fitness", "relative_sd",
               "relative_sd_model")
  if (shared) {
    out[known, columns] <- object$landscape[node[known], columns]
  } else {
    own <- object$samples[[sample]]
    state <- match(labels, own$karyotype)
    inside <- !is.na(state)
    out[inside, columns] <- own[state[inside], columns]
    out$source[inside] <- "sample"
  }
  out
}

#' @export
print.alfak2_sparse_joint_cohort_fit <- function(x, ...) {
  d <- x$diagnostics
  cat("<alfak2_sparse_joint_cohort_fit>\n")
  cat("  samples:", d$n_samples, " karyotypes:", d$n_nodes, "(", d$shared_nodes, "in more than one sample )\n")
  cat("  sharing:", d$sharing, " trend:", d$trend, " between-sample deviation:",
      signif(x$trend$scales[["deviation"]], 3), "\n")
  cat("  baseline birth rate:", d$birth_rate, " convergence:", d$convergence, "\n")
  invisible(x)
}
