# The two multi-time-point methods of fit_sparse_joint() beside the anchored one
# (`series = "pairwise"` and `"state_space"`). Both give every pair of
# consecutive time points its own majors and domain, as two time points would have;
# their states are the union of the pairs' domains, majors first. The
# templates are src/sparse_joint_pairwise_model.h and
# src/sparse_joint_state_space_model.h.

sparse_joint_series_methods <- c("anchored", "pairwise", "state_space")

# Domains of the pairs of consecutive time points and their union. `member` and
# `pair_major` have one column per pair: whether each state is in that pair's
# domain, and whether it is one of the pair's majors. The other fields are
# those of sparse_joint_domain(), for the union.
sparse_joint_pair_domains <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn) {
  n_time <- ncol(counts)
  pairs <- lapply(seq_len(n_time - 1L), function(k) {
    y <- counts[, c(k, k + 1L), drop = FALSE]
    y <- y[rowSums(y) > 0, , drop = FALSE]
    sparse_joint_domain(y, major_frequency, major_count, max_states, min_cn, max_cn)
  })
  major_labels <- unique(unlist(lapply(pairs, function(d) d$labels[seq_len(d$n_major)])))
  labels <- union(major_labels, unlist(lapply(pairs, `[[`, "labels")))
  observed <- rownames(counts)
  in_data <- match(labels, observed)
  y <- matrix(0, length(labels), n_time)
  y[!is.na(in_data), ] <- counts[in_data[!is.na(in_data)], , drop = FALSE]
  outside <- !observed %in% labels
  member <- vapply(pairs, function(d) as.numeric(labels %in% d$labels), numeric(length(labels)))
  pair_major <- vapply(pairs, function(d) as.numeric(labels %in% d$labels[seq_len(d$n_major)]),
                       numeric(length(labels)))
  list(labels = labels, karyotypes = parse_karyotypes(labels), major = parse_karyotypes(major_labels),
       n_major = length(major_labels), y = y, y0 = y[, 1L], y1 = y[, 2L],
       dropped = observed[outside], dropped_cells = sum(counts[outside, ]),
       member = matrix(member, nrow = length(labels)), pair_major = matrix(pair_major, nrow = length(labels)))
}

# Coverage of each pair's own domain: one row per pair and time point, as
# sparse_joint_domain_coverage() gives for two time points, with a `pair` column.
sparse_joint_pair_coverage <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn) {
  rows <- lapply(seq_len(ncol(counts) - 1L), function(k) {
    y <- counts[, c(k, k + 1L), drop = FALSE]
    y <- y[rowSums(y) > 0, , drop = FALSE]
    domain <- sparse_joint_domain(y, major_frequency, major_count, max_states, min_cn, max_cn)
    cbind(data.frame(pair = k), sparse_joint_domain_coverage(y, domain, major_frequency, major_count))
  })
  do.call(rbind, rows)
}

# The parts of a problem that both methods share: domain, kernel, trend
# design and the TMB data of the landscape prior.
sparse_joint_series_method_base <- function(data, settings) {
  domain <- sparse_joint_pair_domains(data$counts, settings$major_frequency, settings$major_count,
                                      settings$max_states, as.integer(settings$min_cn), as.integer(settings$max_cn))
  kernel <- sparse_division_kernel(domain$karyotypes, data$beta)
  total <- rowSums(domain$y)
  centre <- colSums(domain$karyotypes * total) / sum(total)
  design <- sparse_joint_design(domain$karyotypes, centre, domain$major, settings$shell_power)
  times <- alfak2_data_times(data)
  tmb_data <- list(
    edge_from = kernel$from - 1L, edge_to = kernel$to - 1L,
    edge_log_prob = kernel$log_prob, log_self = kernel$log_self,
    design = design$design, coef_group = design$group,
    intercept_mean = log(settings$birth_rate), intercept_sd = settings$birth_rate_log_sd,
    shell_sd = settings$shell_sd,
    log_scale_prior_mean = rep(log(0.01), 4L), log_scale_prior_sd = rep(1.5, 4L),
    balance = as.integer(settings$initial_state == "balance"), balance_sd = settings$balance_sd,
    weight = rep(1, length(domain$labels))
  )
  list(domain = domain, kernel = kernel, total = total, centre = centre, design = design, times = times,
       interval = diff(times), tmb_data = tmb_data)
}

# Pairwise: each pair as two time points, sharing the landscape. A pair starts at
# balance over `settings$history` plus the time since the first time point (the
# population's age), unless `settings$pair_history` is "same".
sparse_joint_pairwise_problem <- function(data, settings) {
  base <- sparse_joint_series_method_base(data, settings)
  domain <- base$domain
  n <- length(domain$labels)
  n_pair <- ncol(domain$y) - 1L
  member <- domain$member
  y_start <- domain$y[, seq_len(n_pair), drop = FALSE] * member
  y_end <- domain$y[, seq_len(n_pair) + 1L, drop = FALSE] * member
  prior_mean <- prior_sd <- matrix(0, n, n_pair)
  for (k in seq_len(n_pair)) {
    inside <- member[, k] > 0
    prior <- sparse_joint_initial_prior(y_start[inside, k], y_end[inside, k], settings$initial_frequency)
    prior_mean[inside, k] <- prior$mean
    prior_sd[inside, k] <- prior$sd
  }
  prior_sd[member == 0] <- 1
  reference <- vapply(seq_len(n_pair), function(k) which.max((y_start[, k] + y_end[, k]) * member[, k]) - 1L,
                      integer(1))
  age <- if (identical(settings$pair_history, "same")) rep(0, n_pair) else base$times[seq_len(n_pair)]
  tmb_data <- c(list(model = "sparse_joint_pairwise", y = domain$y, y_start = y_start, y_end = y_end,
                     interval = base$interval, member = member, pair_major = domain$pair_major,
                     reference = reference, zeta_prior_mean = prior_mean, zeta_prior_sd = prior_sd,
                     pair_history = settings$history + age, balance_offset = matrix(0, n, n_pair)),
                base$tmb_data)
  theta <- log(colSums(y_end) / colSums(y_start))
  parameters <- list(
    zeta = prior_mean,
    g = rep(log(settings$birth_rate), n),
    coef = c(log(settings$birth_rate), rep(0, ncol(base$design$design) - 1L)),
    theta = theta,
    log_scale = tmb_data$log_scale_prior_mean
  )
  free_scale <- c(TRUE, TRUE, FALSE, TRUE)
  # Initial abundances outside a pair's domain are not part of the model.
  zeta_map <- factor(ifelse(member > 0, seq_len(n * n_pair), NA))
  list(domain = domain, kernel = base$kernel, total = base$total, centre = base$centre, design = base$design,
       used = domain$y, tmb_data = tmb_data, parameters = parameters, free_scale = free_scale,
       template = "pairwise", scale_names = c("dosage", "curvature", "interaction", "residual"),
       map = list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA)), zeta = zeta_map),
       random = c("zeta", "g", "coef", "theta"))
}

# State space: a latent abundance per state and time point, linked from time point
# to time point by the closed form plus drift. The drift is a fifth scale with a
# prior centred on 1e-4 per day (a bottleneck of about 2000 cells every 5
# days), or fixed at `settings$drift`.
sparse_joint_state_space_problem <- function(data, settings) {
  base <- sparse_joint_series_method_base(data, settings)
  domain <- base$domain
  n <- length(domain$labels)
  y <- domain$y
  n_time <- ncol(y)
  depth <- colSums(y)
  seen_first <- y[, 1L] > 0
  prior_mean <- ifelse(seen_first, log(pmax(y[, 1L], 1)), log(depth[[1L]]) + settings$initial_frequency[[1L]])
  prior_sd <- ifelse(seen_first, 10, settings$initial_frequency[[2L]])
  # Frequency of each state about each interval: the mean of its share at
  # the two time points, plus half a cell of the two.
  share <- sweep(y, 2L, depth, "/")
  frequency <- (share[, -n_time, drop = FALSE] + share[, -1L, drop = FALSE]) / 2
  frequency <- sweep(frequency, 2L, 1 / (depth[-n_time] + depth[-1L]), "+")
  drift_prior <- c(log(1e-4), 2.5)
  tmb_data <- c(list(model = "sparse_joint_state_space", y = y, interval = base$interval,
                     reference = which.max(base$total) - 1L, zeta_prior_mean = prior_mean, zeta_prior_sd = prior_sd,
                     n_major = domain$n_major, history = settings$history, balance_offset = rep(0, n),
                     drift_weight = 1 / frequency, max_drift_var = 4,
                     transition_offset = matrix(0, n, n_time - 1L)),
                base$tmb_data)
  tmb_data$log_scale_prior_mean <- c(tmb_data$log_scale_prior_mean, drift_prior[[1L]])
  tmb_data$log_scale_prior_sd <- c(tmb_data$log_scale_prior_sd, drift_prior[[2L]])
  later_seen <- y[, -1L, drop = FALSE] > 0
  zeta_later <- ifelse(later_seen, log(pmax(y[, -1L, drop = FALSE], 1)),
                       matrix(log(depth[-1L]) + settings$initial_frequency[[1L]], n, n_time - 1L, byrow = TRUE))
  log_scale <- tmb_data$log_scale_prior_mean
  if (!is.null(settings$drift)) log_scale[[5L]] <- log(settings$drift)
  parameters <- list(
    zeta = prior_mean,
    zeta_later = matrix(zeta_later, n, n_time - 1L),
    g = rep(log(settings$birth_rate), n),
    coef = c(log(settings$birth_rate), rep(0, ncol(base$design$design) - 1L)),
    theta = log(depth[-1L] / depth[-n_time]),
    log_scale = log_scale
  )
  free_scale <- c(TRUE, TRUE, FALSE, TRUE, is.null(settings$drift))
  # The drift lives on another scale than the landscape's: from no drift to
  # a deviation of order one per day.
  bounds <- list(lower = c(rep(log(1e-4), 3L), log(1e-9)), upper = c(rep(log(0.5), 3L), 0))
  bounds <- lapply(bounds, function(b) b[seq_len(sum(free_scale))])
  list(domain = domain, kernel = base$kernel, total = base$total, centre = base$centre, design = base$design,
       used = y, tmb_data = tmb_data, parameters = parameters, free_scale = free_scale, scale_bounds = bounds,
       template = "state_space", scale_names = c("dosage", "curvature", "interaction", "residual", "drift"),
       map = list(log_scale = factor(ifelse(free_scale, cumsum(free_scale), NA))),
       random = c("zeta", "zeta_later", "g", "coef", "theta"))
}

# Expected counts of the fitted model, one column per time point, and the data
# with the counts replaced by them (and the prior residuals that are not
# log-linear set to zero), for the expected information.
sparse_joint_series_method_expected <- function(tmb_data, parameters, map, mode, report) {
  latent <- function(name) {
    value <- parameters[[name]]
    use <- if (is.null(map[[name]])) seq_along(value) else which(!is.na(map[[name]]))
    value[use] <- unname(mode[names(mode) == name])
    value
  }
  expected <- tmb_data
  if (identical(tmb_data$model, "sparse_joint_pairwise")) {
    zeta <- latent("zeta")
    start <- exp(zeta) * tmb_data$member
    end <- exp(report$log_intensity) * tmb_data$member
    expected$y_start <- start
    expected$y_end <- end
    expected$balance_offset <- report$balance_residual
    # A time point's expectation: from the pair that starts there (the first),
    # else from the pair that ends there.
    counts <- cbind(start[, 1L], end)
    counts[tmb_data$member[, 1L] == 0, 1L] <- NA
    counts[, -1L][tmb_data$member == 0] <- NA
  } else {
    counts <- cbind(exp(latent("zeta")), exp(latent("zeta_later")))
    expected$y <- counts
    expected$balance_offset <- as.numeric(report$balance_residual)
    expected$transition_offset <- report$transition_residual
  }
  list(data = expected, counts = counts)
}

# The latent variables of a fit in the shapes of its parameters: entries
# mapped away keep their starting values.
sparse_joint_fill_latent <- function(parameters, map, mode, mode_names, names) {
  stats::setNames(lapply(names, function(name) {
    value <- parameters[[name]]
    estimate <- mode[mode_names == name]
    if (is.null(map[[name]])) {
      if (!is.matrix(value)) return(estimate)
      value[] <- estimate
      return(value)
    }
    value[!is.na(map[[name]])] <- estimate
    value
  }), names)
}

# Position among the latent variables of each state's log abundance at the
# first time point (NA where a state has none): `zeta` of the anchored and
# state-space models, the first pair's column for the pairwise one.
sparse_joint_first_sample_latent <- function(fit, solver) {
  zeta_index <- which(solver$names == "zeta")
  n <- nrow(fit$summary)
  if (!identical(fit$settings$series, "pairwise")) return(zeta_index[seq_len(n)])
  member <- matrix(fit$posterior$zeta, nrow = n) * 0 + NA
  problem <- sparse_joint_problem(fit$data, fit$settings)
  inside <- problem$tmb_data$member > 0
  member[inside] <- zeta_index
  member[, 1L]
}

# Coverage for the pairwise and state-space methods: each pair's own domain
# must hold a cell at each of its two time points. An error of class
# alfak2_domain_coverage_error names the pair otherwise.
sparse_joint_check_pair_coverage <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn) {
  coverage <- sparse_joint_pair_coverage(counts, major_frequency, major_count, max_states, min_cn, max_cn)
  if (any(coverage$cells_in_domain <= 0)) {
    named <- coverage
    named$sample <- coverage$pair
    sparse_joint_coverage_error(named, "A timepoint has no cell inside the model domain.", by_sample = TRUE,
                                unit = "Pair")
  }
  coverage
}

# Time points that the pairwise and state-space methods leave out. Each pair of
# consecutive time points has its own domain. A time point of a pair with no cell in
# that domain carries no information about it: its depth offset is free, as
# for an anchored time point in between (sparse_joint_uncounted_samples()).
# - Such a time point in between is left out, and its neighbours become a pair.
# - The first and the last time point are kept, as in the anchored method, where
#   they define the domain. If one of them has no cell in its pair's domain,
#   its neighbour in between is left out instead, so that it pairs with a
#   time point further in.
# This repeats, from the first uncovered pair on, until every pair is
# covered. At worst the first and the last time point remain: the two-time-point fit,
# which checks its own coverage. Returns the positions left out (in order)
# and, for each, the timepoint it was left out for.
sparse_joint_pair_left_out <- function(counts, major_frequency, major_count, max_states, min_cn, max_cn) {
  if (is.null(colnames(counts))) colnames(counts) <- paste0("t", seq_len(ncol(counts)) - 1L)
  kept <- seq_len(ncol(counts))
  left_out <- integer(0)
  uncovered <- character(0)
  while (length(kept) > 2L) {
    coverage <- sparse_joint_pair_coverage(counts[, kept, drop = FALSE], major_frequency, major_count, max_states,
                                           min_cn, max_cn)
    # Two rows per pair, its earlier and its later time point; position in `kept`.
    position <- coverage$pair + rep(c(0L, 1L), length.out = nrow(coverage))
    empty <- which(coverage$cells_in_domain <= 0)
    if (!length(empty)) break
    pair <- coverage$pair[[empty[[1L]]]]
    bad <- position[empty][coverage$pair[empty] == pair]
    inner <- bad[bad > 1L & bad < length(kept)]
    drop <- if (length(inner)) inner[[1L]] else setdiff(c(pair, pair + 1L), bad)[[1L]]
    left_out <- c(left_out, kept[[drop]])
    uncovered <- c(uncovered, colnames(counts)[kept[[bad[[1L]]]]])
    kept <- kept[-drop]
  }
  order <- order(left_out)
  list(left_out = left_out[order], uncovered = uncovered[order])
}

# The state of each latent variable of the pairwise and state-space methods,
# for naming flat directions (NULL for the others: position among its name).
sparse_joint_latent_states <- function(problem, names) {
  if (!problem$template %in% c("pairwise", "state_space")) return(NULL)
  n <- length(problem$domain$labels)
  index <- as.integer(stats::ave(seq_along(names), names, FUN = seq_along))
  out <- index
  if (problem$template == "pairwise") {
    zeta <- names == "zeta"
    out[zeta] <- ((which(problem$tmb_data$member > 0) - 1L) %% n + 1L)[index[zeta]]
  } else {
    later <- names == "zeta_later"
    out[later] <- (index[later] - 1L) %% n + 1L
  }
  out
}

# Each state's fitted log abundance at the first time point, as stored in the
# fit (NA where a state has none).
sparse_joint_first_sample_zeta <- function(fit) {
  if (!identical(fit$settings$series, "pairwise")) return(fit$posterior$zeta)
  member <- sparse_joint_problem(fit$data, fit$settings)$tmb_data$member[, 1L]
  ifelse(member > 0, as.matrix(fit$posterior$zeta)[, 1L], NA_real_)
}
