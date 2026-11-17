#!/usr/bin/env Rscript
# Development-only paired process attribution. No reserved round-two test data.
# The finite process and deterministic mean start from the identical integer
# population; count samples use common seeds, with exactly identical t0 counts.

process_chromosome_probability <- function(parent, child, beta) {
  difference <- abs(parent - child)
  if (difference > parent) return(0)
  out <- 0
  for (errors in seq.int(difference, parent, by = 2L)) {
    lost <- (errors + parent - child) / 2
    out <- out + dbinom(errors, parent, beta) * dbinom(lost, errors, .5)
  }
  out
}

process_independent_operator <- function(karyotypes, fitness, beta) {
  states <- sort(unique(as.numeric(karyotypes)))
  single <- outer(states, states, Vectorize(function(i, j)
    process_chromosome_probability(i, j, beta)))
  probability <- matrix(1, nrow(karyotypes), nrow(karyotypes))
  for (chr in seq_len(ncol(karyotypes))) {
    index <- match(karyotypes[, chr], states)
    probability <- probability * single[index, index, drop = FALSE]
  }
  # P is a single-daughter substochastic kernel. Out-of-domain daughters die;
  # their probability is never redistributed onto the retained states.
  generator <- sweep(2 * t(probability) - diag(nrow(karyotypes)), 2L, fitness, '*')
  list(generator = generator, probability = probability)
}

process_deterministic_mean <- function(generator, initial_population, dt, tolerance = 1e-14) {
  # Uniformization/Taylor action on a Metzler generator: an independent positive
  # mean solver, rather than the inference code's time-stepping implementation.
  shift <- max(0, -min(diag(generator)))
  B <- generator + diag(shift, nrow(generator))
  term <- total <- as.numeric(initial_population)
  for (order in seq_len(1000L)) {
    term <- as.numeric(B %*% term) * dt / order
    total <- total + term
    if (max(abs(term)) <= tolerance * max(1, max(abs(total)))) break
  }
  if (order == 1000L) stop('Deterministic mean series did not converge')
  exp(-shift * dt) * total
}

process_exact_gillespie <- function(karyotypes, fitness, beta, dt, initial_population,
                                     seed, max_events = 1000000L) {
  labels <- apply(karyotypes, 1L, paste, collapse = '.')
  population <- as.numeric(initial_population)
  set.seed(seed)
  elapsed <- 0; events <- 0L
  repeat {
    rates <- population * fitness
    total_rate <- sum(rates)
    if (total_rate <= 0) break
    next_time <- elapsed + rexp(1, total_rate)
    if (next_time > dt) break
    elapsed <- next_time
    parent <- sample.int(length(rates), 1L, prob = rates)
    x <- karyotypes[parent, ]
    population[parent] <- population[parent] - 1
    errors <- rbinom(1, sum(x), beta)
    if (!errors) population[parent] <- population[parent] + 2 else {
      copies <- rep(seq_along(x), x)
      picked <- copies[sample.int(length(copies), errors)]
      signs <- sample(c(-1L, 1L), errors, replace = TRUE)
      delta <- tabulate(picked[signs > 0], length(x)) - tabulate(picked[signs < 0], length(x))
      for (daughter in list(x + delta, x - delta)) {
        index <- match(paste(daughter, collapse = '.'), labels)
        if (!is.na(index)) population[index] <- population[index] + 1
      }
    }
    events <- events + 1L
    if (events > max_events) stop('Gillespie max_events exceeded')
    if (!sum(population)) stop('Gillespie extinction')
  }
  list(population = population, events = events)
}

process_truth <- function(karyotypes, family, seed) {
  set.seed(seed)
  x <- karyotypes - 2L
  founder <- which(rowSums(abs(x)) == 0L)
  if (family == 'smooth') {
    value <- 1 + as.numeric(x %*% runif(ncol(x), -.16, .16)) - .06 * rowSums(x^2)
  } else if (family == 'rough') {
    centroids <- matrix(runif(16L * ncol(x), -5, 8), ncol = ncol(x))
    raw <- numeric(nrow(x))
    for (j in seq_len(nrow(centroids)))
      raw <- raw + sin(sqrt(rowSums(sweep(karyotypes, 2L, centroids[j, ], '-')^2)) / .22)
    value <- 1 + .18 * (raw - raw[founder]) / sd(raw)
  } else {
    coefficient <- runif(3, .12, .25) * sample(c(-1, 1), 3, replace = TRUE)
    value <- 1 + .06 * x[, 1L] + coefficient[1L] * x[, 1L] * x[, 2L] +
      coefficient[2L] * x[, 2L] * x[, 3L] + coefficient[3L] * x[, 1L] * x[, 2L] * x[, 3L]
  }
  stopifnot(all(is.finite(value)), all(value > 0))
  as.numeric(value)
}

prepare_process_development <- function(output) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(output, 'dataset_manifest.csv')))
    stop('Frozen process dataset already exists; use another directory to regenerate.')
  k <- as.matrix(expand.grid(rep(list(1:3), 3L)))
  storage.mode(k) <- 'integer'
  labels <- apply(k, 1L, paste, collapse = '.')
  rownames(k) <- labels
  distances <- as.matrix(dist(k, method = 'manhattan'))
  edges <- which(distances == 1 & row(distances) < col(distances), arr.ind = TRUE)
  founder <- which(rowSums(abs(k - 2L)) == 0L)
  p0 <- numeric(nrow(k)); p0[founder] <- .85
  p0[distances[founder, ] == 1L] <- .15 / 6L
  config <- list(n_chr = 3L, min_cn = 1L, max_cn = 3L, max_nodes = 27L,
    beta = .005, dt = 1.5, initial_population = 1200L, n0 = 150L, n1 = 150L,
    interval_levels = c(.5, .8, .95), protocol = 'round2_process_development_only',
    gauge = 'same t0-count reference in the two process arms',
    distance_definition = 'nearest union of observed nodes in each matched sampling pair',
    local_shell_depth = 1L,
    joint_control = list(main_sd = .07, curvature_sd = .04, interaction_sd = .04,
                         residual_sd = .125, edge_sd = .25, collapse_gaussian = FALSE))
  manifest <- discrepancy <- list()
  families <- c('smooth', 'rough', 'interaction')
  for (family_id in seq_along(families)) for (landscape in 1:2) {
    family <- families[[family_id]]
    landscape_id <- sprintf('process_%s_L%02d', family, landscape)
    truth_seed <- 260000L + family_id * 1000L + landscape
    truth <- process_truth(k, family, truth_seed)
    set.seed(truth_seed + 100000L)
    initial <- as.numeric(rmultinom(1, config$initial_population, p0))
    operator <- process_independent_operator(k, truth, config$beta)
    expected <- process_deterministic_mean(operator$generator, initial, config$dt)
    expected_half <- process_deterministic_mean(operator$generator, initial, config$dt / 2)
    expected_semigroup <- process_deterministic_mean(operator$generator, expected_half, config$dt / 2)
    stopifnot(max(abs(expected_semigroup - expected)) < 1e-8)
    abm <- process_exact_gillespie(k, truth, config$beta, config$dt, initial,
                                   seed = truth_seed + 200000L)
    pi0 <- initial / sum(initial)
    mean_pi1 <- expected / sum(expected)
    process_pi1 <- abm$population / sum(abm$population)
    discrepancy[[length(discrepancy) + 1L]] <- data.frame(
      family = family, landscape_id = landscape_id, events = abm$events,
      process_frequency_rmse = sqrt(mean((process_pi1 - mean_pi1)^2)),
      sampling_frequency_rmse_expected = sqrt(sum(mean_pi1 * (1 - mean_pi1)) / config$n1 / nrow(k)),
      expected_log_growth = log(sum(expected) / sum(initial)),
      process_log_growth = log(sum(abm$population) / sum(initial)),
      mean_semigroup_max_error = max(abs(expected_semigroup - expected)))
    for (sampling in 1:2) {
      seed <- truth_seed + 300000L + 100L * sampling
      arms <- list()
      for (mode in c('deterministic_mean', 'finite_population')) {
        endpoint <- if (mode == 'deterministic_mean') expected else abm$population
        set.seed(seed)
        y0 <- as.numeric(rmultinom(1, config$n0, pi0))
        y1 <- as.numeric(rmultinom(1, config$n1, endpoint))
        counts <- cbind(t0 = y0, t1 = y1); rownames(counts) <- labels
        arms[[mode]] <- list(counts = counts, observed = which(rowSums(counts) > 0), endpoint = endpoint)
      }
      stopifnot(identical(arms[[1L]]$counts[, 1L], arms[[2L]]$counts[, 1L]))
      union_observed <- sort(unique(c(arms[[1L]]$observed, arms[[2L]]$observed)))
      common_distance <- apply(distances[, union_observed, drop = FALSE], 1L, min)
      gauge <- arms[[1L]]$counts[, 1L] / sum(arms[[1L]]$counts[, 1L])
      pair_id <- sprintf('%s_S%d', landscape_id, sampling)
      for (mode in names(arms)) {
        arm <- arms[[mode]]
        dataset_id <- sprintf('%s_%s_S%d', landscape_id, mode, sampling)
        dataset <- list(dataset_id = dataset_id, pair_id = pair_id, family = family,
          landscape_id = landscape_id, split = 'development', mode = mode,
          sampling_replicate = sampling, sampling_seed = seed, abm_seed = truth_seed + 200000L,
          counts = arm$counts[arm$observed, , drop = FALSE], target_karyotypes = k,
          target_labels = labels, target_distance = common_distance,
          native_target_distance = apply(distances[, arm$observed, drop = FALSE], 1L, min),
          target_edges = edges, gauge_weight = gauge, truth = truth,
          truth_relative = truth - sum(gauge * truth), truth_pi0 = pi0,
          truth_pi1 = arm$endpoint / sum(arm$endpoint), deterministic_mean_pi1 = mean_pi1,
          latent_population_t0 = initial, latent_population_t1 = arm$endpoint,
          truth_seed = truth_seed, config = config)
        path <- file.path(output, paste0(dataset_id, '.rds'))
        saveRDS(dataset, path)
        manifest[[length(manifest) + 1L]] <- data.frame(dataset_id = dataset_id,
          pair_id = pair_id, landscape_id = landscape_id, family = family,
          mode = mode, sampling_replicate = sampling, split = 'development',
          observed_nodes = length(arm$observed), target_nodes = nrow(k),
          max_distance = max(common_distance), dataset_file = normalizePath(path))
      }
    }
  }
  write.csv(do.call(rbind, manifest), file.path(output, 'dataset_manifest.csv'), row.names = FALSE)
  write.csv(do.call(rbind, discrepancy), file.path(output, 'latent_process_discrepancy.csv'), row.names = FALSE)
  saveRDS(config, file.path(output, 'config.rds'))
  invisible(do.call(rbind, manifest))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  output <- if (length(args)) args[[1L]] else '/private/tmp/alfak2_iteration2_20261104_0518/process-data'
  manifest <- prepare_process_development(output)
  cat(sprintf('Frozen %d paired process datasets from six development landscapes in %s\n',
              nrow(manifest), normalizePath(output)))
}
