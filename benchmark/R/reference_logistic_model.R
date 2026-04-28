reference_lattice_size <- function(n_chr, min_cn, max_cn) {
  states_per_chr <- max_cn - min_cn + 1L
  states_per_chr^n_chr
}

encode_reference_karyotypes <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- as.matrix(karyotypes)
  states_per_chr <- max_cn - min_cn + 1L
  powers <- states_per_chr^(seq_len(ncol(karyotypes)) - 1L)
  as.integer(1L + as.vector((karyotypes - min_cn) %*% powers))
}

enumerate_reference_lattice <- function(n_chr = 6, min_cn = 1, max_cn = 8) {
  n_states <- reference_lattice_size(n_chr, min_cn, max_cn)
  if (n_states > 5e6) {
    stop(
      "Requested reference lattice has ", format(n_states, big.mark = ","),
      " states. Use implicit or sampled benchmarking for this scale.",
      call. = FALSE
    )
  }
  states_per_chr <- max_cn - min_cn + 1L
  idx <- seq_len(n_states) - 1L
  mat <- matrix(0L, nrow = n_states, ncol = n_chr)
  for (chr in seq_len(n_chr)) {
    mat[, chr] <- min_cn + (idx %/% states_per_chr^(chr - 1L)) %% states_per_chr
  }
  mat
}

reference_karyotype_labels <- function(karyotypes) {
  apply(karyotypes, 1L, paste, collapse = ".")
}

generate_reference_fitness_landscape <- function(n_chr = 6,
                                                 min_cn = 1,
                                                 max_cn = 8,
                                                 diploid_fitness = 1,
                                                 family = c("structured_epistatic", "smooth_additive", "rugged"),
                                                 seed = 1) {
  family <- match.arg(family)
  set.seed(seed)
  karyotypes <- enumerate_reference_lattice(n_chr, min_cn, max_cn)
  labels <- reference_karyotype_labels(karyotypes)
  diploid <- rep(2L, n_chr)
  if (any(diploid < min_cn | diploid > max_cn)) {
    stop("Diploid state is outside the requested copy-number bounds.", call. = FALSE)
  }
  diploid_label <- paste(diploid, collapse = ".")
  diploid_index <- match(diploid_label, labels)

  centered <- sweep(karyotypes, 2L, diploid, "-")
  additive <- stats::rnorm(n_chr, mean = 0, sd = 0.055)
  quadratic <- stats::runif(n_chr, min = 0.012, max = 0.045)
  log_fitness <- as.vector(centered %*% additive) -
    as.vector((centered^2) %*% quadratic)

  if (family %in% c("structured_epistatic", "rugged")) {
    pair_count <- max(0L, n_chr - 1L)
    pair_effect <- stats::rnorm(pair_count, mean = 0, sd = 0.025)
    for (i in seq_len(pair_count)) {
      log_fitness <- log_fitness + pair_effect[i] * centered[, i] * centered[, i + 1L]
    }
  }

  mean_cn <- rowMeans(karyotypes)
  log_fitness <- log_fitness +
    ifelse(mean_cn >= 3.2, 0.04, 0) -
    ifelse(mean_cn <= 1.4, 0.10, 0) -
    0.015 * abs(mean_cn - 2)

  if (family == "rugged") {
    n_shock <- max(1L, floor(0.01 * length(log_fitness)))
    shock_idx <- sample.int(length(log_fitness), n_shock)
    log_fitness[shock_idx] <- log_fitness[shock_idx] + stats::rnorm(n_shock, 0, 0.18)
  }

  log_fitness <- log_fitness - log_fitness[diploid_index] + log(diploid_fitness)
  fitness <- exp(log_fitness)
  fitness[diploid_index] <- diploid_fitness

  structure(
    list(
      labels = labels,
      karyotypes = karyotypes,
      fitness = fitness,
      log_fitness = log_fitness,
      diploid_label = diploid_label,
      diploid_index = diploid_index,
      n_chr = n_chr,
      min_cn = min_cn,
      max_cn = max_cn,
      family = family,
      seed = seed,
      diploid_fitness = diploid_fitness
    ),
    class = "alfak2_reference_landscape"
  )
}

build_paired_missegregation_neighbors <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- as.matrix(karyotypes)
  n_chr <- ncol(karyotypes)
  event_pairs <- expand.grid(
    gain_chr = seq_len(n_chr),
    loss_chr = seq_len(n_chr),
    KEEP.OUT.ATTRS = FALSE
  )
  event_pairs <- event_pairs[event_pairs$gain_chr != event_pairs$loss_chr, , drop = FALSE]
  neighbor_index <- matrix(NA_integer_, nrow = nrow(karyotypes), ncol = nrow(event_pairs))

  for (event in seq_len(nrow(event_pairs))) {
    child <- karyotypes
    gain_chr <- event_pairs$gain_chr[event]
    loss_chr <- event_pairs$loss_chr[event]
    child[, gain_chr] <- child[, gain_chr] + 1L
    child[, loss_chr] <- child[, loss_chr] - 1L
    valid <- child[, gain_chr] <= max_cn & child[, loss_chr] >= min_cn
    neighbor_index[valid, event] <- encode_reference_karyotypes(child[valid, , drop = FALSE], min_cn, max_cn)
  }

  list(
    neighbor_index = neighbor_index,
    event_pairs = event_pairs
  )
}

step_logistic_missegregation <- function(counts,
                                         fitness,
                                         neighbors,
                                         pm,
                                         K,
                                         step_size,
                                         active_floor = 1e-10) {
  total <- sum(counts)
  if (!is.finite(total) || total <= 0) stop("Population went extinct.", call. = FALSE)
  logistic_factor <- max(0, 1 - total / K)
  if (logistic_factor <= 0) return(counts)

  active <- which(counts > active_floor)
  births <- counts[active] * fitness[active] * logistic_factor * step_size
  p_event <- min(1, length(unique(c(neighbors$event_pairs$gain_chr, neighbors$event_pairs$loss_chr))) * pm)
  p_event <- max(0, p_event)
  event_share <- if (ncol(neighbors$neighbor_index)) p_event / ncol(neighbors$neighbor_index) else 0

  next_counts <- counts
  next_counts[active] <- next_counts[active] + births * (1 - p_event)
  if (event_share > 0) {
    event_births <- births * event_share
    for (event in seq_len(ncol(neighbors$neighbor_index))) {
      child <- neighbors$neighbor_index[active, event]
      valid <- !is.na(child)
      if (!any(valid)) next
      add <- rowsum(event_births[valid], group = child[valid], reorder = FALSE)
      next_counts[as.integer(rownames(add))] <- next_counts[as.integer(rownames(add))] + as.numeric(add)
    }
  }
  next_counts
}

deterministic_census_counts <- function(weights, size) {
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) stop("Cannot census from an empty population.", call. = FALSE)
  expected <- weights / sum(weights) * size
  out <- floor(expected)
  remainder <- as.integer(size - sum(out))
  if (remainder > 0L) {
    frac <- expected - out
    add <- order(frac, decreasing = TRUE)[seq_len(remainder)]
    out[add] <- out[add] + 1L
  }
  as.integer(out)
}

round_reference_time <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be one finite number.", call. = FALSE)
  }
  out <- as.integer(round(x))
  if (out <= 0L) {
    stop("Rounded `", name, "` must be a positive integer.", call. = FALSE)
  }
  out
}

simulate_logistic_missegregation_reference <- function(landscape,
                                                       pm = 5e-4,
                                                       K = 1e9,
                                                       N0 = 100,
                                                       T = NULL,
                                                       dt = NULL,
                                                       step_size = 0.2,
                                                       census_size = K,
                                                       active_floor = 1e-8) {
  if (!inherits(landscape, "alfak2_reference_landscape")) {
    stop("`landscape` must be from generate_reference_fitness_landscape().", call. = FALSE)
  }
  raw_T <- if (is.null(T)) log((K - N0) / N0) / landscape$diploid_fitness else T
  raw_dt <- if (is.null(dt)) 0.8 * raw_T else dt
  T <- round_reference_time(raw_T, "T")
  dt <- round_reference_time(raw_dt, "dt")
  targets <- c(T, T + dt)
  total_time <- max(targets)

  neighbors <- build_paired_missegregation_neighbors(
    landscape$karyotypes,
    min_cn = landscape$min_cn,
    max_cn = landscape$max_cn
  )
  counts <- numeric(length(landscape$fitness))
  counts[landscape$diploid_index] <- N0

  snapshots <- vector("list", length(targets))
  names(snapshots) <- c("T", "T_plus_dt")
  time <- 0
  target_idx <- 1L
  while (time < total_time - 1e-12) {
    next_target <- targets[target_idx]
    h <- min(step_size, next_target - time, total_time - time)
    if (h <= 1e-12) {
      snapshots[[target_idx]] <- counts
      target_idx <- target_idx + 1L
      if (target_idx > length(targets)) break
      next
    }
    counts <- step_logistic_missegregation(
      counts = counts,
      fitness = landscape$fitness,
      neighbors = neighbors,
      pm = pm,
      K = K,
      step_size = h,
      active_floor = active_floor
    )
    time <- time + h
    if (abs(time - next_target) <= 1e-10) {
      snapshots[[target_idx]] <- counts
      target_idx <- target_idx + 1L
      if (target_idx > length(targets)) break
    }
  }
  if (any(vapply(snapshots, is.null, logical(1)))) {
    stop("Failed to record both reference timepoints.", call. = FALSE)
  }

  y0 <- deterministic_census_counts(snapshots[[1]], census_size)
  y1 <- deterministic_census_counts(snapshots[[2]], census_size)
  keep <- (y0 + y1) > 0L
  counts_mat <- cbind(t0 = y0[keep], t1 = y1[keep])
  rownames(counts_mat) <- landscape$labels[keep]

  list(
    counts = counts_mat,
    full_population = list(T = snapshots[[1]], T_plus_dt = snapshots[[2]]),
    time = list(T = T, dt = dt, T_plus_dt = T + dt),
    pm = pm,
    K = K,
    N0 = N0,
    census_size = census_size,
    landscape = landscape,
    sparsity = list(
      occupied_states = nrow(counts_mat),
      fraction_occupied = nrow(counts_mat) / length(landscape$labels),
      total_reference_states = length(landscape$labels)
    )
  )
}

stratified_downsample_counts <- function(counts,
                                         target_depth,
                                         seed = 1,
                                         detection_threshold = 1,
                                         dropout_prob = 0,
                                         strata_probs = c(0, 0.5, 0.9, 1)) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  set.seed(seed)
  out <- matrix(0L, nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
  for (j in seq_len(ncol(counts))) {
    positive <- which(counts[, j] > 0)
    if (!length(positive)) next
    abundance <- counts[positive, j]
    cuts <- unique(stats::quantile(log1p(abundance), probs = strata_probs, names = FALSE, type = 8))
    if (length(cuts) < 2L) {
      strata <- factor(rep(1L, length(positive)))
    } else {
      strata <- cut(log1p(abundance), breaks = cuts, include.lowest = TRUE, labels = FALSE)
    }
    split_idx <- split(positive, strata)
    weights <- sqrt(vapply(split_idx, length, integer(1)))
    alloc <- floor(target_depth * weights / sum(weights))
    remainder <- target_depth - sum(alloc)
    if (remainder > 0L) {
      alloc[seq_len(remainder)] <- alloc[seq_len(remainder)] + 1L
    }
    for (s in seq_along(split_idx)) {
      idx <- split_idx[[s]]
      if (!length(idx) || alloc[s] <= 0L) next
      prob <- counts[idx, j] / sum(counts[idx, j])
      out[idx, j] <- out[idx, j] + as.integer(stats::rmultinom(1L, size = alloc[s], prob = prob))
    }
  }

  keep <- rowSums(out) >= detection_threshold
  if (dropout_prob > 0) keep <- keep & (stats::runif(length(keep)) >= dropout_prob)
  if (!any(keep)) keep[which.max(rowSums(out))] <- TRUE
  out[keep, , drop = FALSE]
}

stratified_downsample_reference_population <- function(reference,
                                                       target_depth,
                                                       seed = 1,
                                                       detection_threshold = 1,
                                                       dropout_prob = 0,
                                                       strata_probs = c(0, 0.5, 0.9, 1)) {
  if (is.null(reference$full_population$T) ||
      is.null(reference$full_population$T_plus_dt) ||
      is.null(reference$landscape$labels)) {
    stop("`reference` must contain full_population at T/T_plus_dt and landscape labels.", call. = FALSE)
  }
  population_weights <- cbind(
    t0 = reference$full_population$T,
    t1 = reference$full_population$T_plus_dt
  )
  rownames(population_weights) <- reference$landscape$labels
  stratified_downsample_counts(
    population_weights,
    target_depth = target_depth,
    seed = seed,
    detection_threshold = detection_threshold,
    dropout_prob = dropout_prob,
    strata_probs = strata_probs
  )
}
