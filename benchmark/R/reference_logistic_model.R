encode_reference_karyotypes <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- as.matrix(karyotypes)
  states_per_chr <- max_cn - min_cn + 1L
  powers <- states_per_chr^(seq_len(ncol(karyotypes)) - 1L)
  as.integer(1L + as.vector((karyotypes - min_cn) %*% powers))
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
  required <- c("labels", "karyotypes", "fitness", "min_cn", "max_cn")
  if (!is.list(landscape) || !all(required %in% names(landscape))) {
    stop("`landscape` must contain labels, karyotypes, fitness, min_cn, and max_cn.", call. = FALSE)
  }
  diploid_cn <- if (!is.null(landscape$diploid_cn)) landscape$diploid_cn else 2L
  diploid_index <- if (!is.null(landscape$diploid_index)) {
    landscape$diploid_index
  } else {
    which(rowSums(abs(as.matrix(landscape$karyotypes) - diploid_cn)) == 0)
  }
  if (length(diploid_index) != 1L) {
    stop("Could not identify exactly one diploid state in `landscape`.", call. = FALSE)
  }
  diploid_fitness <- if (!is.null(landscape$diploid_fitness)) {
    landscape$diploid_fitness
  } else {
    landscape$fitness[diploid_index]
  }
  raw_T <- if (is.null(T)) log((K - N0) / N0) / diploid_fitness else T
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
  counts[diploid_index] <- N0

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
