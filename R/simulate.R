grf_one_step_neighbours <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- coerce_grf_karyotypes(karyotypes)
  p <- ncol(karyotypes)
  out <- vector("list", nrow(karyotypes) * p * 2L)
  idx <- 1L
  for (i in seq_len(nrow(karyotypes))) {
    x <- karyotypes[i, ]
    for (j in seq_len(p)) {
      y <- x
      y[j] <- y[j] - 1L
      if (all(y >= min_cn) && all(y <= max_cn)) {
        out[[idx]] <- y
        idx <- idx + 1L
      }
      y <- x
      y[j] <- y[j] + 1L
      if (all(y >= min_cn) && all(y <= max_cn)) {
        out[[idx]] <- y
        idx <- idx + 1L
      }
    }
  }
  out <- out[seq_len(idx - 1L)]
  if (!length(out)) return(karyotypes[0L, , drop = FALSE])
  mat <- do.call(rbind, out)
  mat[!duplicated(format_grf_labels(mat)), , drop = FALSE]
}

grf_initial_karyotypes <- function(founder, depth, min_cn, max_cn) {
  current <- matrix(founder, nrow = 1L)
  if (depth <= 0L) return(current)
  all_nodes <- current
  frontier <- current
  for (i in seq_len(depth)) {
    frontier <- grf_one_step_neighbours(frontier, min_cn, max_cn)
    if (!nrow(frontier)) break
    all_nodes <- rbind(all_nodes, frontier)
    keep <- !duplicated(format_grf_labels(all_nodes))
    all_nodes <- all_nodes[keep, , drop = FALSE]
    frontier <- all_nodes
  }
  all_nodes
}

grf_multinomial_draw <- function(prob, size, overdispersion = 0) {
  if (length(prob) == 0L) return(integer(0))
  prob <- prob / sum(prob)
  if (overdispersion > 0) {
    alpha <- pmax(1e-10, overdispersion * prob)
    prob <- stats::rgamma(length(prob), shape = alpha, rate = 1)
    prob <- prob / sum(prob)
  }
  as.integer(stats::rmultinom(1L, size = size, prob = prob)[, 1L])
}

grf_population_env <- function(population) {
  env <- new.env(parent = emptyenv(), hash = TRUE)
  for (i in seq_along(population)) {
    assign(names(population)[i], as.numeric(population[i]), envir = env)
  }
  env
}

grf_env_population <- function(env) {
  keys <- ls(envir = env)
  if (!length(keys)) return(stats::setNames(numeric(0), character(0)))
  vals <- vapply(keys, function(k) get(k, envir = env, inherits = FALSE), numeric(1L))
  keep <- vals > 0
  stats::setNames(vals[keep], keys[keep])
}

grf_add_count <- function(env, key, value) {
  if (value == 0) return(invisible(NULL))
  old <- if (exists(key, envir = env, inherits = FALSE)) {
    get(key, envir = env, inherits = FALSE)
  } else {
    0
  }
  assign(key, old + value, envir = env)
  invisible(NULL)
}

grf_viable <- function(karyotype, min_cn, max_cn) {
  all(karyotype >= min_cn) && all(karyotype <= max_cn)
}

grf_missegregate_daughters <- function(karyotype, beta) {
  copies <- rep.int(seq_along(karyotype), karyotype)
  total <- length(copies)
  k <- stats::rbinom(1L, total, beta)
  while (k == 0L && beta > 0) {
    k <- stats::rbinom(1L, total, beta)
  }
  picked <- sample(copies, size = k, replace = FALSE)
  signs <- sample(c(-1L, 1L), size = k, replace = TRUE)
  delta <- integer(length(karyotype))
  for (i in seq_len(k)) delta[picked[i]] <- delta[picked[i]] + signs[i]
  list(karyotype + delta, karyotype - delta)
}

grf_downsample_population <- function(population, size, overdispersion = 0) {
  total <- sum(population)
  if (total <= size) return(population)
  draw <- grf_multinomial_draw(population / total, as.integer(size), overdispersion)
  out <- stats::setNames(as.numeric(draw), names(population))
  out[out > 0]
}

grf_sample_counts <- function(population, size, overdispersion = 0) {
  total <- sum(population)
  if (total <= 0) stop("Population is empty.", call. = FALSE)
  draw <- grf_multinomial_draw(population / total, as.integer(size), overdispersion)
  stats::setNames(draw, names(population))
}

#' Simulate sparse count data from a lazy GRF truth landscape
#'
#' Simulates karyotype evolution without enumerating the full copy-number
#' lattice. The simulator keeps only karyotypes currently present in the
#' population. When mutations create a new karyotype, its fitness is queried from
#' the lazy GRF oracle with `predict_landscape_fitness()`.
#'
#' @param landscape Object returned by `simulate_grf_landscape()`.
#' @param beta Per-chromosome missegregation probability.
#' @param dt Time interval between the two sampled count columns.
#' @param n0,n1 Sampling depths at the start and end of the interval.
#' @param detection_threshold Minimum total observed count retained.
#' @param dropout_prob Random node dropout probability after sampling.
#' @param seed Reproducibility seed.
#' @param init_concentration Concentration around the founder for the initial
#'   local population.
#' @param overdispersion Dirichlet-multinomial concentration; `0` means
#'   multinomial.
#' @param initial_population Number of cells in the simulated population before
#'   the first sample. Defaults to `5 * max(n0, n1)`.
#' @param initial_shell_depth Initial population support around the founder in
#'   one-missegregation steps.
#' @param time_step Internal ABM time step.
#' @param carrying_capacity Optional population cap. Use `Inf` to disable.
#' @param passage_fraction Fraction retained if `carrying_capacity` is exceeded.
#' @param max_unique Maximum number of unique live karyotypes allowed during
#'   simulation.
#'
#' @return A list containing a two-column count matrix and sparsity metadata.
#' @export
simulate_sparse_counts <- function(landscape,
                                   beta = 0.00005,
                                   dt = 1,
                                   n0 = 500,
                                   n1 = 500,
                                   detection_threshold = 1,
                                   dropout_prob = 0,
                                   seed = 1,
                                   init_concentration = 1.5,
                                   overdispersion = 0,
                                   initial_population = NULL,
                                   initial_shell_depth = 1,
                                   time_step = 0.1,
                                   carrying_capacity = Inf,
                                   passage_fraction = 0.2,
                                   max_unique = 50000) {
  if (!inherits(landscape, "alfak2_grf_landscape")) {
    stop("`landscape` must be an `alfak2_grf_landscape` object.", call. = FALSE)
  }
  beta <- validate_grf_scalar(beta, "beta", lower = 0, upper = 1)
  dt <- validate_grf_scalar(dt, "dt", lower = 0, lower_open = TRUE)
  n0 <- validate_grf_integer(n0, "n0", lower = 1L)
  n1 <- validate_grf_integer(n1, "n1", lower = 1L)
  detection_threshold <- validate_grf_integer(detection_threshold, "detection_threshold", lower = 0L)
  dropout_prob <- validate_grf_scalar(dropout_prob, "dropout_prob", lower = 0, upper = 1)
  seed <- validate_grf_integer(seed, "seed", lower = 0L)
  init_concentration <- validate_grf_scalar(init_concentration, "init_concentration", lower = 0)
  overdispersion <- validate_grf_scalar(overdispersion, "overdispersion", lower = 0)
  if (is.null(initial_population)) initial_population <- max(100L, 5L * max(n0, n1))
  initial_population <- validate_grf_integer(initial_population, "initial_population", lower = 1L)
  initial_shell_depth <- validate_grf_integer(initial_shell_depth, "initial_shell_depth", lower = 0L)
  time_step <- validate_grf_scalar(time_step, "time_step", lower = 0, lower_open = TRUE)
  if (!is.numeric(carrying_capacity) || length(carrying_capacity) != 1L ||
      is.na(carrying_capacity) || carrying_capacity < 1) {
    stop("`carrying_capacity` must be one positive number or Inf.", call. = FALSE)
  }
  passage_fraction <- validate_grf_scalar(passage_fraction, "passage_fraction",
                                          lower = 0, upper = 1,
                                          lower_open = TRUE, upper_open = TRUE)
  max_unique <- validate_grf_integer(max_unique, "max_unique", lower = 1L)

  with_grf_seed(seed, {
    init_k <- grf_initial_karyotypes(
      landscape$founder,
      initial_shell_depth,
      landscape$min_cn,
      landscape$max_cn
    )
    d2 <- rowSums(sweep(init_k, 2L, landscape$founder, FUN = "-")^2)
    init_prob <- exp(-init_concentration * d2)
    init_prob <- init_prob / sum(init_prob)
    init_counts <- grf_multinomial_draw(init_prob, initial_population, overdispersion = 0)
    population <- stats::setNames(as.numeric(init_counts), format_grf_labels(init_k))
    population <- population[population > 0]
    visited <- new.env(parent = emptyenv(), hash = TRUE)
    for (label in names(population)) assign(label, TRUE, envir = visited)

    y0 <- grf_sample_counts(population, n0, overdispersion)
    n_steps <- max(1L, ceiling(dt / time_step))
    step <- dt / n_steps

    for (step_id in seq_len(n_steps)) {
      labels <- names(population)
      karyotypes <- coerce_grf_karyotypes(labels, landscape$n_chr)
      fitness <- predict_landscape_fitness(landscape, karyotypes)
      next_env <- grf_population_env(population)

      for (i in seq_along(labels)) {
        count <- as.integer(round(population[i]))
        if (count <= 0L) next
        division_rate <- max(0, fitness[i])
        divisions <- min(count, stats::rpois(1L, count * division_rate * step))
        if (divisions <= 0L) next

        label <- labels[i]
        x <- karyotypes[i, ]
        grf_add_count(next_env, label, -divisions)

        faithful_prob <- (1 - beta)^sum(x)
        faithful <- stats::rbinom(1L, divisions, faithful_prob)
        if (faithful > 0L) grf_add_count(next_env, label, 2L * faithful)

        error_divisions <- divisions - faithful
        if (error_divisions > 0L) {
          for (j in seq_len(error_divisions)) {
            daughters <- grf_missegregate_daughters(x, beta)
            for (daughter in daughters) {
              if (grf_viable(daughter, landscape$min_cn, landscape$max_cn)) {
                daughter_label <- format_grf_labels(matrix(daughter, nrow = 1L))
                grf_add_count(next_env, daughter_label, 1)
                assign(daughter_label, TRUE, envir = visited)
              }
            }
          }
        }
      }

      population <- grf_env_population(next_env)
      if (!length(population)) {
        stop("All karyotypes went extinct during simulation.", call. = FALSE)
      }
      if (length(population) > max_unique) {
        stop("Simulation exceeded `max_unique`; reduce `dt`, `beta`, or population size.",
             call. = FALSE)
      }
      if (is.finite(carrying_capacity) && sum(population) > carrying_capacity) {
        retained <- max(1L, round(carrying_capacity * passage_fraction))
        population <- grf_downsample_population(population, retained, overdispersion = 0)
      }
    }

    y1 <- grf_sample_counts(population, n1, overdispersion)
    all_labels <- union(names(y0), names(y1))
    counts <- matrix(0L, nrow = length(all_labels), ncol = 2L,
                     dimnames = list(all_labels, c("t0", "t1")))
    counts[names(y0), "t0"] <- as.integer(y0)
    counts[names(y1), "t1"] <- as.integer(y1)

    totals <- rowSums(counts)
    keep <- totals >= detection_threshold
    if (dropout_prob > 0) keep <- keep & (stats::runif(length(keep)) >= dropout_prob)
    if (!any(keep)) keep[which.max(totals)] <- TRUE
    counts <- counts[keep, , drop = FALSE]

    obs_labels <- rownames(counts)
    obs_karyotypes <- coerce_grf_karyotypes(obs_labels, landscape$n_chr)
    truth_keep <- predict_landscape_fitness(landscape, obs_karyotypes)

    entropy <- 0.0
    p_obs <- rowSums(counts) / (n0 + n1)
    p_obs <- p_obs[p_obs > 0]
    if (length(p_obs)) entropy <- -sum(p_obs * log(p_obs))

    list(
      counts = counts,
      observed_labels = obs_labels,
      truth_observed = truth_keep,
      landscape = landscape,
      sparsity = list(
        observed_nodes = nrow(counts),
        visited_nodes = length(ls(envir = visited)),
        live_nodes = length(population),
        fraction_observed = nrow(counts) / max(1L, length(ls(envir = visited))),
        effective_support_entropy = entropy,
        sampling_depth = n0 + n1,
        detection_threshold = detection_threshold,
        dropout_prob = dropout_prob
      )
    )
  })
}
