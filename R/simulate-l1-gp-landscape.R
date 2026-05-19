validate_grf_integer <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x)) {
    stop(sprintf("`%s` must be one integer-like finite number.", name), call. = FALSE)
  }
  out <- as.integer(x)
  if (out < lower || out > upper) {
    stop(sprintf("`%s` must be in [%s, %s].", name, lower, upper), call. = FALSE)
  }
  out
}

validate_grf_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("`%s` must be one finite number.", name), call. = FALSE)
  }
  lower_ok <- if (lower_open) x > lower else x >= lower
  upper_ok <- if (upper_open) x < upper else x <= upper
  if (!lower_ok || !upper_ok) {
    lbr <- if (lower_open) "(" else "["
    rbr <- if (upper_open) ")" else "]"
    stop(sprintf("`%s` must be in %s%s, %s%s.", name, lbr, lower, upper, rbr), call. = FALSE)
  }
  as.numeric(x)
}

validate_grf_vector <- function(x, name, n_chr) {
  if (!is.numeric(x) || length(x) != n_chr || any(!is.finite(x)) || any(x != round(x))) {
    stop(sprintf("`%s` must be an integer-like numeric vector of length `n_chr`.", name),
         call. = FALSE)
  }
  as.integer(x)
}

with_grf_seed <- function(seed, expr) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

format_grf_labels <- function(karyotypes) {
  if (exists("format_karyotypes", mode = "function")) {
    return(as.character(format_karyotypes(karyotypes)))
  }
  stop("`format_grf_labels()` requires `format_karyotypes()`.", call. = FALSE)
}

coerce_grf_karyotypes <- function(karyotypes, n_chr = NULL) {
  if (is.character(karyotypes)) {
    if (!exists("parse_karyotypes", mode = "function")) {
      stop("Character karyotypes require `parse_karyotypes()`.", call. = FALSE)
    }
    karyotypes <- parse_karyotypes(karyotypes)
  } else if (is.numeric(karyotypes) && is.null(dim(karyotypes))) {
    karyotypes <- matrix(as.integer(karyotypes), nrow = 1L)
  } else {
    karyotypes <- as.matrix(karyotypes)
  }
  if (!is.numeric(karyotypes) && !is.integer(karyotypes)) {
    stop("`karyotypes` must be a numeric matrix, numeric vector, or character labels.",
         call. = FALSE)
  }
  if (any(!is.finite(karyotypes)) || any(karyotypes != round(karyotypes))) {
    stop("`karyotypes` must contain integer-like finite copy numbers.", call. = FALSE)
  }
  storage.mode(karyotypes) <- "integer"
  if (!is.null(n_chr) && ncol(karyotypes) != n_chr) {
    stop(sprintf("`karyotypes` must have %d columns.", n_chr), call. = FALSE)
  }
  karyotypes
}

raw_grf_values <- function(karyotypes, centroids, lambda) {
  n <- nrow(karyotypes)
  m <- nrow(centroids)
  out <- numeric(n)
  for (i in seq_len(m)) {
    d <- sqrt(rowSums(sweep(karyotypes, 2L, centroids[i, ], FUN = "-")^2))
    out <- out + sin(d / lambda)
  }
  out / (pi * sqrt(m))
}

#' Simulate a paper-style lazy GRF truth landscape
#'
#' Constructs a lightweight Gaussian-random-field fitness oracle following the
#' ALFA-K synthetic-landscape form
#' `sum_i sin(||k - r_i||_2 / lambda) / (pi * sqrt(M))`. The object stores the
#' random centroids and scaling parameters only; karyotype fitness values are
#' computed on demand by `predict_landscape_fitness()`. No complete karyotype
#' lattice is generated.
#'
#' Fitness is anchored at `founder`: `fitness(founder) == founder_fitness`.
#'
#' @param n_chr Number of chromosomes.
#' @param n_centroids Number of random GRF centroids/waves.
#' @param lambda Positive wavelength parameter. Larger values produce smoother
#'   landscapes.
#' @param founder Integer founder/diploid karyotype. Defaults to all
#'   `founder_cn`.
#' @param founder_cn Copy number used to build the default founder.
#' @param founder_fitness Fitness assigned exactly to `founder`.
#' @param scale Multiplicative scale applied after founder-centering the raw
#'   GRF.
#' @param seed Reproducibility seed for centroid generation.
#' @param centroid_min,centroid_max Inclusive integer range used for random
#'   centroids.
#' @param min_cn,max_cn Viability bounds used by lazy simulators.
#' @param cache Whether to attach an environment cache for repeated fitness
#'   queries.
#'
#' @return An `alfak2_grf_landscape` object.
#' @export
simulate_grf_landscape <- function(
  n_chr = 22,
  n_centroids = 30,
  lambda = 0.8,
  founder = NULL,
  founder_cn = 2,
  founder_fitness = 1,
  scale = 1,
  seed = 1,
  centroid_min = -10,
  centroid_max = 20,
  min_cn = 1,
  max_cn = 8,
  cache = TRUE
) {
  n_chr <- validate_grf_integer(n_chr, "n_chr", lower = 1L)
  n_centroids <- validate_grf_integer(n_centroids, "n_centroids", lower = 1L)
  lambda <- validate_grf_scalar(lambda, "lambda", lower = 0, lower_open = TRUE)
  founder_cn <- validate_grf_integer(founder_cn, "founder_cn")
  founder_fitness <- validate_grf_scalar(founder_fitness, "founder_fitness")
  scale <- validate_grf_scalar(scale, "scale")
  seed <- validate_grf_integer(seed, "seed", lower = 0L)
  centroid_min <- validate_grf_integer(centroid_min, "centroid_min")
  centroid_max <- validate_grf_integer(centroid_max, "centroid_max")
  min_cn <- validate_grf_integer(min_cn, "min_cn")
  max_cn <- validate_grf_integer(max_cn, "max_cn")
  if (centroid_min > centroid_max) {
    stop("`centroid_min` must be <= `centroid_max`.", call. = FALSE)
  }
  if (min_cn > max_cn) {
    stop("`min_cn` must be <= `max_cn`.", call. = FALSE)
  }
  if (is.null(founder)) {
    founder <- rep.int(founder_cn, n_chr)
  } else {
    founder <- validate_grf_vector(founder, "founder", n_chr)
  }
  if (any(founder < min_cn | founder > max_cn)) {
    stop("`founder` must lie inside [`min_cn`, `max_cn`].", call. = FALSE)
  }
  if (!is.logical(cache) || length(cache) != 1L || is.na(cache)) {
    stop("`cache` must be TRUE or FALSE.", call. = FALSE)
  }

  centroids <- with_grf_seed(seed, {
    matrix(
      sample(seq.int(centroid_min, centroid_max), n_centroids * n_chr, replace = TRUE),
      nrow = n_centroids,
      ncol = n_chr
    )
  })
  storage.mode(centroids) <- "integer"
  founder_mat <- matrix(founder, nrow = 1L)
  founder_raw <- raw_grf_values(founder_mat, centroids, lambda)[1L]

  out <- list(
    family = "paper_sine_grf",
    n_chr = n_chr,
    n_centroids = n_centroids,
    lambda = lambda,
    centroids = centroids,
    founder = founder,
    founder_label = format_grf_labels(founder_mat),
    founder_fitness = founder_fitness,
    founder_raw = founder_raw,
    scale = scale,
    seed = seed,
    centroid_min = centroid_min,
    centroid_max = centroid_max,
    min_cn = min_cn,
    max_cn = max_cn,
    covariance = NA_character_,
    formula = "founder_fitness + scale * (sum_i sin(||k-r_i||_2/lambda)/(pi*sqrt(M)) - raw(founder))"
  )
  if (isTRUE(cache)) out$cache <- new.env(parent = emptyenv())
  class(out) <- c("alfak2_grf_landscape", "alfak2_landscape")
  out
}

#' Simulate a lazy GRF truth landscape
#'
#' Backward-compatible entry point for the package's synthetic landscape
#' generator. The implementation now follows the paper-style lazy GRF oracle and
#' no longer enumerates the full karyotype lattice.
#'
#' @inheritParams simulate_grf_landscape
#' @param ell Deprecated. Use `lambda`.
#' @param diploid_cn Deprecated alias for `founder_cn`.
#' @param diploid_fitness Deprecated alias for `founder_fitness`.
#' @param lower,upper,use_full_range,include_table Deprecated full-lattice
#'   options that are ignored by the lazy GRF implementation.
#'
#' @return An `alfak2_grf_landscape` object.
#' @export
simulate_l1_gp_landscape <- function(
  n_chr = 22,
  min_cn = 1,
  max_cn = 8,
  diploid_cn = 2,
  diploid_fitness = 1,
  lower = NULL,
  upper = NULL,
  ell = NULL,
  seed = 1,
  use_full_range = NULL,
  include_table = NULL,
  n_centroids = 30,
  lambda = if (is.null(ell)) 0.8 else ell,
  founder = NULL,
  founder_cn = diploid_cn,
  founder_fitness = diploid_fitness,
  scale = 1,
  centroid_min = -10,
  centroid_max = 20,
  cache = TRUE
) {
  simulate_grf_landscape(
    n_chr = n_chr,
    n_centroids = n_centroids,
    lambda = lambda,
    founder = founder,
    founder_cn = founder_cn,
    founder_fitness = founder_fitness,
    scale = scale,
    seed = seed,
    centroid_min = centroid_min,
    centroid_max = centroid_max,
    min_cn = min_cn,
    max_cn = max_cn,
    cache = cache
  )
}

#' Predict fitness from a lazy GRF landscape
#'
#' Computes paper-style GRF fitness values for requested karyotypes. The
#' function is deterministic for a fixed landscape object and karyotype.
#'
#' @param landscape Object returned by `simulate_grf_landscape()`.
#' @param karyotypes Numeric matrix/vector or dot-separated character labels.
#' @param use_cache Whether to reuse and populate `landscape$cache` when present.
#'
#' @return Numeric vector of fitness values.
#' @export
predict_landscape_fitness <- function(landscape, karyotypes, use_cache = TRUE) {
  if (!inherits(landscape, "alfak2_grf_landscape")) {
    stop("`landscape` must be an `alfak2_grf_landscape` object.", call. = FALSE)
  }
  if (!is.logical(use_cache) || length(use_cache) != 1L || is.na(use_cache)) {
    stop("`use_cache` must be TRUE or FALSE.", call. = FALSE)
  }
  mat <- coerce_grf_karyotypes(karyotypes, landscape$n_chr)
  labels <- format_grf_labels(mat)
  out <- rep(NA_real_, nrow(mat))

  cache <- landscape$cache
  can_cache <- isTRUE(use_cache) && is.environment(cache)
  missing <- rep(TRUE, length(out))
  if (can_cache) {
    for (i in seq_along(labels)) {
      if (exists(labels[i], envir = cache, inherits = FALSE)) {
        out[i] <- get(labels[i], envir = cache, inherits = FALSE)
        missing[i] <- FALSE
      }
    }
  }

  if (any(missing)) {
    raw <- raw_grf_values(mat[missing, , drop = FALSE], landscape$centroids, landscape$lambda)
    vals <- landscape$founder_fitness + landscape$scale * (raw - landscape$founder_raw)
    out[missing] <- vals
    if (can_cache) {
      miss_idx <- which(missing)
      for (j in seq_along(miss_idx)) {
        assign(labels[miss_idx[j]], vals[j], envir = cache)
      }
    }
  }
  out
}
