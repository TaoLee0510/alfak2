validate_l1_gp_integer <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x)) {
    stop(sprintf("`%s` must be one integer-like finite number.", name), call. = FALSE)
  }
  out <- as.integer(x)
  if (out < lower || out > upper) {
    stop(sprintf("`%s` must be in [%s, %s].", name, lower, upper), call. = FALSE)
  }
  out
}

validate_l1_gp_scalar <- function(x, name, lower = -Inf, upper = Inf,
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

mode_multiply_l1_gp <- function(arr, mat, mode) {
  dims <- dim(arr)
  nd <- length(dims)
  ord <- c(mode, setdiff(seq_len(nd), mode))
  arr_perm <- aperm(arr, ord)
  mat_arr <- matrix(arr_perm, nrow = dims[mode])
  res <- mat %*% mat_arr
  res_arr <- array(res, dim = dim(arr_perm))
  inv_ord <- match(seq_len(nd), ord)
  aperm(res_arr, inv_ord)
}

format_l1_gp_labels <- function(karyotypes) {
  if (exists("format_karyotypes", mode = "function")) {
    out <- try(format_karyotypes(karyotypes), silent = TRUE)
    if (!inherits(out, "try-error")) return(as.character(out))
  }
  apply(karyotypes, 1L, paste, collapse = ".")
}

#' Simulate an L1-correlated Gaussian random-field truth landscape
#'
#' Generates a complete in silico benchmark fitness landscape directly on the
#' original bounded karyotype copy-number lattice. No UMAP, PCA, t-SNE, or other
#' embedding is used. Correlation decays exponentially with Manhattan distance:
#' `Cov(F(k), F(k')) = exp(-Manhattan(k, k') / ell)`, where `ell` controls
#' smoothness. The implementation avoids constructing the full covariance matrix
#' by using the separable/Kronecker structure of this covariance across
#' chromosome dimensions.
#'
#' The diploid karyotype is anchored exactly to `diploid_fitness`, and all
#' fitness values are scaled to lie inside `[lower, upper]`. The return object is
#' compatible with `simulate_sparse_counts()`.
#'
#' @param n_chr Number of chromosomes.
#' @param min_cn,max_cn Copy-number bounds.
#' @param diploid_cn Diploid copy number used as the anchor state.
#' @param diploid_fitness Fitness value assigned exactly to the diploid state.
#' @param lower,upper Lower and upper bounds for generated fitness values.
#' @param ell Positive Manhattan-distance length scale.
#' @param seed Reproducibility seed.
#' @param use_full_range Fraction of the feasible centered range to use.
#' @param include_table Whether to include a long table with one row per
#'   karyotype.
#'
#' @return An `alfak2_landscape` object with `labels`, `karyotypes`, `fitness`,
#'   metadata, and optionally `table`.
#' @export
simulate_l1_gp_landscape <- function(
  n_chr = 6,
  min_cn = 1,
  max_cn = 8,
  diploid_cn = 2,
  diploid_fitness = 1,
  lower = -5,
  upper = 5,
  ell = 2.5,
  seed = 1,
  use_full_range = 0.95,
  include_table = TRUE
) {
  n_chr <- validate_l1_gp_integer(n_chr, "n_chr", lower = 1L)
  min_cn <- validate_l1_gp_integer(min_cn, "min_cn")
  max_cn <- validate_l1_gp_integer(max_cn, "max_cn")
  diploid_cn <- validate_l1_gp_integer(diploid_cn, "diploid_cn")
  seed <- validate_l1_gp_integer(seed, "seed", lower = 0L)
  if (min_cn > max_cn) stop("`min_cn` must be <= `max_cn`.", call. = FALSE)
  if (diploid_cn < min_cn || diploid_cn > max_cn) {
    stop("`diploid_cn` must satisfy `min_cn <= diploid_cn <= max_cn`.", call. = FALSE)
  }
  lower <- validate_l1_gp_scalar(lower, "lower")
  upper <- validate_l1_gp_scalar(upper, "upper")
  diploid_fitness <- validate_l1_gp_scalar(diploid_fitness, "diploid_fitness")
  if (!(lower < diploid_fitness && diploid_fitness < upper)) {
    stop("`lower < diploid_fitness < upper` must hold.", call. = FALSE)
  }
  ell <- validate_l1_gp_scalar(ell, "ell", lower = 0, lower_open = TRUE)
  use_full_range <- validate_l1_gp_scalar(
    use_full_range,
    "use_full_range",
    lower = 0,
    upper = 1,
    lower_open = TRUE
  )
  if (!is.logical(include_table) || length(include_table) != 1L || is.na(include_table)) {
    stop("`include_table` must be TRUE or FALSE.", call. = FALSE)
  }

  states <- seq.int(min_cn, max_cn)
  m <- length(states)
  n_states <- m^n_chr
  if (!is.finite(n_states) || n_states > .Machine$integer.max) {
    stop("Requested lattice is too large for an explicit R landscape object.", call. = FALSE)
  }

  D1 <- abs(outer(states, states, "-"))
  K1 <- exp(-D1 / ell)
  K1 <- K1 + 1e-10 * diag(m)
  L1 <- t(chol(K1))

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
  arr <- array(stats::rnorm(n_states), dim = rep(m, n_chr))
  for (j in seq_len(n_chr)) {
    arr <- mode_multiply_l1_gp(arr, L1, j)
  }
  raw <- as.vector(arr)

  grid <- expand.grid(
    replicate(n_chr, states, simplify = FALSE),
    KEEP.OUT.ATTRS = FALSE
  )
  names(grid) <- paste0("chr", seq_len(n_chr))
  karyotypes <- as.matrix(grid)
  storage.mode(karyotypes) <- "integer"
  labels <- format_l1_gp_labels(karyotypes)

  if (nrow(karyotypes) != length(raw) || length(labels) != length(raw)) {
    stop("Internal lattice enumeration failed to align with the Gaussian field.", call. = FALSE)
  }
  diploid_idx <- which(rowSums(abs(karyotypes - diploid_cn)) == 0)
  if (length(diploid_idx) != 1L) {
    stop("Could not identify exactly one diploid anchor state.", call. = FALSE)
  }

  residual <- raw - raw[diploid_idx]
  rmax <- max(residual)
  rmin <- min(residual)
  if (!is.finite(rmax) || !is.finite(rmin)) {
    stop("Generated Gaussian field contains non-finite values.", call. = FALSE)
  }
  pos_scale <- if (rmax > 0) (upper - diploid_fitness) / rmax else Inf
  neg_scale <- if (rmin < 0) (diploid_fitness - lower) / abs(rmin) else Inf
  scale <- use_full_range * min(pos_scale, neg_scale)
  if (!is.finite(scale)) {
    fitness <- rep(diploid_fitness, length(raw))
  } else {
    fitness <- diploid_fitness + scale * residual
  }
  fitness[diploid_idx] <- diploid_fitness
  fitness <- pmin(pmax(fitness, lower), upper)
  fitness[diploid_idx] <- diploid_fitness

  out <- list(
    labels = labels,
    karyotypes = karyotypes,
    fitness = as.numeric(fitness),
    family = "l1_correlated_gaussian_random_field",
    seed = seed,
    min_cn = min_cn,
    max_cn = max_cn,
    n_chr = n_chr,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    use_full_range = use_full_range,
    covariance = "exp(-Manhattan_distance / ell)"
  )
  if (isTRUE(include_table)) {
    out$table <- data.frame(
      label = labels,
      as.data.frame(karyotypes, stringsAsFactors = FALSE),
      fitness = as.numeric(fitness),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(out$table)[seq_len(n_chr) + 1L] <- paste0("chr", seq_len(n_chr))
  }
  class(out) <- "alfak2_landscape"
  out
}
