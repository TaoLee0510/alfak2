validate_count_matrix <- function(counts) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || ncol(counts) != 2L) {
    stop("`counts` must be a matrix or data frame with exactly two columns.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
    stop("`counts` must have non-empty karyotype row names.", call. = FALSE)
  }
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("`counts` must contain finite non-negative values.", call. = FALSE)
  }
  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop("`counts` must contain integer counts.", call. = FALSE)
  }
  counts <- matrix(as.integer(round(counts)), nrow = nrow(counts),
                   dimnames = dimnames(counts))
  keep <- rowSums(counts) > 0L
  if (!any(keep)) stop("At least one karyotype must have non-zero counts.", call. = FALSE)
  counts[keep, , drop = FALSE]
}

validate_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be one finite number in [%s, %s].", name, lower, upper), call. = FALSE)
  }
  invisible(x)
}

match_observation_model <- function(x) {
  match.arg(x, c("multinomial", "dirichlet_multinomial"))
}
