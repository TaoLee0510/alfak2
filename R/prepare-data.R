validate_observation_weights <- function(weights, counts, name = "observation_weights") {
  if (is.null(weights)) return(NULL)
  if (is.data.frame(weights)) weights <- as.matrix(weights)
  if (is.null(dim(weights))) {
    if (length(weights) != nrow(counts)) {
      stop("`", name, "` vector length must match count rows.", call. = FALSE)
    }
    weights <- cbind(t0 = weights, t1 = weights)
    rownames(weights) <- rownames(counts)
  }
  if (!is.matrix(weights) || ncol(weights) != 2L) {
    stop("`", name, "` must be a matrix/data frame with two columns or a row vector.", call. = FALSE)
  }
  if (is.null(rownames(weights)) || any(!nzchar(rownames(weights)))) {
    stop("`", name, "` must have karyotype row names.", call. = FALSE)
  }
  missing <- setdiff(rownames(counts), rownames(weights))
  if (length(missing)) {
    stop("`", name, "` is missing weights for: ", paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  }
  weights <- weights[rownames(counts), , drop = FALSE]
  storage.mode(weights) <- "double"
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`", name, "` must contain finite non-negative values.", call. = FALSE)
  }
  colnames(weights) <- colnames(counts)
  weights
}

subset_observation_weights <- function(weights, labels) {
  if (is.null(weights)) return(NULL)
  weights <- validate_observation_weights(weights, matrix(0L, nrow = length(labels), ncol = 2L,
                                                          dimnames = list(labels, c("t0", "t1"))))
  weights
}

#' Prepare two-timepoint karyotype count data
#'
#' @param counts Matrix or data frame with two columns and dot-separated
#'   karyotype labels as row names.
#' @param dt Positive time interval between the two samples.
#' @param beta Missegregation rate used by graph transition kernels.
#' @param metadata Optional list stored on the returned object.
#'
#' @return An `alfak2_data` object.
#' @export
prepare_alfak2_data <- function(counts, dt = 1, beta = 0.00005, metadata = list()) {
  observation_weights <- attr(counts, "observation_weights", exact = TRUE)
  soft_minobs <- attr(counts, "soft_minobs", exact = TRUE)
  counts <- validate_count_matrix(counts)
  colnames(counts) <- c("t0", "t1")
  validate_scalar(dt, "dt", lower = .Machine$double.eps)
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  if (is.null(metadata)) metadata <- list()
  if (!is.list(metadata)) stop("`metadata` must be a list.", call. = FALSE)
  if (!is.null(observation_weights)) {
    metadata$observation_weights <- validate_observation_weights(observation_weights, counts)
  } else if (!is.null(metadata$observation_weights)) {
    metadata$observation_weights <- validate_observation_weights(metadata$observation_weights, counts)
  }
  if (!is.null(soft_minobs)) metadata$soft_minobs <- soft_minobs
  parse_karyotypes(rownames(counts))
  new_alfak2_data(counts = counts, dt = dt, beta = beta, metadata = metadata)
}
