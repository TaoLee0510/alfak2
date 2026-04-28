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
  counts <- validate_count_matrix(counts)
  validate_scalar(dt, "dt", lower = .Machine$double.eps)
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  parse_karyotypes(rownames(counts))
  colnames(counts) <- c("t0", "t1")
  new_alfak2_data(counts = counts, dt = dt, beta = beta, metadata = metadata)
}
