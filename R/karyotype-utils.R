#' Parse karyotype labels
#'
#' Converts dot-separated copy-number labels such as `"2.2.3.1"` to an integer
#' matrix. Parsing is compiled and validates dimensional consistency.
#'
#' @param labels Character vector of karyotype labels.
#' @return Integer matrix with one row per label.
#' @keywords internal
parse_karyotypes <- function(labels) {
  alfak2_parse_karyotypes_cpp(as.character(labels))
}

#' Format karyotype matrix rows
#'
#' @param karyotypes Integer matrix.
#' @return Character vector of dot-separated labels.
#' @keywords internal
format_karyotypes <- function(karyotypes) {
  alfak2_stringify_karyotypes_cpp(as.matrix(karyotypes))
}
