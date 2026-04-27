#' alfak2: compiled hierarchical karyotype fitness inference
#'
#' `alfak2` fits one inference pipeline for sparse two-timepoint karyotype count
#' data: a local hierarchical Bayesian posterior implemented in TMB followed by
#' a graph Gaussian posterior with ordered copy-number epistasis penalties solved
#' in compiled sparse linear algebra.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib alfak2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(".data"))
