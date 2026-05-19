#' ALFA-K2: Adaptive Local-to-Global Fitness Landscapes for Aneuploid Karyotypes
#'
#' `alfak2` implements ALFA-K2 (Adaptive Local-to-Global Fitness Landscapes for
#' Aneuploid Karyotypes). It provides methods to infer fitness landscapes for
#' aneuploid karyotypes from sparse two-timepoint count data and to study
#' karyotype evolution. The package builds bounded karyotype graphs, fits a
#' local hierarchical Bayesian model, and propagates estimates across the graph
#' with a Gaussian-process prior and ordered copy-number epistasis penalties.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib alfak2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(".data"))
