#' Build a bounded karyotype graph
#'
#' @param x An `alfak2_data` object or character vector of observed karyotypes.
#' @param beta Missegregation rate. Defaults to `x$beta` for data objects.
#' @param shell_depth Number of one-step copy-number shells to add.
#' @param min_cn,max_cn Copy-number bounds.
#' @param max_nodes Hard bound on graph expansion.
#'
#' @return An `alfak2_graph` object.
#' @export
build_karyotype_graph <- function(x,
                                  beta = NULL,
                                  shell_depth = 2,
                                  min_cn = 0,
                                  max_cn = 5,
                                  max_nodes = 5000) {
  if (inherits(x, "alfak2_data")) {
    labels <- x$labels
    y0 <- x$counts[, 1]
    y1 <- x$counts[, 2]
    if (is.null(beta)) beta <- x$beta
  } else {
    labels <- as.character(x)
    parse_karyotypes(labels)
    y0 <- rep.int(1L, length(labels))
    y1 <- rep.int(0L, length(labels))
    if (is.null(beta)) beta <- 0.01
  }
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  graph <- alfak2_build_graph_cpp(
    labels = labels,
    y0 = as.integer(y0),
    y1 = as.integer(y1),
    beta = beta,
    shell_depth = as.integer(shell_depth),
    min_cn = as.integer(min_cn),
    max_cn = as.integer(max_cn),
    max_nodes = as.integer(max_nodes)
  )
  new_alfak2_graph(graph)
}
