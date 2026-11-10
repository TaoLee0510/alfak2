#' Build a bounded karyotype graph
#'
#' @param x An `alfak2_data` object or character vector of observed karyotypes.
#' @param beta Missegregation rate. Defaults to `x$beta` for data objects.
#' @param transition_kernel Transition weighting used by graph edges and local
#'   transport. `"exact"` uses chromosome-level missegregation probabilities;
#'   `"linear"` preserves the historical first-order approximation.
#' @param shell_depth Number of one-step copy-number shells to add.
#' @param min_cn,max_cn Copy-number bounds.
#' @param max_nodes Hard bound on graph expansion.
#'
#' @return An `alfak2_graph` object.
#' @export
build_karyotype_graph <- function(x,
                                  beta = NULL,
                                  transition_kernel = c("exact", "linear"),
                                  shell_depth = 2,
                                  min_cn = 0,
                                  max_cn = 5,
                                  max_nodes = 5000) {
  transition_kernel <- match.arg(transition_kernel)
  if (inherits(x, "alfak2_data")) {
    labels <- x$labels
    y0 <- x$counts[, 1]
    y1 <- x$counts[, 2]
    holdout_labels <- if (!is.null(x$metadata$holdout_mode$labels)) {
      as.character(x$metadata$holdout_mode$labels)
    } else {
      character()
    }
    holdout_labels <- holdout_labels[nzchar(holdout_labels)]
    extra_labels <- setdiff(holdout_labels, labels)
    if (length(extra_labels)) {
      labels <- c(labels, extra_labels)
      y0 <- c(y0, stats::setNames(rep(0, length(extra_labels)), extra_labels))
      y1 <- c(y1, stats::setNames(rep(0, length(extra_labels)), extra_labels))
    }
    observation_weights <- x$metadata$observation_weights
    if (!is.null(observation_weights)) {
      weights <- matrix(0, nrow = length(labels), ncol = 2L,
                        dimnames = list(labels, c("t0", "t1")))
      observation_weights <- validate_observation_weights(observation_weights, x$counts)
      weights[rownames(observation_weights), ] <- observation_weights
      seed_mass <- as.numeric(y0) * weights[, 1] + as.numeric(y1) * weights[, 2]
      y0 <- as.integer(is.finite(seed_mass) & seed_mass > 0)
      y1 <- rep.int(0L, length(labels))
      names(y0) <- names(y1) <- labels
    }
    if (is.null(beta)) beta <- x$beta
  } else {
    labels <- as.character(x)
    parse_karyotypes(labels)
    y0 <- rep.int(1L, length(labels))
    y1 <- rep.int(0L, length(labels))
    if (is.null(beta)) beta <- 0.01
  }
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  max_nodes <- validate_grf_integer(max_nodes, "max_nodes", lower = 1L)
  if (length(labels) > max_nodes)
    stop("Observed and explicit target nodes exceed `max_nodes` before graph expansion.", call. = FALSE)
  graph <- alfak2_build_graph_cpp(
    labels = labels,
    y0 = as.integer(y0),
    y1 = as.integer(y1),
    beta = beta,
    transition_kernel = transition_kernel,
    shell_depth = as.integer(shell_depth),
    min_cn = as.integer(min_cn),
    max_cn = as.integer(max_cn),
    max_nodes = as.integer(max_nodes)
  )
  new_alfak2_graph(graph)
}

# A finite query domain can include distant states without enumerating every
# intervening shell of a high-dimensional lattice. Missing intermediate states
# remain computational outflow in division dynamics and are never folded back.
build_explicit_target_graph <- function(data, target_karyotypes, min_cn, max_cn, max_nodes,
                                        transition_kernel = "exact") {
  targets <- if (is.character(target_karyotypes)) target_karyotypes else
    format_karyotypes(as.matrix(target_karyotypes))
  targets <- unique(as.character(targets))
  if (!length(targets)) stop("`target_karyotypes` must contain at least one state.", call. = FALSE)
  parsed <- parse_karyotypes(unique(c(data$labels, targets)))
  if (any(parsed < min_cn | parsed > max_cn))
    stop("Explicit targets must be within the copy-number bounds.", call. = FALSE)
  augmented <- data
  holdout <- augmented$metadata$holdout_mode
  if (is.null(holdout)) holdout <- list()
  holdout$labels <- unique(c(holdout$labels, targets))
  augmented$metadata$holdout_mode <- holdout
  graph <- build_karyotype_graph(augmented, shell_depth = 0, min_cn = min_cn,
    max_cn = max_cn, max_nodes = max_nodes, transition_kernel = transition_kernel)
  weights <- data$metadata$observation_weights
  effective <- if (is.null(weights)) rowSums(data$counts) else
    rowSums(data$counts * validate_observation_weights(weights, data$counts))
  seeds <- parse_karyotypes(data$labels[effective > 0])
  if (!nrow(seeds)) stop("Explicit target graph requires a positive effective-count seed.", call. = FALSE)
  distance <- rep(Inf, length(graph$labels))
  for (i in seq_len(nrow(seeds)))
    distance <- pmin(distance, rowSums(abs(sweep(graph$karyotypes, 2L, seeds[i, ], "-"))))
  graph$support_distance <- as.integer(distance)
  graph$support_tier <- ifelse(distance == 0, "directly_informed",
    ifelse(distance == 1, "local_borrowed", ifelse(distance == 2, "weakly_supported", "graph_borrowed")))
  graph$domain_mode <- "explicit_targets"
  graph$target_labels <- targets
  graph
}
