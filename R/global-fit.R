#' Fit the global graph Gaussian posterior
#'
#' Uses the local posterior as diagonal anchor observations and solves the
#' Gaussian graph posterior in compiled RcppEigen sparse linear algebra.
#'
#' @param local_fit An `alfak2_local_fit` object.
#' @param graph Bounded global `alfak2_graph`.
#' @param lambda_l_grid,lambda_e_grid,sigma_obs_grid Hyperparameter grids for
#'   compiled leave-one-anchor-out tuning.
#' @param eps Diagonal numerical stabilizer.
#'
#' @return A list containing node summaries and tuning diagnostics.
#' @export
fit_graph_posterior <- function(local_fit,
                                graph = NULL,
                                lambda_l_grid = c(0.2, 1, 5),
                                lambda_e_grid = c(0.05, 0.25, 1),
                                sigma_obs_grid = c(0.02, 0.05, 0.1),
                                eps = 1e-5) {
  if (!inherits(local_fit, "alfak2_local_fit")) {
    stop("`local_fit` must be an alfak2_local_fit object.", call. = FALSE)
  }
  if (is.null(graph)) graph <- local_fit$graph
  if (!inherits(graph, "alfak2_graph")) stop("`graph` must be an alfak2_graph object.", call. = FALSE)

  anchor_match <- match(local_fit$summary$karyotype, as.character(graph$labels))
  keep <- which(!is.na(anchor_match) & is.finite(local_fit$summary$fitness_mean))
  if (!length(keep)) stop("No local posterior anchors are present in the global graph.", call. = FALSE)
  anchor_var <- local_fit$summary$fitness_sd[keep]^2
  anchor_var[!is.finite(anchor_var) | anchor_var <= 0] <- stats::median(anchor_var[is.finite(anchor_var) & anchor_var > 0], na.rm = TRUE)
  anchor_var[!is.finite(anchor_var) | anchor_var <= 0] <- 0.25

  res <- alfak2_graph_posterior_cpp(
    karyotypes = graph$karyotypes,
    edge_from = as.integer(graph$edge_from),
    edge_to = as.integer(graph$edge_to),
    edge_weight = as.numeric(graph$edge_weight),
    anchor_index = as.integer(anchor_match[keep]),
    anchor_mean = as.numeric(local_fit$summary$fitness_mean[keep]),
    anchor_var = as.numeric(anchor_var),
    lambda_l_grid = as.numeric(lambda_l_grid),
    lambda_e_grid = as.numeric(lambda_e_grid),
    sigma_obs_grid = as.numeric(sigma_obs_grid),
    eps = eps
  )

  tier <- as.character(graph$support_tier)
  tier[graph$support_distance > 2L] <- "graph_borrowed"
  prior_dominated <- res$sd > stats::quantile(res$sd, 0.9, na.rm = TRUE) &
    !(tier %in% c("directly_informed", "local_borrowed", "weakly_supported"))
  tier[prior_dominated] <- "prior_dominated"
  summary <- data.frame(
    node_id = seq_along(graph$labels),
    karyotype = as.character(graph$labels),
    support_tier = tier,
    support_distance = as.integer(graph$support_distance),
    fitness_mean = as.numeric(res$mean),
    fitness_sd = as.numeric(res$sd),
    conf_low = as.numeric(res$mean) - 1.959963984540054 * as.numeric(res$sd),
    conf_high = as.numeric(res$mean) + 1.959963984540054 * as.numeric(res$sd),
    stringsAsFactors = FALSE
  )
  list(
    graph = graph,
    summary = summary,
    anchors = data.frame(
      node_id = as.integer(anchor_match[keep]),
      karyotype = local_fit$summary$karyotype[keep],
      mean = local_fit$summary$fitness_mean[keep],
      variance = anchor_var,
      stringsAsFactors = FALSE
    ),
    hyperparameters = list(
      lambda_l = res$lambda_l,
      lambda_e = res$lambda_e,
      sigma_obs = res$sigma_obs,
      cv_score = res$cv_score
    ),
    tuning_grid = res$grid,
    diagnostics = list(factorization_status = res$factorization_status)
  )
}
