#' Fit the complete alfak2 pipeline
#'
#' Runs data validation, local graph construction, TMB local posterior fitting,
#' bounded global graph expansion, and compiled graph Gaussian posterior solving.
#'
#' @param counts Two-column count matrix or an `alfak2_data` object.
#' @param dt Time interval if `counts` is not already prepared.
#' @param beta Missegregation rate.
#' @param local_shell_depth Local graph shell depth.
#' @param global_extra_shell Additional bounded shell depth for graph borrowing.
#' @param min_cn,max_cn Copy-number bounds.
#' @param max_nodes Hard graph-size bound.
#' @param lambda_l_grid,lambda_e_grid,sigma_obs_grid Hyperparameter grids passed
#'   to the compiled graph posterior.
#' @param ... Passed to `fit_local_posterior()`.
#'
#' @return An `alfak2_fit` object.
#' @export
fit_alfak2 <- function(counts,
                       dt = 1,
                       beta = 0.00005,
                       local_shell_depth = 2,
                       global_extra_shell = 1,
                       min_cn = 0,
                       max_cn = 5,
                       max_nodes = 5000,
                       lambda_l_grid = c(0.2, 1, 5),
                       lambda_e_grid = c(0.05, 0.25, 1),
                       sigma_obs_grid = c(0.02, 0.05, 0.1),
                       ...) {
  data <- if (inherits(counts, "alfak2_data")) counts else prepare_alfak2_data(counts, dt = dt, beta = beta)
  local_graph <- build_karyotype_graph(
    data,
    shell_depth = local_shell_depth,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  local <- fit_local_posterior(data, local_graph, ...)
  global_graph <- build_karyotype_graph(
    data,
    shell_depth = local_shell_depth + global_extra_shell,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  global <- fit_graph_posterior(
    local,
    global_graph,
    lambda_l_grid = lambda_l_grid,
    lambda_e_grid = lambda_e_grid,
    sigma_obs_grid = sigma_obs_grid
  )
  new_alfak2_fit(list(
    data = data,
    local = local,
    global = global,
    diagnostics = list(
      local = local$diagnostics,
      global = global$diagnostics,
      backend = c(local = "TMB", global = "RcppEigen")
    )
  ))
}
