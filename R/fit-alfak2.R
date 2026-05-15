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
#' @param graph_edge_weight Edge weighting used by the global graph prior.
#' @param anchor_support_tiers Local support tiers used as global graph anchors;
#'   defaults to `"all"` to preserve the legacy public API behavior.
#' @param anchor_covariance_inflation Named variance multipliers for local
#'   covariance statuses.
#' @param anchor_count_reference Optional low-count downweighting reference for
#'   global graph anchors.
#' @param anchor_count_power Exponent for count-based anchor variance inflation.
#' @param input_depth `"raw"` uses the supplied count matrix directly. `"effective"`
#'   converts each timepoint to frequencies and re-counts them at a controlled
#'   effective depth before fitting.
#' @param effective_depth Effective sequencing depth used when
#'   `input_depth = "effective"`. Required for `"cap"` and `"fixed"` modes; ignored
#'   by `"min"` mode.
#' @param effective_depth_mode Effective-depth rule. `"min"` uses the smaller raw
#'   timepoint depth for both columns, `"cap"` caps each raw depth at
#'   `effective_depth`, and `"fixed"` uses `effective_depth` for every timepoint.
#' @param observation_model Optional observation model passed to
#'   `fit_local_posterior()`. If omitted, raw input uses `"multinomial"` and
#'   effective-depth input uses `"dirichlet_multinomial"`.
#' @param dm_concentration Optional Dirichlet-multinomial concentration. If omitted,
#'   effective-depth Dirichlet-multinomial fits use 50 and other fits use 200.
#' @param alfakR_scale Logical; if `TRUE`, append `alfakR`-scale legacy fitness
#'   columns without changing the native `alfak2` posterior estimates.
#' @param n0,nb Initial and final population sizes used to calibrate the
#'   `alfakR` absolute-growth scale. Required when `alfakR_scale = TRUE`.
#' @param correct_efflux Logical; if `TRUE`, apply the `alfakR`
#'   missegregation-efflux viability correction to the legacy-scale columns.
#' @param legacy_weight Weighting scheme used to set the legacy-scale zero point.
#'   `"pi0"` uses the local TMB posterior initial frequencies, `"directly_informed"`
#'   uses observed graph nodes, and `"uniform"` weights all finite summary rows
#'   equally.
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
                       graph_edge_weight = c("mutation", "unit", "normalized"),
                       anchor_support_tiers = "all",
                       anchor_covariance_inflation = c(
                         TMB_sdreport = 1,
                         untrusted_gradient = 4,
                         untrusted_nonconverged = 9,
                         untrusted_sdreport_missing = 4,
                         untrusted_sdreport_nonfinite = 4,
                         unknown = 4
                       ),
                       anchor_count_reference = NULL,
                       anchor_count_power = 1,
                       input_depth = c("raw", "effective"),
                       effective_depth = NULL,
                       effective_depth_mode = c("min", "cap", "fixed"),
                       observation_model = NULL,
                       dm_concentration = NULL,
                       alfakR_scale = FALSE,
                       n0 = NULL,
                       nb = NULL,
                       correct_efflux = FALSE,
                       legacy_weight = c("pi0", "directly_informed", "uniform"),
                       ...) {
  if (!is.logical(alfakR_scale) || length(alfakR_scale) != 1L || is.na(alfakR_scale)) {
    stop("`alfakR_scale` must be TRUE or FALSE.", call. = FALSE)
  }
  modes <- validate_effective_depth_mode(input_depth, effective_depth_mode)
  input_depth <- modes$input_depth
  effective_depth_mode <- modes$effective_depth_mode
  graph_edge_weight <- match.arg(graph_edge_weight)
  obs_controls <- resolve_fit_observation_controls(input_depth, observation_model, dm_concentration)
  legacy_weight <- match.arg(legacy_weight)

  data <- prepare_counts_for_input_depth(
    counts,
    dt = dt,
    beta = beta,
    input_depth = input_depth,
    effective_depth = effective_depth,
    effective_depth_mode = effective_depth_mode
  )
  local_graph <- build_karyotype_graph(
    data,
    shell_depth = local_shell_depth,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  local <- fit_local_posterior(
    data,
    local_graph,
    observation_model = obs_controls$observation_model,
    dm_concentration = obs_controls$dm_concentration,
    ...
  )
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
    sigma_obs_grid = sigma_obs_grid,
    graph_edge_weight = graph_edge_weight,
    anchor_support_tiers = anchor_support_tiers,
    anchor_covariance_inflation = anchor_covariance_inflation,
    anchor_count_reference = anchor_count_reference,
    anchor_count_power = anchor_count_power
  )
  input_depth_diag <- data$metadata$input_depth
  if (is.null(input_depth_diag)) input_depth_diag <- list(input_depth = "raw")
  fit <- list(
    data = data,
    local = local,
    global = global,
    diagnostics = list(
      local = local$diagnostics,
      global = global$diagnostics,
      backend = c(local = "TMB", global = "RcppEigen"),
      input_depth = input_depth_diag,
      observation_model = obs_controls$observation_model,
      dm_concentration = obs_controls$dm_concentration,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = anchor_support_tiers,
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power
    )
  )
  if (isTRUE(alfakR_scale)) {
    if (is.null(n0) || is.null(nb)) {
      stop("`n0` and `nb` are required when `alfakR_scale = TRUE`.", call. = FALSE)
    }
    fit <- add_alfakR_scale_to_fit(
      fit,
      n0 = n0,
      nb = nb,
      correct_efflux = correct_efflux,
      legacy_weight = legacy_weight
    )
  } else {
    fit$diagnostics$alfakR_scale <- list(enabled = FALSE)
  }
  new_alfak2_fit(fit)
}
