#' Fit the complete alfak2 pipeline
#'
#' Runs data validation, local graph construction, TMB local posterior fitting,
#' bounded global graph expansion, and compiled graph Gaussian posterior solving.
#'
#' @param counts Two-column count matrix or an `alfak2_data` object.
#' @param dt Time interval if `counts` is not already prepared.
#' @param beta Missegregation rate.
#' @param transition_kernel Transition weighting used by graph edges and local
#'   transport. `"exact"` uses chromosome-level missegregation probabilities;
#'   `"linear"` preserves the historical first-order approximation.
#' @param local_shell_depth Local graph shell depth.
#' @param global_extra_shell Additional bounded shell depth for graph borrowing.
#' @param min_cn,max_cn Copy-number bounds.
#' @param max_nodes Hard graph-size bound.
#' @param lambda_l_grid,lambda_e_grid,sigma_obs_grid Hyperparameter grids passed
#'   to the compiled graph posterior.
#' @param graph_edge_weight Edge weighting used by the global graph prior.
#' @param extrapolation_method Second-layer extrapolation method. The default
#'   `"graph_gaussian_baseline"` preserves the historical global graph posterior
#'   path. The legacy alias `"graph_gaussian"` is still accepted.
#' @param max_prediction_distance Maximum support distance reported as reliable
#'   for near-field extrapolation methods.
#' @param extrapolation_control Optional named list of method-specific controls
#'   passed to `fit_extrapolation_layer()`.
#' @param anchor_support_tiers Local support tiers used as global graph anchors;
#'   defaults to `"all"` to preserve the legacy public API behavior.
#' @param anchor_exclude Optional karyotype labels to exclude from the global
#'   anchor set, primarily for neighbor-holdout validation.
#' @param anchor_covariance_inflation Named variance multipliers for local
#'   covariance statuses.
#' @param anchor_count_reference Optional low-count downweighting reference for
#'   global graph anchors.
#' @param anchor_count_power Exponent for count-based anchor variance inflation.
#' @param anchor_min_effective_count Minimum effective count required for local
#'   rows to be used as global anchors. Set to `NULL` to preserve legacy
#'   all-finite-row anchoring.
#' @param input_depth `"raw"` uses the supplied count matrix directly. `"effective"`
#'   converts each timepoint to frequencies and re-counts them at a controlled
#'   effective depth before fitting.
#' @param effective_depth Effective sequencing depth used when
#'   `input_depth = "effective"`. Required for `"cap"` and `"fixed"` modes; ignored
#'   by `"min"` mode.
#' @param effective_depth_mode Effective-depth rule. `"min"` uses the smaller raw
#'   timepoint depth for both columns, `"cap"` caps each raw depth at
#'   `effective_depth`, and `"fixed"` uses `effective_depth` for every timepoint.
#' @param effective_depth_rounding Rounding rule for effective-depth expected
#'   counts. `"hash"` uses deterministic label-based tie breaking,
#'   `"largest_remainder"` preserves row-order tie breaking, and `"stochastic"`
#'   samples remainders with `effective_depth_seed`.
#' @param effective_depth_seed Optional seed used only by stochastic
#'   effective-depth rounding.
#' @param observation_model Optional observation model passed to
#'   `fit_local_posterior()`. If omitted, raw input uses `"multinomial"` and
#'   effective-depth input uses `"dirichlet_multinomial"`.
#' @param dm_concentration Optional Dirichlet-multinomial concentration value or
#'   grid. If omitted, effective-depth Dirichlet-multinomial fits use 50 and other
#'   fits use 200.
#' @param observation_weight_mode How observation weights enter the
#'   Dirichlet-multinomial likelihood. See `fit_local_posterior()`.
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
#' @param pipeline `"two_layer"` retains separate local and extrapolation fits;
#'   `"joint_count"` fits a single count likelihood and structured landscape
#'   prior on the entire bounded target graph. The joint model is an explicitly
#'   selectable candidate pending independent all-distance validation.
#'   `"sparse_joint"` fits the sparse joint count model, which is built for two
#'   sparse samples: see [fit_sparse_joint()], whose result it returns. Of the
#'   other arguments only `counts`, `dt`, `beta` and `sparse_joint_control`
#'   then apply.
#' @param sparse_joint_control Named list of further arguments for
#'   [fit_sparse_joint()] when `pipeline = "sparse_joint"`, such as
#'   `birth_rate`, `min_cn` and `max_cn`.
#' @param joint_control Named prior controls passed to the joint landscape.
#'   Experimental `interaction_support = "observed_nullspace"` separates pair
#'   directions that vanish at observed design rows; this is not a likelihood
#'   identifiability test. `unidentified_interaction_sd` controls their prior.
#'   Optional `prior_components` is a list of prior-control lists, with
#'   `model_prior` supplying component weights before conditional evidence
#'   updating. It cannot be combined with `model_scales`.
#' @param target_karyotypes Optional explicit finite query states, as labels or a
#'   copy-number matrix. They define the global domain together with observed
#'   nodes, permitting distant queries without complete shell enumeration.
#' @param forward_tolerance Joint division-model step-halving tolerance on both
#'   maximum frequency change and log total growth change.
#' @param forward_max_refinements Maximum joint forward step halvings and refits.
#'
#' @return An `alfak2_fit` object, or an `alfak2_sparse_joint_fit` object for
#'   `pipeline = "sparse_joint"`.
#' @export
fit_alfak2 <- function(counts,
                       dt = 1,
                       beta = 0.00005,
                       transition_kernel = c("exact", "linear"),
                       local_shell_depth = 2,
                       global_extra_shell = 1,
                       min_cn = 0,
                       max_cn = 5,
                       max_nodes = 5000,
                       lambda_l_grid = c(0.2, 1, 5),
                       lambda_e_grid = c(0.05, 0.25, 1),
                       sigma_obs_grid = c(0.02, 0.05, 0.1),
                       graph_edge_weight = c("mutation", "unit", "normalized"),
                       extrapolation_method = c(
                         "graph_gaussian_baseline",
                         "edge_effect_empirical_bayes",
                         "edge_effect_interaction_path_ensemble",
                         "kronecker_or_graph_trend_filtering",
                         "local_NNGP_or_GPnn",
                         "delta_tree_ensemble",
                         "tabpfn_nearfield_feature_model",
                         "truncated_nearfield_gmrf",
                         "local_polynomial_stencil"
                       ),
                       max_prediction_distance = 2,
                       extrapolation_control = list(),
                       anchor_support_tiers = "all",
                       anchor_exclude = character(),
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
                       anchor_min_effective_count = 0,
                       input_depth = c("raw", "effective"),
                       effective_depth = NULL,
                       effective_depth_mode = c("min", "cap", "fixed"),
                       effective_depth_rounding = c("hash", "largest_remainder", "stochastic"),
                       effective_depth_seed = NULL,
                       observation_model = NULL,
                       dm_concentration = NULL,
                       observation_weight_mode = c("likelihood", "fractional_count"),
                       alfakR_scale = FALSE,
                       n0 = NULL,
                       nb = NULL,
                       correct_efflux = FALSE,
                       legacy_weight = c("pi0", "directly_informed", "uniform"),
                       pipeline = c("two_layer", "joint_count", "sparse_joint"),
                       joint_control = list(),
                       sparse_joint_control = list(),
                       target_karyotypes = NULL,
                       forward_tolerance = 1e-5,
                       forward_max_refinements = 3L,
                       ...) {
  pipeline <- match.arg(pipeline)
  if (identical(pipeline, "sparse_joint")) {
    if (!is.list(sparse_joint_control) || (length(sparse_joint_control) &&
        (is.null(names(sparse_joint_control)) || any(!nzchar(names(sparse_joint_control)))))) {
      stop("`sparse_joint_control` must be a named list.", call. = FALSE)
    }
    # A model of its own, with its own domain, result and prediction method.
    return(do.call(fit_sparse_joint, c(list(counts, dt = dt, beta = beta), sparse_joint_control)))
  }
  if (!is.list(joint_control) || (length(joint_control) &&
      (is.null(names(joint_control)) || any(!nzchar(names(joint_control))))))
    stop("`joint_control` must be a named list.", call. = FALSE)
  if (!is.logical(alfakR_scale) || length(alfakR_scale) != 1L || is.na(alfakR_scale)) {
    stop("`alfakR_scale` must be TRUE or FALSE.", call. = FALSE)
  }
  modes <- validate_effective_depth_mode(input_depth, effective_depth_mode)
  input_depth <- modes$input_depth
  effective_depth_mode <- modes$effective_depth_mode
  graph_edge_weight <- match.arg(graph_edge_weight)
  extrapolation_method <- match_extrapolation_method(extrapolation_method)
  max_prediction_distance <- validate_prediction_distance(max_prediction_distance)
  if (!is.list(extrapolation_control) ||
      (length(extrapolation_control) && (is.null(names(extrapolation_control)) || any(!nzchar(names(extrapolation_control)))))) {
    stop("`extrapolation_control` must be a named list.", call. = FALSE)
  }
  transition_kernel <- match.arg(transition_kernel)
  effective_depth_rounding <- match.arg(effective_depth_rounding)
  observation_weight_mode <- match_observation_weight_mode(observation_weight_mode)
  obs_controls <- resolve_fit_observation_controls(input_depth, observation_model, dm_concentration)
  legacy_weight <- match.arg(legacy_weight)

  data <- prepare_counts_for_input_depth(
    counts,
    dt = dt,
    beta = beta,
    input_depth = input_depth,
    effective_depth = effective_depth,
    effective_depth_mode = effective_depth_mode,
    effective_depth_rounding = effective_depth_rounding,
    effective_depth_seed = effective_depth_seed
  )
  local_graph <- if (identical(pipeline, "two_layer")) build_karyotype_graph(
    data,
    transition_kernel = transition_kernel,
    shell_depth = local_shell_depth,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  ) else NULL
  global_graph <- if (!is.null(target_karyotypes)) {
    build_explicit_target_graph(data, target_karyotypes, min_cn, max_cn,
      max_nodes, transition_kernel)
  } else build_karyotype_graph(
    data,
    transition_kernel = transition_kernel,
    shell_depth = local_shell_depth + global_extra_shell,
    min_cn = min_cn,
    max_cn = max_cn,
    max_nodes = max_nodes
  )
  if (identical(pipeline, "joint_count")) {
    joint <- fit_joint_landscape(data, global_graph, joint_control = joint_control,
      forward_tolerance = forward_tolerance, forward_max_refinements = forward_max_refinements,
      observation_model = obs_controls$observation_model,
      dm_concentration = obs_controls$dm_concentration,
      observation_weight_mode = observation_weight_mode, ...)
    local <- joint$local
    global <- joint$global
  } else {
    local <- fit_local_posterior(data, local_graph,
      observation_model = obs_controls$observation_model,
      dm_concentration = obs_controls$dm_concentration,
      observation_weight_mode = observation_weight_mode, ...)
    global_args <- c(
    list(
      local_fit = local,
      graph = global_graph,
      method = extrapolation_method,
      max_prediction_distance = max_prediction_distance,
      lambda_l_grid = lambda_l_grid,
      lambda_e_grid = lambda_e_grid,
      sigma_obs_grid = sigma_obs_grid,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = anchor_support_tiers,
      anchor_exclude = anchor_exclude,
      anchor_covariance_inflation = anchor_covariance_inflation,
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power,
      anchor_min_effective_count = anchor_min_effective_count
    ),
    extrapolation_control
  )
  global <- do.call(fit_extrapolation_layer, global_args)
  }
  input_depth_diag <- data$metadata$input_depth
  if (is.null(input_depth_diag)) input_depth_diag <- list(input_depth = "raw")
  fit <- list(
    data = data,
    local = local,
    global = global,
    diagnostics = list(
      local = local$diagnostics,
      global = global$diagnostics,
      backend = c(local = "TMB", global = if (identical(pipeline, "joint_count")) "TMB" else "RcppEigen"),
      pipeline = pipeline,
      input_depth = input_depth_diag,
      observation_model = obs_controls$observation_model,
      observation_weight_mode = observation_weight_mode,
      dm_concentration = local$diagnostics$dm_concentration,
      dm_concentration_grid = obs_controls$dm_concentration,
      transition_kernel = transition_kernel,
      graph_edge_weight = graph_edge_weight,
      extrapolation_method = if (identical(pipeline, "joint_count")) "joint_count_landscape" else extrapolation_method,
      max_prediction_distance = max_prediction_distance,
      anchor_support_tiers = anchor_support_tiers,
      anchor_exclude = as.character(anchor_exclude),
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power,
      anchor_min_effective_count = anchor_min_effective_count
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
