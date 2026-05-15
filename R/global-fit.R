graph_edge_weights <- function(edge_weight,
                               mode = c("mutation", "unit", "normalized")) {
  mode <- match.arg(mode)
  edge_weight <- as.numeric(edge_weight)
  if (identical(mode, "mutation")) return(edge_weight)
  positive <- is.finite(edge_weight) & edge_weight > 0
  if (identical(mode, "unit")) {
    out <- edge_weight
    out[positive] <- 1
    return(out)
  }
  scale <- stats::median(edge_weight[positive], na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) return(edge_weight)
  edge_weight / scale
}

anchor_covariance_multiplier <- function(status,
                                         inflation = c(
                                           TMB_sdreport = 1,
                                           untrusted_gradient = 4,
                                           untrusted_nonconverged = 9,
                                           untrusted_sdreport_missing = 4,
                                           untrusted_sdreport_nonfinite = 4,
                                           unknown = 4
                                         )) {
  if (is.null(inflation) || !length(inflation)) return(rep(1, length(status)))
  inflation_names <- names(inflation)
  inflation <- as.numeric(inflation)
  names(inflation) <- inflation_names
  if (is.null(names(inflation)) || any(!nzchar(names(inflation)))) {
    stop("`anchor_covariance_inflation` must be a named numeric vector.", call. = FALSE)
  }
  status <- as.character(status)
  status[is.na(status) | !nzchar(status)] <- "unknown"
  out <- inflation[match(status, names(inflation))]
  fallback <- inflation[["unknown"]]
  if (!is.finite(fallback) || fallback <= 0) fallback <- 4
  out[!is.finite(out) | out <= 0] <- fallback
  as.numeric(out)
}

count_anchor_multiplier <- function(count_total,
                                    count_reference = NULL,
                                    power = 1) {
  if (is.null(count_reference)) return(rep(1, length(count_total)))
  validate_scalar(as.numeric(count_reference), "anchor_count_reference", lower = .Machine$double.eps)
  validate_scalar(as.numeric(power), "anchor_count_power", lower = 0)
  count_total <- as.numeric(count_total)
  out <- rep(1, length(count_total))
  ok <- is.finite(count_total) & count_total > 0
  out[ok] <- pmax(1, (as.numeric(count_reference) / count_total[ok]) ^ as.numeric(power))
  out[!ok] <- as.numeric(count_reference) ^ as.numeric(power)
  out
}

resolve_anchor_tiers <- function(anchor_support_tiers) {
  if (is.null(anchor_support_tiers) || identical(anchor_support_tiers, "all")) {
    return(NULL)
  }
  anchor_support_tiers <- as.character(anchor_support_tiers)
  anchor_support_tiers[nzchar(anchor_support_tiers)]
}

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
#' @param graph_edge_weight Edge weighting used by the Gaussian graph prior.
#'   `"mutation"` preserves the transition weights implied by `beta`, `"unit"`
#'   gives every graph edge equal smoothing weight, and `"normalized"` preserves
#'   relative mutation weights while removing the global `beta` scale.
#' @param anchor_support_tiers Local support tiers to use as Gaussian anchors.
#'   Defaults to `"all"` for the legacy behavior of anchoring every finite local
#'   posterior row.
#' @param anchor_covariance_inflation Named variance multipliers for local
#'   covariance statuses.
#' @param anchor_count_reference Optional count total below which observed
#'   anchors are downweighted by increasing their variance.
#' @param anchor_count_power Exponent used by count-based variance inflation.
#' @param anchor_exclude Optional karyotype labels to exclude from the anchor
#'   set, used for neighbor-holdout validation.
#'
#' @return A list containing node summaries and tuning diagnostics.
#' @export
fit_graph_posterior <- function(local_fit,
                                graph = NULL,
                                lambda_l_grid = c(0.2, 1, 5),
                                lambda_e_grid = c(0.05, 0.25, 1),
                                sigma_obs_grid = c(0.02, 0.05, 0.1),
                                eps = 1e-5,
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
                                anchor_exclude = character()) {
  if (!inherits(local_fit, "alfak2_local_fit")) {
    stop("`local_fit` must be an alfak2_local_fit object.", call. = FALSE)
  }
  if (is.null(graph)) graph <- local_fit$graph
  if (!inherits(graph, "alfak2_graph")) stop("`graph` must be an alfak2_graph object.", call. = FALSE)
  graph_edge_weight <- match.arg(graph_edge_weight)
  anchor_tiers <- resolve_anchor_tiers(anchor_support_tiers)

  anchor_match <- match(local_fit$summary$karyotype, as.character(graph$labels))
  tier_ok <- rep(TRUE, nrow(local_fit$summary))
  if (!is.null(anchor_tiers)) {
    tier_ok <- as.character(local_fit$summary$support_tier) %in% anchor_tiers
  }
  if (!is.null(anchor_exclude) && length(anchor_exclude)) {
    tier_ok <- tier_ok & !(as.character(local_fit$summary$karyotype) %in% as.character(anchor_exclude))
  }
  keep <- which(!is.na(anchor_match) & is.finite(local_fit$summary$fitness_mean) & tier_ok)
  if (!length(keep)) stop("No local posterior anchors are present in the global graph.", call. = FALSE)
  anchor_var_base <- local_fit$summary$fitness_sd[keep]^2
  anchor_var <- anchor_var_base
  anchor_var[!is.finite(anchor_var) | anchor_var <= 0] <- stats::median(anchor_var[is.finite(anchor_var) & anchor_var > 0], na.rm = TRUE)
  anchor_var[!is.finite(anchor_var) | anchor_var <= 0] <- 0.25
  covariance_status <- if ("covariance_status" %in% names(local_fit$summary)) {
    as.character(local_fit$summary$covariance_status[keep])
  } else {
    status <- local_fit$diagnostics$covariance_status
    if (is.null(status) || !length(status) || is.na(status[[1L]]) || !nzchar(as.character(status[[1L]]))) {
      status <- "unknown"
    }
    rep(as.character(status[[1L]]), length(keep))
  }
  covariance_mult <- anchor_covariance_multiplier(covariance_status, anchor_covariance_inflation)
  count_mult <- if ("count_total" %in% names(local_fit$summary)) {
    count_anchor_multiplier(local_fit$summary$count_total[keep], anchor_count_reference, anchor_count_power)
  } else {
    rep(1, length(keep))
  }
  anchor_var_multiplier <- covariance_mult * count_mult
  anchor_var <- anchor_var * anchor_var_multiplier
  edge_weight <- graph_edge_weights(graph$edge_weight, graph_edge_weight)

  res <- alfak2_graph_posterior_cpp(
    karyotypes = graph$karyotypes,
    edge_from = as.integer(graph$edge_from),
    edge_to = as.integer(graph$edge_to),
    edge_weight = as.numeric(edge_weight),
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
      variance_base = anchor_var_base,
      variance = anchor_var,
      variance_multiplier = anchor_var_multiplier,
      covariance_status = covariance_status,
      count_total = if ("count_total" %in% names(local_fit$summary)) local_fit$summary$count_total[keep] else NA_integer_,
      stringsAsFactors = FALSE
    ),
    hyperparameters = list(
      lambda_l = res$lambda_l,
      lambda_e = res$lambda_e,
      sigma_obs = res$sigma_obs,
      cv_score = res$cv_score,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = if (is.null(anchor_tiers)) "all" else paste(anchor_tiers, collapse = ","),
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = anchor_count_power,
      anchor_excluded = length(anchor_exclude)
    ),
    tuning_grid = res$grid,
    diagnostics = list(
      factorization_status = res$factorization_status,
      graph_edge_weight = graph_edge_weight,
      anchor_support_tiers = if (is.null(anchor_tiers)) "all" else anchor_tiers,
      anchor_excluded = as.character(anchor_exclude),
      anchor_variance_multiplier_summary = summary(anchor_var_multiplier)
    )
  )
}
