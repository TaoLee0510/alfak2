extract_tmb_vector <- function(sdrep, name, n) {
  if (inherits(sdrep, "try-error") || is.null(sdrep)) return(NULL)
  tab <- try(suppressWarnings(as.data.frame(summary(sdrep, "report"))), silent = TRUE)
  if (inherits(tab, "try-error") || !nrow(tab)) return(NULL)
  rn <- rownames(tab)
  idx <- which(rn == name | grepl(paste0("^", name, "(\\.|\\[)"), rn))
  if (length(idx) != n) return(NULL)
  list(mean = tab$Estimate[idx], sd = tab$`Std. Error`[idx])
}

tmb_gradient_norm <- function(obj, par) {
  grad <- try(obj$gr(par), silent = TRUE)
  if (inherits(grad, "try-error")) return(NA_real_)
  max(abs(grad))
}

run_local_tmb_attempt <- function(obj, par, control, n) {
  opt <- stats::nlminb(par, obj$fn, obj$gr, control = control)
  sdrep <- try(suppressWarnings(TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = FALSE)),
               silent = TRUE)
  list(
    opt = opt,
    sdrep = if (inherits(sdrep, "try-error")) NULL else sdrep,
    frep = extract_tmb_vector(sdrep, "f", n),
    gradient_norm = tmb_gradient_norm(obj, opt$par)
  )
}

assess_local_covariance <- function(attempt, gradient_tolerance) {
  if (!isTRUE(attempt$opt$convergence == 0)) return("untrusted_nonconverged")
  if (!is.finite(attempt$gradient_norm) || attempt$gradient_norm > gradient_tolerance) {
    return("untrusted_gradient")
  }
  if (is.null(attempt$frep)) return("untrusted_sdreport_missing")
  if (any(!is.finite(attempt$frep$sd) | attempt$frep$sd <= 0)) {
    return("untrusted_sdreport_nonfinite")
  }
  "TMB_sdreport"
}

fallback_local_fitness_sd <- function(report, support_distance) {
  sigma_anchor <- if (!inherits(report, "try-error") && !is.null(report$sigma_anchor)) {
    as.numeric(report$sigma_anchor)[1]
  } else {
    NA_real_
  }
  sigma_neighbor <- if (!inherits(report, "try-error") && !is.null(report$sigma_neighbor)) {
    as.numeric(report$sigma_neighbor)[1]
  } else {
    NA_real_
  }
  if (!is.finite(sigma_anchor) || sigma_anchor <= 0) sigma_anchor <- 1
  if (!is.finite(sigma_neighbor) || sigma_neighbor <= 0) sigma_neighbor <- sigma_anchor

  out <- ifelse(support_distance == 0L, sigma_anchor, sigma_neighbor)
  out[support_distance >= 2L] <- out[support_distance >= 2L] * 1.75
  pmax(as.numeric(out), .Machine$double.eps)
}

local_attempt_diagnostics <- function(attempt, covariance_status) {
  list(
    convergence = attempt$opt$convergence,
    message = attempt$opt$message,
    objective = attempt$opt$objective,
    gradient_norm = attempt$gradient_norm,
    covariance_status = covariance_status
  )
}

#' Fit the local hierarchical posterior
#'
#' The likelihood and hierarchical prior objective are implemented in TMB. This
#' R function validates inputs, builds initial values, runs `nlminb`, and
#' collects TMB-reported posterior marginals.
#'
#' @param data An `alfak2_data` object.
#' @param graph Optional `alfak2_graph`; defaults to a local depth-2 graph.
#' @param observation_model `"multinomial"` or `"dirichlet_multinomial"`.
#' @param dm_concentration Dirichlet-multinomial concentration when requested.
#' @param control List passed to `stats::nlminb`.
#' @param gradient_tolerance Maximum absolute gradient tolerated before TMB
#'   covariance is marked untrusted.
#' @param retry_on_untrusted_covariance Whether to retry with `retry_control`
#'   when the first optimizer pass does not yield trusted covariance.
#' @param retry_control Optional larger `stats::nlminb` control list used for the
#'   retry pass.
#'
#' @return An `alfak2_local_fit` object.
#' @export
fit_local_posterior <- function(data,
                                graph = NULL,
                                observation_model = c("multinomial", "dirichlet_multinomial"),
                                dm_concentration = 200,
                                control = list(eval.max = 500, iter.max = 500),
                                gradient_tolerance = 1e-3,
                                retry_on_untrusted_covariance = TRUE,
                                retry_control = list(eval.max = 2000, iter.max = 2000)) {
  if (!inherits(data, "alfak2_data")) {
    stop("`data` must be an alfak2_data object.", call. = FALSE)
  }
  observation_model <- match_observation_model(observation_model)
  if (!is.numeric(gradient_tolerance) || length(gradient_tolerance) != 1L ||
      !is.finite(gradient_tolerance) || gradient_tolerance <= 0) {
    stop("`gradient_tolerance` must be one positive finite number.", call. = FALSE)
  }
  if (!is.logical(retry_on_untrusted_covariance) ||
      length(retry_on_untrusted_covariance) != 1L ||
      is.na(retry_on_untrusted_covariance)) {
    stop("`retry_on_untrusted_covariance` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(graph)) graph <- build_karyotype_graph(data, shell_depth = 2)
  if (!inherits(graph, "alfak2_graph")) {
    stop("`graph` must be an alfak2_graph object.", call. = FALSE)
  }

  n <- length(graph$labels)
  observed_index <- match(data$labels, as.character(graph$labels))
  if (anyNA(observed_index)) {
    missing <- data$labels[is.na(observed_index)]
    stop("Graph is missing observed karyotypes: ", paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  }
  y0 <- numeric(n)
  y1 <- numeric(n)
  y0[observed_index] <- data$counts[, 1]
  y1[observed_index] <- data$counts[, 2]
  observation_weights <- data$metadata$observation_weights
  if (is.null(observation_weights)) {
    observation_weights <- matrix(1, nrow = nrow(data$counts), ncol = 2L,
                                  dimnames = list(rownames(data$counts), c("t0", "t1")))
  } else {
    observation_weights <- validate_observation_weights(observation_weights, data$counts)
  }
  obs_weight0 <- rep(1, n)
  obs_weight1 <- rep(1, n)
  obs_weight0[observed_index] <- observation_weights[, 1]
  obs_weight1[observed_index] <- observation_weights[, 2]
  use_observation_weights <- any(abs(c(obs_weight0[observed_index], obs_weight1[observed_index]) - 1) > 1e-12)
  likelihood_model <- if (isTRUE(use_observation_weights) && identical(observation_model, "dirichlet_multinomial")) {
    "weighted_dirichlet_multinomial"
  } else if (isTRUE(use_observation_weights)) {
    "weighted_multinomial"
  } else {
    observation_model
  }

  y0_init <- y0 * obs_weight0
  y1_init <- y1 * obs_weight1
  p0 <- (y0_init + 0.5) / sum(y0_init + 0.5)
  p1 <- (y1_init + 0.5) / sum(y1_init + 0.5)
  f0 <- log(p1) - log(p0)
  f0 <- f0 - mean(f0)
  n_context <- length(graph$context_label)
  n_group <- max(unlist(graph$context_group0), 0L) + 1L
  parameters <- list(
    eta = log(p0),
    f = f0 / max(data$dt, .Machine$double.eps),
    delta_context = rep(0, n_context),
    mu_group = rep(0, n_group),
    log_sigma_neighbor = log(0.35),
    log_sigma_anchor = log(0.6),
    log_tau_group = rep(log(0.2), n_group)
  )
  tmb_data <- list(
    y0 = as.numeric(y0),
    y1 = as.numeric(y1),
    n_nodes = as.integer(n),
    trans_from = as.integer(unlist(graph$transition_from0)),
    trans_to = as.integer(unlist(graph$transition_to0)),
    trans_weight = as.numeric(unlist(graph$transition_weight)),
    support_distance = as.integer(graph$support_distance),
    parent_from = as.integer(unlist(graph$parent_from0)),
    parent_to = as.integer(unlist(graph$parent_to0)),
    parent_weight = as.numeric(unlist(graph$parent_weight)),
    parent_context = as.integer(unlist(graph$parent_context0)),
    context_group = as.integer(unlist(graph$context_group0)),
    dt = as.numeric(data$dt),
    obs_weight0 = as.numeric(obs_weight0),
    obs_weight1 = as.numeric(obs_weight1),
    use_observation_weights = as.integer(use_observation_weights),
    anchor_prior_scale = 1.0,
    mu_prior_scale = 0.5,
    scale_prior_scale = 1.0,
    observation_model = as.integer(observation_model == "dirichlet_multinomial"),
    dm_concentration = as.numeric(dm_concentration)
  )

  obj <- TMB::MakeADFun(tmb_data, parameters, DLL = "alfak2", silent = TRUE)
  attempt <- run_local_tmb_attempt(obj, obj$par, control, n)
  covariance_status <- assess_local_covariance(attempt, gradient_tolerance)
  attempts <- list(initial = local_attempt_diagnostics(attempt, covariance_status))
  retry_attempted <- FALSE
  retry_reason <- NULL
  if (!identical(covariance_status, "TMB_sdreport") &&
      isTRUE(retry_on_untrusted_covariance) &&
      !is.null(retry_control)) {
    retry_attempted <- TRUE
    retry_reason <- covariance_status
    merged_retry_control <- utils::modifyList(control, retry_control)
    retry_attempt <- run_local_tmb_attempt(obj, attempt$opt$par, merged_retry_control, n)
    retry_status <- assess_local_covariance(retry_attempt, gradient_tolerance)
    attempts$retry <- local_attempt_diagnostics(retry_attempt, retry_status)
    attempt <- retry_attempt
    covariance_status <- retry_status
  }
  opt <- attempt$opt
  plist <- obj$env$parList(opt$par)
  f_mean <- as.numeric(plist$f)
  report <- try(obj$report(opt$par), silent = TRUE)
  covariance_fallback <- !identical(covariance_status, "TMB_sdreport")
  f_sd <- if (!covariance_fallback && !is.null(attempt$frep)) {
    as.numeric(attempt$frep$sd)
  } else {
    fallback_local_fitness_sd(report, graph$support_distance)
  }
  bad_sd <- !is.finite(f_sd) | f_sd <= 0
  if (any(bad_sd)) {
    replacement <- fallback_local_fitness_sd(report, graph$support_distance)
    f_sd[bad_sd] <- replacement[bad_sd]
    covariance_fallback <- TRUE
    if (identical(covariance_status, "TMB_sdreport")) {
      covariance_status <- "untrusted_sdreport_nonfinite"
    }
  }

  summary <- data.frame(
    node_id = seq_len(n),
    karyotype = as.character(graph$labels),
    support_tier = as.character(graph$support_tier),
    support_distance = as.integer(graph$support_distance),
    count_t0 = as.integer(y0),
    count_t1 = as.integer(y1),
    count_total = as.integer(y0 + y1),
    observation_weight_t0 = as.numeric(obs_weight0),
    observation_weight_t1 = as.numeric(obs_weight1),
    effective_count_total = as.numeric(y0 * obs_weight0 + y1 * obs_weight1),
    is_observed = as.logical((y0 + y1) > 0),
    fitness_mean = f_mean,
    fitness_sd = f_sd,
    conf_low = f_mean - 1.959963984540054 * f_sd,
    conf_high = f_mean + 1.959963984540054 * f_sd,
    covariance_status = covariance_status,
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    convergence = opt$convergence,
    message = opt$message,
    objective = opt$objective,
    gradient_norm = attempt$gradient_norm,
    covariance_status = covariance_status,
    covariance_fallback = covariance_fallback,
    fitness_sd_source = if (isTRUE(covariance_fallback)) "fallback_prior_scale" else "TMB_sdreport",
    retry_attempted = retry_attempted,
    retry_reason = retry_reason,
    attempts = attempts,
    observation_model = observation_model,
    likelihood_model = likelihood_model,
    use_observation_weights = use_observation_weights,
    observation_weight_summary = summary(c(obs_weight0[observed_index], obs_weight1[observed_index]))
  )
  new_alfak2_local_fit(list(
    data = data,
    graph = graph,
    summary = summary,
    parameter_mode = plist,
    sdreport = attempt$sdrep,
    posterior_predictive = if (inherits(report, "try-error")) NULL else report[c("pi0", "pi1")],
    diagnostics = diagnostics
  ))
}
