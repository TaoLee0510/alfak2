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
  if (inherits(grad, "try-error")) {
    return(list(value = NA_real_, error = conditionMessage(attr(grad, "condition"))))
  }
  list(value = max(abs(grad)), error = NULL)
}

run_local_tmb_attempt <- function(obj, par, control, n) {
  opt <- stats::nlminb(par, obj$fn, obj$gr, control = control)
  sdrep <- try(suppressWarnings(TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = FALSE)),
               silent = TRUE)
  gradient <- tmb_gradient_norm(obj, opt$par)
  list(
    opt = opt,
    sdrep = if (inherits(sdrep, "try-error")) NULL else sdrep,
    sdreport_error = if (inherits(sdrep, "try-error")) conditionMessage(attr(sdrep, "condition")) else NULL,
    frep = extract_tmb_vector(sdrep, "f", n),
    gradient_norm = gradient$value,
    gradient_error = gradient$error
  )
}

assess_local_covariance <- function(attempt, gradient_tolerance) {
  if (!isTRUE(attempt$opt$convergence == 0)) return("untrusted_nonconverged")
  if (is.finite(gradient_tolerance) &&
      (!is.finite(attempt$gradient_norm) || attempt$gradient_norm > gradient_tolerance)) {
    return("untrusted_gradient")
  }
  if (is.null(attempt$frep)) return("untrusted_sdreport_missing")
  if (any(!is.finite(attempt$frep$sd) | attempt$frep$sd <= 0)) {
    return("untrusted_sdreport_nonfinite")
  }
  "TMB_sdreport"
}

local_attempt_diagnostics <- function(attempt, covariance_status) {
  list(
    convergence = attempt$opt$convergence,
    message = attempt$opt$message,
    objective = attempt$opt$objective,
    gradient_norm = attempt$gradient_norm,
    gradient_error = attempt$gradient_error,
    sdreport_error = attempt$sdreport_error,
    covariance_status = covariance_status
  )
}

support_mass_by_tier <- function(values, support_tier) {
  values <- as.numeric(values)
  support_tier <- as.character(support_tier)
  values[!is.finite(values)] <- 0
  out <- rowsum(values, group = support_tier, reorder = FALSE)[, 1]
  total <- sum(out)
  if (is.finite(total) && total > 0) out <- out / total
  out
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
#' @param dm_concentration Dirichlet-multinomial concentration value or grid when
#'   requested. Multiple values are fit explicitly and the finite objective minimum
#'   is selected.
#' @param observation_weight_mode How observation weights enter the
#'   Dirichlet-multinomial likelihood. `"likelihood"` uses a weighted likelihood
#'   score that preserves the original count scale; `"fractional_count"` uses the
#'   historical fractional-count pseudo-likelihood.
#' @param control List passed to `stats::nlminb`.
#' @param gradient_tolerance Maximum absolute gradient tolerated before TMB
#'   covariance is marked untrusted. Set to `Inf` to rely on optimizer
#'   convergence and finite `sdreport` uncertainty.
#' @param retry_on_untrusted_covariance Whether to retry with `retry_control`
#'   when the first optimizer pass does not yield trusted covariance.
#' @param retry_control Optional larger `stats::nlminb` control list used for the
#'   retry pass.
#' @param eta_prior_sd Prior standard deviation for directly observed initial
#'   log-abundances.
#' @param eta_borrowed_prior_mean Prior mean for unobserved graph-node initial
#'   log-abundances.
#' @param eta_borrowed_prior_sd Prior standard deviation for unobserved graph-node
#'   initial log-abundances.
#' @param eta_distance_penalty Additional negative prior-mean shift per support
#'   shell beyond the first borrowed shell.
#'
#' @return An `alfak2_local_fit` object.
#' @export
fit_local_posterior <- function(data,
                                graph = NULL,
                                observation_model = c("multinomial", "dirichlet_multinomial"),
                                dm_concentration = 200,
                                observation_weight_mode = c("likelihood", "fractional_count"),
                                control = list(eval.max = 500, iter.max = 500),
                                gradient_tolerance = Inf,
                                retry_on_untrusted_covariance = TRUE,
                                retry_control = list(eval.max = 2000, iter.max = 2000),
                                eta_prior_sd = 5,
                                eta_borrowed_prior_mean = -6,
                                eta_borrowed_prior_sd = 1.5,
                                eta_distance_penalty = 0.75) {
  if (!inherits(data, "alfak2_data")) {
    stop("`data` must be an alfak2_data object.", call. = FALSE)
  }
  observation_model <- match_observation_model(observation_model)
  observation_weight_mode <- match_observation_weight_mode(observation_weight_mode)
  dm_concentration_grid <- validate_positive_grid(dm_concentration, "dm_concentration")
  if (!is.numeric(gradient_tolerance) || length(gradient_tolerance) != 1L ||
      is.na(gradient_tolerance) || gradient_tolerance <= 0) {
    stop("`gradient_tolerance` must be one positive number or Inf.", call. = FALSE)
  }
  if (!is.logical(retry_on_untrusted_covariance) ||
      length(retry_on_untrusted_covariance) != 1L ||
      is.na(retry_on_untrusted_covariance)) {
    stop("`retry_on_untrusted_covariance` must be TRUE or FALSE.", call. = FALSE)
  }
  validate_scalar(eta_prior_sd, "eta_prior_sd", lower = .Machine$double.eps)
  validate_scalar(eta_borrowed_prior_mean, "eta_borrowed_prior_mean")
  validate_scalar(eta_borrowed_prior_sd, "eta_borrowed_prior_sd", lower = .Machine$double.eps)
  validate_scalar(eta_distance_penalty, "eta_distance_penalty", lower = 0)
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
    paste0("weighted_dirichlet_multinomial_", observation_weight_mode)
  } else if (isTRUE(use_observation_weights)) {
    "weighted_multinomial"
  } else {
    observation_model
  }

  y0_init <- y0 * obs_weight0
  y1_init <- y1 * obs_weight1
  effective_count_total <- y0_init + y1_init
  borrowed_eta <- effective_count_total <= 0 & graph$support_distance > 0L
  eta_prior_mean <- rep(0, n)
  eta_prior_sd_vec <- rep(as.numeric(eta_prior_sd), n)
  eta_prior_mean[borrowed_eta] <- as.numeric(eta_borrowed_prior_mean) -
    as.numeric(eta_distance_penalty) * pmax(0, as.integer(graph$support_distance[borrowed_eta]) - 1L)
  eta_prior_sd_vec[borrowed_eta] <- as.numeric(eta_borrowed_prior_sd)
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
    eta_prior_mean = as.numeric(eta_prior_mean),
    eta_prior_sd = as.numeric(eta_prior_sd_vec),
    anchor_prior_scale = 1.0,
    mu_prior_scale = 0.5,
    scale_prior_scale = 1.0,
    observation_model = as.integer(observation_model == "dirichlet_multinomial"),
    observation_weight_mode = as.integer(observation_weight_mode == "likelihood"),
    dm_concentration = as.numeric(dm_concentration_grid[1])
  )

  concentration_results <- lapply(dm_concentration_grid, function(phi) {
    candidate_data <- tmb_data
    candidate_data$dm_concentration <- as.numeric(phi)
    obj <- tryCatch(
      TMB::MakeADFun(candidate_data, parameters, DLL = "alfak2", silent = TRUE),
      error = function(e) {
        alfak2_abort(
          "Local TMB objective construction failed.",
          diagnostics = list(
            stage = "local_make_adfun",
            dm_concentration = as.numeric(phi),
            error = conditionMessage(e)
          )
        )
      }
    )
    attempt <- tryCatch(
      run_local_tmb_attempt(obj, obj$par, control, n),
      error = function(e) {
        alfak2_abort(
          "Local TMB optimizer failed.",
          diagnostics = list(
            stage = "local_optimizer",
            dm_concentration = as.numeric(phi),
            error = conditionMessage(e)
          )
        )
      }
    )
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
      retry_attempt <- tryCatch(
        run_local_tmb_attempt(obj, attempt$opt$par, merged_retry_control, n),
        error = function(e) {
          alfak2_abort(
            "Local TMB retry optimizer failed.",
            diagnostics = list(
              stage = "local_optimizer_retry",
              dm_concentration = as.numeric(phi),
              initial = attempts$initial,
              error = conditionMessage(e)
            )
          )
        }
      )
      retry_status <- assess_local_covariance(retry_attempt, gradient_tolerance)
      attempts$retry <- local_attempt_diagnostics(retry_attempt, retry_status)
      attempt <- retry_attempt
      covariance_status <- retry_status
    }
    if (!identical(covariance_status, "TMB_sdreport")) {
      alfak2_abort(
        "Local TMB covariance is untrusted; aborting without substituting prior-scale uncertainty.",
        diagnostics = list(
          stage = "local_covariance",
          dm_concentration = as.numeric(phi),
          covariance_status = covariance_status,
          retry_attempted = retry_attempted,
          retry_reason = retry_reason,
          attempts = attempts
        )
      )
    }
    report <- try(obj$report(attempt$opt$par), silent = TRUE)
    if (inherits(report, "try-error")) {
      alfak2_abort(
        "Local TMB report extraction failed.",
        diagnostics = list(
          stage = "local_report",
          dm_concentration = as.numeric(phi),
          error = conditionMessage(attr(report, "condition")),
          attempts = attempts
        )
      )
    }
    if (is.null(report$pi0) || is.null(report$pi1)) {
      alfak2_abort(
        "Local TMB report is missing posterior predictive probabilities.",
        diagnostics = list(
          stage = "local_report",
          dm_concentration = as.numeric(phi),
          report_names = names(report),
          attempts = attempts
        )
      )
    }
    f_sd <- as.numeric(attempt$frep$sd)
    bad_sd <- !is.finite(f_sd) | f_sd <= 0
    if (any(bad_sd)) {
      alfak2_abort(
        "Local TMB covariance returned non-finite fitness uncertainty.",
        diagnostics = list(
          stage = "local_fitness_sd",
          dm_concentration = as.numeric(phi),
          bad_indices = which(bad_sd),
          attempts = attempts
        )
      )
    }
    list(
      dm_concentration = as.numeric(phi),
      obj = obj,
      attempt = attempt,
      covariance_status = covariance_status,
      attempts = attempts,
      retry_attempted = retry_attempted,
      retry_reason = retry_reason,
      report = report,
      objective = as.numeric(attempt$opt$objective)
    )
  })
  objectives <- vapply(concentration_results, `[[`, numeric(1), "objective")
  if (any(!is.finite(objectives))) {
    alfak2_abort(
      "Local TMB concentration grid produced non-finite objectives.",
      diagnostics = list(
        stage = "local_dm_concentration_grid",
        dm_concentration_grid = as.numeric(dm_concentration_grid),
        objectives = objectives
      )
    )
  }
  selected_idx <- which.min(objectives)
  selected <- concentration_results[[selected_idx]]
  obj <- selected$obj
  attempt <- selected$attempt
  covariance_status <- selected$covariance_status
  attempts <- selected$attempts
  retry_attempted <- selected$retry_attempted
  retry_reason <- selected$retry_reason
  dm_concentration_selected <- selected$dm_concentration
  opt <- attempt$opt
  plist <- obj$env$parList(opt$par)
  f_mean <- as.numeric(plist$f)
  report <- selected$report
  pi0_report <- as.numeric(report$pi0)
  pi1_report <- as.numeric(report$pi1)
  f_sd <- as.numeric(attempt$frep$sd)

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
    effective_count_total = as.numeric(effective_count_total),
    eta_prior_mean = as.numeric(eta_prior_mean),
    eta_prior_sd = as.numeric(eta_prior_sd_vec),
    pi0 = as.numeric(pi0_report),
    pi1 = as.numeric(pi1_report),
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
    fitness_sd_source = "TMB_sdreport",
    retry_attempted = retry_attempted,
    retry_reason = retry_reason,
    attempts = attempts,
    observation_model = observation_model,
    observation_weight_mode = observation_weight_mode,
    likelihood_model = likelihood_model,
    dm_concentration = dm_concentration_selected,
    dm_concentration_grid = as.numeric(dm_concentration_grid),
    dm_concentration_objective = stats::setNames(objectives, as.character(dm_concentration_grid)),
    use_observation_weights = use_observation_weights,
    observation_weight_summary = summary(c(obs_weight0[observed_index], obs_weight1[observed_index])),
    eta_prior = list(
      direct_sd = as.numeric(eta_prior_sd),
      borrowed_mean = as.numeric(eta_borrowed_prior_mean),
      borrowed_sd = as.numeric(eta_borrowed_prior_sd),
      distance_penalty = as.numeric(eta_distance_penalty),
      n_borrowed_shrunk = sum(borrowed_eta)
    ),
    pi0_support_mass = support_mass_by_tier(pi0_report, graph$support_tier),
    pi1_support_mass = support_mass_by_tier(pi1_report, graph$support_tier)
  )
  new_alfak2_local_fit(list(
    data = data,
    graph = graph,
    summary = summary,
    parameter_mode = plist,
    sdreport = attempt$sdrep,
    posterior_predictive = report[c("pi0", "pi1")],
    diagnostics = diagnostics
  ))
}
