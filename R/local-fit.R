extract_tmb_vector <- function(sdrep, name, n) {
  if (inherits(sdrep, "try-error") || is.null(sdrep)) return(NULL)
  tab <- try(suppressWarnings(as.data.frame(summary(sdrep, "report"))), silent = TRUE)
  if (inherits(tab, "try-error") || !nrow(tab)) return(NULL)
  rn <- rownames(tab)
  idx <- which(rn == name | grepl(paste0("^", name, "(\\.|\\[)"), rn))
  if (length(idx) != n) return(NULL)
  list(mean = tab$Estimate[idx], sd = tab$`Std. Error`[idx], index = idx)
}

local_tmb_parameter_map <- function(parameters, fixed_sigma_anchor,
                                    fixed_sigma_neighbor, fixed_tau_group) {
  out <- list()
  fixed <- c(log_sigma_anchor = is.finite(fixed_sigma_anchor),
             log_sigma_neighbor = is.finite(fixed_sigma_neighbor),
             log_tau_group = is.finite(fixed_tau_group))
  for (nm in names(fixed)[fixed]) {
    out[[nm]] <- factor(rep(NA_integer_, length(parameters[[nm]])))
  }
  out
}

local_tmb_parameter_list <- function(obj, par) {
  # Avoid repeating an inner Newton solve when the fitted objective is already
  # at this fixed vector (notably the zero-dimensional structured-prior case).
  fixed_index <- setdiff(seq_along(obj$env$last.par), obj$env$random)
  if (!identical(as.numeric(obj$env$last.par[fixed_index]), as.numeric(par))) obj$fn(par)
  # parList fills its random entries from last.par and expects only the fixed
  # vector as its argument, including when the fixed vector has length zero.
  obj$env$parList(par)
}

extract_local_fitness_covariance <- function(sdrep, frep, covariance_status) {
  if (!identical(covariance_status, "TMB_sdreport") || is.null(frep) ||
      is.null(sdrep$cov) || is.null(frep$index)) return(NULL)
  if (max(frep$index) > nrow(sdrep$cov)) return(NULL)
  out <- as.matrix(sdrep$cov[frep$index, frep$index, drop = FALSE])
  if (any(!is.finite(out))) return(NULL)
  (out + t(out)) / 2
}

tmb_gradient_norm <- function(obj, par) {
  if (!length(par)) return(list(value = 0, error = NULL))
  grad <- try(obj$gr(par), silent = TRUE)
  if (inherits(grad, "try-error")) {
    return(list(value = NA_real_, error = conditionMessage(attr(grad, "condition"))))
  }
  list(value = if (length(grad)) max(abs(grad)) else 0, error = NULL)
}

local_conditional_gradient <- function(obj) {
  if (!length(obj$env$random)) return(list(value = 0, gradient = numeric(), error = NULL))
  grad <- try(as.numeric(obj$env$f(obj$env$last.par, order = 1)), silent = TRUE)
  if (inherits(grad, "try-error")) {
    return(list(value = NA_real_, gradient = NULL,
                error = conditionMessage(attr(grad, "condition"))))
  }
  list(value = max(abs(grad[obj$env$random])), gradient = grad, error = NULL)
}

local_conditional_warmstart <- function(data, parameters, parameter_map, control) {
  started <- proc.time()[["elapsed"]]
  diagnostics <- list(attempted = TRUE, accepted = FALSE, status = "failed")
  answer <- tryCatch({
    objective <- TMB::MakeADFun(data, parameters, map = parameter_map,
                                random = NULL, DLL = "alfak2", silent = TRUE)
    initial_objective <- objective$fn(objective$par)
    opt <- stats::nlminb(objective$par, objective$fn, objective$gr,
                         control = control)
    final_objective <- objective$fn(opt$par)
    gradient <- as.numeric(objective$gr(opt$par))
    accepted <- all(is.finite(opt$par)) && all(is.finite(gradient)) && is.finite(final_objective) &&
      final_objective <= initial_objective
    list(parameters = if (accepted) objective$env$parList(opt$par) else parameters,
      diagnostics = list(attempted = TRUE, accepted = accepted,
        status = if (accepted) "accepted_seed" else "rejected_seed",
        convergence = opt$convergence, message = opt$message,
        initial_objective = initial_objective, objective = final_objective,
        gradient_norm = if (length(gradient)) max(abs(gradient)) else 0,
        iterations = opt$iterations, evaluations = opt$evaluations))
  }, error = function(e) list(parameters = parameters,
    diagnostics = c(diagnostics, list(error = conditionMessage(e)))))
  answer$diagnostics$elapsed_seconds <- proc.time()[["elapsed"]] - started
  answer
}

run_local_tmb_attempt <- function(obj, par, control, n, require_latent = FALSE) {
  opt <- if (length(par)) {
    stats::nlminb(par, obj$fn, obj$gr, control = control)
  } else {
    list(par = par, objective = obj$fn(par), convergence = 0L,
         message = "All hierarchical scales fixed; conditional random effects optimized.")
  }
  sdrep <- try(suppressWarnings(TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = FALSE)),
               silent = TRUE)
  gradient <- tmb_gradient_norm(obj, opt$par)
  conditional_gradient <- local_conditional_gradient(obj)
  list(
    opt = opt,
    sdrep = if (inherits(sdrep, "try-error")) NULL else sdrep,
    sdreport_error = if (inherits(sdrep, "try-error")) conditionMessage(attr(sdrep, "condition")) else NULL,
    frep = {
      out <- extract_tmb_vector(sdrep, "f_report", n)
      if (is.null(out)) out <- extract_tmb_vector(sdrep, "f", n)
      out
    },
    latent_frep = extract_tmb_vector(sdrep, "f_latent", n),
    requires_latent = require_latent,
    gradient_norm = gradient$value,
    gradient_error = gradient$error,
    has_random_effects = length(obj$env$random) > 0L,
    conditional_gradient_norm = conditional_gradient$value,
    conditional_gradient_error = conditional_gradient$error,
    pd_hessian = !inherits(sdrep, "try-error") && isTRUE(sdrep$pdHess),
    scale_boundary = names(opt$par)[grepl("^log_(sigma|tau)", names(opt$par)) &
                                      opt$par < log(1e-4)]
  )
}

assess_local_covariance <- function(attempt, gradient_tolerance) {
  if (!isTRUE(attempt$opt$convergence == 0)) return("untrusted_nonconverged")
  if (is.finite(gradient_tolerance) &&
      (!is.finite(attempt$gradient_norm) || attempt$gradient_norm > gradient_tolerance)) {
    return("untrusted_gradient")
  }
  if (isTRUE(attempt$has_random_effects) && is.finite(gradient_tolerance) &&
      (!is.finite(attempt$conditional_gradient_norm) ||
       attempt$conditional_gradient_norm > gradient_tolerance)) {
    return("untrusted_inner_gradient")
  }
  if (is.null(attempt$frep)) return("untrusted_sdreport_missing")
  if (length(attempt$scale_boundary)) return("untrusted_scale_boundary")
  if (!isTRUE(attempt$pd_hessian)) return("untrusted_hessian")
  if (any(!is.finite(attempt$frep$sd) | attempt$frep$sd <= 0)) {
    return("untrusted_sdreport_nonfinite")
  }
  if (isTRUE(attempt$requires_latent) && is.null(attempt$latent_frep)) {
    return("untrusted_sdreport_missing")
  }
  if (isTRUE(attempt$requires_latent) &&
      any(!is.finite(attempt$latent_frep$sd) | attempt$latent_frep$sd <= 0)) {
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
    gradient_error = attempt$gradient_error,
    conditional_gradient_norm = attempt$conditional_gradient_norm,
    conditional_gradient_error = attempt$conditional_gradient_error,
    sdreport_error = attempt$sdreport_error,
    pd_hessian = attempt$pd_hessian,
    scale_boundary = attempt$scale_boundary,
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

resolve_support_tier_multiplier <- function(graph, support_tier_f_sd_multiplier = NULL) {
  out <- rep(1, length(graph$labels))
  if (is.null(support_tier_f_sd_multiplier)) return(out)
  if (is.null(names(support_tier_f_sd_multiplier)) || any(!nzchar(names(support_tier_f_sd_multiplier)))) {
    stop("`support_tier_f_sd_multiplier` must be a named numeric vector.", call. = FALSE)
  }
  multiplier_names <- names(support_tier_f_sd_multiplier)
  if (anyDuplicated(multiplier_names)) {
    stop("`support_tier_f_sd_multiplier` names must be unique.", call. = FALSE)
  }
  support_tier_f_sd_multiplier <- stats::setNames(as.numeric(support_tier_f_sd_multiplier), multiplier_names)
  if (any(!is.finite(support_tier_f_sd_multiplier) | support_tier_f_sd_multiplier <= 0)) {
    stop("`support_tier_f_sd_multiplier` values must be positive finite numbers.", call. = FALSE)
  }
  idx <- match(as.character(graph$support_tier), names(support_tier_f_sd_multiplier))
  out[!is.na(idx)] <- support_tier_f_sd_multiplier[idx[!is.na(idx)]]
  out
}

resolve_residual_prior_sd <- function(graph, borrowed_residual_sd = NULL, weakly_supported_residual_sd = NULL) {
  out <- rep(0, length(graph$labels))
  if (!is.null(borrowed_residual_sd)) {
    validate_scalar(as.numeric(borrowed_residual_sd), "borrowed_residual_sd", lower = .Machine$double.eps)
    out[graph$support_distance > 0L] <- as.numeric(borrowed_residual_sd)
  }
  if (!is.null(weakly_supported_residual_sd)) {
    validate_scalar(as.numeric(weakly_supported_residual_sd), "weakly_supported_residual_sd", lower = .Machine$double.eps)
    out[as.character(graph$support_tier) == "weakly_supported"] <- as.numeric(weakly_supported_residual_sd)
  }
  out
}

resolve_f_center_weights <- function(local_centering, local_centering_weight_mode,
                                     graph, y0, y1, effective_count_total) {
  local_centering <- match.arg(local_centering, c("none", "direct_weighted_mean", "reference_direct", "observed_weighted_mean"))
  local_centering_weight_mode <- match.arg(local_centering_weight_mode, c("effective_count", "count_total", "uniform"))
  n <- length(graph$labels)
  out <- rep(0, n)
  if (identical(local_centering, "none")) return(out)
  base_weight <- switch(
    local_centering_weight_mode,
    effective_count = as.numeric(effective_count_total),
    count_total = as.numeric(y0 + y1),
    uniform = rep(1, n)
  )
  base_weight[!is.finite(base_weight) | base_weight < 0] <- 0
  support_tier <- as.character(graph$support_tier)
  direct <- support_tier == "directly_informed" | as.integer(graph$support_distance) == 0L
  observed <- (y0 + y1) > 0
  if (identical(local_centering, "direct_weighted_mean")) {
    out[direct] <- base_weight[direct]
  } else if (identical(local_centering, "observed_weighted_mean")) {
    out[observed] <- base_weight[observed]
  } else if (identical(local_centering, "reference_direct")) {
    candidates <- which(direct)
    if (length(candidates)) {
      ref <- candidates[which.max(base_weight[candidates])]
      out[ref] <- 1
    }
  }
  if (!any(out > 0)) out[] <- 0
  out
}

safe_rowsum_local <- function(x, group, n) {
  out <- numeric(n)
  if (!length(x) || !length(group)) return(out)
  tab <- rowsum(as.numeric(x), group = as.integer(group), reorder = FALSE)
  out[as.integer(rownames(tab))] <- as.numeric(tab[, 1])
  out
}

local_gradient_diagnostics <- function(obj, opt, graph, y0, y1, effective_count_total,
                                       eta_prior_mean, eta_prior_sd_vec, report,
                                       summary = NULL, top_n = 20) {
  grad <- try(obj$gr(opt$par), silent = TRUE)
  if (inherits(grad, "try-error")) {
    return(list(error = conditionMessage(attr(grad, "condition"))))
  }
  par <- local_tmb_parameter_list(obj, opt$par)
  # Marginal gradients only contain free hyperparameters, rather than the full
  # vector of conditional random-effect modes. Keep the diagnostic blocks aligned.
  gpl <- lapply(par, function(x) rep(0, length(x)))
  for (nm in unique(names(opt$par))) {
    gpl[[nm]] <- as.numeric(grad)[names(opt$par) == nm]
  }
  conditional_gradient <- local_conditional_gradient(obj)
  random_blocks <- unique(names(obj$env$last.par)[obj$env$random])
  for (nm in random_blocks) {
    gpl[[nm]] <- if (is.null(conditional_gradient$gradient)) {
      rep(NA_real_, length(par[[nm]]))
    } else conditional_gradient$gradient[names(obj$env$last.par) == nm]
  }
  max_abs <- function(x) if (length(x) && any(is.finite(x))) max(abs(as.numeric(x)), na.rm = TRUE) else NA_real_
  block <- data.frame(
    grad_eta_max_abs = max_abs(gpl$eta),
    grad_f_max_abs = max_abs(gpl$f),
    grad_delta_context_max_abs = max_abs(gpl$delta_context),
    grad_mu_group_max_abs = max_abs(gpl$mu_group),
    grad_log_sigma_neighbor_abs = max_abs(gpl$log_sigma_neighbor),
    grad_log_sigma_anchor_abs = max_abs(gpl$log_sigma_anchor),
    grad_log_tau_group_max_abs = max_abs(gpl$log_tau_group),
    conditional_gradient_norm = conditional_gradient$value,
    global_gradient_norm = max_abs(grad),
    stringsAsFactors = FALSE
  )
  block_names <- c("eta", "f", "delta_context", "mu_group", "log_sigma_neighbor", "log_sigma_anchor", "log_tau_group")
  block_values <- c(
    block$grad_eta_max_abs,
    block$grad_f_max_abs,
    block$grad_delta_context_max_abs,
    block$grad_mu_group_max_abs,
    block$grad_log_sigma_neighbor_abs,
    block$grad_log_sigma_anchor_abs,
    block$grad_log_tau_group_max_abs
  )
  block$max_gradient_block_name <- block_names[which.max(replace(block_values, !is.finite(block_values), -Inf))]

  support_tier <- as.character(graph$support_tier)
  support_distance <- as.integer(graph$support_distance)
  support_scope <- ifelse(support_tier == "directly_informed", "direct",
                          ifelse(support_tier == "local_borrowed", "local_borrowed",
                                 ifelse(support_tier == "weakly_supported", "weakly_supported", "other")))
  grad_f <- as.numeric(gpl$f)
  by_tier <- do.call(rbind, lapply(split(seq_along(grad_f), interaction(support_tier, support_scope, drop = TRUE)), function(ii) {
    data.frame(
      support_tier = support_tier[ii[[1]]],
      support_scope = support_scope[ii[[1]]],
      n_nodes = length(ii),
      grad_f_max_abs = max_abs(grad_f[ii]),
      grad_f_median_abs = if (any(is.finite(grad_f[ii]))) stats::median(abs(grad_f[ii]), na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  block$grad_f_direct_max_abs <- max_abs(grad_f[support_scope == "direct"])
  block$grad_f_local_borrowed_max_abs <- max_abs(grad_f[support_scope == "local_borrowed"])
  block$grad_f_weakly_supported_max_abs <- max_abs(grad_f[support_scope == "weakly_supported"])
  block$grad_f_graph_borrowed_max_abs <- max_abs(grad_f[support_scope == "other"])

  n <- length(graph$labels)
  parent_count <- tabulate(as.integer(unlist(graph$parent_to0)) + 1L, nbins = n)
  child_count <- tabulate(as.integer(unlist(graph$parent_from0)) + 1L, nbins = n)
  parent_weight_sum <- safe_rowsum_local(as.numeric(unlist(graph$parent_weight)), as.integer(unlist(graph$parent_to0)) + 1L, n)
  child_weight_sum <- safe_rowsum_local(as.numeric(unlist(graph$parent_weight)), as.integer(unlist(graph$parent_from0)) + 1L, n)
  n_transition_in <- tabulate(as.integer(unlist(graph$transition_to0)) + 1L, nbins = n)
  n_transition_out <- tabulate(as.integer(unlist(graph$transition_from0)) + 1L, nbins = n)
  ord <- order(abs(grad_f), decreasing = TRUE)
  ord <- ord[seq_len(min(top_n, length(ord)))]
  fitness_mean <- if (!is.null(summary) && "fitness_mean" %in% names(summary)) as.numeric(summary$fitness_mean) else as.numeric(report$f_report)
  fitness_sd <- if (!is.null(summary) && "fitness_sd" %in% names(summary)) as.numeric(summary$fitness_sd) else rep(NA_real_, n)
  top <- data.frame(
    node_id = ord,
    karyotype = as.character(graph$labels)[ord],
    support_tier = support_tier[ord],
    support_distance = support_distance[ord],
    count_t0 = as.numeric(y0)[ord],
    count_t1 = as.numeric(y1)[ord],
    count_total = as.numeric(y0 + y1)[ord],
    effective_count_total = as.numeric(effective_count_total)[ord],
    pi0 = as.numeric(report$pi0)[ord],
    pi1 = as.numeric(report$pi1)[ord],
    fitness_mean = fitness_mean[ord],
    fitness_sd = fitness_sd[ord],
    grad_f = grad_f[ord],
    eta = as.numeric(par$eta)[ord],
    grad_eta = as.numeric(gpl$eta)[ord],
    eta_prior_mean = as.numeric(eta_prior_mean)[ord],
    eta_prior_sd = as.numeric(eta_prior_sd_vec)[ord],
    parent_count = parent_count[ord],
    parent_weight_sum = parent_weight_sum[ord],
    child_count = child_count[ord],
    child_weight_sum = child_weight_sum[ord],
    n_transition_in = n_transition_in[ord],
    n_transition_out = n_transition_out[ord],
    is_observed = as.logical((y0 + y1) > 0)[ord],
    stringsAsFactors = FALSE
  )
  list(
    gradient = as.numeric(grad),
    parameter_block_gradient = gpl,
    gradient_block_summary = block,
    grad_f_by_support_tier = by_tier,
    top_gradient_nodes = top
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
#' @param dm_concentration Dirichlet-multinomial concentration value or grid when
#'   requested. Multiple values are fit explicitly and the finite objective minimum
#'   is selected.
#' @param observation_weight_mode How observation weights enter the
#'   Dirichlet-multinomial likelihood. `"likelihood"` uses a weighted likelihood
#'   score that preserves the original count scale; `"fractional_count"` uses the
#'   historical fractional-count pseudo-likelihood.
#' @param control List passed to `stats::nlminb`.
#' @param inner_control TMB inner Newton controls for Laplace random-effect
#'   optimization. The default preserves a strict `tol = 1e-10`; reducing
#'   `maxit` does not bypass conditional-gradient validity checks.
#' @param conditional_warmstart `"auto"` uses first-derivative conditional
#'   optimization to initialize structured Laplace fits with at least 50 nodes.
#'   `"always"` allows smaller fixed-scale problems; `"none"` disables it.
#'   Hierarchical scales must all be fixed for this initialization to be used.
#' @param inference_method `"laplace"` marginalizes initial abundances, fitness,
#'   and context effects while estimating hierarchical scales. `"joint_map"`
#'   retains the legacy joint optimization for controlled comparisons.
#' @param forward_model Evolution model passed to the compiled likelihood.
#'   `"division_process"` uses repeated division-driven transport;
#'   `"selection_transport"` retains the historical endpoint transport model.
#' @param forward_step Maximum numerical time step for the division-process
#'   solver. Decrease this value to check numerical convergence.
#' @param population_sizes Optional positive numeric pair `c(N0, Nb)` containing
#'   actual initial and final population sizes, not the numbers of sampled cells.
#' @param population_log_growth_sd Standard deviation of the observation error
#'   for `log(Nb/N0)`, used only when `population_sizes` is supplied.
#' @param oracle_initial_probability Diagnostic-only named initial probability
#'   vector on every graph node, including true zeros. When supplied, initial
#'   frequencies are fixed exactly rather than estimated from pseudo-counts.
#'   Results use additional oracle information and are not a fair algorithm
#'   comparison against ordinary count-only fits.
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
#' @param local_parameterization Local TMB fitness parameterization. `"f"` uses
#'   per-time fitness directly. `"g_equivalent"` optimizes `g = dt * f` but
#'   reports fitness on the original per-time scale.
#' @param return_optimizer_diagnostics Whether to attach gradient block and top
#'   node diagnostics to `diagnostics`.
#' @param return_tmb_objects Whether to attach TMB optimizer objects to the
#'   returned fit. Defaults to `FALSE` because these objects are large.
#' @param support_tier_f_sd_multiplier Optional named multiplier for local
#'   fitness prior standard deviations by support tier.
#' @param borrowed_residual_sd Optional extra residual prior standard deviation
#'   around parent-propagated means for borrowed nodes.
#' @param weakly_supported_residual_sd Optional extra residual prior standard
#'   deviation for weakly supported nodes.
#' @param local_centering Optional strong fitness centering penalty for local
#'   identifiability probes.
#' @param local_centering_weight Penalty weight used when `local_centering` is
#'   not `"none"`.
#' @param local_centering_weight_mode Weighting scheme for centering penalties.
#' @param fixed_sigma_anchor,fixed_sigma_neighbor,fixed_tau_group Optional fixed
#'   local scale values. `NA` keeps the corresponding scale estimated.
#' @param sdreport_mode TMB ADREPORT selection used for covariance diagnostics.
#'   Defaults to the historical all-fitness plus probability reports.
#' @param initial_jitter_sd Optional standard deviation for benchmark/debug
#'   multistart jitter applied to eta, fitness, and context-effect initial values.
#' @param initial_seed Optional seed used when `initial_jitter_sd > 0`.
#' @param structured_prior Internal joint-landscape prior specification. Ordinary
#'   local fits leave this as `NULL`.
#'
#' @return An `alfak2_local_fit` object.
#' @export
fit_local_posterior <- function(data,
                                graph = NULL,
                                observation_model = c("multinomial", "dirichlet_multinomial"),
                                dm_concentration = 200,
                                observation_weight_mode = c("likelihood", "fractional_count"),
                                control = list(eval.max = 500, iter.max = 500),
                                inner_control = list(maxit = 1000, tol = 1e-10),
                                conditional_warmstart = c("auto", "none", "always"),
                                inference_method = c("laplace", "joint_map"),
                                forward_model = c("selection_transport", "division_process"),
                                forward_step = 0.1,
                                population_sizes = NULL,
                                population_log_growth_sd = 0.1,
                                oracle_initial_probability = NULL,
                                gradient_tolerance = 1e-3,
                                retry_on_untrusted_covariance = TRUE,
                                retry_control = list(eval.max = 2000, iter.max = 2000),
                                eta_prior_sd = 5,
                                eta_borrowed_prior_mean = -6,
                                eta_borrowed_prior_sd = 1.5,
                                eta_distance_penalty = 0.75,
                                local_parameterization = c("f", "g_equivalent"),
                                return_optimizer_diagnostics = FALSE,
                                return_tmb_objects = FALSE,
                                support_tier_f_sd_multiplier = NULL,
                                borrowed_residual_sd = NULL,
                                weakly_supported_residual_sd = NULL,
                                local_centering = c("none", "direct_weighted_mean", "reference_direct", "observed_weighted_mean"),
                                local_centering_weight = 0,
                                local_centering_weight_mode = c("effective_count", "count_total", "uniform"),
                                fixed_sigma_anchor = NA_real_,
                                fixed_sigma_neighbor = NA_real_,
                                fixed_tau_group = NA_real_,
                                sdreport_mode = c("all_f_current", "none", "direct_f_only", "direct_and_local_borrowed_f", "pi0_pi1_only", "fitness_only"),
                                initial_jitter_sd = 0,
                                initial_seed = NULL,
                                structured_prior = NULL) {
  if (!inherits(data, "alfak2_data")) {
    stop("`data` must be an alfak2_data object.", call. = FALSE)
  }
  observation_model <- match_observation_model(observation_model)
  observation_weight_mode <- match_observation_weight_mode(observation_weight_mode)
  inference_method <- match.arg(inference_method)
  conditional_warmstart <- match.arg(conditional_warmstart)
  if (!is.list(inner_control) ||
      (length(inner_control) && (is.null(names(inner_control)) || any(!nzchar(names(inner_control)))))) {
    stop("`inner_control` must be a named list of TMB inner optimizer controls.", call. = FALSE)
  }
  inner_control <- utils::modifyList(list(maxit = 1000, tol = 1e-10), inner_control)
  validate_scalar(inner_control$maxit, "inner_control$maxit", lower = 1)
  validate_scalar(inner_control$tol, "inner_control$tol", lower = .Machine$double.eps)
  forward_model <- match.arg(forward_model)
  validate_scalar(population_log_growth_sd, "population_log_growth_sd", lower = .Machine$double.eps)
  observed_log_growth <- 0
  log_growth_sd <- -1
  if (!is.null(population_sizes)) {
    if (!is.numeric(population_sizes) || length(population_sizes) != 2L ||
        any(!is.finite(population_sizes)) || any(population_sizes <= 0)) {
      stop("`population_sizes` must be a positive finite numeric pair c(N0, Nb) of actual population sizes.", call. = FALSE)
    }
    population_sizes <- stats::setNames(as.numeric(population_sizes), c("N0", "Nb"))
    observed_log_growth <- log(population_sizes[2L]) - log(population_sizes[1L])
    log_growth_sd <- as.numeric(population_log_growth_sd)
  }
  validate_scalar(forward_step, "forward_step", lower = .Machine$double.eps)
  forward_steps <- max(4, ceiling(data$dt / forward_step))
  if (!is.finite(forward_steps) || forward_steps > .Machine$integer.max) {
    stop("`forward_step` requests too many numerical steps.", call. = FALSE)
  }
  local_parameterization <- match.arg(local_parameterization)
  local_centering <- match.arg(local_centering)
  local_centering_weight_mode <- match.arg(local_centering_weight_mode)
  sdreport_mode <- match.arg(sdreport_mode)
  validate_scalar(as.numeric(local_centering_weight), "local_centering_weight", lower = 0)
  if (identical(forward_model, "division_process") &&
      !identical(local_centering, "none") && local_centering_weight > 0)
    stop("Relative-fitness centering is not a birth-rate constraint; use population growth information with division dynamics.", call. = FALSE)
  if (identical(forward_model, "division_process") &&
      !identical(local_centering, "none") && local_centering_weight > 0) {
    stop("Fitness centering penalties are only supported for `forward_model = \"selection_transport\"`; division-process birth rates have a physical baseline.", call. = FALSE)
  }
  validate_scalar(as.numeric(initial_jitter_sd), "initial_jitter_sd", lower = 0)
  initial_jitter_sd <- as.numeric(initial_jitter_sd)
  if (!is.null(initial_seed)) {
    if (!is.numeric(initial_seed) || length(initial_seed) != 1L || !is.finite(initial_seed)) {
      stop("`initial_seed` must be NULL or one finite numeric seed.", call. = FALSE)
    }
    initial_seed <- as.integer(initial_seed)
  }
  validate_optional_positive <- function(x, name) {
    if (is.null(x) || length(x) != 1L || is.na(x)) return(NA_real_)
    validate_scalar(as.numeric(x), name, lower = .Machine$double.eps)
    as.numeric(x)
  }
  fixed_sigma_anchor <- validate_optional_positive(fixed_sigma_anchor, "fixed_sigma_anchor")
  fixed_sigma_neighbor <- validate_optional_positive(fixed_sigma_neighbor, "fixed_sigma_neighbor")
  fixed_tau_group <- validate_optional_positive(fixed_tau_group, "fixed_tau_group")
  sdreport_mode_code <- match(sdreport_mode, c("none", "direct_f_only", "direct_and_local_borrowed_f", "all_f_current", "pi0_pi1_only", "fitness_only")) - 1L
  if (!is.logical(return_optimizer_diagnostics) || length(return_optimizer_diagnostics) != 1L || is.na(return_optimizer_diagnostics)) {
    stop("`return_optimizer_diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_tmb_objects) || length(return_tmb_objects) != 1L || is.na(return_tmb_objects)) {
    stop("`return_tmb_objects` must be TRUE or FALSE.", call. = FALSE)
  }
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
  use_oracle_initial_probability <- !is.null(oracle_initial_probability)
  oracle_pi0 <- numeric(n)
  if (use_oracle_initial_probability) {
    if (!is.numeric(oracle_initial_probability) || !is.null(dim(oracle_initial_probability)) ||
        length(oracle_initial_probability) != n || is.null(names(oracle_initial_probability)) ||
        anyDuplicated(names(oracle_initial_probability)) ||
        !setequal(names(oracle_initial_probability), as.character(graph$labels)) ||
        any(!is.finite(oracle_initial_probability)) || any(oracle_initial_probability < 0) ||
        abs(sum(oracle_initial_probability) - 1) > 1e-8) {
      stop("`oracle_initial_probability` must be a finite non-negative named probability vector summing to one on every graph node.", call. = FALSE)
    }
    oracle_pi0 <- as.numeric(oracle_initial_probability[as.character(graph$labels)])
    oracle_pi0 <- oracle_pi0 / sum(oracle_pi0)
  }
  structured <- !is.null(structured_prior)
  structure_design <- matrix(numeric(), nrow = n, ncol = 0L)
  structure_offset <- numeric(n)
  structure_precision <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n),
                                              x = rep(1, n), dims = c(n, n))
  structure_logdet <- 0
  structure_coef_sd <- numeric()
  if (structured) {
    required <- c("design", "offset", "precision", "logdet", "coef_sd")
    if (!is.list(structured_prior) || !all(required %in% names(structured_prior))) {
      stop("`structured_prior` must provide design, offset, precision, logdet, and coef_sd.", call. = FALSE)
    }
    structure_design <- as.matrix(structured_prior$design)
    structure_offset <- as.numeric(structured_prior$offset)
    structure_precision <- methods::as(Matrix::Matrix(structured_prior$precision, sparse = TRUE), "generalMatrix")
    structure_logdet <- as.numeric(structured_prior$logdet)
    structure_coef_sd <- as.numeric(structured_prior$coef_sd)
    if (nrow(structure_design) != n || length(structure_offset) != n ||
        !identical(dim(structure_precision), c(n, n)) ||
        length(structure_logdet) != 1L ||
        length(structure_coef_sd) != ncol(structure_design) ||
        any(!is.finite(c(structure_design, structure_offset, structure_precision@x,
                        structure_logdet, structure_coef_sd))) ||
        any(structure_coef_sd <= 0)) {
      stop("`structured_prior` dimensions or values are invalid.", call. = FALSE)
    }
  }
  observed_index <- match(data$labels, as.character(graph$labels))
  if (anyNA(observed_index)) {
    missing <- data$labels[is.na(observed_index)]
    stop("Graph is missing observed karyotypes: ", paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  }
  y0 <- numeric(n)
  y1 <- numeric(n)
  y0[observed_index] <- data$counts[, 1]
  y1[observed_index] <- data$counts[, 2]
  if (use_oracle_initial_probability && any(y0 > 0 & oracle_pi0 == 0)) {
    stop("Oracle initial probabilities contradict positive time-zero counts.", call. = FALSE)
  }
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
  f_prior_sd_multiplier <- resolve_support_tier_multiplier(graph, support_tier_f_sd_multiplier)
  residual_prior_sd <- resolve_residual_prior_sd(graph, borrowed_residual_sd, weakly_supported_residual_sd)
  f_center_weights <- resolve_f_center_weights(local_centering, local_centering_weight_mode, graph, y0, y1, effective_count_total)
  p0 <- (y0_init + 0.5) / sum(y0_init + 0.5)
  p1 <- (y1_init + 0.5) / sum(y1_init + 0.5)
  f0 <- log(p1) - log(p0)
  f0 <- f0 - mean(f0)
  f0_scale <- f0 / max(data$dt, .Machine$double.eps)
  if (identical(forward_model, "division_process")) {
    if (!is.null(population_sizes) && observed_log_growth > 0) {
      # A unit-aware starting point only; the actual likelihood also accounts
      # for daughter loss and transport rather than assuming pure net growth.
      baseline_birth <- as.numeric(observed_log_growth) / data$dt
      f0_scale <- log(pmax(0.05 * baseline_birth, baseline_birth + f0_scale))
    } else {
      f0_scale <- pmax(-2, pmin(2, f0_scale))
    }
    f0 <- data$dt * f0_scale
  }
  if (structured) {
    f0_scale <- f0_scale - structure_offset
    f0 <- data$dt * f0_scale
  }
  n_context <- length(graph$context_label)
  n_group <- max(unlist(graph$context_group0), 0L) + 1L
  parameters <- list(
    eta = log(p0),
    f = if (identical(local_parameterization, "g_equivalent")) f0 else f0_scale,
    delta_context = rep(0, n_context),
    mu_group = rep(0, n_group),
    log_sigma_neighbor = log(if (is.finite(fixed_sigma_neighbor)) pmax(fixed_sigma_neighbor - 1e-5, 1e-8) else 0.35),
    log_sigma_anchor = log(if (is.finite(fixed_sigma_anchor)) pmax(fixed_sigma_anchor - 1e-5, 1e-8) else 0.6),
    log_tau_group = rep(log(if (is.finite(fixed_tau_group)) pmax(fixed_tau_group - 1e-5, 1e-8) else 0.2), n_group),
    structured_coefficients = rep(0, ncol(structure_design))
  )
  if (initial_jitter_sd > 0) {
    if (!is.null(initial_seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(initial_seed)
    }
    parameters$eta <- parameters$eta + stats::rnorm(length(parameters$eta), 0, initial_jitter_sd)
    parameters$f <- parameters$f + stats::rnorm(length(parameters$f), 0, initial_jitter_sd)
    parameters$delta_context <- parameters$delta_context + stats::rnorm(length(parameters$delta_context), 0, initial_jitter_sd)
  }
  if (use_oracle_initial_probability) parameters$eta[] <- 0
  division_transition <- if (identical(forward_model, "division_process")) {
    division_transition_data(graph, beta = data$beta)
  } else NULL
  tmb_data <- list(
    y0 = as.numeric(y0),
    y1 = as.numeric(y1),
    n_nodes = as.integer(n),
    trans_from = if (is.null(division_transition)) as.integer(unlist(graph$transition_from0)) else as.integer(division_transition$from0),
    trans_to = if (is.null(division_transition)) as.integer(unlist(graph$transition_to0)) else as.integer(division_transition$to0),
    trans_weight = if (is.null(division_transition)) as.numeric(unlist(graph$transition_weight)) else as.numeric(division_transition$weight),
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
    dm_concentration = as.numeric(dm_concentration_grid[1]),
    local_parameterization = as.integer(identical(local_parameterization, "g_equivalent")),
    f_prior_sd_multiplier = as.numeric(f_prior_sd_multiplier),
    residual_prior_sd = as.numeric(residual_prior_sd),
    f_center_weights = as.numeric(f_center_weights),
    f_centering_weight = as.numeric(local_centering_weight),
    fixed_sigma_neighbor = if (is.finite(fixed_sigma_neighbor)) as.numeric(fixed_sigma_neighbor) else -1,
    fixed_sigma_anchor = if (is.finite(fixed_sigma_anchor)) as.numeric(fixed_sigma_anchor) else -1,
    fixed_tau_group = if (is.finite(fixed_tau_group)) as.numeric(fixed_tau_group) else -1,
    sdreport_mode = as.integer(sdreport_mode_code)
  )
  tmb_data$forward_model <- as.integer(match(forward_model, c("division_process", "selection_transport")) - 1L)
  tmb_data$forward_steps <- as.integer(forward_steps)
  tmb_data$prior_model <- as.integer(structured)
  tmb_data$structure_design <- structure_design
  tmb_data$structure_offset <- structure_offset
  tmb_data$structure_precision <- structure_precision
  tmb_data$structure_logdet <- structure_logdet
  tmb_data$structure_coef_sd <- structure_coef_sd
  tmb_data$observed_log_growth <- as.numeric(observed_log_growth)
  tmb_data$log_growth_sd <- log_growth_sd
  tmb_data$oracle_pi0 <- oracle_pi0
  tmb_data$use_oracle_pi0 <- as.integer(use_oracle_initial_probability)
  parameter_map <- local_tmb_parameter_map(parameters, fixed_sigma_anchor,
                                            fixed_sigma_neighbor, fixed_tau_group)
  random_effects <- if (identical(inference_method, "laplace")) {
    if (structured) c("eta", "f", "structured_coefficients") else c("eta", "f", "delta_context", "mu_group")
  } else NULL
  if (structured) {
    for (nm in c("delta_context", "mu_group", "log_sigma_anchor", "log_sigma_neighbor", "log_tau_group")) {
      parameter_map[[nm]] <- factor(rep(NA_integer_, length(parameters[[nm]])))
    }
  }
  if (use_oracle_initial_probability) parameter_map$eta <- factor(rep(NA_integer_, n))
  random_effects <- random_effects[vapply(parameters[random_effects], length, integer(1)) > 0L]
  random_effects <- random_effects[!vapply(random_effects, function(nm)
    !is.null(parameter_map[[nm]]) && all(is.na(parameter_map[[nm]])), logical(1))]
  scale_names <- c("log_sigma_anchor", "log_sigma_neighbor", "log_tau_group")
  all_scales_fixed <- all(vapply(scale_names, function(nm)
    !is.null(parameter_map[[nm]]) && all(is.na(parameter_map[[nm]])), logical(1)))
  if (identical(conditional_warmstart, "always") && !all_scales_fixed) {
    stop("`conditional_warmstart = \"always\"` requires every hierarchical scale to be fixed.", call. = FALSE)
  }
  use_warmstart <- identical(inference_method, "laplace") && all_scales_fixed &&
    (identical(conditional_warmstart, "always") ||
     (identical(conditional_warmstart, "auto") && structured && n >= 50L))

  concentration_results <- lapply(dm_concentration_grid, function(phi) {
    candidate_data <- tmb_data
    candidate_data$dm_concentration <- as.numeric(phi)
    warmstart <- if (use_warmstart) {
      local_conditional_warmstart(candidate_data, parameters, parameter_map, control)
    } else list(parameters = parameters, diagnostics = list(
      attempted = FALSE, accepted = FALSE, status = "disabled_or_ineligible", elapsed_seconds = 0))
    obj <- tryCatch(
      TMB::MakeADFun(candidate_data, warmstart$parameters, map = parameter_map,
                     random = random_effects, DLL = "alfak2", silent = TRUE,
                     inner.control = inner_control),
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
      run_local_tmb_attempt(obj, obj$par, control, n,
                            require_latent = identical(forward_model, "division_process")),
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
        run_local_tmb_attempt(obj, attempt$opt$par, merged_retry_control, n,
                              require_latent = identical(forward_model, "division_process")),
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
    covariance_fallback <- !identical(covariance_status, "TMB_sdreport")
    warning_log_path <- NA_character_
    warning_diagnostics <- NULL
    if (isTRUE(covariance_fallback)) {
      warning_diagnostics <- list(
        stage = "local_covariance",
        dm_concentration = as.numeric(phi),
        covariance_status = covariance_status,
        retry_attempted = retry_attempted,
        retry_reason = retry_reason,
        attempts = attempts
      )
      warning_log_path <- alfak2_warn(
        "Local TMB covariance is untrusted; using prior-scale uncertainty fallback.",
        diagnostics = warning_diagnostics
      )
    }
    report <- try({
      local_tmb_parameter_list(obj, attempt$opt$par)
      obj$report()
    }, silent = TRUE)
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
    fitness_sd_source <- if (isTRUE(covariance_fallback)) "fallback_prior_scale" else "TMB_sdreport"
    f_sd <- if (isTRUE(covariance_fallback)) {
      fallback_local_fitness_sd(report, graph$support_distance)
    } else {
      as.numeric(attempt$frep$sd)
    }
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
      objective = as.numeric(attempt$opt$objective),
      f_sd = f_sd,
      covariance_fallback = covariance_fallback,
      fitness_sd_source = fitness_sd_source,
      warning_log_path = warning_log_path,
      warning_diagnostics = warning_diagnostics,
      warmstart = warmstart$diagnostics
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
  covariance_fallback <- selected$covariance_fallback
  fitness_sd_source <- selected$fitness_sd_source
  warning_log_path <- selected$warning_log_path
  warning_diagnostics <- selected$warning_diagnostics
  attempts <- selected$attempts
  retry_attempted <- selected$retry_attempted
  retry_reason <- selected$retry_reason
  dm_concentration_selected <- selected$dm_concentration
  opt <- attempt$opt
  plist <- local_tmb_parameter_list(obj, opt$par)
  report <- selected$report
  f_mean <- if (!is.null(report$f_report)) as.numeric(report$f_report) else as.numeric(plist$f)
  plist_raw <- plist
  plist$f_raw <- plist_raw$f
  plist$f <- f_mean
  pi0_report <- as.numeric(report$pi0)
  pi1_report <- as.numeric(report$pi1)
  f_sd <- selected$f_sd
  fitness_covariance <- extract_local_fitness_covariance(attempt$sdrep, attempt$frep,
                                                        covariance_status)
  latent_fitness_covariance <- if (identical(forward_model, "division_process")) {
    extract_local_fitness_covariance(attempt$sdrep, attempt$latent_frep, covariance_status)
  } else NULL
  latent_fitness_mean <- if (!is.null(attempt$latent_frep)) attempt$latent_frep$mean else NULL
  conf_low <- f_mean - 1.959963984540054 * f_sd
  conf_high <- f_mean + 1.959963984540054 * f_sd
  point_estimate_type <- if (identical(inference_method, "laplace")) "conditional_mode_at_marginal_scale_mode" else "joint_map"
  if (identical(forward_model, "division_process") && !covariance_fallback) {
    latent_sd <- attempt$latent_frep$sd
    conf_low <- exp(latent_fitness_mean - 1.959963984540054 * latent_sd)
    conf_high <- exp(latent_fitness_mean + 1.959963984540054 * latent_sd)
    if (identical(inference_method, "laplace") && !is.null(latent_fitness_covariance)) {
      rate_mean <- exp(latent_fitness_mean + 0.5 * diag(latent_fitness_covariance))
      rate_cov <- tcrossprod(rate_mean) * expm1(latent_fitness_covariance)
      if (all(is.finite(rate_mean)) && all(is.finite(rate_cov))) {
        f_mean <- as.numeric(rate_mean)
        fitness_covariance <- rate_cov
        f_sd <- sqrt(pmax(0, diag(rate_cov)))
        point_estimate_type <- "laplace_lognormal_mean"
      } else {
        covariance_status <- "untrusted_lognormal_moments"
        covariance_fallback <- TRUE
        fitness_covariance <- NULL
        fitness_sd_source <- "fallback_prior_scale"
        f_sd <- fallback_local_fitness_sd(report, graph$support_distance)
        warning_log_path <- alfak2_warn("Lognormal posterior moments are non-finite; fitness intervals are withheld.",
                                        diagnostics = list(stage = "local_lognormal_moments"))
      }
    }
  }
  if (isTRUE(covariance_fallback)) conf_low[] <- conf_high[] <- NA_real_
  if (!is.null(fitness_covariance)) {
    dimnames(fitness_covariance) <- list(as.character(graph$labels), as.character(graph$labels))
  }
  if (!is.null(latent_fitness_covariance)) {
    dimnames(latent_fitness_covariance) <- list(as.character(graph$labels), as.character(graph$labels))
    names(latent_fitness_mean) <- as.character(graph$labels)
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
    effective_count_total = as.numeric(effective_count_total),
    eta_prior_mean = as.numeric(eta_prior_mean),
    eta_prior_sd = as.numeric(eta_prior_sd_vec),
    pi0 = as.numeric(pi0_report),
    pi1 = as.numeric(pi1_report),
    is_observed = as.logical((y0 + y1) > 0),
    fitness_mean = f_mean,
    fitness_sd = f_sd,
    fitness_sd_source = fitness_sd_source,
    conf_low = conf_low,
    conf_high = conf_high,
    interval_valid = !covariance_fallback,
    covariance_status = covariance_status,
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    convergence = opt$convergence,
    message = opt$message,
    objective = opt$objective,
    gradient_norm = attempt$gradient_norm,
    gradient_tolerance = gradient_tolerance,
    conditional_gradient_norm = attempt$conditional_gradient_norm,
    conditional_gradient_error = attempt$conditional_gradient_error,
    pd_hessian = attempt$pd_hessian,
    scale_boundary = attempt$scale_boundary,
    covariance_status = covariance_status,
    covariance_fallback = covariance_fallback,
    fitness_sd_source = fitness_sd_source,
    warning_log_path = warning_log_path,
    warning_diagnostics = warning_diagnostics,
    retry_attempted = retry_attempted,
    retry_reason = retry_reason,
    attempts = attempts,
    observation_model = observation_model,
    observation_weight_mode = observation_weight_mode,
    likelihood_model = likelihood_model,
    local_parameterization = local_parameterization,
    inference_method = inference_method,
    forward_model = forward_model,
    forward_step = forward_step,
    forward_steps = as.integer(forward_steps),
    division_transport = if (is.null(division_transition)) NULL else division_transition[c("retained", "viable", "leakage")],
    fitness_interpretation = if (identical(forward_model, "division_process")) {
      "birth_rate"
    } else if (is.null(population_sizes)) "relative_net_growth_rate" else "net_growth_rate",
    absolute_level = if (!is.null(population_sizes)) {
      "informed_by_population_growth"
    } else if (identical(forward_model, "division_process")) "prior_identified_without_external_growth" else "prior_centered_relative_scale",
    population_growth = list(population_sizes = population_sizes,
      observed_log_growth = if (is.null(population_sizes)) NA_real_ else as.numeric(observed_log_growth),
      log_growth_sd = if (is.null(population_sizes)) NA_real_ else log_growth_sd),
    use_oracle_initial_probability = use_oracle_initial_probability,
    oracle_information = if (use_oracle_initial_probability) "exact_initial_frequencies_for_diagnostic_only" else "none",
    random_effects = random_effects,
    structured_prior = structured,
    point_estimate_type = point_estimate_type,
    probability_estimate_type = "conditional_mode",
    uncertainty_method = if (identical(point_estimate_type, "laplace_lognormal_mean")) {
      "laplace_log_gaussian_moments"
    } else if (identical(inference_method, "laplace")) "laplace_delta_method" else "joint_hessian_delta_method",
    hierarchical_scales = if (structured) NULL else list(
      sigma_anchor = as.numeric(report$sigma_anchor),
      sigma_neighbor = as.numeric(report$sigma_neighbor),
      tau_group = if (is.finite(fixed_tau_group)) rep(fixed_tau_group, n_group) else 1e-5 + exp(plist$log_tau_group)
    ),
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
    f_prior = list(
      support_tier_sd_multiplier = support_tier_f_sd_multiplier,
      borrowed_residual_sd = borrowed_residual_sd,
      weakly_supported_residual_sd = weakly_supported_residual_sd
    ),
    local_centering = list(
      mode = local_centering,
      weight = as.numeric(local_centering_weight),
      weight_mode = local_centering_weight_mode,
      n_weighted_nodes = sum(f_center_weights > 0)
    ),
    fixed_scale = list(
      sigma_anchor = fixed_sigma_anchor,
      sigma_neighbor = fixed_sigma_neighbor,
      tau_group = fixed_tau_group
    ),
    sdreport_mode = sdreport_mode,
    inner_control = inner_control,
    conditional_warmstart = c(list(mode = conditional_warmstart), selected$warmstart),
    initial_jitter = list(
      sd = initial_jitter_sd,
      seed = initial_seed
    ),
    pi0_support_mass = support_mass_by_tier(pi0_report, graph$support_tier),
    pi1_support_mass = support_mass_by_tier(pi1_report, graph$support_tier)
  )
  if (isTRUE(return_optimizer_diagnostics)) {
    diagnostics$optimizer <- local_gradient_diagnostics(
      obj = obj,
      opt = opt,
      graph = graph,
      y0 = y0,
      y1 = y1,
      effective_count_total = effective_count_total,
      eta_prior_mean = eta_prior_mean,
      eta_prior_sd_vec = eta_prior_sd_vec,
      report = report,
      summary = summary
    )
  }
  out <- list(
    data = data,
    graph = graph,
    summary = summary,
    parameter_mode = plist,
    fitness_covariance = fitness_covariance,
    latent_fitness_mean = latent_fitness_mean,
    latent_fitness_covariance = latent_fitness_covariance,
    sdreport = attempt$sdrep,
    posterior_predictive = report[c("pi0", "pi1", "log_total_growth")],
    diagnostics = diagnostics
  )
  if (isTRUE(return_tmb_objects)) {
    out$optimizer <- list(
      obj = obj,
      opt = opt,
      par = opt$par,
      gradient = try(obj$gr(opt$par), silent = TRUE)
    )
  }
  new_alfak2_local_fit(out)
}
