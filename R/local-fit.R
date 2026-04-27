extract_tmb_vector <- function(sdrep, name, n) {
  if (inherits(sdrep, "try-error") || is.null(sdrep)) return(NULL)
  tab <- try(suppressWarnings(as.data.frame(summary(sdrep, "report"))), silent = TRUE)
  if (inherits(tab, "try-error") || !nrow(tab)) return(NULL)
  rn <- rownames(tab)
  idx <- which(rn == name | grepl(paste0("^", name, "(\\.|\\[)"), rn))
  if (length(idx) != n) return(NULL)
  list(mean = tab$Estimate[idx], sd = tab$`Std. Error`[idx])
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
#'
#' @return An `alfak2_local_fit` object.
#' @export
fit_local_posterior <- function(data,
                                graph = NULL,
                                observation_model = c("multinomial", "dirichlet_multinomial"),
                                dm_concentration = 200,
                                control = list(eval.max = 500, iter.max = 500)) {
  if (!inherits(data, "alfak2_data")) {
    stop("`data` must be an alfak2_data object.", call. = FALSE)
  }
  observation_model <- match_observation_model(observation_model)
  if (is.null(graph)) graph <- build_karyotype_graph(data, shell_depth = 2)
  if (!inherits(graph, "alfak2_graph")) {
    stop("`graph` must be an alfak2_graph object.", call. = FALSE)
  }

  n <- length(graph$labels)
  y0 <- numeric(n)
  y1 <- numeric(n)
  y0[graph$observed_index] <- data$counts[, 1]
  y1[graph$observed_index] <- data$counts[, 2]

  p0 <- (y0 + 0.5) / sum(y0 + 0.5)
  p1 <- (y1 + 0.5) / sum(y1 + 0.5)
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
    anchor_prior_scale = 1.0,
    mu_prior_scale = 0.5,
    scale_prior_scale = 1.0,
    observation_model = as.integer(observation_model == "dirichlet_multinomial"),
    dm_concentration = as.numeric(dm_concentration)
  )

  obj <- TMB::MakeADFun(tmb_data, parameters, DLL = "alfak2", silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = control)
  sdrep <- try(suppressWarnings(TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = FALSE)), silent = TRUE)
  frep <- extract_tmb_vector(sdrep, "f", n)
  if (inherits(sdrep, "try-error") || is.null(frep)) {
    stop("TMB sdreport did not return local fitness uncertainty.", call. = FALSE)
  }
  plist <- obj$env$parList(opt$par)
  f_mean <- as.numeric(plist$f)
  f_sd <- as.numeric(frep$sd)
  f_sd[!is.finite(f_sd)] <- NA_real_
  report <- try(obj$report(opt$par), silent = TRUE)

  summary <- data.frame(
    node_id = seq_len(n),
    karyotype = as.character(graph$labels),
    support_tier = as.character(graph$support_tier),
    support_distance = as.integer(graph$support_distance),
    fitness_mean = f_mean,
    fitness_sd = f_sd,
    conf_low = f_mean - 1.959963984540054 * f_sd,
    conf_high = f_mean + 1.959963984540054 * f_sd,
    stringsAsFactors = FALSE
  )
  grad <- try(obj$gr(opt$par), silent = TRUE)
  diagnostics <- list(
    convergence = opt$convergence,
    message = opt$message,
    objective = opt$objective,
    gradient_norm = if (inherits(grad, "try-error")) NA_real_ else max(abs(grad)),
    covariance_status = "TMB_sdreport",
    observation_model = observation_model
  )
  new_alfak2_local_fit(list(
    data = data,
    graph = graph,
    summary = summary,
    parameter_mode = plist,
    sdreport = if (inherits(sdrep, "try-error")) NULL else sdrep,
    posterior_predictive = if (inherits(report, "try-error")) NULL else report[c("pi0", "pi1")],
    diagnostics = diagnostics
  ))
}
