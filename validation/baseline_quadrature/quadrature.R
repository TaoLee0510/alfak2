baseline_component_prior <- function(prior, value) {
  column <- which(colnames(prior$design) == "intercept")
  if (length(column) != 1L || any(abs(prior$design[, column] - 1) > 1e-12))
    stop("One proper independent intercept prior is required")
  output <- prior
  output$design <- prior$design[, -column, drop = FALSE]
  output$coef_sd <- prior$coef_sd[-column]
  output$offset <- prior$offset + value
  output$diagnostics$conditional_intercept <- value
  output$diagnostics$prior_birth_mean <- NA_real_
  output
}

baseline_proposal <- function(case, fit) {
  column <- which(colnames(case$prior$design) == "intercept")
  sd0 <- case$prior$coef_sd[column]
  cross <- rep(sd0^2, case$n)
  coefficient <- as.numeric(solve(case$covariance, cross))
  mean <- sum(coefficient * (fit$latent_fitness_mean - case$prior$offset))
  conditional_var <- sd0^2 - sum(cross * coefficient)
  variance <- conditional_var + as.numeric(crossprod(coefficient,
    fit$latent_fitness_covariance %*% coefficient))
  if (!is.finite(variance) || variance <= 0) stop("Invalid intercept proposal")
  list(mean = mean, sd = sqrt(variance), prior_sd = sd0,
    conditional_prior_variance = conditional_var)
}

baseline_combine <- function(components, weight) {
  mean <- Reduce(`+`, Map(function(x, w) w * x$mean, components, weight))
  latent_mean <- Reduce(`+`, Map(function(x, w) w * x$latent_mean, components, weight))
  covariance <- Reduce(`+`, Map(function(x, w)
    w * (x$covariance + tcrossprod(x$mean - mean)), components, weight))
  latent_covariance <- Reduce(`+`, Map(function(x, w)
    w * (x$latent_covariance + tcrossprod(x$latent_mean - latent_mean)), components, weight))
  list(mean = mean, covariance = covariance, latent_mean = latent_mean,
    latent_covariance = latent_covariance, components = components, weights = weight)
}

baseline_integrate <- function(case, original_fit, runtime, orders = c(9L, 17L, 33L, 65L),
                                cache_directory = NULL, cache_key = NULL) {
  if (!is.null(cache_directory) && is.null(cache_key)) stop("Conditional caches require a source/case fingerprint")
  proposal <- baseline_proposal(case, original_fit)
  history <- list(); previous <- NULL; conditional_status <- list()
  A <- diag(case$n) - matrix(case$weight, case$n, case$n, byrow = TRUE)
  for (order in orders) {
    rule <- oracle_normal_quadrature(order)
    values <- proposal$mean + proposal$sd * rule$nodes
    pieces <- list(); log_weights <- numeric(); kept <- integer()
    for (j in seq_along(values)) {
      if (rule$weights[j] == 0) next
      value <- values[j]
      cache <- if (is.null(cache_directory)) NULL else file.path(cache_directory,
        sprintf("%s_order%d_node%d.rds", case$id, order, j))
      node_key <- list(source_case = cache_key, intercept = value)
      if (!is.null(cache) && file.exists(cache)) {
        cached <- readRDS(cache)
        if (!identical(cached$key, node_key)) stop("Conditional-node cache fingerprint mismatch")
        fit <- cached$fit
      } else {
        data <- runtime$pkg$prepare_alfak2_data(case$counts, dt = case$specification$dt, beta = case$beta)
        fit <- runtime$pkg$fit_local_posterior(data, graph = case$graph,
          structured_prior = baseline_component_prior(case$prior, value),
          forward_model = "division_process", forward_step = .1,
          eta_prior_sd = 2, eta_borrowed_prior_mean = 0, eta_borrowed_prior_sd = 2,
          eta_distance_penalty = 0,
          population_sizes = if (case$growth) c(1000, 1000 * exp(case$specification$growth)) else NULL,
          population_log_growth_sd = .1, sdreport_mode = "fitness_only",
          conditional_warmstart = "none")
        if (!is.null(cache)) saveRDS(list(key = node_key, fit = fit), cache)
      }
      conditional_status[[length(conditional_status) + 1L]] <- data.frame(order = order,
        node = j, intercept = value, objective = fit$diagnostics$objective,
        covariance_status = fit$diagnostics$covariance_status)
      if (!identical(fit$diagnostics$covariance_status, "TMB_sdreport"))
        return(list(status = "untrusted_conditional_node", proposal = proposal,
          history = if (length(history)) do.call(rbind, history) else data.frame(),
          conditional_status = do.call(rbind, conditional_status)))
      pieces[[length(pieces) + 1L]] <- list(mean = fit$summary$fitness_mean,
        covariance = fit$fitness_covariance, latent_mean = fit$latent_fitness_mean,
        latent_covariance = fit$latent_fitness_covariance, intercept = value)
      log_weights <- c(log_weights, log(rule$weights[j]) - fit$diagnostics$objective +
        dnorm(value, 0, proposal$prior_sd, log = TRUE) -
        dnorm(value, proposal$mean, proposal$sd, log = TRUE))
      kept <- c(kept, j)
    }
    weight <- exp(log_weights - max(log_weights)); weight <- weight / sum(weight)
    current <- baseline_combine(pieces, weight)
    current$log_evidence <- max(log_weights) + log(sum(exp(log_weights - max(log_weights))))
    mean_change <- covariance_change <- relative_change <- Inf
    if (!is.null(previous)) {
      mean_change <- max(abs(current$mean - previous$mean) / pmax(1, abs(current$mean)))
      covariance_change <- norm(current$covariance - previous$covariance, "F") / norm(current$covariance, "F")
      relative_change <- norm(A %*% (current$covariance - previous$covariance) %*% t(A), "F") /
        norm(A %*% current$covariance %*% t(A), "F")
    }
    passed <- mean_change < 1e-3 && covariance_change < .005 && relative_change < .005
    history[[as.character(order)]] <- data.frame(order = order, mean_change = mean_change,
      covariance_relative_change = covariance_change, contrast_covariance_relative_change = relative_change,
      numerical_convergence = passed, min_weight = min(weight), max_weight = max(weight))
    if (passed) break
    previous <- current
  }
  c(current, list(status = if (passed) "passed_quadrature_refinement" else "quadrature_not_converged",
    proposal = proposal, history = do.call(rbind, history),
    conditional_status = do.call(rbind, conditional_status)))
}
