# Build a proper, dimension-scaled prior on the fixed target domain.
# The optional rotation conditions on the observed support selected from counts;
# its design nullspace is not a count-likelihood identifiability result.
build_joint_landscape_prior <- function(graph, forward_model = "division_process",
                                        baseline_sd = 0.5, main_sd = 0.07,
                                        curvature_sd = 0.04, interaction_sd = 0.04,
                                        residual_sd = 0.125, edge_sd = 0.25,
                                        include_interactions = TRUE,
                                        interaction_support = c("isotropic", "observed_nullspace"),
                                        unidentified_interaction_sd = 0.2,
                                        collapse_gaussian = FALSE,
                                        covariance_probe_seed = 71109L,
                                        covariance_probes = 64L) {
  for (name in c("baseline_sd", "main_sd", "curvature_sd", "interaction_sd", "residual_sd", "edge_sd"))
    validate_scalar(get(name), name, lower = .Machine$double.eps)
  interaction_support <- match.arg(interaction_support)
  validate_scalar(unidentified_interaction_sd, "unidentified_interaction_sd", lower = .Machine$double.eps)
  x <- graph$karyotypes
  n <- nrow(x); p <- ncol(x)
  if (!is.logical(collapse_gaussian) || length(collapse_gaussian) != 1L || is.na(collapse_gaussian))
    stop("`collapse_gaussian` must be TRUE or FALSE.", call. = FALSE)
  if (isTRUE(collapse_gaussian) && n > 32L)
    stop("Dense Gaussian marginalization is limited to 32 nodes; use the explicit Gaussian representation for larger domains.", call. = FALSE)
  reference <- pmin(graph$max_cn, pmax(graph$min_cn, rep(2, p)))
  dx <- sweep(x, 2L, reference, "-")
  design <- cbind(intercept = 1, dx, dx^2)
  colnames(design) <- c("intercept", paste0("dosage_", seq_len(p)), paste0("curvature_", seq_len(p)))
  coef_sd <- c(baseline_sd, rep(main_sd / sqrt(p), p), rep(curvature_sd / sqrt(p), p))
  interaction_diagnostics <- list(mode = interaction_support, rank = 0L, n_directions = 0L,
                                  n_unidentified = 0L, observed_nullspace_max_error = 0)
  if (isTRUE(include_interactions) && p > 1L) {
    pairs <- utils::combn(p, 2L)
    interactions <- vapply(seq_len(ncol(pairs)), function(j)
      dx[, pairs[1L, j]] * dx[, pairs[2L, j]], numeric(n))
    if (!is.matrix(interactions)) interactions <- matrix(interactions, nrow = n)
    colnames(interactions) <- paste0("interaction_", pairs[1L, ], "_", pairs[2L, ])
    if (identical(interaction_support, "observed_nullspace")) {
      if (length(graph$support_distance) != n || any(!is.finite(graph$support_distance)))
        stop("Observed-nullspace interaction priors require actual observed support distances.", call. = FALSE)
      observed <- which(graph$support_distance == 0L)
      if (!length(observed)) stop("No observed states define the interaction support.", call. = FALSE)
      decomposition <- svd(interactions[observed, , drop = FALSE], nu = 0L, nv = ncol(pairs))
      tolerance <- max(1, decomposition$d) * max(length(observed), ncol(pairs)) * 1e-10
      rank <- sum(decomposition$d > tolerance)
      identified <- if (rank) interactions %*% decomposition$v[, seq_len(rank), drop = FALSE] else matrix(numeric(), n, 0L)
      null <- if (rank < ncol(pairs)) interactions %*% decomposition$v[, (rank + 1L):ncol(pairs), drop = FALSE] else matrix(numeric(), n, 0L)
      error <- if (ncol(null)) max(abs(null[observed, , drop = FALSE])) else 0
      if (error > 1e-8) stop("Interaction nullspace numerical check failed.", call. = FALSE)
      if (ncol(null)) null[observed, ] <- 0
      if (ncol(identified)) colnames(identified) <- paste0("identified_interaction_", seq_len(ncol(identified)))
      if (ncol(null)) colnames(null) <- paste0("unidentified_interaction_", seq_len(ncol(null)))
      interactions <- cbind(identified, null)
      interaction_sd_vector <- c(rep(interaction_sd / sqrt(ncol(pairs)), rank),
        rep(unidentified_interaction_sd / sqrt(ncol(pairs)), ncol(pairs) - rank))
      interaction_diagnostics <- list(mode = interaction_support, rank = rank,
        n_directions = ncol(pairs), n_unidentified = ncol(pairs) - rank,
        observed_nullspace_max_error = error,
        interpretation = "null at observed design rows only; normalization and mutation flow can still inform these directions")
    } else {
      interaction_sd_vector <- rep(interaction_sd / sqrt(ncol(pairs)), ncol(pairs))
      interaction_diagnostics$n_directions <- ncol(pairs)
    }
    design <- cbind(design, interactions)
    coef_sd <- c(coef_sd, interaction_sd_vector)
  }
  # Prior topology is independent of beta: absence of mutation is not evidence
  # that neighboring dosage effects are statistically unrelated.
  from <- to <- integer()
  labels <- as.character(graph$labels)
  for (chr in seq_len(p)) {
    neighbor <- x
    neighbor[, chr] <- neighbor[, chr] + 1L
    j <- match(format_karyotypes(neighbor), labels)
    keep <- which(!is.na(j))
    from <- c(from, keep); to <- c(to, j[keep])
  }
  degree <- tabulate(c(from, to), nbins = n)
  precision <- Matrix::Diagonal(n, x = rep(1 / residual_sd^2, n))
  if (length(from)) {
    edge_scale <- 1 / (edge_sd^2 * max(1, mean(degree)))
    precision <- precision + Matrix::sparseMatrix(
      i = c(from, to, from, to), j = c(from, to, to, from),
      x = c(rep(edge_scale, 2L * length(from)), rep(-edge_scale, 2L * length(from))),
      dims = c(n, n))
  }
  precision <- methods::as(precision, "generalMatrix")
  logdet <- as.numeric(Matrix::determinant(precision, logarithm = TRUE)$modulus)
  if (n <= 500L) {
    residual_covariance <- solve(as.matrix(precision))
    residual_variance <- diag(residual_covariance)
    variance_method <- "exact"
    probe_mcse <- 0
  } else {
    probes <- with_grf_seed(covariance_probe_seed,
      matrix(sample(c(-1, 1), n * covariance_probes, replace = TRUE), nrow = n))
    products <- probes * as.matrix(Matrix::solve(precision, probes))
    residual_variance <- pmin(residual_sd^2, pmax(0, rowMeans(products)))
    probe_mcse <- max(apply(products, 1L, stats::sd) / sqrt(covariance_probes))
    variance_method <- "fixed_seed_hutchinson"
  }
  variance <- rowSums(sweep(design, 2L, coef_sd, "*")^2) + residual_variance
  # The correction keeps marginal E[b(x)]=1 throughout the target domain for
  # the Gaussian log-birth prior, rather than inflating distant birth means.
  offset <- if (identical(forward_model, "division_process")) -0.5 * variance else rep(0, n)
  n_coefficients <- ncol(design)
  representation <- "explicit_gaussian_coefficients_and_graph_residual"
  if (isTRUE(collapse_gaussian) && n <= 500L) {
    # The likelihood sees only X*a+u. Integrate its proper Gaussian components
    # analytically before Laplace approximation, eliminating unused directions
    # while preserving exactly the same prior on every node's latent fitness.
    covariance <- residual_covariance + tcrossprod(sweep(design, 2L, coef_sd, "*"))
    precision <- methods::as(Matrix::Matrix(solve(covariance), sparse = TRUE), "generalMatrix")
    logdet <- -as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
    design <- matrix(numeric(), nrow = n, ncol = 0L)
    coef_sd <- numeric()
    representation <- "analytically_marginalized_gaussian_landscape"
  }
  list(design = design, offset = offset, precision = precision, logdet = logdet,
       coef_sd = coef_sd,
       diagnostics = list(reference = reference, n_coefficients = n_coefficients,
         gaussian_representation = representation,
         interaction = interaction_diagnostics,
         prior_birth_mean = if (identical(forward_model, "division_process") &&
                                  identical(variance_method, "exact")) 1 else NA_real_,
         prior_birth_mean_correction = if (identical(variance_method, "exact")) "exact" else "stochastic_diagonal_approximation",
         residual_variance_method = variance_method, residual_variance_mcse_max = probe_mcse,
         scales = c(baseline = baseline_sd, main = main_sd, curvature = curvature_sd,
                    interaction = interaction_sd, residual = residual_sd, edge = edge_sd)))
}

fit_joint_landscape <- function(data, graph, joint_control = list(),
                                 forward_tolerance = 1e-5,
                                 forward_max_refinements = 3L, ...) {
  args <- list(...)
  if (!is.null(joint_control$model_scales) || !is.null(joint_control$prior_components)) {
    args$forward_tolerance <- forward_tolerance
    args$forward_max_refinements <- forward_max_refinements
    return(fit_joint_landscape_ensemble(data, graph, joint_control, args))
  }
  validate_scalar(forward_tolerance, "forward_tolerance", lower = .Machine$double.eps)
  if (!is.null(args$population_sizes)) {
    growth_sd <- if (is.null(args$population_log_growth_sd)) 0.1 else args$population_log_growth_sd
    forward_tolerance <- min(forward_tolerance, 0.1 * growth_sd)
  }
  forward_max_refinements <- validate_grf_integer(forward_max_refinements, "forward_max_refinements", lower = 0L)
  forward_model <- if (is.null(args$forward_model)) "division_process" else args$forward_model[1L]
  args$forward_model <- forward_model
  if (is.null(args$sdreport_mode)) args$sdreport_mode <- "fitness_only"
  prior <- do.call(build_joint_landscape_prior,
                  c(list(graph = graph, forward_model = forward_model), joint_control))
  refinements <- 0L
  repeat {
    posterior <- do.call(fit_local_posterior,
      c(list(data = data, graph = graph, structured_prior = prior), args))
    if (!identical(forward_model, "division_process")) break
    accuracy <- division_forward_accuracy(posterior, graph, data$dt, forward_tolerance)
    accuracy$refinements <- refinements
    posterior$diagnostics$forward_accuracy <- accuracy
    if (identical(accuracy$status, "passed_step_halving")) break
    if (refinements >= forward_max_refinements) {
      posterior$diagnostics$covariance_status <- "untrusted_forward_accuracy"
      posterior$summary$interval_valid <- FALSE
      posterior$summary$conf_low <- posterior$summary$conf_high <- NA_real_
      posterior$summary$covariance_status <- "untrusted_forward_accuracy"
      break
    }
    args$forward_step <- (if (is.null(args$forward_step)) 0.1 else args$forward_step) / 2
    refinements <- refinements + 1L
  }
  summary <- posterior$summary
  summary$extrapolation_method <- "joint_count_landscape"
  summary$prediction_status <- ifelse(summary$interval_valid, "joint_posterior", "untrusted_joint_posterior")
  observed <- summary$effective_count_total > 0
  global <- list(graph = graph, summary = summary,
    fitness_covariance = posterior$fitness_covariance,
    anchors = summary[observed, , drop = FALSE], hyperparameters = prior$diagnostics$scales,
    tuning_grid = data.frame(),
    diagnostics = c(posterior$diagnostics, list(
      method = "joint_count_landscape", prior = prior$diagnostics,
      prior_handling = "single_joint_count_likelihood_and_structured_prior",
      uncertainty_includes = "latent_and_fixed_prior_scale_laplace_uncertainty",
      scales_selected_on = "user_or_independent_training_data",
      n_predictions = sum(is.finite(summary$fitness_mean)), n_nodes = nrow(summary))))
  class(global) <- "alfak2_global_fit"
  list(local = posterior, global = global)
}

draw_joint_component <- function(component, n_draws) {
  log_scale <- identical(component$forward_model, "division_process")
  mu <- if (log_scale) component$latent_mean else component$mean
  covariance <- if (log_scale) component$latent_covariance else component$covariance
  if (is.null(covariance) || any(!is.finite(covariance)))
    stop("A joint posterior component has no trusted covariance.", call. = FALSE)
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (min(decomposition$values) < -1e-7)
    stop("A joint posterior component covariance is indefinite.", call. = FALSE)
  transform <- sweep(decomposition$vectors, 2L, sqrt(pmax(0, decomposition$values)), "*")
  draws <- sweep(transform %*% matrix(stats::rnorm(length(mu) * n_draws), nrow = length(mu)), 1L, mu, "+")
  if (log_scale) exp(draws) else draws
}

fit_joint_landscape_ensemble <- function(data, graph, joint_control, args) {
  if (!is.null(args$inference_method) && !identical(args$inference_method[1L], "laplace"))
    stop("Joint scale evidence averaging requires `inference_method = \"laplace\"`.", call. = FALSE)
  if (!is.null(args$dm_concentration) && length(args$dm_concentration) != 1L)
    stop("Joint scale evidence averaging requires one shared Dirichlet-multinomial concentration.", call. = FALSE)
  observation_weights <- data$metadata$observation_weights
  if (!is.null(observation_weights) && any(abs(observation_weights - 1) > 1e-12))
    stop("Joint evidence averaging requires an unweighted, normalized count likelihood.", call. = FALSE)
  grouped <- !is.null(joint_control$prior_components)
  if (grouped && !is.null(joint_control$model_scales))
    stop("Choose prior_components or model_scales, not both.", call. = FALSE)
  component_controls <- joint_control$prior_components
  if (grouped && (!is.list(component_controls) || !length(component_controls) ||
      any(!vapply(component_controls, is.list, logical(1)))))
    stop("prior_components must be a nonempty list of named prior-control lists.", call. = FALSE)
  scales <- if (grouped) rep(NA_real_, length(component_controls)) else
    validate_positive_grid(joint_control$model_scales, "joint_control$model_scales")
  labels <- if (grouped && !is.null(names(component_controls))) names(component_controls) else
    paste0("component_", seq_along(scales))
  probabilities <- joint_control$model_prior
  if (is.null(probabilities)) probabilities <- rep(1 / length(scales), length(scales))
  if (length(probabilities) != length(scales) || any(!is.finite(probabilities)) || any(probabilities <= 0))
    stop("`model_prior` must give a positive finite probability for every model scale.", call. = FALSE)
  probabilities <- probabilities / sum(probabilities)
  joint_control$model_scales <- joint_control$model_prior <- joint_control$prior_components <- NULL
  defaults <- list(main_sd = 0.07, curvature_sd = 0.04, interaction_sd = 0.04,
                   residual_sd = 0.125, edge_sd = 0.25)
  actual_controls <- lapply(seq_along(scales), function(k) {
    control <- if (grouped) utils::modifyList(joint_control, component_controls[[k]]) else joint_control
    if (!grouped) for (name in names(defaults))
      control[[name]] <- scales[k] * if (is.null(control[[name]])) defaults[[name]] else control[[name]]
    if (any(c("prior_components", "model_scales", "model_prior") %in% names(control)))
      stop("Nested evidence ensembles are not supported.", call. = FALSE)
    control
  })
  components <- lapply(actual_controls, function(control) {
    do.call(fit_joint_landscape, c(list(data = data, graph = graph, joint_control = control), args))
  })
  trusted <- vapply(components, function(x)
    identical(x$local$diagnostics$covariance_status, "TMB_sdreport"), logical(1))
  if (!all(trusted))
    stop("All joint ensemble components must pass posterior covariance checks; failed components are not silently discarded.", call. = FALSE)
  log_evidence <- -vapply(components, function(x) x$local$diagnostics$objective, numeric(1))
  log_weight <- log(probabilities) + log_evidence
  weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
  means <- do.call(cbind, lapply(components, function(x) x$local$summary$fitness_mean))
  mean <- as.numeric(means %*% weight)
  covariance <- matrix(0, nrow(means), nrow(means))
  descriptors <- lapply(seq_along(components), function(k) {
    component <- components[[k]]$local
    delta <- means[, k] - mean
    covariance <<- covariance + weight[k] * (component$fitness_covariance + tcrossprod(delta))
    list(weight = weight[k], mean = means[, k], covariance = component$fitness_covariance,
         latent_mean = component$latent_fitness_mean,
         latent_covariance = component$latent_fitness_covariance,
         forward_model = component$diagnostics$forward_model,
         log_evidence = log_evidence[k], scale = scales[k], label = labels[k], prior_control = actual_controls[[k]])
  })
  dimnames(covariance) <- list(graph$labels, graph$labels)
  # Raw marginal intervals have an exact mixture CDF under the component
  # Gaussian/lognormal approximations. Reference contrasts require joint draws.
  quantile_at <- function(i, probability) {
    log_scale <- identical(descriptors[[1L]]$forward_model, "division_process")
    mu <- vapply(descriptors, function(d) if (log_scale) d$latent_mean[i] else d$mean[i], numeric(1))
    sd <- vapply(descriptors, function(d) sqrt(if (log_scale) d$latent_covariance[i, i] else d$covariance[i, i]), numeric(1))
    component_quantile <- if (log_scale) exp(mu + sd * stats::qnorm(probability)) else
      mu + sd * stats::qnorm(probability)
    lower <- min(component_quantile); upper <- max(component_quantile)
    if (upper - lower <= 1e-12) return((lower + upper) / 2)
    cdf <- function(value) {
      z <- if (log_scale) log(value) else value
      sum(weight * stats::pnorm(z, mu, sd)) - probability
    }
    stats::uniroot(cdf, c(lower, upper), tol = 1e-9)$root
  }
  interval <- t(vapply(seq_len(nrow(means)), function(i)
    c(quantile_at(i, 0.025), quantile_at(i, 0.975)), numeric(2)))
  result <- components[[1L]]
  for (layer in c("local", "global")) {
    result[[layer]]$summary$fitness_mean <- mean
    result[[layer]]$summary$fitness_sd <- sqrt(pmax(0, diag(covariance)))
    result[[layer]]$summary$conf_low <- interval[, 1L]
    result[[layer]]$summary$conf_high <- interval[, 2L]
    result[[layer]]$fitness_covariance <- covariance
    result[[layer]]$posterior_components <- descriptors
    result[[layer]]$diagnostics$point_estimate_type <- "conditional_graph_Laplace_evidence_mixture_mean"
    result[[layer]]$diagnostics$uncertainty_method <- "component_mixture_including_between_scale_uncertainty"
    result[[layer]]$diagnostics$model_scales <- scales
    result[[layer]]$diagnostics$component_labels <- labels
    result[[layer]]$diagnostics$component_controls <- actual_controls
    result[[layer]]$diagnostics$model_weights <- weight
    result[[layer]]$diagnostics$model_log_evidence <- log_evidence
    result[[layer]]$summary$fitness_sd_source <- "total_component_mixture_covariance"
  }
  result$local$latent_fitness_mean <- result$local$latent_fitness_covariance <- NULL
  result$local$sdreport <- NULL
  result$local$parameter_mode <- NULL
  result$local$diagnostics$objective <- -(max(log_weight) + log(sum(exp(log_weight - max(log_weight)))))
  result$global$diagnostics$objective <- result$local$diagnostics$objective
  result$local$posterior_predictive <- list(
    pi0 = Reduce(`+`, lapply(seq_along(components), function(k) weight[k] * components[[k]]$local$posterior_predictive$pi0)),
    pi1 = Reduce(`+`, lapply(seq_along(components), function(k) weight[k] * components[[k]]$local$posterior_predictive$pi1)),
    log_total_growth = {
      growth <- vapply(components, function(x) as.numeric(x$local$posterior_predictive$log_total_growth), numeric(1))
      center <- max(growth)
      center + log(sum(weight * exp(growth - center)))
    }, averaging = "weighted_component_mode_frequencies_and_total_growth")
  result$local$summary$pi0 <- result$local$posterior_predictive$pi0
  result$local$summary$pi1 <- result$local$posterior_predictive$pi1
  result$global$summary$pi0 <- result$local$summary$pi0
  result$global$summary$pi1 <- result$local$summary$pi1
  result$global$anchors <- result$global$summary[result$global$summary$effective_count_total > 0, , drop = FALSE]
  result$global$hyperparameters <- list(model_scales = scales, component_controls = actual_controls,
    component_labels = labels, model_prior = probabilities, model_weights = weight)
  result
}
