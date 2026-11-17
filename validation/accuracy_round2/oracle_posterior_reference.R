#!/usr/bin/env Rscript
# Fixed-data posterior integration, not SBC and not a coverage experiment.
oracle_normal_quadrature <- function(order) {
  jacobi <- matrix(0, order, order)
  jacobi[cbind(seq_len(order - 1L), 2:order)] <- sqrt(seq_len(order - 1L))
  jacobi <- jacobi + t(jacobi)
  e <- eigen(jacobi, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1L, ]^2)
}

oracle_weighted_quantile <- function(value, weight, probabilities = c(0.025, 0.975)) {
  index <- order(value)
  cumulative <- cumsum(weight[index]) / sum(weight)
  vapply(probabilities, function(p) value[index[which(cumulative >= p)[1L]]], numeric(1))
}

oracle_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

oracle_division_reference <- function(prior_mean, prior_covariance, counts,
    dt, order = 72L, eta_sd = 2, eta_step = 0.025, eta_limit = 5,
    known_pi0 = NULL, observed_growth = NULL, growth_sd = 0.1,
    proposal_mean = prior_mean, proposal_covariance = prior_covariance) {
  quadrature <- oracle_normal_quadrature(order)
  z <- as.matrix(expand.grid(quadrature$nodes, quadrature$nodes))
  latent <- sweep(z %*% chol(proposal_covariance), 2L, proposal_mean, "+")
  rates <- exp(latent)
  prior_weight <- as.numeric(outer(quadrature$weights, quadrature$weights))
  log_normal_density <- function(x, mean, covariance) {
    factor <- chol(covariance)
    standardized <- sweep(x, 2L, mean, "-") %*% solve(factor)
    -0.5 * rowSums(standardized^2) - sum(log(diag(factor))) - log(2 * pi)
  }
  # Importance quadrature is centered on a broad Laplace-guided proposal, but
  # every weight uses the original normalized prior and exact beta-zero count
  # likelihood. It does not replace the target by the Gaussian proposal.
  log_density_ratio <- log_normal_density(latent, prior_mean, prior_covariance) -
    log_normal_density(latent, proposal_mean, proposal_covariance)
  if (is.null(known_pi0)) {
    eta <- seq(-eta_limit, eta_limit, by = eta_step)
    trapezoid <- rep(eta_step, length(eta)); trapezoid[c(1L, length(eta))] <- eta_step / 2
    log_eta_weight <- stats::dnorm(eta, 0, sqrt(2) * eta_sd, log = TRUE) + log(trapezoid)
  } else {
    eta <- stats::qlogis(known_pi0[1L])
    log_eta_weight <- 0
  }
  p0 <- stats::plogis(eta)
  log_weight0 <- log_eta_weight + stats::dbinom(counts[1L, 1L], sum(counts[, 1L]), p0, log = TRUE)
  likelihood <- vapply(seq_len(nrow(rates)), function(i) {
    log_weight <- log_weight0 + stats::dbinom(counts[1L, 2L], sum(counts[, 2L]),
      stats::plogis(eta + dt * (rates[i, 1L] - rates[i, 2L])), log = TRUE)
    if (!is.null(observed_growth)) {
      left <- log(p0) + dt * rates[i, 1L]
      right <- log1p(-p0) + dt * rates[i, 2L]
      shift <- pmax(left, right)
      growth <- shift + log(exp(left - shift) + exp(right - shift))
      log_weight <- log_weight + stats::dnorm(growth, observed_growth, growth_sd, log = TRUE)
    }
    oracle_logsumexp(log_weight)
  }, numeric(1))
  log_weight <- log(prior_weight) + log_density_ratio + likelihood
  weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
  mean <- as.numeric(crossprod(weight, rates))
  centered <- sweep(rates, 2L, mean, "-")
  covariance <- crossprod(centered * sqrt(weight))
  contrast <- 0.5 * (rates[, 1L] - rates[, 2L])
  contrast_mean <- sum(weight * contrast)
  contrast_sd <- sqrt(sum(weight * (contrast - contrast_mean)^2))
  list(mean = mean, covariance = covariance, latent_mean = as.numeric(crossprod(weight, latent)),
       contrast_mean = contrast_mean, contrast_sd = contrast_sd,
       contrast_quantile = oracle_weighted_quantile(contrast, weight),
       rate_quantile = t(vapply(seq_len(2L), function(j)
         oracle_weighted_quantile(rates[, j], weight), numeric(2))),
       order = order, eta_step = eta_step, eta_limit = eta_limit)
}

oracle_prior_sampling_check <- function(prior_mean, prior_covariance, counts, dt,
    eta_sd = 2, known_pi0 = NULL, observed_growth = NULL, growth_sd = 0.1,
    n = 500000L, seed = 9040534L) {
  set.seed(seed)
  latent <- sweep(matrix(stats::rnorm(2L * n), nrow = n) %*% chol(prior_covariance),
                  2L, prior_mean, "+")
  birth <- exp(latent)
  eta <- if (is.null(known_pi0)) stats::rnorm(n, sd = sqrt(2) * eta_sd) else rep(stats::qlogis(known_pi0[1L]), n)
  p0 <- stats::plogis(eta)
  log_weight <- stats::dbinom(counts[1L, 1L], sum(counts[, 1L]), p0, log = TRUE) +
    stats::dbinom(counts[1L, 2L], sum(counts[, 2L]),
      stats::plogis(eta + dt * (birth[, 1L] - birth[, 2L])), log = TRUE)
  if (!is.null(observed_growth)) {
    left <- log(p0) + dt * birth[, 1L]; right <- log1p(-p0) + dt * birth[, 2L]
    shift <- pmax(left, right)
    growth <- shift + log(exp(left - shift) + exp(right - shift))
    log_weight <- log_weight + stats::dnorm(growth, observed_growth, growth_sd, log = TRUE)
  }
  weight <- exp(log_weight - max(log_weight)); weight <- weight / sum(weight)
  values <- cbind(birth, 0.5 * (birth[, 1L] - birth[, 2L]))
  mean <- as.numeric(crossprod(weight, values))
  centered <- sweep(values, 2L, mean, "-")
  variance <- as.numeric(crossprod(weight, centered^2))
  mean_mcse <- sqrt(colSums(centered^2 * weight^2))
  variance_mcse <- sqrt(colSums(sweep(centered^2, 2L, variance, "-")^2 * weight^2))
  list(mean = mean, sd = sqrt(variance), mean_mcse = mean_mcse,
       sd_mcse = variance_mcse / (2 * sqrt(variance)), ess = 1 / sum(weight^2), n = n)
}

oracle_posterior_reference_main <- function(args = commandArgs(TRUE)) {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
  repo <- normalizePath(file.path(dirname(script), "..", ".."))
  source(file.path(dirname(script), "oracle_diagnostics.R"))
  run_root <- oracle_read_argument(args, "run-root", "/private/tmp/alfak2_iteration2_20261104_0518")
  output <- oracle_read_argument(args, "output", file.path(run_root, "oracle-reference"))
  frozen <- oracle_read_argument(args, "v1-source", file.path(run_root, "v1-source"))
  dll <- oracle_read_argument(args, "dll", file.path(run_root, "current-lib/alfak2/libs/alfak2.so"))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  requireNamespace("Rcpp"); requireNamespace("TMB")
  binary <- dyn.load(normalizePath(dll))
  package <- new.env(parent = globalenv())
  for (name in names(getDLLRegisteredRoutines(binary)$`.Call`))
    assign(name, getDLLRegisteredRoutines(binary)$`.Call`[[name]], package)
  files <- sort(list.files(file.path(frozen, "R"), "[.]R$", full.names = TRUE))
  for (file in files) sys.source(file, package)
  sys.source(file.path(repo, "R/local-fit.R"), package)
  # The graph and all priors are specified before the fixed count observations.
  graph <- package$build_karyotype_graph(c("1", "2"), beta = 0, shell_depth = 0,
                                        min_cn = 1, max_cn = 2)
  prior <- package$build_joint_landscape_prior(graph, forward_model = "division_process",
    baseline_sd = 0.5, main_sd = 0.07, curvature_sd = 0.04,
    interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25,
    collapse_gaussian = FALSE)
  covariance <- solve(as.matrix(prior$precision)) +
    tcrossprod(sweep(prior$design, 2L, prior$coef_sd, "*"))
  counts <- matrix(c(20, 11, 20, 29), 2L, byrow = TRUE,
    dimnames = list(c("1", "2"), c("t0", "t1")))
  dt <- 1.5; eta_sd <- 2
  population_growth <- log(0.5 * exp(dt * 0.7) + 0.5 * exp(dt * 1.3))
  population <- c(1000, 1000 * exp(population_growth))
  data <- package$prepare_alfak2_data(counts, dt = dt, beta = 0)
  arms <- list(count_only = list(), known_pi0 = list(oracle_initial_probability = c("1" = 0.5, "2" = 0.5)),
    population_growth = list(population_sizes = population, population_log_growth_sd = 0.1),
    known_pi0_and_growth = list(oracle_initial_probability = c("1" = 0.5, "2" = 0.5),
      population_sizes = population, population_log_growth_sd = 0.1))
  rows <- convergence <- results <- checks <- list()
  for (arm in names(arms)) {
    fit <- do.call(package$fit_local_posterior, c(list(data = data, graph = graph,
      structured_prior = prior, forward_model = "division_process", eta_prior_sd = eta_sd,
      sdreport_mode = "fitness_only", conditional_warmstart = "none"), arms[[arm]]))
    if (!identical(fit$diagnostics$covariance_status, "TMB_sdreport")) stop("Untrusted reference comparison")
    reference <- lapply(c(40L, 64L, 96L, 160L), function(order)
      oracle_division_reference(prior$offset, covariance, counts, dt, order = order,
        eta_sd = eta_sd, eta_step = if (order >= 96L) 0.0125 else 0.025,
        eta_limit = if (order >= 96L) 6 else 5,
        known_pi0 = if (grepl("known_pi0", arm, fixed = TRUE)) c(0.5, 0.5) else NULL,
        observed_growth = if (grepl("growth", arm, fixed = TRUE)) population_growth else NULL,
        proposal_mean = fit$latent_fitness_mean,
        proposal_covariance = 2 * fit$latent_fitness_covariance))
    final <- reference[[length(reference)]]
    check <- oracle_prior_sampling_check(prior$offset, covariance, counts, dt,
      eta_sd = eta_sd,
      known_pi0 = if (grepl("known_pi0", arm, fixed = TRUE)) c(0.5, 0.5) else NULL,
      observed_growth = if (grepl("growth", arm, fixed = TRUE)) population_growth else NULL)
    checks[[arm]] <- data.frame(arm = arm, quantity = c("birth_1", "birth_2", "relative_contrast_1"),
      quadrature_mean = c(final$mean, final$contrast_mean), sampling_mean = check$mean,
      sampling_mean_mcse = check$mean_mcse, quadrature_sd = c(sqrt(diag(final$covariance)), final$contrast_sd),
      sampling_sd = check$sd, sampling_sd_mcse = check$sd_mcse, prior_samples = check$n,
      importance_ess = check$ess)
    contrast_vector <- c(0.5, -0.5)
    laplace_contrast_mean <- sum(contrast_vector * fit$summary$fitness_mean)
    laplace_contrast_sd <- sqrt(as.numeric(crossprod(contrast_vector,
      fit$fitness_covariance %*% contrast_vector)))
    rows[[arm]] <- data.frame(arm = arm, quantity = c("birth_1", "birth_2", "relative_contrast_1"),
      reference_mean = c(final$mean, final$contrast_mean),
      laplace_mean = c(fit$summary$fitness_mean, laplace_contrast_mean),
      reference_sd = c(sqrt(diag(final$covariance)), final$contrast_sd),
      laplace_sd = c(fit$summary$fitness_sd, laplace_contrast_sd))
    for (j in seq_len(length(reference) - 1L)) convergence[[paste(arm, j)]] <- data.frame(arm = arm,
      coarse_order = reference[[j]]$order, fine_order = final$order,
      max_mean_difference = max(abs(reference[[j]]$mean - final$mean)),
      max_covariance_difference = max(abs(reference[[j]]$covariance - final$covariance)),
      contrast_mean_difference = abs(reference[[j]]$contrast_mean - final$contrast_mean),
      contrast_sd_difference = abs(reference[[j]]$contrast_sd - final$contrast_sd),
      max_contrast_quantile_difference = max(abs(reference[[j]]$contrast_quantile - final$contrast_quantile)))
    results[[arm]] <- list(reference = reference, laplace_summary = fit$summary,
      laplace_covariance = fit$fitness_covariance, diagnostics = fit$diagnostics,
      independent_prior_sampling_check = check)
    cat("posterior reference", arm, "completed\n")
  }
  comparison <- do.call(rbind, rows)
  comparison$mean_difference <- comparison$laplace_mean - comparison$reference_mean
  comparison$sd_ratio <- comparison$laplace_sd / comparison$reference_sd
  utils::write.csv(comparison, file.path(output, "oracle_reference_comparison.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, convergence), file.path(output, "oracle_reference_convergence.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, checks), file.path(output, "oracle_reference_sampling_check.csv"), row.names = FALSE)
  saveRDS(list(results = results, prior = prior, aggregate_prior_covariance = covariance,
    counts = counts, dt = dt, eta_prior_sd = eta_sd,
    limitation = "Fixed-data numerical posterior reference; not SBC or repeated-sampling coverage.",
    interval_reference_status = "not_certified_discontinuous_Gauss_Hermite_quantiles_are_not_used_for_coverage_claims",
    source_hashes = tools::md5sum(c(files, file.path(repo, "R/local-fit.R"), script, dll))),
    file.path(output, "oracle_reference_details.rds"))
}

if (sys.nframe() == 0L) oracle_posterior_reference_main()
