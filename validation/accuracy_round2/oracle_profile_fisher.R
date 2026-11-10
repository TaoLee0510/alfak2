#!/usr/bin/env Rscript
# Local expected-model information planning, not a global-identifiability test.
oracle_fisher_propagator <- function(probability, log_birth, dt, tolerance = 1e-14) {
  n <- length(log_birth)
  generator <- sweep(2 * t(probability) - diag(n), 2L, exp(log_birth), "*")
  # Independent positive Taylor/uniformization matrix exponential. This does
  # not use the inference solver's splitting scheme or its AD derivatives.
  shift <- max(0, -min(diag(generator)))
  positive <- generator + diag(shift, n)
  term <- total <- diag(n)
  for (order in seq_len(1000L)) {
    term <- positive %*% term * dt / order
    total <- total + term
    if (max(abs(term)) <= tolerance * max(1, max(abs(total)))) break
  }
  if (order == 1000L) stop("Independent Fisher propagator failed to converge")
  exp(-shift * dt) * total
}

oracle_fisher_probability <- function(pi0, log_birth, probability, dt) {
  propagator <- oracle_fisher_propagator(probability, log_birth, dt)
  mass <- as.numeric(propagator %*% pi0)
  list(pi0 = pi0, pi1 = mass / sum(mass), propagator = propagator, total = sum(mass))
}

oracle_fisher_jacobian <- function(pi0, log_birth, probability, dt, epsilon = 1e-5) {
  n <- length(pi0)
  prediction <- oracle_fisher_probability(pi0, log_birth, probability, dt)
  j0 <- diag(pi0) - tcrossprod(pi0)
  j1_eta <- (diag(n) - tcrossprod(prediction$pi1, rep(1, n))) %*%
    prediction$propagator %*% j0 / prediction$total
  j1_birth <- vapply(seq_len(n), function(j) {
    delta <- numeric(n); delta[j] <- epsilon
    plus <- oracle_fisher_probability(pi0, log_birth + delta, probability, dt)$pi1
    minus <- oracle_fisher_probability(pi0, log_birth - delta, probability, dt)$pi1
    (plus - minus) / (2 * epsilon)
  }, numeric(n))
  list(j0_eta = j0, j1_eta = j1_eta, j1_birth = j1_birth,
       pi0 = pi0, pi1 = prediction$pi1, prediction = prediction)
}

oracle_fisher_design <- function(jacobian, depths, profile_eta = TRUE, tolerance = 1e-10) {
  n <- length(jacobian$pi0)
  standardize <- function(j, p, count) {
    weight <- numeric(length(p)); positive <- p > 0
    weight[positive] <- sqrt(count / p[positive])
    j * weight
  }
  eta <- rbind(standardize(jacobian$j0_eta, jacobian$pi0, depths[1L]),
               standardize(jacobian$j1_eta, jacobian$pi1, depths[2L]))
  birth <- rbind(matrix(0, n, n),
                 standardize(jacobian$j1_birth, jacobian$pi1, depths[2L]))
  rank <- 0L
  if (profile_eta) {
    decomposition <- svd(eta, nu = min(dim(eta)), nv = 0L)
    cutoff <- max(dim(eta)) * max(1, decomposition$d) * tolerance
    rank <- sum(decomposition$d > cutoff)
    if (rank) {
      basis <- decomposition$u[, seq_len(rank), drop = FALSE]
      birth <- birth - basis %*% crossprod(basis, birth)
    }
  }
  list(design = birth, fisher = crossprod(birth), eta_rank = rank)
}

oracle_fisher_spectrum <- function(design, prior_covariance, log_birth, depth_scale = 1) {
  factor <- t(chol(prior_covariance))
  white_design <- sqrt(depth_scale) * design %*% factor
  fisher <- crossprod(white_design)
  eig <- eigen(fisher, symmetric = TRUE)
  if (min(eig$values) < -1e-8) stop("Expected Fisher construction is not PSD")
  values <- pmax(0, eig$values)
  # A common additive birth shift is d(log b)=1/b, rather than a common
  # multiplicative birth shift d(log b)=1. These coincide only for equal b.
  absolute <- as.numeric(solve(factor, exp(-log_birth)))
  absolute <- absolute / sqrt(sum(absolute^2))
  common_log <- as.numeric(solve(factor, rep(1, length(log_birth))))
  common_log <- common_log / sqrt(sum(common_log^2))
  shape_basis <- qr.Q(qr(matrix(absolute, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  shape <- eigen(crossprod(white_design %*% shape_basis), symmetric = TRUE)
  shape_values <- pmax(0, shape$values)
  cutoff <- max(1, max(values)) * length(values) * 1e-9
  null <- values < cutoff
  weak <- values < 1
  physical_weakest <- exp(log_birth) * as.numeric(factor %*% eig$vectors[, length(values)])
  alignment <- abs(sum(physical_weakest)) /
    (sqrt(length(physical_weakest)) * sqrt(sum(physical_weakest^2)))
  list(fisher = fisher, eigenvalues = values, shape_eigenvalues = shape_values,
    factor = factor, eigenvectors = eig$vectors,
    shape_weakest = as.numeric(factor %*% shape_basis %*% shape$vectors[, length(shape_values)]),
    effective_directions = sum(values / (1 + values)),
    shape_effective_directions = sum(shape_values / (1 + shape_values)),
    rank = sum(!null), rank_cutoff = cutoff,
    fixed_shape_absolute_rayleigh = as.numeric(crossprod(absolute, fisher %*% absolute)),
    common_log_birth_rayleigh = as.numeric(crossprod(common_log, fisher %*% common_log)),
    absolute_projection_on_data_null = sum(crossprod(eig$vectors[, null, drop = FALSE], absolute)^2),
    absolute_projection_on_prior_dominated = sum(crossprod(eig$vectors[, weak, drop = FALSE], absolute)^2),
    weakest_physical_absolute_alignment = alignment)
}

oracle_profile_fisher_main <- function(args = commandArgs(TRUE)) {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
  repo <- normalizePath(file.path(dirname(script), "..", ".."))
  source(file.path(dirname(script), "oracle_diagnostics.R"))
  root <- oracle_read_argument(args, "run-root", "/private/tmp/alfak2_iteration2_20261104_0518")
  output <- oracle_read_argument(args, "output", file.path(root, "oracle-fisher"))
  frozen <- oracle_read_argument(args, "v1-source", file.path(root, "v1-source"))
  dll <- oracle_read_argument(args, "dll", file.path(root, "current-lib/alfak2/libs/alfak2.so"))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  requireNamespace("Rcpp"); requireNamespace("TMB")
  binary <- dyn.load(normalizePath(dll)); package <- new.env(parent = globalenv())
  for (name in names(getDLLRegisteredRoutines(binary)$`.Call`))
    assign(name, getDLLRegisteredRoutines(binary)$`.Call`[[name]], package)
  sources <- sort(list.files(file.path(frozen, "R"), "[.]R$", full.names = TRUE))
  for (file in sources) sys.source(file, package)
  sys.source(file.path(repo, "R/local-fit.R"), package)
  process <- new.env(parent = globalenv())
  process_file <- file.path(dirname(script), "process_prepare_data.R")
  sys.source(process_file, process)
  manifest <- utils::read.csv(file.path(root, "process-data/dataset_manifest.csv"), stringsAsFactors = FALSE)
  manifest <- manifest[manifest$mode == "finite_population" & manifest$sampling_replicate == 1L, ]
  manifest <- manifest[!duplicated(manifest$family), ]
  summaries <- spectra <- stability <- directions <- list()
  for (case in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[case])
    fit <- package$fit_alfak2(dataset$counts, dt = dataset$config$dt, beta = dataset$config$beta,
      min_cn = 1, max_cn = 3, max_nodes = 27, target_karyotypes = dataset$target_karyotypes,
      pipeline = "joint_count", forward_model = "division_process", sdreport_mode = "fitness_only",
      joint_control = list(baseline_sd = 0.5, main_sd = 0.07, curvature_sd = 0.04,
        interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25, collapse_gaussian = FALSE))
    if (!identical(fit$local$diagnostics$covariance_status, "TMB_sdreport")) stop("Untrusted Fisher expansion point")
    graph <- fit$global$graph
    prior <- package$build_joint_landscape_prior(graph, baseline_sd = 0.5, main_sd = 0.07,
      curvature_sd = 0.04, interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25)
    covariance <- solve(as.matrix(prior$precision)) + tcrossprod(sweep(prior$design, 2L, prior$coef_sd, "*"))
    pi0 <- fit$local$posterior_predictive$pi0
    log_birth <- as.numeric(fit$local$latent_fitness_mean)
    birth <- exp(log_birth)
    probability <- process$process_independent_operator(graph$karyotypes, birth, dataset$config$beta)$probability
    depth <- colSums(dataset$counts)
    same <- oracle_fisher_jacobian(pi0, log_birth, probability, dataset$config$dt, epsilon = 1e-5)
    half <- oracle_fisher_jacobian(pi0, log_birth, probability, dataset$config$dt, epsilon = 5e-6)
    design <- oracle_fisher_design(same, depth)
    fine_design <- oracle_fisher_design(half, depth)
    eta <- as.numeric(fit$local$parameter_mode$eta)
    softmax <- function(x) {value <- exp(x - max(x)); value / sum(value)}
    j0_fd <- vapply(seq_along(pi0), function(j) {
      delta <- numeric(length(pi0)); delta[j] <- 5e-6
      (softmax(eta + delta) - softmax(eta - delta)) / 1e-5
    }, numeric(length(pi0)))
    j1_eta_fd <- vapply(seq_along(pi0), function(j) {
      delta <- numeric(length(pi0)); delta[j] <- 5e-6
      plus <- same$prediction$propagator %*% softmax(eta + delta)
      minus <- same$prediction$propagator %*% softmax(eta - delta)
      as.numeric(plus / sum(plus) - minus / sum(minus)) / 1e-5
    }, numeric(length(pi0)))
    stability[[case]] <- data.frame(family = dataset$family, dataset_id = dataset$dataset_id,
      pi0_eta_jacobian_max_error = max(abs(j0_fd - half$j0_eta)),
      pi1_eta_jacobian_max_error = max(abs(j1_eta_fd - half$j1_eta)),
      logbirth_jacobian_epsilon_halving_max_change = max(abs(same$j1_birth - half$j1_birth)),
      profiled_fisher_relative_change = norm(design$fisher - fine_design$fisher, "F") /
        max(1e-14, norm(fine_design$fisher, "F")),
      mode_cpp_vs_independent_exact_pi1 = max(abs(half$pi1 - fit$local$posterior_predictive$pi1)))
    truth_pi0 <- dataset$truth_pi0[match(graph$labels, dataset$target_labels)]
    oracle <- oracle_fisher_jacobian(truth_pi0, log_birth, probability, dataset$config$dt, epsilon = 5e-6)
    scenarios <- list(unknown_pi0_profile = fine_design,
      fixed_pi0_same_mode = oracle_fisher_design(half, depth, profile_eta = FALSE),
      fixed_true_pi0_same_birth = oracle_fisher_design(oracle, depth, profile_eta = FALSE))
    increment <- scenarios$fixed_pi0_same_mode$fisher - fine_design$fisher
    stability[[case]]$known_minus_profile_min_eigenvalue <- min(eigen((increment + t(increment)) / 2,
      symmetric = TRUE, only.values = TRUE)$values)
    for (scenario in names(scenarios)) for (planned_depth in c(150, 1000, 5000)) {
      decomposition <- oracle_fisher_spectrum(scenarios[[scenario]]$design, covariance, log_birth,
                                             depth_scale = planned_depth / depth[1L])
      summaries[[length(summaries) + 1L]] <- data.frame(family = dataset$family,
        dataset_id = dataset$dataset_id, scenario = scenario, planned_depth = planned_depth,
        n_nodes = length(pi0), eta_nuisance_rank = scenarios[[scenario]]$eta_rank,
        data_rank = decomposition$rank, rank_cutoff = decomposition$rank_cutoff,
        effective_learned_directions = decomposition$effective_directions,
        shape_effective_directions = decomposition$shape_effective_directions,
        maximum_prior_relative_information = max(decomposition$eigenvalues),
        fixed_shape_absolute_rayleigh = decomposition$fixed_shape_absolute_rayleigh,
        common_log_birth_rayleigh = decomposition$common_log_birth_rayleigh,
        absolute_projection_on_data_null = decomposition$absolute_projection_on_data_null,
        absolute_projection_on_prior_dominated = decomposition$absolute_projection_on_prior_dominated,
        weakest_physical_absolute_alignment = decomposition$weakest_physical_absolute_alignment)
      for (space in c("full_logbirth", "shape_fixed_absolute")) {
        values <- if (space == "full_logbirth") decomposition$eigenvalues else decomposition$shape_eigenvalues
        spectra[[length(spectra) + 1L]] <- data.frame(family = dataset$family,
          dataset_id = dataset$dataset_id, scenario = scenario, planned_depth = planned_depth,
          space = space, index = seq_along(values), eigenvalue = values,
          local_gaussian_learned_fraction = values / (1 + values))
      }
      if (planned_depth == 150L) {
        for (j in seq.int(max(1, length(pi0) - 2L), length(pi0))) {
          log_direction <- as.numeric(decomposition$factor %*% decomposition$eigenvectors[, j])
          physical <- birth * log_direction
          directions[[length(directions) + 1L]] <- data.frame(family = dataset$family,
            scenario = scenario, direction = paste0("weak_full_", j),
            eigenvalue = decomposition$eigenvalues[j], karyotype = graph$labels,
            logbirth_loading = log_direction, birth_loading = physical,
            pi0_mode = pi0, pi1_mode = half$pi1, native_distance = graph$support_distance)
        }
      }
    }
    saveRDS(list(dataset_id = dataset$dataset_id, graph = graph, pi0_mode = pi0,
      logbirth_mode = log_birth, prior_covariance = covariance, daughter_probability = probability,
      jacobian = half, scenarios = scenarios, fit_diagnostics = fit$local$diagnostics),
      file.path(output, paste0("oracle_fisher_", dataset$family, ".rds")))
    cat("Fisher representative", dataset$family, "completed\n")
  }
  utils::write.csv(do.call(rbind, summaries), file.path(output, "oracle_fisher_summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, spectra), file.path(output, "oracle_fisher_spectra.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, stability), file.path(output, "oracle_fisher_numerical_checks.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, directions), file.path(output, "oracle_fisher_weak_directions.csv"), row.names = FALSE)
  saveRDS(list(source_hashes = tools::md5sum(c(sources, file.path(repo, "R/local-fit.R"),
    process_file, script, dll)), data_hashes = tools::md5sum(manifest$dataset_file),
    scope = "local_expected_multinomial_Fisher_at_three_v1_development_modes",
    limitations = c("unregularized eta profiling, not actual posterior precision",
      "prior-relative effective directions depend on the chosen prior",
      "depth scales expected Fisher only; no new data or accuracy claim",
      "fixed_true_pi0 changes the expansion point; same_mode fixes only nuisance information",
      "no population growth information or finite-process variance included")),
    file.path(output, "oracle_fisher_provenance.rds"))
}

if (sys.nframe() == 0L) oracle_profile_fisher_main()
