local_accuracy_counts <- function() {
  matrix(c(80, 72, 24, 18, 18, 24, 20, 22, 18, 17, 22, 21, 18, 20),
         nrow = 7L, byrow = TRUE,
         dimnames = list(c("2.2.2", "1.2.2", "3.2.2", "2.1.2", "2.3.2", "2.2.1", "2.2.3"),
                         c("t0", "t1")))
}

local_accuracy_problem <- function() {
  data <- prepare_alfak2_data(local_accuracy_counts(), dt = 1, beta = 0.00005)
  graph <- build_karyotype_graph(data, shell_depth = 1, min_cn = 1, max_cn = 3,
                                 max_nodes = 50)
  list(data = data, graph = graph)
}

test_that("fitness integration coordinates preserve marginal evidence", {
  counts <- matrix(c(50, 30, 50, 70), 2L, byrow = TRUE,
    dimnames = list(c("1", "2"), c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 2, beta = 0)
  graph <- build_karyotype_graph(data, shell_depth = 0, min_cn = 1, max_cn = 2)
  args <- list(data = data, graph = graph, forward_model = "selection_transport",
    fixed_sigma_anchor = 0.6, fixed_sigma_neighbor = 0.35, fixed_tau_group = 0.2)
  fitness <- do.call(fit_local_posterior, c(args, list(local_parameterization = "f")))
  growth <- do.call(fit_local_posterior, c(args, list(local_parameterization = "g_equivalent")))
  expect_equal(fitness$summary$fitness_mean, growth$summary$fitness_mean, tolerance = 1e-7)
  expect_equal(fitness$fitness_covariance, growth$fitness_covariance, tolerance = 1e-7)
  expect_equal(fitness$diagnostics$objective, growth$diagnostics$objective, tolerance = 1e-7)
})

test_that("support-tier prior multipliers retain their named alignment", {
  graph <- list(labels = letters[1:4],
                support_tier = c("directly_informed", "local_borrowed", "weakly_supported", "graph_borrowed"))
  expect_equal(resolve_support_tier_multiplier(graph,
                c(local_borrowed = 3, directly_informed = 2)), c(2, 3, 1, 1))
  expect_error(resolve_support_tier_multiplier(graph,
                c(local_borrowed = 2, local_borrowed = 3)), "unique")
})

test_that("positive marginal errors cannot certify an indefinite Hessian", {
  attempt <- list(opt = list(convergence = 0L), gradient_norm = 1e-6,
                  frep = list(sd = c(0.1, 0.2)), pd_hessian = FALSE)
  expect_equal(assess_local_covariance(attempt, 1e-3), "untrusted_hessian")
  attempt$pd_hessian <- TRUE
  attempt$gradient_norm <- 0.1
  expect_equal(assess_local_covariance(attempt, 1e-3), "untrusted_gradient")
  attempt$gradient_norm <- 1e-6
  attempt$has_random_effects <- TRUE
  attempt$conditional_gradient_norm <- 0.1
  expect_equal(assess_local_covariance(attempt, 1e-3), "untrusted_inner_gradient")
})

test_that("Laplace inference avoids the legacy lower-bound scale collapse", {
  problem <- local_accuracy_problem()
  fit <- fit_local_posterior(problem$data, problem$graph,
                             forward_model = "selection_transport")
  direct <- fit$summary$support_distance == 0L
  expect_equal(fit$diagnostics$inference_method, "laplace")
  expect_equal(fit$diagnostics$covariance_status, "TMB_sdreport")
  expect_true(fit$diagnostics$pd_hessian)
  expect_lt(fit$diagnostics$gradient_norm, 1e-3)
  # This fixture previously returned every observed SD at 1.1667e-5, a value
  # set by the artificial scale floor and the number of observed nodes.
  expect_gt(min(fit$summary$fitness_sd[direct]), 1e-3)
  expect_true(all(fit$summary$interval_valid))
  expect_equal(dim(fit$fitness_covariance), rep(nrow(fit$summary), 2L))
  expect_equal(diag(fit$fitness_covariance), fit$summary$fitness_sd^2,
               ignore_attr = TRUE, tolerance = 1e-7)
})

test_that("fixed scales do not leave unused Hessian columns", {
  problem <- local_accuracy_problem()
  fit <- fit_local_posterior(problem$data, problem$graph,
                             forward_model = "selection_transport",
                             fixed_sigma_anchor = 0.6,
                             fixed_sigma_neighbor = 0.35,
                             fixed_tau_group = 0.2,
                             return_tmb_objects = TRUE)
  expect_length(fit$optimizer$par, 0L)
  expect_equal(fit$diagnostics$covariance_status, "TMB_sdreport")
  expect_true(fit$diagnostics$pd_hessian)
  expect_true(all(is.finite(fit$summary$fitness_sd)))
})

test_that("Laplace contrast uncertainty agrees with independent posterior quadrature", {
  counts <- matrix(c(50, 30, 50, 70), 2L, byrow = TRUE,
                   dimnames = list(c("1", "2"), c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 1, beta = 0)
  graph <- build_karyotype_graph(data, shell_depth = 0, min_cn = 1, max_cn = 2)
  fit <- fit_local_posterior(data, graph, forward_model = "selection_transport",
                             fixed_sigma_anchor = 0.6, fixed_sigma_neighbor = 0.35,
                             fixed_tau_group = 0.2)
  # With beta=0, the two-node likelihood depends only on the initial log odds
  # and the fitness contrast. Their exact prior variances are 2*5^2 and 2*.6^2.
  # Integrate this two-dimensional posterior independently of TMB.
  eta <- seq(-1.5, 1.5, by = 0.01)
  contrast <- seq(-2.5, 1, by = 0.01)
  log_weight <- outer(eta, contrast, function(e, f)
    stats::dbinom(30, 100, stats::plogis(e + f), log = TRUE))
  log_weight <- log_weight + stats::dbinom(50, 100, stats::plogis(eta), log = TRUE) +
    stats::dnorm(eta, 0, sqrt(2) * 5, log = TRUE)
  log_weight <- sweep(log_weight, 2L,
    stats::dnorm(contrast, 0, sqrt(2) * 0.6, log = TRUE), "+")
  marginal <- colSums(exp(log_weight - max(log_weight)))
  marginal <- marginal / sum(marginal)
  exact_mean <- sum(contrast * marginal)
  exact_sd <- sqrt(sum((contrast - exact_mean)^2 * marginal))
  C <- fit$fitness_covariance
  laplace_sd <- sqrt(C[1, 1] + C[2, 2] - 2 * C[1, 2])
  expect_equal(diff(rev(fit$summary$fitness_mean)), exact_mean, tolerance = 0.02)
  expect_equal(laplace_sd, exact_sd, tolerance = 0.02)
})

test_that("division-process summaries preserve positive support and joint uncertainty", {
  problem <- local_accuracy_problem()
  fit <- fit_local_posterior(problem$data, problem$graph, forward_model = "division_process")
  expect_equal(fit$diagnostics$point_estimate_type, "laplace_lognormal_mean")
  expect_equal(fit$diagnostics$fitness_interpretation, "birth_rate")
  expect_true(all(fit$summary$fitness_mean > 0))
  expect_true(all(fit$summary$conf_low > 0))
  expect_true(all(fit$summary$conf_high > fit$summary$conf_low))
  expect_equal(diag(fit$fitness_covariance), fit$summary$fitness_sd^2,
               ignore_attr = TRUE, tolerance = 1e-7)
  expect_gt(min(eigen(fit$fitness_covariance, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  expect_true(all(fit$diagnostics$division_transport$retained <= 1 + 1e-10))
  expect_equal(fit$diagnostics$forward_steps, 10L)
  fitness_only <- fit_local_posterior(problem$data, problem$graph,
    forward_model = "division_process", sdreport_mode = "fitness_only")
  expect_equal(fitness_only$summary$fitness_mean, fit$summary$fitness_mean, tolerance = 1e-7)
  expect_equal(fitness_only$fitness_covariance, fit$fitness_covariance, tolerance = 1e-7)
  expect_equal(length(fitness_only$sdreport$value), 2L * nrow(fit$summary))
})

test_that("joint structured residual prior matches the independent local Gaussian prior", {
  counts <- matrix(c(50, 30, 50, 70), 2L, byrow = TRUE,
                   dimnames = list(c("1", "2"), c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 1, beta = 0)
  graph <- build_karyotype_graph(data, shell_depth = 0, min_cn = 1, max_cn = 2)
  reference <- fit_local_posterior(data, graph, forward_model = "selection_transport",
                                   fixed_sigma_anchor = 0.6, fixed_sigma_neighbor = 0.35,
                                   fixed_tau_group = 0.2)
  structure <- list(design = matrix(numeric(), 2L, 0L), offset = numeric(2),
                    precision = Matrix::Diagonal(2L, 1 / 0.6^2),
                    logdet = 2 * log(1 / 0.6^2), coef_sd = numeric())
  joint <- fit_local_posterior(data, graph, forward_model = "selection_transport",
                               structured_prior = structure, return_tmb_objects = TRUE)
  expect_equal(joint$diagnostics$covariance_status, "TMB_sdreport")
  expect_length(joint$optimizer$par, 0L)
  expect_equal(joint$summary$fitness_mean, reference$summary$fitness_mean, tolerance = 1e-6)
  expect_equal(joint$fitness_covariance, reference$fitness_covariance, tolerance = 1e-6)
})

test_that("actual population growth identifies an otherwise unobserved birth baseline", {
  counts <- matrix(c(100, 100), 1L, dimnames = list("2", c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 2, beta = 0)
  graph <- build_karyotype_graph(data, shell_depth = 0, min_cn = 1, max_cn = 3)
  controls <- list(data = data, graph = graph, forward_model = "division_process",
                   fixed_sigma_anchor = 1, fixed_sigma_neighbor = 0.35,
                   fixed_tau_group = 0.2)
  unobserved <- do.call(fit_local_posterior, controls)
  observed <- do.call(fit_local_posterior, c(controls,
    list(population_sizes = c(1000, 1000 * exp(2 * 0.7)), population_log_growth_sd = 0.02)))
  expect_equal(unobserved$diagnostics$absolute_level, "prior_identified_without_external_growth")
  expect_equal(observed$diagnostics$absolute_level, "informed_by_population_growth")
  expect_equal(observed$summary$fitness_mean, 0.7, tolerance = 0.01)
  expect_gt(abs(unobserved$summary$fitness_mean - 0.7), 0.1)
  expect_equal(observed$posterior_predictive$log_total_growth, 1.4,
               ignore_attr = TRUE, tolerance = 0.01)
  expect_equal(observed$diagnostics$population_growth$observed_log_growth, 1.4)
  expect_error(do.call(fit_local_posterior, c(controls, list(population_sizes = c(0, 100)))), "positive finite")
  expect_error(do.call(fit_local_posterior, c(controls, list(population_sizes = 100))), "pair")
})

test_that("conditional warm starts preserve the target posterior and final precision", {
  problem <- local_accuracy_problem()
  controls <- list(data = problem$data, graph = problem$graph,
    forward_model = "division_process", fixed_sigma_anchor = 0.6,
    fixed_sigma_neighbor = 0.35, fixed_tau_group = 0.2, return_tmb_objects = TRUE)
  cold <- do.call(fit_local_posterior, c(controls, list(conditional_warmstart = "none")))
  warm <- do.call(fit_local_posterior, c(controls, list(conditional_warmstart = "always")))
  expect_true(warm$diagnostics$conditional_warmstart$accepted)
  expect_gte(warm$diagnostics$conditional_warmstart$elapsed_seconds, 0)
  expect_equal(warm$summary$fitness_mean, cold$summary$fitness_mean, tolerance = 1e-6)
  expect_equal(warm$fitness_covariance, cold$fitness_covariance, tolerance = 1e-6)
  joint_objective <- function(fit) {
    obj <- fit$optimizer$obj
    as.numeric(obj$env$f(obj$env$last.par))
  }
  expect_equal(joint_objective(warm), joint_objective(cold), tolerance = 1e-7)
  expect_lt(warm$diagnostics$conditional_gradient_norm, 1e-3)
  expect_lt(cold$diagnostics$conditional_gradient_norm, 1e-3)
  expect_error(fit_local_posterior(problem$data, problem$graph,
    conditional_warmstart = "always"), "every hierarchical scale")
})

test_that("legacy joint MAP remains explicit and untrusted intervals are withheld", {
  problem <- local_accuracy_problem()
  fit <- suppressWarnings(fit_local_posterior(problem$data, problem$graph,
                             inference_method = "joint_map",
                             forward_model = "selection_transport",
                             control = list(eval.max = 1, iter.max = 1),
                             retry_on_untrusted_covariance = FALSE))
  expect_equal(fit$diagnostics$inference_method, "joint_map")
  expect_false(any(fit$summary$interval_valid))
  expect_true(all(is.na(fit$summary$conf_low)))
  expect_true(all(is.na(fit$summary$conf_high)))
  expect_null(fit$fitness_covariance)
})
