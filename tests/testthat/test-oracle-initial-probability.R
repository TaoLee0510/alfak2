oracle_probability_problem <- function() {
  counts <- matrix(c(100, 90, 0, 10), 2L, byrow = TRUE,
    dimnames = list(c("1", "2"), c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 1, beta = 0.05)
  graph <- build_karyotype_graph(data, shell_depth = 0, min_cn = 1, max_cn = 2)
  list(data = data, graph = graph)
}

test_that("oracle initial probabilities are fixed exactly including structural zeros", {
  problem <- oracle_probability_problem()
  probability <- c("2" = 0, "1" = 1)
  fit <- fit_local_posterior(problem$data, problem$graph,
    forward_model = "division_process", oracle_initial_probability = probability,
    fixed_sigma_anchor = 0.6, fixed_sigma_neighbor = 0.35, fixed_tau_group = 0.2,
    return_tmb_objects = TRUE)
  expected <- as.numeric(probability[fit$summary$karyotype])
  expect_identical(fit$summary$pi0, expected)
  expect_identical(fit$data$counts, problem$data$counts)
  expect_true(fit$diagnostics$use_oracle_initial_probability)
  expect_false("eta" %in% fit$diagnostics$random_effects)
  expect_true(all(fit$parameter_mode$eta == 0))
  expect_equal(fit$diagnostics$covariance_status, "TMB_sdreport")
  expect_lt(fit$diagnostics$conditional_gradient_norm, 1e-3)
})

test_that("oracle inputs cannot hide a contradictory or incomplete initial distribution", {
  problem <- oracle_probability_problem()
  fit <- function(probability) fit_local_posterior(problem$data, problem$graph,
    oracle_initial_probability = probability)
  expect_error(fit(c(1, 0)), "named probability")
  expect_error(fit(c("1" = 0.5, "2" = 0)), "summing to one")
  expect_error(fit(c("1" = 1.1, "2" = -0.1)), "non-negative")
  expect_error(fit(c("1" = 0, "2" = 1)), "contradict")
})
