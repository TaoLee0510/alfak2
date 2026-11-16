test_that("one missegregation rate fits as before, with no rate profile", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4, history = 60)
  expect_null(fit$diagnostics$beta_profile)
  expect_equal(fit$data$beta, sim$beta)
})

test_that("several missegregation rates are candidates and the lowest objective is kept", {
  sim <- simulate_sparse_joint_counts()
  rates <- sim$beta * c(1 / 4, 1, 4)
  fit_at <- function(beta) fit_sparse_joint(sim$counts, dt = sim$dt, beta = beta, max_cn = 4, history = 60)
  single <- lapply(rates, fit_at)
  objective <- vapply(single, function(f) f$diagnostics$objective, numeric(1))
  fit <- fit_at(rates)
  profile <- fit$diagnostics$beta_profile
  expect_equal(profile$beta, rates)
  expect_equal(profile$objective, objective)
  best <- which.min(objective)
  expect_equal(fit$data$beta, rates[[best]])
  expect_identical(fit$summary, single[[best]]$summary)
  expect_identical(fit$diagnostics$objective, single[[best]]$diagnostics$objective)
  # The history is still chosen inside each candidate.
  expect_equal(profile$history, rep(60, 3))
})

test_that("beta = \"auto\" fits the grid around the paper's rate", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = "auto", max_cn = 4, history = 60)
  expect_equal(fit$diagnostics$beta_profile$beta, 5e-5 * c(1 / 16, 1 / 4, 1, 4, 16, 64))
  expect_equal(fit$data$beta, fit$diagnostics$beta_profile$beta[[which.min(fit$diagnostics$beta_profile$objective)]])
})

test_that("the missegregation rate of an alfak2_data object is replaced by each candidate", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_alfak2_data(sim$counts, dt = sim$dt, beta = sim$beta)
  rates <- sim$beta * c(1, 4)
  fit <- fit_sparse_joint(data, beta = rates, max_cn = 4, history = 60)
  direct <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = rates, max_cn = 4, history = 60)
  expect_equal(fit$diagnostics$beta_profile$objective, direct$diagnostics$beta_profile$objective)
  # One value keeps the object's own rate, as before.
  fixed <- fit_sparse_joint(data, beta = 1e-3, max_cn = 4, history = 60)
  expect_equal(fixed$data$beta, sim$beta)
})

test_that("the pairwise method takes the rate the anchored model chooses", {
  sim <- simulate_series_counts()
  rates <- sim$beta * c(1 / 4, 4)
  fit_with <- function(series, beta) {
    suppressWarnings(fit_sparse_joint(sim$counts, times = sim$times, beta = beta, max_cn = 4, history = 60,
                                      series = series))
  }
  anchored <- fit_with("anchored", rates)
  pairwise <- fit_with("pairwise", rates)
  profile <- pairwise$diagnostics$beta_profile
  expect_equal(profile$objective, anchored$diagnostics$beta_profile$objective)
  expect_equal(unique(profile$fitted_with), "anchored")
  expect_equal(pairwise$data$beta, anchored$data$beta)
  expect_identical(pairwise$summary, fit_with("pairwise", anchored$data$beta)$summary)
})

test_that("invalid missegregation rates are refused, and cohorts take one rate", {
  sim <- simulate_sparse_joint_counts()
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = "best"), "missegregation rates")
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = c(1e-4, -1)), "missegregation rates")
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = 2), "missegregation rates")
  cohort <- simulate_cohort_counts()
  expect_error(fit_sparse_joint_cohort(cohort$samples, dt = cohort$dt, beta = "auto"), "one missegregation rate")
})
