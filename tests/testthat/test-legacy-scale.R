test_that("alfakR scale columns preserve native estimates and match finite-time growth", {
  counts <- matrix(
    c(30, 20,
      10, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  fit <- fit_alfak2(pipeline = "two_layer", 
    counts,
    dt = 2,
    beta = 0.01,
    local_shell_depth = 1,
    global_extra_shell = 0,
    min_cn = 1,
    max_cn = 3,
    max_nodes = 100,
    lambda_l_grid = 1,
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05,
    alfakR_scale = TRUE,
    n0 = 100,
    nb = 1000,
    control = list(eval.max = 120, iter.max = 120)
  )

  s <- summarize_alfak2(fit, layer = "local")
  expect_true(all(c(
    "fitness_mean_alfakR_scale",
    "fitness_sd_alfakR_scale",
    "conf_low_alfakR_scale",
    "conf_high_alfakR_scale",
    "efflux_viability"
  ) %in% names(s)))
  expect_true(all(is.finite(s$fitness_mean_alfakR_scale)))
  expect_true(all(is.finite(s$fitness_mean)))

  pi0 <- as.numeric(fit$local$posterior_predictive$pi0)
  pi0 <- pi0 / sum(pi0)
  g0 <- log(1000 / 100) / 2
  expect_equal(log(sum(pi0 * exp(2 * s$fitness_mean_alfakR_scale))) / 2, g0, tolerance = 1e-8)
  C <- attr(s, "fitness_covariance_alfakR_scale")
  expect_equal(diag(C), s$fitness_sd_alfakR_scale^2, tolerance = 1e-12)
})

test_that("alfakR scale can apply efflux viability correction", {
  counts <- matrix(
    c(30, 20,
      10, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  fit <- fit_alfak2(pipeline = "two_layer", 
    counts,
    dt = 2,
    beta = 0.01,
    local_shell_depth = 1,
    global_extra_shell = 0,
    min_cn = 1,
    max_cn = 3,
    max_nodes = 100,
    lambda_l_grid = 1,
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05,
    alfakR_scale = TRUE,
    n0 = 100,
    nb = 1000,
    correct_efflux = TRUE,
    control = list(eval.max = 120, iter.max = 120)
  )

  s <- summarize_alfak2(fit, layer = "local")
  pi0 <- as.numeric(fit$local$posterior_predictive$pi0)
  pi0 <- pi0 / sum(pi0)
  g0 <- log(1000 / 100) / 2
  expect_true(all(s$efflux_viability > 0))
  expect_equal(log(sum(pi0 * exp(2 * s$efflux_viability * s$fitness_mean_alfakR_scale))) / 2,
               g0, tolerance = 1e-8)
})

test_that("finite-time growth calibration removes the nonlinear centering bias", {
  summary <- data.frame(karyotype = c("1", "2"), fitness_mean = c(-0.5, 0.5), fitness_sd = c(0.1, 0.1))
  result <- add_alfakR_scale_to_summary(summary, matrix(1:2, ncol = 1), c(0.5, 0.5),
    dt = 1, beta = 0, n0 = 100, nb = 100 * exp(0.7))
  expect_equal(log(mean(exp(result$fitness_mean_alfakR_scale))), 0.7, tolerance = 1e-12)
  expect_true(all(result$fitness_sd_alfakR_scale > 0))
})

test_that("alfakR scale requires calibration population sizes", {
  counts <- matrix(
    c(30, 20,
      10, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  expect_error(
    fit_alfak2(pipeline = "two_layer", 
      counts,
      local_shell_depth = 1,
      global_extra_shell = 0,
      min_cn = 1,
      max_cn = 3,
      max_nodes = 100,
      alfakR_scale = TRUE
    ),
    "n0"
  )
})
