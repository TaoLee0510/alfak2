stable_counts_input <- function() {
  matrix(
    c(80, 72,
      24, 18,
      18, 24,
      20, 22,
      18, 17,
      22, 21,
      18, 20),
    nrow = 7,
    byrow = TRUE,
    dimnames = list(
      c("2.2.2", "1.2.2", "3.2.2", "2.1.2", "2.3.2", "2.2.1", "2.2.3"),
      c("t0", "t1")
    )
  )
}

sparse_stress_input <- function() {
  land <- simulate_grf_landscape(
    n_chr = 3,
    n_centroids = 8,
    min_cn = 1,
    max_cn = 3,
    lambda = 0.8,
    scale = 0.5,
    seed = 1
  )
  sim <- simulate_sparse_counts(
    land,
    beta = 0.01,
    dt = 1,
    n0 = 200,
    n1 = 200,
    seed = 2,
    initial_population = 600,
    time_step = 0.25,
    detection_threshold = 1
  )
  list(landscape = land, sim = sim)
}

test_that("local TMB fit smoke test returns finite intervals and tiers", {
  dat <- prepare_alfak2_data(stable_counts_input(), dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  expect_s3_class(fit, "alfak2_local_fit")
  expect_true(all(is.finite(fit$summary$fitness_mean)))
  expect_true(all(is.finite(fit$summary$fitness_sd)))
  expect_true(all(is.finite(fit$summary$conf_low)))
  expect_true("directly_informed" %in% fit$summary$support_tier)
  expect_true(fit$diagnostics$covariance_status %in% c(
    "TMB_sdreport",
    "untrusted_nonconverged",
    "untrusted_gradient",
    "untrusted_sdreport_missing",
    "untrusted_sdreport_nonfinite"
  ))
})

test_that("graph posterior and end-to-end fit smoke tests work", {
  dat <- prepare_alfak2_data(stable_counts_input(), dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  local <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  gp <- fit_graph_posterior(local, graph, lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05)
  expect_true(all(is.finite(gp$summary$fitness_mean)))
  expect_true(all(is.finite(gp$summary$fitness_sd)))

  full <- fit_alfak2(dat, local_shell_depth = 1, global_extra_shell = 0,
                     min_cn = 1, max_cn = 3, max_nodes = 200,
                     control = list(eval.max = 120, iter.max = 120))
  expect_s3_class(full, "alfak2_fit")
  expect_gt(nrow(summarize_alfak2(full)), 0)
})

test_that("untrusted local covariance uses finite fallback uncertainty", {
  dat <- prepare_alfak2_data(stable_counts_input(), dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(
    dat,
    graph,
    control = list(eval.max = 1, iter.max = 1),
    retry_on_untrusted_covariance = FALSE
  )
  expect_true(fit$diagnostics$covariance_fallback)
  expect_equal(fit$diagnostics$fitness_sd_source, "fallback_prior_scale")
  expect_true(all(is.finite(fit$summary$fitness_sd)))
  expect_true(all(is.finite(fit$summary$conf_low)))
})

test_that("sparse stochastic stress input does not leak non-finite local intervals", {
  x <- sparse_stress_input()
  dat <- prepare_alfak2_data(x$sim$counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  expect_true(all(is.finite(fit$summary$fitness_mean)))
  expect_true(all(is.finite(fit$summary$fitness_sd)))
  expect_true(all(is.finite(fit$summary$conf_low)))
})
