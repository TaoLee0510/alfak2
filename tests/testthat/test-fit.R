tiny_input <- function() {
  land <- simulate_toy_landscape(n_chr = 3, min_cn = 1, max_cn = 3, seed = 5)
  sim <- simulate_sparse_counts(land, n0 = 80, n1 = 80, seed = 6)
  list(landscape = land, sim = sim)
}

test_that("local TMB fit smoke test returns finite intervals and tiers", {
  x <- tiny_input()
  dat <- prepare_alfak2_data(x$sim$counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  expect_s3_class(fit, "alfak2_local_fit")
  expect_true(all(is.finite(fit$summary$fitness_mean)))
  expect_true(all(is.finite(fit$summary$fitness_sd)))
  expect_true(all(is.finite(fit$summary$conf_low)))
  expect_true("directly_informed" %in% fit$summary$support_tier)
})

test_that("graph posterior and end-to-end fit smoke tests work", {
  x <- tiny_input()
  dat <- prepare_alfak2_data(x$sim$counts, dt = 1, beta = 0.01)
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
