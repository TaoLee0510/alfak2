test_that("experimental native NN/Kriging backend returns direct, NN, and kriging nodes", {
  counts <- matrix(
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
  fit <- suppressWarnings(fit_alfak2_nn_kriging_experimental(
    counts,
    dt = 1,
    beta = 0.01,
    local_shell_depth = 0,
    global_extra_shell = 2,
    min_cn = 1,
    max_cn = 3,
    max_nodes = 500,
    graph_edge_weight = "unit",
    nn_prior = "empirical",
    nboot = 2,
    kriging_max_anchors = 20,
    control = list(eval.max = 120, iter.max = 120),
    retry_control = list(eval.max = 160, iter.max = 160)
  ))

  expect_s3_class(fit, "alfak2_native_nn_kriging_fit")
  s <- summarize_alfak2(fit)
  expect_true(all(c("direct", "nn", "kriging") %in% s$support_scope))
  expect_true(all(is.finite(s$fitness_mean)))
  expect_true(all(is.finite(s$fitness_sd)))
  expect_equal(nrow(fit$posterior_samples), 2)
  expect_true(ncol(fit$posterior_samples) == nrow(s))
  expect_true(nrow(fit$nn$diagnostics) > 0)
  expect_true(fit$kriging$diagnostics$n_anchors > 1)
  expect_equal(fit$kriging$diagnostics$engine, "cpp")
})
