test_that("calibration gate reports failure and valid-shape states", {
  calibration_script <- file.path(getwd(), "benchmark", "scr", "run_grf_alfak2_parameter_calibration.R")
  if (!file.exists(calibration_script)) {
    calibration_script <- file.path(getwd(), "..", "..", "benchmark", "scr", "run_grf_alfak2_parameter_calibration.R")
  }
  source(normalizePath(calibration_script, mustWork = TRUE))
  ranked_fail <- data.frame(
    recommended_status = c("amplitude_collapse", "wrong_direction", "delta_untrusted"),
    stringsAsFactors = FALSE
  )
  fail <- calibration_gate_summary(ranked_fail)
  expect_equal(fail$recommended_status, "no_valid_shape_configuration")
  expect_equal(fail$n_valid_shape_configs, 0L)

  ranked_valid <- data.frame(
    recommended_status = c("wrong_direction", "valid_shape_config"),
    stringsAsFactors = FALSE
  )
  valid <- calibration_gate_summary(ranked_valid)
  expect_equal(valid$recommended_status, "valid_shape_config")
  expect_equal(valid$n_valid_shape_configs, 1L)
})

test_that("global compute_sd false preserves posterior means", {
  counts <- matrix(
    c(80, 72, 24, 18, 18, 24, 20, 22, 18, 17, 22, 21, 18, 20),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(
      c("2.2.2", "2.2.3", "2.3.2", "3.2.2", "2.2.1", "2.1.2", "1.2.2"),
      c("0", "90")
    )
  )
  dat <- prepare_alfak2_data(counts, dt = 90)
  graph <- build_karyotype_graph(dat, shell_depth = 1, max_nodes = 1000)
  local <- fit_local_posterior(
    dat,
    graph,
    control = list(eval.max = 80, iter.max = 80),
    retry_on_untrusted_covariance = FALSE
  )
  gp_sd <- fit_graph_posterior(local, graph, lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05, compute_sd = TRUE)
  gp_mean <- fit_graph_posterior(local, graph, lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05, compute_sd = FALSE)
  expect_lt(max(abs(gp_sd$summary$fitness_mean - gp_mean$summary$fitness_mean)), 1e-10)
  expect_true(all(is.na(gp_mean$summary$fitness_sd)))
  expect_true(all(is.na(gp_mean$summary$conf_low)))
  expect_false(gp_mean$diagnostics$compute_sd)
})

test_that("local optimizer diagnostics are optional and include f block", {
  counts <- matrix(
    c(80, 72, 24, 18, 18, 24, 20, 22, 18, 17, 22, 21, 18, 20),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(
      c("2.2.2", "2.2.3", "2.3.2", "3.2.2", "2.2.1", "2.1.2", "1.2.2"),
      c("0", "90")
    )
  )
  dat <- prepare_alfak2_data(counts, dt = 90)
  graph <- build_karyotype_graph(dat, shell_depth = 1, max_nodes = 1000)
  plain <- fit_local_posterior(
    dat,
    graph,
    control = list(eval.max = 60, iter.max = 60),
    retry_on_untrusted_covariance = FALSE
  )
  expect_null(plain$diagnostics$optimizer)

  diag <- fit_local_posterior(
    dat,
    graph,
    control = list(eval.max = 60, iter.max = 60),
    retry_on_untrusted_covariance = FALSE,
    return_optimizer_diagnostics = TRUE,
    local_parameterization = "g_equivalent",
    local_centering = "direct_weighted_mean",
    local_centering_weight = 10,
    fixed_sigma_neighbor = 0.1
  )
  expect_s3_class(diag, "alfak2_local_fit")
  expect_true(all(c("grad_f_max_abs", "max_gradient_block_name") %in% names(diag$diagnostics$optimizer$gradient_block_summary)))
  expect_true(is.finite(diag$diagnostics$optimizer$gradient_block_summary$grad_f_max_abs))
})
