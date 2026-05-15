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
  expect_true(all(c("count_t0", "count_t1", "count_total", "covariance_status") %in% names(fit$summary)))
  expect_true(all(c("observation_weight_t0", "observation_weight_t1", "effective_count_total",
                    "eta_prior_mean", "eta_prior_sd", "pi0", "pi1") %in% names(fit$summary)))
  expect_true("directly_informed" %in% fit$summary$support_tier)
  expect_lt(sum(fit$summary$pi0[fit$summary$support_distance > 0]), 0.05)
  expect_lt(sum(fit$summary$pi1[fit$summary$support_distance > 0]), 0.1)
  expect_gt(fit$diagnostics$eta_prior$n_borrowed_shrunk, 0)
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
  expect_true(all(c("variance_base", "variance_multiplier", "covariance_status", "count_total") %in% names(gp$anchors)))
  expect_equal(nrow(gp$anchors), sum(local$summary$effective_count_total > 0))
  expect_gt(gp$diagnostics$anchor_count_excluded, 0)

  full <- fit_alfak2(dat, local_shell_depth = 1, global_extra_shell = 0,
                     min_cn = 1, max_cn = 3, max_nodes = 200,
                     graph_edge_weight = "unit",
                     anchor_count_reference = 20,
                     anchor_exclude = "1.2.2",
                     control = list(eval.max = 120, iter.max = 120))
  expect_s3_class(full, "alfak2_fit")
  expect_gt(nrow(summarize_alfak2(full)), 0)
  expect_equal(full$diagnostics$graph_edge_weight, "unit")
  expect_equal(full$diagnostics$transition_kernel, "exact")
  expect_equal(full$diagnostics$anchor_support_tiers, "all")
  expect_equal(full$diagnostics$anchor_exclude, "1.2.2")
  expect_false("1.2.2" %in% as.character(full$global$anchors$karyotype))
})

test_that("observation weights downweight low-support rows in the local likelihood", {
  counts <- stable_counts_input()
  weights <- matrix(1, nrow = nrow(counts), ncol = 2L, dimnames = list(rownames(counts), c("t0", "t1")))
  weights["1.2.2", ] <- 0.25
  attr(counts, "observation_weights") <- weights
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  row <- fit$summary[fit$summary$karyotype == "1.2.2", , drop = FALSE]
  expect_equal(row$observation_weight_t0, 0.25)
  expect_equal(row$effective_count_total, sum(counts["1.2.2", ]) * 0.25)
  expect_true(fit$diagnostics$use_observation_weights)
  expect_equal(fit$diagnostics$likelihood_model, "weighted_multinomial")

  gp <- fit_graph_posterior(local_fit = fit, graph = graph,
                            lambda_l_grid = 1, lambda_e_grid = 0.1,
                            sigma_obs_grid = 0.05,
                            anchor_count_reference = 20)
  anchor <- gp$anchors[gp$anchors$karyotype == "1.2.2", , drop = FALSE]
  expect_equal(anchor$anchor_count_for_weight, row$effective_count_total)
})

test_that("zero-effective-count local rows are not global anchors by default", {
  counts <- stable_counts_input()
  weights <- matrix(1, nrow = nrow(counts), ncol = 2L, dimnames = list(rownames(counts), c("t0", "t1")))
  weights["1.2.2", ] <- 0
  attr(counts, "observation_weights") <- weights
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  gp <- fit_graph_posterior(fit, graph, lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05)

  expect_equal(fit$summary$effective_count_total[fit$summary$karyotype == "1.2.2"], 0)
  expect_false("1.2.2" %in% as.character(gp$anchors$karyotype))
})

test_that("dirichlet-multinomial uses explicit weighted likelihood when observation weights are present", {
  counts <- stable_counts_input()
  weights <- matrix(1, nrow = nrow(counts), ncol = 2L, dimnames = list(rownames(counts), c("t0", "t1")))
  weights["1.2.2", ] <- 0.25
  attr(counts, "observation_weights") <- weights
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  fit <- fit_local_posterior(
    dat,
    graph,
    observation_model = "dirichlet_multinomial",
    dm_concentration = 50,
    control = list(eval.max = 120, iter.max = 120)
  )
  expect_equal(fit$diagnostics$observation_model, "dirichlet_multinomial")
  expect_equal(fit$diagnostics$likelihood_model, "weighted_dirichlet_multinomial")
})

test_that("global tuning records insufficient-anchor fallback", {
  counts <- stable_counts_input()[1:2, , drop = FALSE]
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 0, min_cn = 1, max_cn = 3, max_nodes = 20)
  local <- fit_local_posterior(dat, graph, control = list(eval.max = 120, iter.max = 120))
  gp <- fit_graph_posterior(
    local,
    graph,
    lambda_l_grid = c(0.2, 1, 5),
    lambda_e_grid = c(0.05, 0.25, 1),
    sigma_obs_grid = c(0.02, 0.05, 0.1)
  )

  expect_equal(gp$diagnostics$cv_status, "insufficient_anchors")
  expect_equal(gp$hyperparameters$lambda_l, 1)
  expect_equal(gp$hyperparameters$lambda_e, 0.25)
  expect_equal(gp$hyperparameters$sigma_obs, 0.05)
  expect_true(all(is.na(gp$tuning_grid$score)))
})

test_that("global tuning skips heldout anchors without component support", {
  graph <- build_karyotype_graph(
    c("1.1", "1.2", "1.3", "3.3"),
    beta = 0.01,
    transition_kernel = "linear",
    shell_depth = 0,
    min_cn = 1,
    max_cn = 3
  )
  local <- structure(
    list(
      summary = data.frame(
        karyotype = as.character(graph$labels),
        support_tier = as.character(graph$support_tier),
        fitness_mean = c(0, 0.1, 0.2, 1),
        fitness_sd = rep(0.1, 4),
        count_total = rep(10, 4),
        effective_count_total = rep(10, 4),
        covariance_status = "TMB_sdreport",
        stringsAsFactors = FALSE
      ),
      diagnostics = list(covariance_status = "TMB_sdreport"),
      graph = graph
    ),
    class = "alfak2_local_fit"
  )
  gp <- fit_graph_posterior(
    local,
    graph,
    lambda_l_grid = c(0.2, 1),
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05
  )

  expect_equal(gp$diagnostics$cv_status, "partial_components")
  expect_equal(gp$diagnostics$cv_evaluated, 3)
  expect_equal(gp$diagnostics$cv_skipped, 1)
  expect_true(all(is.finite(gp$tuning_grid$score)))
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
