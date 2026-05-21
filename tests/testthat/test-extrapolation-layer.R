second_layer_counts_input <- function() {
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

test_that("unified extrapolation method validation works", {
  dat <- prepare_alfak2_data(second_layer_counts_input(), dt = 1, beta = 0.01)
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3, max_nodes = 200)
  local <- suppressWarnings(fit_local_posterior(dat, graph, control = list(eval.max = 100, iter.max = 100)))
  expect_error(fit_extrapolation_layer(local, graph, method = "not_a_method"), "arg")
})

test_that("all extrapolation methods return the common global summary schema", {
  counts <- second_layer_counts_input()
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  local_graph <- build_karyotype_graph(dat, shell_depth = 2, min_cn = 0, max_cn = 4, max_nodes = 500)
  global_graph <- build_karyotype_graph(dat, shell_depth = 3, min_cn = 0, max_cn = 4, max_nodes = 500)
  local <- suppressWarnings(fit_local_posterior(dat, local_graph, control = list(eval.max = 120, iter.max = 120)))
  methods <- alfak2:::alfak2_extrapolation_methods()
  required <- c(
    "node_id", "karyotype", "support_tier", "support_distance",
    "fitness_mean", "fitness_sd", "conf_low", "conf_high",
    "extrapolation_method", "prediction_status"
  )
  for (method in methods) {
    fit <- fit_extrapolation_layer(
      local,
      global_graph,
      method = method,
      lambda_l_grid = 1,
      lambda_e_grid = 0.1,
      sigma_obs_grid = 0.05,
      graph_edge_weight = "unit"
    )
    expect_s3_class(fit, "alfak2_global_fit")
    expect_true(all(c("graph", "summary", "anchors", "hyperparameters", "tuning_grid", "diagnostics") %in% names(fit)))
    expect_true(all(required %in% names(fit$summary)))
    expect_true(any(is.finite(fit$summary$fitness_mean[fit$summary$support_distance <= 2])))
    expect_equal(unique(fit$summary$extrapolation_method), method)
    expect_equal(fit$diagnostics$extrapolation_method, method)
    expect_true(all(c("max_prediction_distance", "n_anchors", "n_predicted", "n_out_of_scope",
                      "convergence_status", "runtime_seconds", "warnings") %in% names(fit$diagnostics)))
  }
})

test_that("nearfield methods do not report reliable farfield means", {
  counts <- second_layer_counts_input()
  dat <- prepare_alfak2_data(counts, dt = 1, beta = 0.01)
  local_graph <- build_karyotype_graph(dat, shell_depth = 2, min_cn = 0, max_cn = 4, max_nodes = 500)
  global_graph <- build_karyotype_graph(dat, shell_depth = 3, min_cn = 0, max_cn = 4, max_nodes = 500)
  local <- suppressWarnings(fit_local_posterior(dat, local_graph, control = list(eval.max = 120, iter.max = 120)))
  methods <- setdiff(alfak2:::alfak2_extrapolation_methods(), "graph_gaussian_baseline")
  for (method in methods) {
    fit <- fit_extrapolation_layer(
      local,
      global_graph,
      method = method,
      lambda_l_grid = 1,
      lambda_e_grid = 0.1,
      sigma_obs_grid = 0.05,
      graph_edge_weight = "unit"
    )
    far <- fit$summary$support_distance > 2
    expect_true(any(far))
    expect_true(all(fit$summary$prediction_status[far] == "out_of_scope"))
    expect_true(all(!is.finite(fit$summary$fitness_mean[far])))
  }
})

test_that("second-layer benchmark metrics match hand-computed examples", {
  eval_nodes <- data.frame(
    karyotype = c("a", "b", "c"),
    support_distance = c(0L, 1L, 2L),
    truth = c(0, 1, 3),
    pred = c(0, 2, 4),
    pred_sd = c(0.5, 0.5, 1),
    stringsAsFactors = FALSE
  )
  eval_edges <- data.frame(
    parent_karyotype = c("a", "b"),
    child_karyotype = c("b", "c"),
    child_distance = c(1L, 2L),
    truth_gradient = c(1, 2),
    pred_gradient = c(2, 2),
    stringsAsFactors = FALSE
  )
  vals <- second_layer_metric_values(eval_nodes, eval_edges, shell = "all_nearfield")
  expect_equal(vals[["rmse"]], 1)
  expect_equal(vals[["mae"]], 1)
  expect_equal(vals[["bias"]], 1)
  expect_equal(vals[["relative_rmse"]], 1 / stats::sd(c(1, 3)))
  expect_equal(vals[["centered_rmse"]], 0)
  expect_equal(vals[["edge_gradient_rmse"]], sqrt(mean(c(1, 0)^2)))
  expect_equal(vals[["sign_accuracy"]], 1)
  expect_equal(vals[["beneficial_sign_accuracy"]], 1)
  expect_equal(vals[["deleterious_sign_accuracy"]], NA_real_)
  expect_equal(vals[["top_k_overlap_count"]], 1)
  expect_equal(vals[["coverage"]], 1)
  expect_true(is.finite(vals[["interval_coverage_95_closeness"]]))
})

test_that("full second-layer run index preserves requested combination counts", {
  idx <- second_layer_build_run_index("full")
  expect_equal(sum(idx$package == "alfak2"), 405)
  expect_equal(sum(idx$package == "alfakR"), 225)
  expect_equal(nrow(idx), 630)
  quick <- second_layer_build_run_index("quick")
  expect_equal(sum(quick$package == "alfak2"), 27)
  expect_equal(sum(quick$package == "alfakR"), 5)
  expect_equal(nrow(quick), 32)
  weighted <- idx$NN_prior_slot[idx$package == "alfakR" & grepl("weighted", idx$NN_prior_slot)]
  expect_true("empirical_censored_weighted_slot4" %in% weighted)
  expect_true("empirical_censored_weighted_slot5" %in% weighted)
})

test_that("balanced ranking and pareto helpers run on hand-built metrics", {
  env <- new.env(parent = globalenv())
  script <- "benchmark/run_full_second_layer_comparison.R"
  if (!file.exists(script)) script <- file.path("..", "..", "benchmark", "run_full_second_layer_comparison.R")
  sys.source(script, envir = env)
  metrics <- data.frame(
    run_id = rep(c("r1", "r2"), each = 17),
    package = rep("alfak2", 34),
    input_mode = rep("full", 34),
    extrapolation_method = rep(c("graph_gaussian_baseline", "edge_effect_empirical_bayes"), each = 17),
    minobs = NA_integer_,
    NN_prior_slot = NA_character_,
    grf_lambda = 0.6,
    landscape_id = "L1",
    landscape_rep = 1L,
    shell = "all_nearfield",
    prediction_scale = "raw",
    metric = rep(c(
      "rmse", "mae", "relative_rmse", "bias_abs", "q90_absolute_error", "uncalibrated_r2",
      "edge_gradient_rmse", "edge_gradient_spearman", "centered_rmse", "spearman", "sign_accuracy",
      "top_k_overlap_fraction", "interval_coverage_95_closeness", "standardized_rmse", "coverage",
      "runtime_seconds", "failure_rate"
    ), 2),
    value = c(c(2,2,2,2,2,0.2, 2,0.2,2,0.2,0.2,0.2, 0.2,2,0.8, 2,0),
              c(1,1,1,1,1,0.8, 1,0.8,1,0.8,0.8,0.8, 0.1,1,1, 1,0)),
    failure_status = "ok",
    stringsAsFactors = FALSE
  )
  ranks <- env$build_landscape_rankings(metrics)
  expect_true("balanced_weighted_rank" %in% names(ranks))
  expect_true(any(is.finite(ranks$balanced_weighted_rank)))
  pf <- env$build_pareto_front(metrics)
  expect_true("is_pareto_optimal" %in% names(pf))
})
