global_accuracy_chain <- function() {
  list(karyotypes = matrix(as.integer(1:3), ncol = 1L),
       edge_from = c(1L, 2L, 2L, 3L), edge_to = c(2L, 1L, 3L, 2L),
       edge_weight = rep(1, 4), anchor_index = 1:3,
       anchor_mean = c(-1, 0, 1), anchor_var = rep(0.01, 3),
       lambda_l_grid = 1, lambda_e_grid = 0.1,
       sigma_obs_grid = c(0.02, 0.05, 0.1, 10), eps = 1e-5)
}

global_accuracy_prior <- function() {
  L <- matrix(c(1, -1, 0, -1, 2, -1, 0, -1, 1), 3)
  2 * L + 0.1 * tcrossprod(c(1, -2, 1)) + diag(1e-5, 3)
}

test_that("graph CV is a proper predictive score including latent variance", {
  args <- global_accuracy_chain()
  out <- do.call(alfak2_graph_posterior_cpp, args)
  K <- solve(global_accuracy_prior())
  expected <- vapply(args$sigma_obs_grid, function(sigma) {
    C <- K + diag(0.01 + sigma^2, 3)
    mean(vapply(1:3, function(hold) {
      train <- setdiff(1:3, hold)
      w <- solve(C[train, train], C[train, hold])
      mu <- sum(w * args$anchor_mean[train])
      v <- C[hold, hold] - sum(C[hold, train] * w)
      -dnorm(args$anchor_mean[hold], mu, sqrt(v), log = TRUE)
    }, numeric(1)))
  }, numeric(1))
  expect_equal(out$grid$score, expected, tolerance = 1e-7)
  expect_lt(out$sigma_obs, 10)
  expect_identical(out$cv_score_type, "gaussian_predictive_nlpd")
  expect_identical(out$cv_scope, "anchor_loo_conditional_on_local_fit")
})

test_that("correlated anchor posterior and CV match exact Gaussian conditioning", {
  args <- global_accuracy_chain()
  covariance <- matrix(c(.01, .006, .002, .006, .01, .004, .002, .004, .01), 3)
  args$anchor_covariance <- covariance
  args$compute_covariance <- TRUE
  out <- do.call(alfak2_graph_posterior_cpp, args)
  noise <- covariance + diag(out$sigma_obs^2, 3)
  expected_cov <- solve(global_accuracy_prior() + solve(noise))
  expected_mean <- as.numeric(expected_cov %*% solve(noise, args$anchor_mean))
  expect_equal(out$mean, expected_mean, tolerance = 1e-8)
  expect_equal(out$covariance, expected_cov, tolerance = 1e-8)
  expect_equal(out$sd^2, diag(expected_cov), tolerance = 1e-8)
  K <- solve(global_accuracy_prior())
  expected_scores <- vapply(args$sigma_obs_grid, function(sigma) {
    C <- K + covariance + diag(sigma^2, 3)
    mean(vapply(1:3, function(hold) {
      train <- setdiff(1:3, hold)
      w <- solve(C[train, train], C[train, hold])
      v <- C[hold, hold] - sum(C[hold, train] * w)
      -dnorm(args$anchor_mean[hold], sum(w * args$anchor_mean[train]), sqrt(v), log = TRUE)
    }, numeric(1)))
  }, numeric(1))
  expect_equal(out$grid$score, expected_scores, tolerance = 1e-7)
})

test_that("local polynomial propagates noisy anchors even with a perfect fit", {
  graph <- list(karyotypes = matrix(1:4, ncol = 1))
  anchors <- data.frame(node_id = 1:3, karyotype = as.character(1:3),
                        mean = 1:3, variance = rep(.01, 3),
                        support_tier = "directly_informed")
  fit <- fit_local_polynomial_target(NULL, 4, graph, anchors, ridge = 0)
  expect_equal(unname(fit[["mean"]]), 4, tolerance = 1e-10)
  expect_equal(unname(fit[["variance"]]), .19, tolerance = 1e-10)
  covariance <- matrix(.005, 3, 3) + diag(.005, 3)
  dimnames(covariance) <- list(anchors$karyotype, anchors$karyotype)
  attr(anchors, "fitness_covariance") <- covariance
  correlated <- fit_local_polynomial_target(NULL, 4, graph, anchors, ridge = 0)
  expect_equal(unname(correlated[["variance"]]), .10, tolerance = 1e-10)
})

global_accuracy_local <- function(graph, means, variances = rep(.01, length(means)), counts = NULL) {
  if (is.null(counts)) counts <- as.numeric(graph$support_distance == 0L)
  structure(list(graph = graph, summary = data.frame(
    karyotype = as.character(graph$labels), support_tier = as.character(graph$support_tier),
    support_distance = graph$support_distance, fitness_mean = means,
    fitness_sd = sqrt(variances), covariance_status = "TMB_sdreport",
    effective_count_total = counts, count_total = counts
  )), class = "alfak2_local_fit")
}

test_that("graph fitting uses aligned joint local covariance and returns contrasts", {
  graph <- build_karyotype_graph(as.character(1:3), shell_depth = 0, min_cn = 1, max_cn = 3)
  local <- global_accuracy_local(graph, c(-1, 0, 1))
  covariance <- matrix(.004, 3, 3) + diag(.006, 3)
  dimnames(covariance) <- list(graph$labels, graph$labels)
  local$fitness_covariance <- covariance[3:1, 3:1]
  fit <- fit_graph_posterior(local, graph, lambda_l_grid = 1, lambda_e_grid = .1,
                             sigma_obs_grid = .05, graph_edge_weight = "unit", compute_covariance = TRUE)
  expect_identical(fit$diagnostics$anchor_covariance_mode, "joint")
  expect_equal(fit$anchor_covariance, covariance)
  expect_equal(diag(fit$fitness_covariance), setNames(fit$summary$fitness_sd^2, graph$labels))
  expect_true(all(eigen(fit$fitness_covariance, symmetric = TRUE)$values > 0))
  local$fitness_covariance[1, 2] <- 100
  expect_error(fit_graph_posterior(local, graph), "symmetric")
})

test_that("truncated graph retains the requested nearfield covariance", {
  graph <- build_karyotype_graph(c("2", "3"), shell_depth = 3, min_cn = 0, max_cn = 7)
  local <- global_accuracy_local(graph, seq_along(graph$labels) / 10)
  fit <- fit_extrapolation_layer(local, graph, method = "truncated_nearfield_gmrf",
    lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05,
    graph_edge_weight = "unit", compute_covariance = TRUE)
  labels <- fit$fitness_covariance_labels
  expect_true(length(labels) > 0)
  index <- match(labels, fit$summary$karyotype)
  expect_equal(diag(fit$fitness_covariance),
    setNames(fit$summary$fitness_sd[index]^2, labels), tolerance = 1e-8)
  expect_true(all(fit$summary$support_distance[index] <= 2))
})

test_that("edge effects update finite borrowed rows with an explicit compatibility mode", {
  graph <- build_karyotype_graph(c("2", "3"), shell_depth = 4, min_cn = 0, max_cn = 8)
  means <- rep(100, length(graph$labels))
  means[match(c("2", "3"), graph$labels)] <- c(0, 1)
  local <- global_accuracy_local(graph, means)
  fit <- fit_extrapolation_layer(local, graph, method = "edge_effect_empirical_bayes",
                                graph_edge_weight = "unit", ridge_delta = 0)
  old <- fit_extrapolation_layer(local, graph, method = "edge_effect_empirical_bayes",
                                graph_edge_weight = "unit", ridge_delta = 0,
                                preserve_local_predictions = TRUE)
  target <- match("4", graph$labels)
  expect_equal(fit$summary$fitness_mean[target], 2, tolerance = 1e-8)
  expect_equal(old$summary$fitness_mean[target], 100)
  far <- graph$support_distance > 2
  expect_true(any(far))
  expect_true(all(is.na(fit$summary$fitness_mean[far])))
  expect_true(all(fit$summary$prediction_status[far] == "out_of_scope"))
})

test_that("tree diagnostics separate node holdout error from training residuals", {
  graph <- build_karyotype_graph(as.character(1:5), shell_depth = 0, min_cn = 1, max_cn = 5)
  local <- global_accuracy_local(graph, c(0, .2, 3, -.5, 1))
  fit <- fit_extrapolation_layer(local, graph, method = "delta_tree_ensemble",
                                 force_ridge = TRUE, ridge_lambda = .01)
  expect_identical(fit$diagnostics$cv_scope, "node_loo_conditional_on_local_fit")
  expect_equal(fit$diagnostics$cv_evaluated, 5)
  expect_true(is.finite(fit$diagnostics$training_delta_rmse))
  expect_gt(abs(fit$diagnostics$cv_score - fit$diagnostics$training_delta_rmse), 1e-4)
  disabled <- fit_extrapolation_layer(local, graph, method = "delta_tree_ensemble",
                                      force_ridge = TRUE, max_cv_anchors = 0)
  expect_identical(disabled$diagnostics$cv_status, "disabled")
  expect_true(is.na(disabled$diagnostics$cv_score))
})
