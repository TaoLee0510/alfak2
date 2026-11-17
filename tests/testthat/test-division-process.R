test_that("division dynamics have the correct zero-time and mutation-free limits", {
  graph <- build_karyotype_graph(c("1", "2", "3"), beta = 0,
    shell_depth = 0, min_cn = 1, max_cn = 3)
  op <- division_transition_data(graph)
  pi0 <- c(0.2, 0.3, 0.5)
  birth <- c(0.5, 1, 1.5)
  expect_equal(division_forward_r(pi0, birth, op, 0)$pi1, pi0)
  out <- division_forward_r(pi0, birth, op, 1)
  expected <- pi0 * exp(birth)
  expect_equal(out$pi1, expected / sum(expected), tolerance = 1e-12)
  expect_equal(out$log_total_growth, log(sum(expected)), tolerance = 1e-12)
})

test_that("death and computational graph outflow are not renormalized", {
  graph <- build_karyotype_graph("1", beta = 0.2,
    shell_depth = 0, min_cn = 1, max_cn = 1)
  op <- division_transition_data(graph)
  expect_equal(op$retained, 0.8, ignore_attr = TRUE)
  expect_equal(op$biological_loss, 0.2)
  expect_equal(op$computational_outflow, 0)
  out <- division_forward_r(1, 2, op, 1)
  expect_equal(out$log_total_growth, 1.2, tolerance = 1e-12)
  graph$min_cn <- 0L; graph$max_cn <- 3L
  op <- division_transition_data(graph)
  expect_equal(op$biological_loss, 0)
  expect_equal(op$computational_outflow, 0.2)
})

test_that("full single-daughter kernel includes simultaneous chromosome changes", {
  graph <- build_karyotype_graph("2.2", beta = 0.05,
    shell_depth = 2, min_cn = 1, max_cn = 3, max_nodes = 20)
  op <- division_transition_data(graph)
  parent <- match("2.2", graph$labels) - 1L
  child <- match("3.3", graph$labels) - 1L
  edge <- which(op$from0 == parent & op$to0 == child)
  expect_length(edge, 1L)
  expect_equal(op$weight[edge], (0.05 * 0.95)^2)
  expect_equal(op$retained, op$viable, tolerance = 1e-10)
})

test_that("sparse splitting converges to an independent matrix exponential", {
  graph <- build_karyotype_graph("2.2", beta = 0.05,
    shell_depth = 2, min_cn = 1, max_cn = 3, max_nodes = 20)
  op <- division_transition_data(graph)
  n <- length(graph$labels)
  pi0 <- rep(1 / n, n); birth <- seq(0.5, 1.5, length.out = n)
  P <- matrix(0, n, n)
  P[cbind(op$from0 + 1L, op$to0 + 1L)] <- op$weight
  generator <- 2 * t(P) %*% diag(birth) - diag(birth)
  decomposition <- eigen(generator)
  expected <- Re(decomposition$vectors %*%
    (exp(decomposition$values) * solve(decomposition$vectors, pi0)))
  expected_frequency <- expected / sum(expected)
  coarse <- division_forward_r(pi0, birth, op, 1, steps = 4L)
  fine <- division_forward_r(pi0, birth, op, 1, steps = 16L)
  expect_lt(max(abs(fine$pi1 - expected_frequency)), 1e-6)
  expect_lt(max(abs(fine$pi1 - expected_frequency)),
            max(abs(coarse$pi1 - expected_frequency)) / 10)
  expect_equal(fine$log_total_growth, log(sum(expected)), tolerance = 1e-5)
})

test_that("joint landscape updates every node with a proper single prior", {
  counts <- matrix(c(80, 72, 24, 18, 18, 24), ncol = 2, byrow = TRUE,
    dimnames = list(c("2.2", "1.2", "3.2"), c("t0", "t1")))
  fit <- fit_alfak2(counts, dt = 1, beta = 0.01, pipeline = "joint_count",
    local_shell_depth = 1, global_extra_shell = 1,
    min_cn = 1, max_cn = 3, max_nodes = 20)
  expect_equal(nrow(fit$global$summary), 9L)
  expect_true(all(fit$global$summary$interval_valid))
  expect_true(all(fit$global$summary$prediction_status == "joint_posterior"))
  expect_equal(fit$global$diagnostics$prior_handling,
    "single_joint_count_likelihood_and_structured_prior")
  expect_equal(dim(fit$global$fitness_covariance), c(9L, 9L))
  prior <- build_joint_landscape_prior(fit$global$graph)
  covariance <- solve(as.matrix(prior$precision))
  total_variance <- rowSums(sweep(prior$design, 2L, prior$coef_sd, "*")^2) + diag(covariance)
  expect_equal(exp(prior$offset + total_variance / 2), rep(1, 9L), tolerance = 1e-12)
})

test_that("analytic Gaussian marginalization preserves posterior and evidence", {
  counts <- matrix(c(60, 50, 20, 30), 2L, byrow = TRUE,
    dimnames = list(c("2.2", "1.2"), c("t0", "t1")))
  args <- list(counts = counts, pipeline = "joint_count", forward_model = "division_process",
    min_cn = 1, max_cn = 3, max_nodes = 20, local_shell_depth = 1, global_extra_shell = 1)
  collapsed <- do.call(fit_alfak2, c(args, list(joint_control = list(collapse_gaussian = TRUE))))
  expanded <- do.call(fit_alfak2, c(args, list(joint_control = list(collapse_gaussian = FALSE))))
  expect_equal(collapsed$global$summary$fitness_mean, expanded$global$summary$fitness_mean, tolerance = 1e-7)
  expect_equal(collapsed$global$fitness_covariance, expanded$global$fitness_covariance, tolerance = 1e-7)
  expect_equal(collapsed$local$diagnostics$objective, expanded$local$diagnostics$objective, tolerance = 1e-7)
})

test_that("explicit high-dimensional targets preserve true support distances", {
  founder <- rep(2L, 22L)
  targets <- matrix(rep(founder, 6L), nrow = 6L, byrow = TRUE)
  for (i in 2:6) targets[i, seq_len(i - 1L)] <- 3L
  counts <- matrix(c(150L, 150L), nrow = 1L,
    dimnames = list(paste(founder, collapse = "."), c("t0", "t1")))
  fit <- fit_alfak2(counts, dt = 1, beta = 5e-5, pipeline = "joint_count",
    forward_model = "selection_transport", target_karyotypes = targets,
    min_cn = 1, max_cn = 3, max_nodes = 10)
  expect_equal(fit$global$graph$domain_mode, "explicit_targets")
  expect_equal(sort(fit$global$summary$support_distance), 0:5)
  expect_equal(nrow(fit$global$summary), 6L)
  expect_true(all(fit$global$summary$interval_valid))
  expect_error(build_karyotype_graph(c("1", "2", "3"), shell_depth = 0, max_nodes = 2), "exceed")
})

test_that("joint evidence mixtures propagate between-model variance", {
  counts <- matrix(c(40, 35, 20, 25), 2L, byrow = TRUE,
    dimnames = list(c("1", "2"), c("t0", "t1")))
  fit <- fit_alfak2(counts, pipeline = "joint_count", forward_model = "selection_transport",
    target_karyotypes = as.character(1:3), min_cn = 1, max_cn = 3, max_nodes = 5,
    joint_control = list(model_scales = c(0.05, 0.2, 0.5)))
  components <- fit$global$posterior_components
  weights <- vapply(components, `[[`, numeric(1), "weight")
  expect_equal(sum(weights), 1, tolerance = 1e-12)
  expected_mean <- Reduce(`+`, lapply(components, function(x) x$weight * x$mean))
  expected_variance <- Reduce(`+`, lapply(components, function(x)
    x$weight * (diag(x$covariance) + (x$mean - expected_mean)^2)))
  expect_equal(fit$global$summary$fitness_mean, expected_mean)
  expect_equal(diag(fit$global$fitness_covariance), expected_variance, ignore_attr = TRUE)
  expect_equal(fit$global$anchors$fitness_mean,
    fit$global$summary$fitness_mean[fit$global$summary$effective_count_total > 0])
  expect_error(fit_alfak2(counts, pipeline = "joint_count", inference_method = "joint_map",
    joint_control = list(model_scales = c(.2, .5))), "laplace")
})
