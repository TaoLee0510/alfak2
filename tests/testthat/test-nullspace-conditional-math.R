nullspace_math_covariance <- function(prior) {
  solve(as.matrix(prior$precision)) +
    tcrossprod(sweep(prior$design, 2L, prior$coef_sd, "*"))
}

nullspace_math_graph <- function(seeds) {
  build_karyotype_graph(seeds, beta = 0, shell_depth = 3,
    min_cn = 1, max_cn = 3, max_nodes = 27)
}

nullspace_math_prior <- function(graph, hidden_sd = 0.2) {
  build_joint_landscape_prior(graph, interaction_support = "observed_nullspace",
    unidentified_interaction_sd = hidden_sd)
}

test_that("partial interaction nullspaces are invariant to support and coordinate order", {
  graph <- nullspace_math_graph(c("2.2.2", "3.3.2"))
  prior <- nullspace_math_prior(graph)
  covariance <- nullspace_math_covariance(prior)
  expect_equal(prior$diagnostics$interaction$rank, 1L)
  expect_equal(prior$diagnostics$interaction$n_unidentified, 2L)

  # Reversing seeds changes graph traversal order, but cannot change its prior.
  reordered <- nullspace_math_graph(c("3.3.2", "2.2.2"))
  other <- nullspace_math_prior(reordered)
  index <- match(graph$labels, reordered$labels)
  expect_equal(nullspace_math_covariance(other)[index, index], covariance,
               ignore_attr = TRUE, tolerance = 1e-11)
  expect_equal(other$offset[index], prior$offset, tolerance = 1e-12)

  # 1.1.2 and 3.3.2 have the same pair-feature row, so adding that observed
  # pattern must not create an extra identifiable interaction direction.
  duplicate_pattern <- nullspace_math_graph(c("2.2.2", "3.3.2", "1.1.2"))
  other <- nullspace_math_prior(duplicate_pattern)
  index <- match(graph$labels, duplicate_pattern$labels)
  expect_equal(other$diagnostics$interaction$rank, 1L)
  expect_equal(nullspace_math_covariance(other)[index, index], covariance,
               ignore_attr = TRUE, tolerance = 1e-11)

  permutation <- c(3L, 1L, 2L)
  permuted_seeds <- apply(graph$karyotypes[graph$support_distance == 0L, permutation, drop = FALSE],
                          1L, paste, collapse = ".")
  permuted <- nullspace_math_graph(permuted_seeds)
  other <- nullspace_math_prior(permuted)
  permuted_labels <- apply(graph$karyotypes[, permutation, drop = FALSE], 1L, paste, collapse = ".")
  index <- match(permuted_labels, permuted$labels)
  expect_equal(nullspace_math_covariance(other)[index, index], covariance,
               ignore_attr = TRUE, tolerance = 1e-11)
  expect_equal(other$offset[index], prior$offset, tolerance = 1e-12)
})

test_that("hidden log variance and offsets have the correct conditional covariance law", {
  graph <- nullspace_math_graph(c("2.2.2", "3.3.2"))
  small <- nullspace_math_prior(graph, 0.04)
  wide <- nullspace_math_prior(graph, 0.2)
  observed <- which(graph$support_distance == 0L)
  hidden <- which(graph$support_distance > 0L)
  conditioning_value <- c(-0.2, 0.1)
  conditional <- function(prior) {
    covariance <- nullspace_math_covariance(prior)
    regression <- covariance[hidden, observed, drop = FALSE] %*%
      solve(covariance[observed, observed, drop = FALSE])
    mean <- as.numeric(prior$offset[hidden] + regression %*%
      (conditioning_value - prior$offset[observed]))
    variance <- covariance[hidden, hidden, drop = FALSE] - regression %*%
      covariance[observed, hidden, drop = FALSE]
    rate_mean <- exp(mean + 0.5 * diag(variance))
    list(mean = mean, covariance = variance, rate_mean = rate_mean,
      rate_covariance = tcrossprod(rate_mean) * expm1(variance))
  }
  a <- conditional(small); b <- conditional(wide)
  increment <- nullspace_math_covariance(wide) - nullspace_math_covariance(small)
  expect_equal(b$mean - a$mean, -0.5 * diag(increment)[hidden], tolerance = 1e-11)
  expect_gt(max(abs(b$mean - a$mean)), 1e-4)
  expect_equal(b$covariance - a$covariance, increment[hidden, hidden, drop = FALSE], tolerance = 1e-11)
  expect_gt(min(eigen(b$covariance - a$covariance, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
  # Independent mean-one lognormal multipliers preserve the physical mean and
  # add a positive-semidefinite covariance, despite lowering the log means.
  expect_equal(b$rate_mean, a$rate_mean, tolerance = 1e-11)
  expect_gt(min(eigen(b$rate_covariance - a$rate_covariance,
    symmetric = TRUE, only.values = TRUE)$values), -1e-10)
})

test_that("genuinely uninformative beta-zero counts preserve prior mixture evidence", {
  counts <- matrix(c(100, 100), 1L, dimnames = list("2.2", c("t0", "t1")))
  data <- prepare_alfak2_data(counts, dt = 1, beta = 0)
  graph <- build_karyotype_graph(data, shell_depth = 2, min_cn = 1, max_cn = 3, max_nodes = 9)
  pi0 <- stats::setNames(as.numeric(graph$labels == "2.2"), graph$labels)
  base <- list(interaction_support = "observed_nullspace", unidentified_interaction_sd = 0.2)
  pieces <- list(weak = list(residual_sd = 0.03, edge_sd = 0.06),
                 rough = list(residual_sd = 0.3, edge_sd = 0.6))
  control <- c(base, list(prior_components = pieces, model_prior = c(0.75, 0.25)))
  fit <- fit_joint_landscape(data, graph, joint_control = control,
    forward_model = "division_process", oracle_initial_probability = pi0)
  evidence <- fit$local$diagnostics$model_log_evidence
  expect_equal(evidence[1L], evidence[2L], ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(fit$local$diagnostics$model_weights, c(0.75, 0.25), tolerance = 1e-8)
  expect_equal(fit$global$summary$fitness_mean, rep(1, nrow(graph$karyotypes)), tolerance = 1e-8)

  prior_rate_covariance <- function(piece) {
    prior <- do.call(build_joint_landscape_prior, c(list(graph = graph), base, piece))
    expm1(nullspace_math_covariance(prior))
  }
  expected <- 0.75 * prior_rate_covariance(pieces$weak) + 0.25 * prior_rate_covariance(pieces$rough)
  expect_equal(fit$global$fitness_covariance, expected, ignore_attr = TRUE, tolerance = 1e-8)

  reversed <- c(base, list(prior_components = rev(pieces), model_prior = c(0.25, 0.75)))
  other <- fit_joint_landscape(data, graph, joint_control = reversed,
    forward_model = "division_process", oracle_initial_probability = pi0)
  expect_equal(other$global$summary$fitness_mean, fit$global$summary$fitness_mean, tolerance = 1e-9)
  expect_equal(other$global$fitness_covariance, fit$global$fitness_covariance, tolerance = 1e-9)
  expect_equal(other$global$summary$conf_low, fit$global$summary$conf_low, tolerance = 1e-8)
  expect_equal(other$global$summary$conf_high, fit$global$summary$conf_high, tolerance = 1e-8)
})
