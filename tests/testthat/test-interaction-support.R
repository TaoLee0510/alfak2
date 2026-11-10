interaction_prior_covariance <- function(prior) {
  tcrossprod(sweep(prior$design, 2L, prior$coef_sd, "*")) + solve(as.matrix(prior$precision))
}

test_that("support rotation with equal scales preserves the entire prior", {
  g <- build_karyotype_graph(c("2.2.2", "3.2.2", "2.3.2"), beta = .005,
    shell_depth = 3, min_cn = 1, max_cn = 3, max_nodes = 30)
  original <- build_joint_landscape_prior(g)
  rotated <- build_joint_landscape_prior(g, interaction_support = "observed_nullspace",
    unidentified_interaction_sd = .04)
  expect_equal(interaction_prior_covariance(original), interaction_prior_covariance(rotated), tolerance = 1e-12)
  expect_equal(original$offset, rotated$offset, tolerance = 1e-12)
})

test_that("unidentified prior variance leaves observed states untouched", {
  g <- build_karyotype_graph("2.2.2", beta = .005,
    shell_depth = 3, min_cn = 1, max_cn = 3, max_nodes = 30)
  a <- build_joint_landscape_prior(g, interaction_support = "observed_nullspace", unidentified_interaction_sd = .04)
  b <- build_joint_landscape_prior(g, interaction_support = "observed_nullspace", unidentified_interaction_sd = .2)
  delta <- interaction_prior_covariance(b) - interaction_prior_covariance(a)
  observed <- which(g$support_distance == 0)
  expect_equal(delta[observed, , drop = FALSE], matrix(0, length(observed), nrow(delta)), tolerance = 1e-12)
  expect_gt(sum(diag(delta)), 0)
  expect_gt(min(eigen(delta, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
  expect_equal(a$offset[observed], b$offset[observed], tolerance = 1e-12)
  expect_equal(b$diagnostics$interaction$n_unidentified, 3L)
  # In exact Gaussian anchor conditioning, the offset change is essential:
  # latent means change, but lognormal conditional birth means stay the same.
  conditional_birth <- function(prior) {
    K <- interaction_prior_covariance(prior)
    cross <- K[, observed, drop = FALSE] %*% solve(K[observed, observed, drop = FALSE])
    mean <- prior$offset + as.numeric(cross %*% (log(1.2) - prior$offset[observed]))
    variance <- diag(K - cross %*% K[observed, , drop = FALSE])
    exp(mean + .5 * variance)
  }
  expect_equal(conditional_birth(a), conditional_birth(b), tolerance = 1e-12)
})

test_that("full observed interaction rank leaves no artificial null directions", {
  labels <- apply(expand.grid(1:3, 1:3), 1L, paste, collapse = ".")
  g <- build_karyotype_graph(labels, shell_depth = 0, min_cn = 1, max_cn = 3)
  prior <- build_joint_landscape_prior(g, interaction_support = "observed_nullspace")
  expect_equal(prior$diagnostics$interaction$rank, 1L)
  expect_equal(prior$diagnostics$interaction$n_unidentified, 0L)
})

test_that("separately controlled prior components are a genuine joint mixture", {
  y <- matrix(c(40, 35, 20, 25), 2L, byrow = TRUE,
    dimnames = list(c("2.2", "3.2"), c("t0", "t1")))
  control <- list(interaction_support = "observed_nullspace", unidentified_interaction_sd = .2,
    prior_components = list(weak = list(residual_sd = .03, edge_sd = .06),
                            rough = list(residual_sd = .3, edge_sd = .6)), model_prior = c(.75, .25))
  fit <- fit_alfak2(y, pipeline = "joint_count", forward_model = "selection_transport",
    min_cn = 1, max_cn = 3, max_nodes = 20, local_shell_depth = 1, global_extra_shell = 1,
    joint_control = control)
  expect_equal(fit$global$diagnostics$component_labels, c("weak", "rough"))
  w <- fit$global$diagnostics$model_weights
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(fit$global$summary$interval_valid))
  expect_equal(diag(fit$global$fitness_covariance),
    setNames(fit$global$summary$fitness_sd^2, fit$global$summary$karyotype), tolerance = 1e-8)
})
