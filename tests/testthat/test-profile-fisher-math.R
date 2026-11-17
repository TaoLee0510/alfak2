profile_fisher_functions <- function() {
  # The Fisher diagnostic lives with the validation scripts, which are not
  # shipped in the built package (see .Rbuildignore).
  fisher_script <- testthat::test_path("..", "..", "validation", "accuracy_round2",
                                       "oracle_profile_fisher.R")
  skip_if_not(file.exists(fisher_script), "validation scripts are not available")
  env <- new.env(parent = globalenv())
  sys.source(fisher_script, env)
  env
}

test_that("beta-zero information has the correct additive-birth null direction", {
  f <- profile_fisher_functions()
  p0 <- c(0.4, 0.6); birth <- c(0.7, 1.3); dt <- 1.5
  jac <- f$oracle_fisher_jacobian(p0, log(birth), diag(2), dt)
  exact <- p0 * exp(dt * birth); exact <- exact / sum(exact)
  exact_jac <- dt * (diag(exact) - tcrossprod(exact)) %*% diag(birth)
  expect_equal(jac$pi1, exact, tolerance = 1e-12)
  expect_equal(jac$j1_birth, exact_jac, tolerance = 1e-8)
  expect_equal(as.numeric(jac$j0_eta %*% c(1, 1)), c(0, 0), tolerance = 1e-12)
  expect_equal(as.numeric(jac$j1_eta %*% c(1, 1)), c(0, 0), tolerance = 1e-12)
  expect_equal(as.numeric(jac$j1_birth %*% (1 / birth)), c(0, 0), tolerance = 1e-8)
  profiled <- f$oracle_fisher_design(jac, c(150, 150))
  fixed <- f$oracle_fisher_design(jac, c(150, 150), profile_eta = FALSE)
  expect_equal(profiled$eta_rank, 1L)
  expect_gt(min(eigen(profiled$fisher, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
  expect_gt(min(eigen(fixed$fisher - profiled$fisher, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
})

test_that("mutation changes the null direction without removing the composition rank bound", {
  f <- profile_fisher_functions()
  p <- matrix(c(0.9, 0.1, 0.2, 0.8), 2L, byrow = TRUE)
  birth <- c(0.7, 1.3)
  jac <- f$oracle_fisher_jacobian(c(0.4, 0.6), log(birth), p, dt = 1.5)
  information <- f$oracle_fisher_design(jac, c(150, 150), profile_eta = FALSE)
  values <- eigen(information$fisher, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(sum(values > max(values) * 1e-8), 1L)
  expect_gt(max(abs(jac$j1_birth %*% (1 / birth))), 1e-4)
  half <- f$oracle_fisher_jacobian(c(0.4, 0.6), log(birth), p, dt = 1.5, epsilon = 5e-6)
  expect_equal(half$j1_birth, jac$j1_birth, tolerance = 1e-8)
})

test_that("prior-whitened information scales linearly while learned fractions saturate", {
  f <- profile_fisher_functions()
  design <- cbind(c(1, 2, 3), c(-1, 1, 0))
  prior <- matrix(c(0.3, 0.1, 0.1, 0.2), 2L)
  base <- f$oracle_fisher_spectrum(design, prior, log(c(0.7, 1.3)))
  deep <- f$oracle_fisher_spectrum(design, prior, log(c(0.7, 1.3)), depth_scale = 1000 / 150)
  deepest <- f$oracle_fisher_spectrum(design, prior, log(c(0.7, 1.3)), depth_scale = 5000 / 150)
  expect_equal(deep$fisher, (1000 / 150) * base$fisher, tolerance = 1e-12)
  expect_equal(deep$eigenvalues, (1000 / 150) * base$eigenvalues, tolerance = 1e-12)
  expect_gt(deep$effective_directions, base$effective_directions)
  expect_gt(deepest$effective_directions, deep$effective_directions)
  expect_lt(deepest$effective_directions, 2)
})
