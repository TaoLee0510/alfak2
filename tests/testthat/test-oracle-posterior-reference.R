test_that("posterior reference integration recovers analytic lognormal prior moments", {
  # The reference integrator lives with the validation scripts, which are not
  # shipped in the built package (see .Rbuildignore).
  reference_script <- testthat::test_path("..", "..", "validation", "accuracy_round2",
                                          "oracle_posterior_reference.R")
  skip_if_not(file.exists(reference_script), "validation scripts are not available")
  reference <- new.env(parent = globalenv())
  sys.source(reference_script, reference)
  mu <- c(0.1, -0.2)
  covariance <- matrix(c(0.3, 0.1, 0.1, 0.2), 2L)
  result <- reference$oracle_division_reference(mu, covariance,
    counts = matrix(0, 2L, 2L), dt = 1, order = 24L, known_pi0 = c(0.5, 0.5))
  exact_mean <- exp(mu + 0.5 * diag(covariance))
  exact_covariance <- tcrossprod(exact_mean) * expm1(covariance)
  expect_equal(result$mean, exact_mean, tolerance = 1e-10)
  expect_equal(result$covariance, exact_covariance, tolerance = 1e-9)
  expect_equal(result$contrast_mean, 0.5 * diff(rev(exact_mean)), tolerance = 1e-10)
  expect_equal(result$contrast_sd,
    sqrt(as.numeric(crossprod(c(0.5, -0.5), exact_covariance %*% c(0.5, -0.5)))),
    tolerance = 1e-9)
})
