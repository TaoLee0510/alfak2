#!/usr/bin/env Rscript
repo <- normalizePath(".")
source("validation/posterior_reference_4to8/reference.R")
source("validation/accuracy_round2/oracle_posterior_reference.R")
runtime <- reference_load(repo)
out <- file.path(Sys.getenv("ALFAK2_REFERENCE_OUTPUT", "validation/posterior_reference_4to8/results"), "validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(testthat)
test_results <- list()
reporter <- ListReporter$new()
with_reporter(reporter, {
  test_that("independent target agrees with TMB and its Gaussian prior is unchanged", {
    for (n in c(4L, 8L)) for (beta in c(0, .005)) {
      case <- reference_case(n, "sparse", beta, TRUE, runtime)
      fit <- reference_fit(case, runtime)
      expect_identical(fit$diagnostics$covariance_status, "TMB_sdreport")
      collapsed <- reference_fit(case, runtime, collapsed = TRUE)
      expect_equal(fit$latent_fitness_mean, collapsed$latent_fitness_mean, tolerance = 1e-6)
      expect_equal(fit$latent_fitness_covariance, collapsed$latent_fitness_covariance, tolerance = 1e-6)
      expect_equal(crossprod(case$specification$eta_basis) / 4, diag(n - 1L), tolerance = 1e-12)
      expect_equal(unname(tcrossprod(case$specification$eta_basis) / 4),
        diag(n) - matrix(1 / n, n, n), tolerance = 1e-12)
      expect_equal(diag(case$covariance), -2 * case$prior$offset, tolerance = 1e-12)
      native <- runtime$pkg$division_transition_data(case$graph)
      P <- matrix(0, n, n); P[cbind(native$from0 + 1L, native$to0 + 1L)] <- native$weight
      expect_equal(case$specification$probability, P, tolerance = 1e-12)
      obj <- fit$optimizer$obj
      for (point in 1:5) {
        set.seed(710000 + 100 * n + 10 * point)
        z <- rnorm(case$specification$q, sd = .75)
        g <- as.numeric(case$specification$offset + case$specification$lower %*% z[1:n])
        eta <- as.numeric(case$specification$eta_basis %*% z[n + 1:(n - 1)])
        full <- obj$env$last.par.best
        full[names(full) == "f"] <- g - case$prior$offset
        full[names(full) == "eta"] <- eta
        full[names(full) == "structured_coefficients"] <- 0
        report <- obj$report(full)
        expected_ll <- sum(case$counts[, 1] * log(report$pi0 + 1e-14)) +
          sum(case$counts[, 2] * log(report$pi1 + 1e-14)) +
          dnorm(report$log_total_growth, 1.5, .1, log = TRUE)
        independent <- runtime$cpp$ref_likelihood_cpp(z, case$specification)
        expect_equal(independent$pi1, as.numeric(report$pi1), tolerance = 1e-11)
        expect_equal(independent$log_growth, as.numeric(report$log_total_growth), tolerance = 1e-11)
        expect_equal(independent$log_likelihood, as.numeric(expected_ll), tolerance = 1e-9)
      }
    }
  })
  test_that("elliptical slice sampling recovers analytic Gaussian-prior moments", {
    case <- reference_case(4L, "balanced", 0, FALSE, runtime)
    case$specification$y0[] <- case$specification$y1[] <- 0
    chains <- reference_chains(case, runtime, retained = 20000L, warmup = 1000L, seed = 710300L)
    d <- reference_diagnostics(chains, case)
    expected <- exp(case$prior$offset + .5 * diag(case$covariance))
    expected_cov <- tcrossprod(expected) * expm1(case$covariance)
    m <- d$summary[startsWith(d$summary$variable, "birth_"), ]
    expect_true(all(abs(m$mean - expected) < 6 * m$mcse_mean))
    expect_true(all(abs(d$covariance$birth - expected_cov) < 6 * d$covariance_mcse$birth))
    expect_lt(d$convergence$max_rhat, 1.01)
    test_results$prior <<- list(summary = d$summary, convergence = d$convergence)
  })
  test_that("MCMC agrees with independent two-node quadrature", {
    old <- readRDS(Sys.getenv("ALFAK2_REFERENCE_TWONODE",
      "validation/posterior_reference_4to8/inputs/two_node_reference.rds"))
    for (growth in c(FALSE, TRUE)) {
      arm <- if (growth) "population_growth" else "count_only"
      sp <- list(n = 2L, q = 3L, steps = 15L, dt = old$dt,
        offset = old$prior$offset, lower = t(chol(old$aggregate_prior_covariance)),
        eta_basis = matrix(c(1, -1) * old$eta_prior_sd / sqrt(2), 2L, 1L),
        known_pi0 = numeric(), probability = diag(2), beta_zero = TRUE,
        y0 = old$counts[, 1], y1 = old$counts[, 2],
        growth = log(.5 * exp(old$dt * .7) + .5 * exp(old$dt * 1.3)),
        growth_sd = if (growth) .1 else 0)
      case <- list(id = paste0("quadrature_", arm), n = 2L, specification = sp,
        weight = c(.5, .5), edges = matrix(c(1L, 2L), 1L, 2L))
      chains <- reference_chains(case, runtime, retained = 30000L, seed = 710400L + 1000L * growth)
      d <- reference_diagnostics(chains, case)
      exact <- tail(old$results[[arm]]$reference, 1)[[1]]
      m <- d$summary[d$summary$variable %in% c("birth_1", "birth_2", "relative_1"), ]
      expected_mean <- c(exact$mean, exact$contrast_mean)
      expected_sd <- c(sqrt(diag(exact$covariance)), exact$contrast_sd)
      expect_true(all(abs(m$mean - expected_mean) <= 6 * m$mcse_mean + 3e-5))
      expect_true(all(abs(m$sd - expected_sd) <= 6 * m$mcse_sd + 3e-5))
      expect_true(all(abs(d$covariance$birth - exact$covariance) <= 6 * d$covariance_mcse$birth + 3e-5))
      expect_true(d$convergence$passed)
      test_results[[arm]] <<- list(summary = d$summary, convergence = d$convergence)
    }
  })
})
saveRDS(test_results, file.path(out, "sampler_reference_validation.rds"))
results <- as.data.frame(reporter$get_results())
results$result <- NULL
utils::write.csv(results, file.path(out, "validation_tests.csv"), row.names = FALSE)
print(colSums(results[, c("nb", "failed", "error", "warning", "passed")]))
if (any(results$failed > 0 | results$error | results$warning > 0)) stop("Reference validation did not pass")
cat("Independent target, gauge, Gaussian-prior and quadrature checks passed.\n")
