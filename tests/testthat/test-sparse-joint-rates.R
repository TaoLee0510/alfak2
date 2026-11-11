test_that("a range of birth rates is fitted at log-spaced rates, with point predictions at the lowest", {
  sim <- simulate_sparse_joint_counts()
  fit_at <- function(point) {
    fit_sparse_joint_rates(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                           birth_rate_range = c(0.25, 4), point = point)
  }
  fit <- fit_at("centre")
  expect_s3_class(fit, "alfak2_sparse_joint_rate_fit")
  expect_equal(fit$rates, c(0.25, 1, 4))
  expect_equal(fit$central, 2L)
  expect_equal(fit$point, 2L)
  expect_equal(vapply(fit$fits, function(f) f$diagnostics$birth_rate, numeric(1)), fit$rates)
  single <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                             birth_rate = 1)
  targets <- c(single$summary$karyotype, "2.2.3.3")
  prediction <- predict(fit, targets)
  at_centre <- predict(single, targets)
  expect_equal(prediction$relative_fitness, at_centre$relative_fitness)
  expect_equal(prediction$per_division, at_centre$relative_fitness)
  # the interval covers the interval at every rate
  z <- stats::qnorm(0.975)
  expect_true(all(prediction$relative_low <= at_centre$relative_fitness - z * at_centre$relative_sd + 1e-12))
  expect_true(all(prediction$relative_high >= at_centre$relative_fitness + z * at_centre$relative_sd - 1e-12))
  expect_true(all(prediction$relative_sd >= at_centre$relative_sd - 1e-12))
  expect_true(all(prediction$per_division_low <= prediction$per_division &
                    prediction$per_division <= prediction$per_division_high))
  # by default, point predictions come from the lowest rate; the interval is the same
  lowest <- fit_at("lowest")
  expect_equal(lowest$point, 1L)
  expect_identical(fit_sparse_joint_rates(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                          max_cn = 4, birth_rate_range = c(0.25, 4))$point, 1L)
  at_lowest <- predict(lowest, targets)
  expect_equal(at_lowest$relative_fitness, predict(lowest$fits[[1L]], targets)$relative_fitness)
  expect_equal(at_lowest[c("relative_low", "relative_high")], prediction[c("relative_low", "relative_high")])
  expect_true(all(at_lowest$relative_fitness - z * at_lowest$relative_sd <= at_lowest$relative_low + 1e-12))
  expect_true(all(at_lowest$relative_fitness + z * at_lowest$relative_sd >= at_lowest$relative_high - 1e-12))
})

test_that("the assumed rate scales the deficit of neighbours and leaves their ranking", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint_rates(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                                birth_rate_range = c(0.5, 2))
  summaries <- lapply(fit$fits, function(f) f$summary)
  minor <- !summaries[[1L]]$major
  deficit <- vapply(summaries, function(s) mean(s$relative_fitness[minor]), numeric(1))
  # per unit of time the mean deficit grows with the assumed rate ...
  expect_true(all(diff(deficit) < 0))
  expect_gt(deficit[1L] / deficit[3L], 0.2)
  expect_lt(deficit[1L] / deficit[3L], 0.5)
  # ... per division it changes much less
  per_division <- deficit / fit$rates
  expect_lt(max(per_division) / min(per_division), 1.5)
  prediction <- predict(fit, summaries[[1L]]$karyotype)
  expect_gt(attr(prediction, "rank_agreement"), 0.8)
})

test_that("invalid ranges are rejected", {
  sim <- simulate_sparse_joint_counts()
  expect_error(fit_sparse_joint_rates(sim$counts, dt = sim$dt, birth_rate_range = c(1, 1)), "birth_rate_range")
  expect_error(fit_sparse_joint_rates(sim$counts, dt = sim$dt, birth_rate_range = 1), "birth_rate_range")
  expect_error(fit_sparse_joint_rates(sim$counts, dt = sim$dt, birth_rate = 1, birth_rate_range = c(0.5, 2)),
               "birth_rate_range")
})
