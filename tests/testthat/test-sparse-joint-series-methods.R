test_that("two samples give the two-sample fit with every series method", {
  sim <- simulate_series_counts()
  pair <- sim$counts[, c(1L, 4L)]
  pair <- pair[rowSums(pair) > 0, ]
  base <- fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60)
  for (method in c("pairwise", "state_space")) {
    fit <- fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60, series = method)
    expect_identical(fit$summary, base$summary)
    expect_identical(fit$diagnostics$objective, base$diagnostics$objective)
    expect_null(fit$settings$series)
  }
})

test_that("pair domains give each pair of samples its own majors", {
  counts <- rbind("2.2.2" = c(30, 10, 0), "3.2.2" = c(0, 25, 30), "1.1.1" = c(0, 0, 2))
  domains <- sparse_joint_pair_domains(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(domains$labels[seq_len(domains$n_major)], c("2.2.2", "3.2.2"))
  # Both are majors of both pairs; 1.1.1 has too few cells to be one and is
  # three steps from them, outside every domain.
  expect_equal(domains$pair_major[1:2, ], cbind(c(1, 1), c(1, 1)))
  expect_identical(domains$dropped, "1.1.1")
  expect_equal(nrow(domains$member), length(domains$labels))
  expect_true(all(domains$member[seq_len(domains$n_major), ] == 1))
  expect_equal(domains$y[match("2.2.2", domains$labels), ], c(30, 10, 0))
})

test_that("the pairwise method fits each pair as two samples with one landscape", {
  sim <- simulate_series_counts()
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, initial_state = "free",
                          series = "pairwise")
  expect_identical(fit$settings$series, "pairwise")
  expect_identical(fit$diagnostics$series, "pairwise")
  expect_equal(fit$diagnostics$times, sim$times)
  expect_equal(ncol(as.matrix(fit$posterior$zeta)), 3L)
  # The pairs as a cohort that shares every birth rate and the trend. The
  # cohort centres the trend on the mean of the pairs' centres, the pairwise
  # method on all samples' counts; the fits agree closely (and exactly with
  # the same centre, validation/sparse_joint/series_methods).
  pairs <- lapply(1:3, function(k) {
    y <- sim$counts[, c(k, k + 1L)]
    y[rowSums(y) > 0, , drop = FALSE]
  })
  names(pairs) <- paste0("p", 1:3)
  cohort <- fit_sparse_joint_cohort(pairs, dt = 20, beta = sim$beta, max_cn = 4, initial_state = "free",
                                    sharing = "full", trend = "shared")
  index <- match(fit$summary$karyotype, cohort$landscape$karyotype)
  expect_false(anyNA(index))
  expect_equal(fit$diagnostics$objective, cohort$diagnostics$objective, tolerance = 1e-4)
  expect_lt(max(abs(fit$summary$fitness_mean - cohort$landscape$fitness_mean[index])), 1e-3)
  # Balance at the start of each pair over the history alone (internal
  # option of the evaluation) rather than plus the population's age.
  old <- options(alfak2.pair_history = "same")
  same <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                           series = "pairwise")
  options(old)
  aged <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                           series = "pairwise")
  expect_identical(same$diagnostics$pair_history, "same")
  expect_identical(aged$diagnostics$pair_history, "age")
  expect_equal(sparse_joint_problem(same$data, same$settings)$tmb_data$pair_history, rep(60, 3))
  expect_equal(sparse_joint_problem(aged$data, aged$settings)$tmb_data$pair_history, 60 + c(0, 20, 40))
  # A rebuilt posterior reproduces the fit's own.
  rebuilt <- sparse_joint_rebuild(fit)
  expect_equal(rebuilt$solver$mode, unname(unlist(lapply(rebuilt$problem$random, function(name) {
    value <- fit$posterior[[name]]
    if (name == "zeta") value[rebuilt$problem$tmb_data$member > 0] else value
  }))), tolerance = 1e-6)
})

test_that("the state-space method follows the anchored fit when there is no drift", {
  sim <- simulate_series_counts()
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                          series = "state_space")
  expect_identical(fit$diagnostics$series, "state_space")
  expect_true(fit$diagnostics$drift > 1e-9 && fit$diagnostics$drift < 1)
  expect_identical(names(fit$trend$scales), c("dosage", "curvature", "interaction", "residual", "drift"))
  expect_equal(dim(fit$posterior$zeta_later), c(nrow(fit$summary), 3L))
  rebuilt <- sparse_joint_rebuild(fit)
  expect_equal(rebuilt$solver$mode, unname(unlist(fit$posterior[rebuilt$problem$random])), tolerance = 1e-6)
  # Nearly no drift: the closed form from sample to sample, as the anchored
  # model with every sample for every karyotype, up to the second order in
  # the mutational flux.
  old <- options(alfak2.state_space_drift = 1e-8)
  on.exit(options(old), add = TRUE)
  tight <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            series = "state_space")
  expect_equal(tight$diagnostics$drift, 1e-8)
  anchored <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                               intermediate = "all")
  index <- match(anchored$summary$karyotype, tight$summary$karyotype)
  expect_false(anyNA(index))
  expect_gt(stats::cor(anchored$summary$relative_fitness, tight$summary$relative_fitness[index]), 0.99)
  major <- anchored$summary$major
  expect_lt(max(abs(anchored$summary$relative_fitness[major] - tight$summary$relative_fitness[index][major])), 1e-3)
})

test_that("predictions and alfakR outputs work for every series method", {
  sim <- simulate_series_counts()
  for (method in c("pairwise", "state_space")) {
    fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            series = method)
    prediction <- predict(fit, c(fit$summary$karyotype[1:2], "1.1.1.1"))
    expect_true(all(is.finite(prediction$relative_fitness)))
    expect_true(all(prediction$relative_sd > 0))
    out <- as_alfakR(fit, minobs = 20, nboot = 5, n0 = 1e5, nb = 1e7)
    expect_true(all(out$landscape$k[out$landscape$fq] %in% rownames(sim$counts)))
    expect_true(is.finite(out$xval$R2R))
  }
})

test_that("a pair with a sample outside its own domain stops the fit", {
  counts <- rbind("2.2.2" = c(0, 20, 15), "3.3.2" = c(1, 0, 0), "1.1.2" = c(1, 0, 0), "2.4.2" = c(1, 0, 0),
                  "4.2.3" = c(1, 0, 0))
  expect_error(fit_sparse_joint(counts, dt = c(10, 10), beta = 0.001, series = "pairwise"),
               class = "alfak2_domain_coverage_error")
})
