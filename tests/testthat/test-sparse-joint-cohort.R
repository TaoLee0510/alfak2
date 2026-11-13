test_that("a cohort of one sample with full sharing is the single-sample fit", {
  sim <- simulate_cohort_counts()
  single <- fit_sparse_joint(sim$samples$A, dt = sim$dt[1L], beta = sim$beta, major_frequency = 0.05, max_cn = 4, initial_state = "free")
  cohort <- fit_sparse_joint_cohort(sim$samples["A"], dt = sim$dt[1L], beta = sim$beta, sharing = "full",
                                    major_frequency = 0.05, max_cn = 4, initial_state = "free")
  expect_s3_class(cohort, "alfak2_sparse_joint_cohort_fit")
  expect_equal(cohort$diagnostics$objective, single$diagnostics$objective, tolerance = 1e-8)
  expect_equal(cohort$samples$A$karyotype, single$summary$karyotype)
  columns <- c("fitness_mean", "fitness_sd", "relative_fitness", "relative_sd")
  expect_equal(cohort$samples$A[columns], single$summary[columns], tolerance = 1e-6)
  expect_equal(cohort$landscape[columns], single$summary[columns], tolerance = 1e-6)
  expect_equal(unname(cohort$trend$scales[["deviation"]]), 0)
  targets <- c(single$summary$karyotype[2L], "4.4.3.3")
  expect_equal(predict(cohort, targets, sample = "A")[columns], predict(single, targets)[columns], tolerance = 1e-6)
})

test_that("samples share the karyotypes they have in common", {
  sim <- simulate_cohort_counts()
  fit <- fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4, initial_state = "free")
  expect_true(fit$diagnostics$finite_uncertainty)
  expect_equal(fit$diagnostics$n_samples, 2L)
  expect_equal(names(fit$samples), c("A", "B"))
  landscape <- fit$landscape
  expect_equal(nrow(landscape), length(union(fit$samples$A$karyotype, fit$samples$B$karyotype)))
  expect_equal(fit$diagnostics$n_states, nrow(fit$samples$A) + nrow(fit$samples$B))
  expect_gt(fit$diagnostics$shared_nodes, 0)
  shared <- landscape[landscape$karyotype == "2.3.2.2", ]
  expect_equal(shared$n_samples, 2L)
  expect_equal(shared$n_major, 1L)
  # Alone, sample A takes "2.3.2.2" for one more unfit neighbour of its
  # dominant clone. Sample B shows it to be the fitter of the two, and in the
  # cohort sample A's own value follows.
  gap <- function(summary) {
    summary$fitness_mean[summary$karyotype == "2.3.2.2"] - summary$fitness_mean[summary$karyotype == "2.2.2.2"]
  }
  truth <- stats::setNames(sim$net, sim$labels)
  true_gap <- truth[["2.3.2.2"]] - truth[["2.2.2.2"]]
  alone <- fit_sparse_joint(sim$samples$A, dt = sim$dt[1L], beta = sim$beta, major_frequency = 0.05, max_cn = 4, initial_state = "free")
  expect_lt(gap(alone$summary), 0)
  expect_gt(gap(fit$samples$A), gap(alone$summary) + 0.04)
  full <- fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, sharing = "full",
                                  major_frequency = 0.05, max_cn = 4, initial_state = "free")
  expect_lt(abs(gap(full$samples$A) - true_gap), 0.015)
  expect_equal(full$samples$A$fitness_mean[full$samples$A$karyotype == "2.3.2.2"],
               full$landscape$fitness_mean[full$landscape$karyotype == "2.3.2.2"])
  expect_gt(stats::cor(landscape$fitness_mean, truth[landscape$karyotype], method = "spearman"), 0.4)
  expect_gt(fit$trend$scales[["deviation"]], 0)
})

test_that("prediction takes the best-informed source", {
  sim <- simulate_cohort_counts()
  samples <- sim$samples
  samples$A <- rbind(samples$A, "4.4.4.4" = c(1L, 0L))
  fit <- fit_sparse_joint_cohort(samples, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4)
  only_b <- setdiff(fit$samples$B$karyotype, fit$samples$A$karyotype)[1L]
  targets <- c("2.2.2.2", only_b, "4.4.4.4", "4.4.3.3")
  in_a <- predict(fit, targets, sample = "A")
  expect_equal(in_a$source, c("sample", "landscape", "observed_level", "trend"))
  expect_equal(in_a$fitness_mean[1L], fit$samples$A$fitness_mean[fit$samples$A$karyotype == "2.2.2.2"])
  node <- fit$landscape[fit$landscape$karyotype == only_b, ]
  expect_equal(in_a$fitness_mean[2L], node$fitness_mean)
  expect_gt(in_a$fitness_sd[2L], node$fitness_sd)
  shared <- predict(fit, targets)
  expect_equal(shared$source, c("landscape", "landscape", "trend", "trend"))
  expect_equal(shared$fitness_sd[2L], node$fitness_sd)
  expect_equal(predict(fit, targets, sample = 1L), in_a)
  expect_true(all(is.finite(in_a$fitness_mean)) && all(in_a$fitness_sd > 0))
  expect_error(predict(fit, targets, sample = "C"), "sample")
})

test_that("invalid cohort inputs are rejected", {
  sim <- simulate_cohort_counts()
  expect_error(fit_sparse_joint_cohort(sim$samples$A, dt = sim$dt[1L]), "list")
  expect_error(fit_sparse_joint_cohort(sim$samples, beta = sim$beta), "dt")
  expect_error(fit_sparse_joint_cohort(sim$samples, dt = c(1, 2, 3), beta = sim$beta), "dt")
  expect_error(fit_sparse_joint_cohort(stats::setNames(sim$samples, c("A", "A")), dt = sim$dt), "names")
  expect_error(fit_sparse_joint_cohort(list(A = sim$samples$A, B = prepare_alfak2_data(sim$samples$B, dt = 60, beta = 0.01)),
                                       dt = 10, beta = sim$beta), "missegregation")
})

test_that("the trend is per sample, shared, or chosen by marginal likelihood", {
  sim <- simulate_cohort_counts()
  fit <- function(trend) {
    fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                            trend = trend)
  }
  sample <- fit("sample")
  shared <- fit("shared")
  expect_equal(sample$diagnostics$trend, "sample")
  expect_equal(shared$diagnostics$trend, "shared")
  # per sample: each sample's trend is centred on its own mean karyotype
  expect_false(isTRUE(all.equal(sample$trend$samples$A$centre, sample$trend$samples$B$centre)))
  expect_false(isTRUE(all.equal(sample$trend$samples$A$coefficients, sample$trend$samples$B$coefficients)))
  # shared: one set of coefficients and one centre
  expect_equal(shared$trend$samples$A$coefficients, shared$trend$samples$B$coefficients)
  expect_equal(shared$trend$samples$A$centre, shared$trend$samples$B$centre)
  # a karyotype far from both samples is predicted from that sample's trend
  far <- "1.4.1.4"
  in_a <- predict(sample, far, sample = "A")
  in_b <- predict(sample, far, sample = "B")
  expect_equal(c(in_a$source, in_b$source), c("trend", "trend"))
  expect_false(isTRUE(all.equal(in_a$fitness_mean, in_b$fitness_mean)))
  expect_true(is.finite(predict(sample, far)$relative_sd))
  # selection keeps the fit with the lower negative log marginal likelihood
  chosen <- fit("select")
  objective <- c(shared = shared$diagnostics$objective, sample = sample$diagnostics$objective)
  expect_equal(chosen$diagnostics$trend_objective, objective)
  expect_equal(chosen$diagnostics$trend, names(which.min(objective)))
  expect_equal(chosen$landscape, if (which.min(objective) == 1L) shared$landscape else sample$landscape)
})

test_that("with balance, each sample's history is chosen as for a single sample", {
  sim <- simulate_cohort_counts()
  cohort <- fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, trend = "sample",
                                    major_frequency = 0.05, max_cn = 4, initial_state = "balance")
  single <- vapply(names(sim$samples), function(s) {
    fit_sparse_joint(sim$samples[[s]], dt = sim$dt[match(s, names(sim$samples))], beta = sim$beta,
                     major_frequency = 0.05, max_cn = 4)$settings$history
  }, numeric(1))
  expect_equal(cohort$diagnostics$history, single)
  expect_equal(names(cohort$diagnostics$history_profiles), names(sim$samples))
  given <- fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, trend = "sample",
                                   major_frequency = 0.05, max_cn = 4, initial_state = "balance",
                                   history = unname(single))
  expect_equal(given$diagnostics$objective, cohort$diagnostics$objective)
  expect_null(given$diagnostics$history_profiles)
  # Free by default: no history to choose.
  expect_null(fit_sparse_joint_cohort(sim$samples, dt = sim$dt, beta = sim$beta, trend = "sample",
                                      major_frequency = 0.05, max_cn = 4)$diagnostics$history_profiles)
})
