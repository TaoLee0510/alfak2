test_that("two samples give the two-sample fit with every series method", {
  sim <- simulate_series_counts()
  pair <- sim$counts[, c(1L, 4L)]
  pair <- pair[rowSums(pair) > 0, ]
  base <- fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60)
  for (method in c("pairwise", "state_space", "auto")) {
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

test_that("a first sample outside every pair's domain stops the fit", {
  # The first sample is kept; its neighbour is left out, and the first and
  # the last sample alone still have no cell of the first in the domain.
  counts <- rbind("2.2.2" = c(0, 20, 15), "3.3.2" = c(1, 0, 0), "1.1.2" = c(1, 0, 0), "2.4.2" = c(1, 0, 0),
                  "4.2.3" = c(1, 0, 0))
  expect_identical(sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)$left_out, 2L)
  expect_error(suppressWarnings(fit_sparse_joint(counts, dt = c(10, 10), beta = 0.001, series = "pairwise")),
               class = "alfak2_domain_coverage_error")
})

# The simulated series with sample `k` replaced by `cells`.
replace_series_sample <- function(sim, k, cells) {
  counts <- sim$counts
  extra <- setdiff(names(cells), rownames(counts))
  counts <- rbind(counts, matrix(0, length(extra), ncol(counts), dimnames = list(extra, colnames(counts))))
  counts[, k] <- 0
  counts[names(cells), k] <- cells
  counts[rowSums(counts) > 0, , drop = FALSE]
}

test_that("a sample in between with no cell in its pair's domain is left out", {
  sim <- simulate_series_counts()
  # Sample t2 holds single cells far from every major: in pair (t1, t2) it
  # has no cell in the domain. It is left out, and t1 and t3 form a pair.
  far <- c("4.4.4.4" = 1, "1.1.4.4" = 1, "4.1.1.4" = 1)
  counts <- replace_series_sample(sim, 3L, far)
  pairs <- sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(pairs$left_out, 3L)
  expect_identical(pairs$uncovered, "t2")
  kept <- counts[, -3L]
  kept <- kept[rowSums(kept) > 0, ]
  for (method in c("pairwise", "state_space")) {
    expect_warning(fit <- fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                           series = method),
                   class = "alfak2_samples_left_out_warning")
    direct <- fit_sparse_joint(kept, times = sim$times[-3L], beta = sim$beta, max_cn = 4, history = 60,
                               series = method)
    expect_identical(fit$summary, direct$summary)
    expect_identical(fit$diagnostics$dropped_samples$sample, 3L)
    expect_identical(fit$diagnostics$dropped_samples$uncovered, "t2")
    expect_equal(fit$input$counts, counts)
  }
})

test_that("the automatic choice keeps the anchored fit without a sweep", {
  sim <- simulate_series_counts()
  expect_identical(sparse_joint_auto_series(0.5, NA), "anchored")
  expect_identical(sparse_joint_auto_series(0.1, 0.9), "pairwise")
  expect_identical(sparse_joint_auto_series(0.1, 0.5), "pairwise")
  expect_identical(sparse_joint_auto_series(0.1, 0.3), "state_space")
  expect_identical(sparse_joint_auto_series(0.1, NA), "pairwise")
  auto <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60, series = "auto")
  anchored <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60)
  choice <- auto$diagnostics$series_choice
  expect_identical(choice$series, "anchored")
  expect_gte(choice$ends_major, sparse_joint_auto_thresholds[["ends_major"]])
  expect_true(is.na(choice$weight))
  expect_identical(auto$summary, anchored$summary)
  expect_error(fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, series = "auto",
                                intermediate = "all"), "intermediate")
})

test_that("the automatic choice fits a sweep with the method it picks", {
  sim <- simulate_series_counts()
  # t2 has swept to a far clone: no cell on the first and last samples'
  # majors, so the anchored method leaves it out and auto looks at drift.
  counts <- replace_series_sample(sim, 3L, c("4.4.4.4" = 300, "4.4.4.3" = 20))
  fit_with <- function(series) {
    suppressWarnings(fit_sparse_joint(counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                      series = series))
  }
  auto <- fit_with("auto")
  choice <- auto$diagnostics$series_choice
  expect_identical(choice$ends_major, 0)
  # Its weight is the anchored fit's own.
  expect_identical(choice$weight, fit_with("anchored")$diagnostics$intermediate_weight)
  expect_identical(choice$series, sparse_joint_auto_series(choice$ends_major, choice$weight))
  expect_identical(auto$summary, fit_with(choice$series)$summary)
})

test_that("a last sample with no cell in its pair's domain pairs with a sample further in", {
  sim <- simulate_series_counts()
  # t2 has swept to a far clone; t3 holds two single cells of the earlier
  # majors. In (t2, t3) the last sample has no cell in the domain, so t2 is
  # left out and t3 pairs with t1, where those karyotypes are majors.
  counts <- replace_series_sample(sim, 3L, c("4.4.4.4" = 300, "4.4.4.3" = 20))
  counts <- replace_series_sample(list(counts = counts), 4L, c("2.2.2.2" = 1, "2.3.2.2" = 1))
  pairs <- sparse_joint_pair_left_out(counts, 0.01, 3, 5000, 1L, 4L)
  expect_identical(pairs$left_out, 3L)
  expect_identical(pairs$uncovered, "t3")
  # Every pair covered: nothing left out.
  expect_identical(sparse_joint_pair_left_out(sim$counts, 0.01, 3, 5000, 1L, 4L)$left_out, integer(0))
})
