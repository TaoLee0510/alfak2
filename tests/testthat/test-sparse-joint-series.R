series_settings <- function(...) {
  utils::modifyList(list(major_frequency = 0.05, major_count = 3, max_states = Inf, min_cn = 1L, max_cn = 4L,
                         shell_power = 0.85, initial_frequency = c(log(1e-5), 3), birth_rate = 1,
                         birth_rate_log_sd = 0.05, shell_sd = 0.5, initial_state = "balance", history = 100,
                         balance_sd = 1), list(...))
}

test_that("count data take two or more timepoints", {
  counts <- cbind(c(5, 1, 0), c(3, 2, 1), c(1, 4, 2))
  rownames(counts) <- c("2.2", "2.3", "3.3")
  two <- prepare_alfak2_data(counts[, 1:2], dt = 10)
  expect_null(two$times)
  expect_equal(two$dt, 10)
  expect_equal(colnames(two$counts), c("t0", "t1"))
  expect_equal(alfak2_data_times(two), c(0, 10))
  three <- prepare_alfak2_data(counts, dt = c(10, 20))
  expect_equal(colnames(three$counts), c("t0", "t1", "t2"))
  expect_equal(three$times, c(0, 10, 30))
  expect_equal(three$dt, 30)
  expect_equal(alfak2_data_times(three), c(0, 10, 30))
  expect_equal(prepare_alfak2_data(counts, times = c(5, 15, 35))$times, c(0, 10, 30))
  expect_equal(prepare_alfak2_data(counts[, 1:2], times = c(3, 7))$dt, 4)
  expect_output(print(three), "times: ")
  expect_error(prepare_alfak2_data(counts, dt = 10), "intervals between consecutive")
  expect_error(prepare_alfak2_data(counts, dt = c(10, 20), times = c(0, 1, 2)), "not both")
  expect_error(prepare_alfak2_data(counts, times = c(0, 2, 1)), "strictly increasing")
  expect_error(prepare_alfak2_data(counts, times = c(0, 2)), "3 finite times")
  expect_error(prepare_alfak2_data(counts[, 1, drop = FALSE]), "at least two")
})

test_that("the series template gives the two-timepoint objective for two timepoints", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_alfak2_data(sim$counts, dt = sim$dt, beta = sim$beta)
  pair <- sparse_joint_problem(data, series_settings())
  series <- sparse_joint_problem(data, series_settings(template = "series"))
  expect_equal(pair$template, "pair")
  expect_equal(series$template, "series")
  expect_equal(series$tmb_data$tau, sim$dt)
  set.seed(3)
  n <- length(pair$domain$labels)
  at_pair <- pair$parameters
  at_pair$zeta <- at_pair$zeta + stats::rnorm(n, 0, 0.5)
  at_pair$g <- stats::rnorm(n, 0, 0.1)
  at_pair$theta1 <- 0.2
  at_series <- at_pair[c("zeta", "g", "coef")]
  at_series$theta <- 0.2
  at_series$log_scale <- at_pair$log_scale
  a <- TMB::MakeADFun(pair$tmb_data, at_pair, DLL = "alfak2", silent = TRUE)
  b <- TMB::MakeADFun(series$tmb_data, at_series, DLL = "alfak2", silent = TRUE)
  expect_equal(b$fn(), a$fn(), tolerance = 1e-12)
  expect_equal(b$gr(), a$gr(), tolerance = 1e-10)
  expect_equal(b$report()$log_intensity[, 1L], as.numeric(a$report()$log_intensity1), tolerance = 1e-12)
})

test_that("two-timepoint fits are the same with either template", {
  sim <- simulate_sparse_joint_counts()
  fit <- function() fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                                     history = 100)
  pair <- fit()
  old <- options(alfak2.sparse_joint_template = "series")
  series <- tryCatch(fit(), finally = options(old))
  expect_null(pair$settings$template)
  expect_equal(series$settings$template, "series")
  expect_equal(names(pair$posterior), c("zeta", "g", "coef", "theta1", "log_scale"))
  expect_equal(names(series$posterior), c("zeta", "g", "coef", "theta", "log_scale"))
  expect_equal(series$diagnostics$objective, pair$diagnostics$objective, tolerance = 1e-8)
  expect_equal(series$summary, pair$summary, tolerance = 1e-6)
  expect_equal(predict(series, "2.2.3.3"), predict(pair, "2.2.3.3"), tolerance = 1e-6)
  # The rebuilt posterior of a series fit is its own.
  expect_equal(sparse_joint_rebuild(series, posterior = FALSE)$g, series$posterior$g, tolerance = 1e-6)
})

test_that("compiled series model evaluates the closed form at every sample time", {
  sim <- simulate_series_counts()
  data <- prepare_alfak2_data(sim$counts, beta = sim$beta, times = sim$times)
  problem <- sparse_joint_problem(data, series_settings(initial_state = "free"))
  expect_equal(problem$template, "series")
  expect_equal(problem$tmb_data$tau, c(20, 40, 60))
  expect_equal(problem$random, c("zeta", "g", "coef", "theta"))
  n <- length(problem$domain$labels)
  set.seed(5)
  parameters <- problem$parameters
  parameters$zeta <- stats::rnorm(n, -1, 2)
  parameters$g <- stats::rnorm(n, 0, 0.1)
  parameters$theta <- c(0.3, -0.2, 0.1)
  objective <- TMB::MakeADFun(problem$tmb_data, parameters, DLL = "alfak2", silent = TRUE)
  report <- objective$report()
  net <- (2 * exp(problem$kernel$log_self) - 1) * exp(parameters$g)
  reference <- problem$tmb_data$reference + 1L
  for (k in 1:3) {
    tau <- problem$tmb_data$tau[[k]]
    expected <- closed_form_forward_r(parameters$zeta, parameters$g, problem$kernel, tau, order = 1L) -
      net[reference] * tau + parameters$theta[[k]]
    expect_equal(report$log_intensity[, k], expected, tolerance = 1e-8)
  }
  expect_true(is.finite(objective$fn()))
  expect_true(all(is.finite(objective$gr())))
})

test_that("a series fit uses every sample", {
  sim <- simulate_series_counts()
  fit <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, birth_rate = 1, major_frequency = 0.05,
                          max_cn = 4)
  expect_s3_class(fit, "alfak2_sparse_joint_fit")
  expect_true(fit$diagnostics$finite_uncertainty)
  expect_equal(fit$diagnostics$times, c(0, 20, 40, 60))
  expect_equal(fit$diagnostics$n_major, 3L)
  summary <- fit$summary
  expect_true(all(c(paste0("count_t", 0:3), paste0("expected_t", 0:3)) %in% names(summary)))
  expect_equal(length(fit$posterior$theta), 3L)
  coverage <- fit$diagnostics$domain_coverage
  expect_equal(coverage$timepoint, paste0("t", 0:3))
  expect_equal(coverage$cells_in_domain, unname(colSums(sim$counts)))
  # The history is at least the time from the first to the last sample.
  expect_equal(fit$diagnostics$history_profile$history, sparse_joint_history_candidates(60, 1))
  # The majors' expected counts follow the observed ones at every timepoint.
  for (k in 0:3) {
    observed <- summary[[paste0("count_t", k)]][summary$major]
    expected <- summary[[paste0("expected_t", k)]][summary$major]
    expect_equal(expected / sum(expected), observed / sum(observed), tolerance = 0.1)
  }
  fitness <- stats::setNames(summary$fitness_mean, summary$karyotype)
  expect_gt(fitness[["2.3.2.2"]], fitness[["2.2.2.2"]])
  expect_gt(fitness[["2.2.2.2"]], fitness[["2.2.1.2"]])
  truth <- stats::setNames(sim$net, sim$labels)[summary$karyotype]
  expect_lt(mean(fitness[!summary$major]), min(fitness[summary$major]))
  expect_gt(stats::cor(fitness, truth, method = "spearman"), 0.5)
  prediction <- predict(fit, c("2.2.2.2", "2.2.3.3"))
  expect_equal(prediction$in_domain, c(TRUE, FALSE))
  expect_true(all(is.finite(prediction$relative_fitness)) && all(prediction$relative_sd > 0))
  expect_equal(sparse_joint_rebuild(fit, posterior = FALSE)$g, fit$posterior$g, tolerance = 1e-6)
  # The main entry point takes the intervals or the times.
  entry <- fit_alfak2(sim$counts, dt = c(20, 20, 20), beta = sim$beta,
                      sparse_joint_control = list(birth_rate = 1, major_frequency = 0.05, max_cn = 4,
                                                  history = fit$settings$history))
  expect_equal(entry$summary$fitness_mean, fit$summary$fitness_mean, tolerance = 1e-8)
})

test_that("a karyotype first seen later is centred on its frequency then", {
  y <- cbind(c(10, 0, 0, 0), c(5, 4, 0, 0), c(5, 2, 6, 0))
  prior <- sparse_joint_initial_prior_matrix(y, c(log(1e-5), 3))
  expect_equal(prior$mean, c(log(10), log(4 * 10 / 9), log(6 * 10 / 13), log(10) + log(1e-5)))
  expect_equal(prior$sd, c(10, 3, 3, 3))
  two <- sparse_joint_initial_prior(y[, 1L], y[, 2L], c(log(1e-5), 3))
  expect_equal(two, sparse_joint_initial_prior_matrix(y[, 1:2], c(log(1e-5), 3)))
})

test_that("a timepoint with no cell next to a major names the timepoints the majors come from", {
  counts <- rbind("2.2.2" = c(0, 20, 15), "3.3.2" = c(1, 0, 0), "1.1.2" = c(1, 0, 0), "2.4.2" = c(1, 0, 0),
                  "4.2.3" = c(1, 0, 0))
  err <- tryCatch(fit_sparse_joint(counts, dt = c(10, 10), beta = 0.001), alfak2_domain_coverage_error = function(e) e)
  expect_s3_class(err, "alfak2_domain_coverage_error")
  expect_match(conditionMessage(err), "^A timepoint has no cell inside the model domain\\. ")
  expect_match(conditionMessage(err), "the majors come from t1, t2, it shares no karyotype with another timepoint",
               fixed = TRUE)
  expect_equal(err$coverage$shared, c(0, 1, 1))
})

test_that("functions that take two timepoints say so", {
  sim <- simulate_series_counts()
  expect_error(fit_alfak2(sim$counts, dt = c(20, 20, 20), beta = sim$beta, pipeline = "two_layer"),
               "two_layer pipeline fits two timepoints")
  expect_error(fit_sparse_joint_cohort(list(a = sim$counts, b = sim$counts[, 1:2]), dt = 20),
               "two timepoints per sample; a has 4")
  empty <- cbind(sim$counts, 0)
  expect_error(fit_sparse_joint(empty, dt = c(20, 20, 20, 20), beta = sim$beta, max_cn = 4),
               "Every timepoint needs at least one cell")
})

test_that("intermediate samples can count for the majors only", {
  sim <- simulate_series_counts()
  all <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                          intermediate = "all")
  majors <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                             intermediate = "majors", intermediate_weight = 1)
  expect_identical(majors$settings$intermediate, "majors")
  expect_null(majors$diagnostics$intermediate_dispersion)
  expect_identical(all$settings$intermediate, "all")
  problem <- sparse_joint_problem(majors$data, majors$settings)
  n_major <- problem$domain$n_major
  weight <- problem$tmb_data$sample_weight
  expect_equal(dim(weight), c(length(problem$domain$labels), 3L))
  # Majors count in every sample; the others in the last one only.
  expect_true(all(weight[seq_len(n_major), ] == 1))
  expect_true(all(weight[-seq_len(n_major), 1:2] == 0))
  expect_true(all(weight[, 3] == 1))
  expect_true(all(problem$used[-seq_len(n_major), 2:3] == 0))
  expect_equal(problem$used[, c(1, 4)], problem$domain$y[, c(1, 4)])
  expect_false(isTRUE(all.equal(majors$summary$relative_fitness, all$summary$relative_fitness)))
  # A rebuilt posterior reproduces the fit's own.
  rebuilt <- sparse_joint_rebuild(majors)
  expect_equal(rebuilt$solver$mode, unname(unlist(majors$posterior[problem$random])), tolerance = 1e-6)
  # The majors come from the first and the last sample.
  ends <- sim$counts[, c(1, 4)]
  expect_setequal(problem$domain$labels[seq_len(n_major)],
                  sparse_joint_major(ends, 0.01, 3, 5000, 1L, 4L))
  # Estimated weight: the inverse of the majors' dispersion in a first fit.
  auto <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                           intermediate = "majors")
  dispersion <- auto$diagnostics$intermediate_dispersion
  expect_true(is.finite(dispersion))
  expect_equal(auto$diagnostics$intermediate_weight, 1 / max(1, dispersion))
  first <- fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                            intermediate = "majors", intermediate_weight = 1)
  expect_equal(dispersion, sparse_joint_intermediate_dispersion(first))
  expect_error(fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                intermediate = "majors", intermediate_weight = 2), "intermediate_weight")
  # Two samples have nothing in between: the option changes nothing.
  pair <- sim$counts[, c(1, 4)]
  expect_equal(fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60, intermediate = "majors")$summary,
               fit_sparse_joint(pair, dt = 60, beta = sim$beta, max_cn = 4, history = 60)$summary)
})
