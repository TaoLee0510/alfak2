labelled_counts <- function(...) {
  rows <- list(...)
  counts <- do.call(rbind, rows)
  dimnames(counts) <- list(names(rows), c("t0", "t1"))
  counts
}

# Counts drawn from the model's own mean dynamics: three well-sampled clones
# with known birth rates and their one-step neighbours, all less fit. The
# missegregation rate is high enough for neighbours to be seen in a few cells,
# so fits use `major_frequency = 0.05` to anchor on the three clones only.
simulate_sparse_joint_counts <- function(seed = 21, depth = 4000, dt = 60, beta = 5e-4) {
  set.seed(seed)
  major <- rbind(c(2L, 2L, 2L, 2L), c(2L, 3L, 2L, 2L), c(2L, 2L, 1L, 2L))
  labels <- union(format_karyotypes(major), format_karyotypes(one_step_neighbours(major, 1L, 4L)))
  karyotypes <- parse_karyotypes(labels)
  n <- length(labels)
  birth <- c(1.00, 1.03, 0.97, stats::runif(n - 3L, 0.80, 0.94))
  frequency0 <- c(0.60, 0.10, 0.30, rep(0, n - 3L))
  kernel <- sparse_division_kernel(karyotypes, beta)
  log_n1 <- exact_forward_reference(log(frequency0), log(birth), kernel, dt)
  frequency1 <- exp(log_n1 - max(log_n1))
  frequency1 <- frequency1 / sum(frequency1)
  counts <- cbind(stats::rmultinom(1L, depth, frequency0)[, 1L], stats::rmultinom(1L, depth, frequency1)[, 1L])
  dimnames(counts) <- list(labels, c("t0", "t1"))
  list(counts = counts[rowSums(counts) > 0, , drop = FALSE], labels = labels, birth = birth,
       net = (2 * exp(kernel$log_self) - 1) * birth, dt = dt, beta = beta)
}

test_that("the domain is the major karyotypes and all their one-step neighbours", {
  counts <- labelled_counts("2.2.2" = c(40, 30), "2.3.2" = c(4, 8), "3.2.2" = c(9, 0),
                            "2.2.3" = c(0, 1), "4.4.2" = c(1, 0))
  domain <- sparse_joint_domain(counts, major_frequency = 0.1, major_count = 3, max_states = Inf,
                                min_cn = 1L, max_cn = 4L)
  expect_equal(domain$n_major, 3L)
  expect_equal(domain$labels[1:3], c("2.2.2", "2.3.2", "3.2.2"))
  expected <- union(c("2.2.2", "2.3.2", "3.2.2"),
                    format_karyotypes(one_step_neighbours(parse_karyotypes(c("2.2.2", "2.3.2", "3.2.2")), 1L, 4L)))
  expect_setequal(domain$labels, expected)
  expect_true(all(domain$karyotypes >= 1L & domain$karyotypes <= 4L))
  # a neighbour seen in one cell is in the domain with its count
  expect_equal(domain$y1[domain$labels == "2.2.3"], 1)
  # an observed karyotype far from every major is left out, with its cell
  expect_false("4.4.2" %in% domain$labels)
  expect_equal(domain$dropped, "4.4.2")
  expect_equal(domain$dropped_cells, 1)
  expect_equal(sum(domain$y0) + sum(domain$y1), sum(counts) - 1)
  distance <- sparse_joint_major_distance(domain$karyotypes, domain$major)
  expect_equal(distance[1:3], c(0, 0, 0))
  expect_true(all(distance[-(1:3)] == 1))
})

test_that("majors need both a share of a sample and a number of cells", {
  counts <- labelled_counts("2.2.2" = c(30, 20), "2.3.2" = c(0, 3), "3.2.2" = c(1, 0), "2.2.1" = c(0, 1))
  # shares of either sample: 0.97, 0.125, 0.032, 0.042
  expect_equal(sparse_joint_major(counts, 0.01, major_count = 3), c("2.2.2", "2.3.2"))
  expect_equal(sparse_joint_major(counts, 0.01, major_count = 4), "2.2.2")
  expect_equal(sparse_joint_major(counts, 0.2, major_count = 3), "2.2.2")
  # a single cell holds 3% of a sample here but does not make a major, so its
  # neighbours are not brought into the domain
  domain <- sparse_joint_domain(counts, 0.01, 3, Inf, 1L, 4L)
  expect_equal(domain$n_major, 2L)
  expect_false("4.2.2" %in% domain$labels)
  expect_true("3.2.2" %in% domain$labels)
  # with no karyotype reaching the count, the most abundant one anchors
  thin <- labelled_counts("2.2.2" = c(2, 0), "2.3.2" = c(0, 1), "3.2.2" = c(1, 0))
  expect_equal(sparse_joint_major(thin, 0.01, 3), "2.2.2")
  # a budget raises the share until the domain fits
  one <- sparse_joint_domain(counts, 0.01, 3, 8, 1L, 4L)
  expect_gt(length(domain$labels), 8)
  expect_equal(one$n_major, 1L)
  expect_lte(length(one$labels), 8)
})

test_that("compiled forward model equals the R closed form", {
  sim <- simulate_sparse_joint_counts()
  domain <- sparse_joint_domain(sim$counts, 0.05, 3, Inf, 1L, 4L)
  n <- length(domain$labels)
  kernel <- sparse_division_kernel(domain$karyotypes, sim$beta)
  design <- sparse_joint_design(domain$karyotypes, rep(2, 4), domain$major, 0.85)
  set.seed(5)
  zeta <- stats::rnorm(n, -1, 2)
  g <- stats::rnorm(n, 0, 0.1)
  theta1 <- 0.3
  tmb_data <- list(
    model = "sparse_joint", y0 = domain$y0, y1 = domain$y1, dt = sim$dt,
    edge_from = kernel$from - 1L, edge_to = kernel$to - 1L, edge_log_prob = kernel$log_prob,
    log_self = kernel$log_self, design = design$design, coef_group = design$group,
    reference = 0L,
    zeta_prior_mean = rep(0, n), zeta_prior_sd = rep(3, n), intercept_mean = 0, intercept_sd = 1,
    shell_sd = 0.5, log_scale_prior_mean = rep(log(0.05), 4), log_scale_prior_sd = rep(1.5, 4)
  )
  parameters <- list(zeta = zeta, g = g, coef = rep(0, ncol(design$design)), theta1 = theta1,
                     log_scale = rep(log(0.05), 4))
  objective <- TMB::MakeADFun(tmb_data, parameters, DLL = "alfak2", silent = TRUE)
  report <- objective$report()
  net <- (2 * exp(kernel$log_self) - 1) * exp(g)
  expect_equal(as.numeric(report$net), net, tolerance = 1e-12)
  expected <- closed_form_forward_r(zeta, g, kernel, sim$dt, order = 1L) - net[1L] * sim$dt + theta1
  expect_equal(as.numeric(report$log_intensity1), expected, tolerance = 1e-8)
  expect_true(is.finite(objective$fn()))
  expect_true(all(is.finite(objective$gr())))
})

test_that("fit recovers the order of clones and the deficit of their neighbours", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, birth_rate = 1, major_frequency = 0.05, max_cn = 4)
  expect_s3_class(fit, "alfak2_sparse_joint_fit")
  expect_true(fit$diagnostics$finite_uncertainty)
  expect_equal(fit$diagnostics$n_major, 3L)
  expect_equal(fit$diagnostics$dropped_karyotypes, 0L)
  summary <- fit$summary
  expect_true(all(is.finite(summary$fitness_mean)) && all(summary$fitness_sd > 0))
  expect_true(all(summary$conf_low < summary$fitness_mean & summary$fitness_mean < summary$conf_high))
  fitness <- stats::setNames(summary$fitness_mean, summary$karyotype)
  # rising clone above the dominant clone above the declining clone
  expect_gt(fitness[["2.3.2.2"]], fitness[["2.2.2.2"]])
  expect_gt(fitness[["2.2.2.2"]], fitness[["2.2.1.2"]])
  truth <- stats::setNames(sim$net, sim$labels)[summary$karyotype]
  expect_equal(unname(fitness[["2.3.2.2"]] - fitness[["2.2.2.2"]]),
               unname(truth[["2.3.2.2"]] - truth[["2.2.2.2"]]), tolerance = 0.5)
  # neighbours are less fit than the clones they surround
  expect_lt(mean(fitness[!summary$major]), min(fitness[summary$major]))
  shell <- fit$trend$coefficients[length(fit$trend$coefficients)]
  expect_lt(shell, 0)
  expect_gt(stats::cor(fitness, truth, method = "spearman"), 0.5)
  # the supplied baseline is respected
  expect_equal(fit$diagnostics$fitted_birth_rate, 1, tolerance = 0.2)
})

test_that("prediction returns the posterior inside the domain and the trend outside it", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, birth_rate = 1, major_frequency = 0.05, max_cn = 4)
  inside <- fit$summary$karyotype[c(1L, 10L)]
  two_steps <- c("2.2.3.3", "3.3.3.2")
  prediction <- predict(fit, c(inside, two_steps))
  expect_equal(prediction$in_domain, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(prediction$fitness_mean[1:2], fit$summary$fitness_mean[c(1L, 10L)])
  expect_equal(prediction$fitness_sd[1:2], fit$summary$fitness_sd[c(1L, 10L)])
  expect_equal(prediction$shell[3:4], c(2, 2))
  expect_true(all(is.finite(prediction$fitness_mean)) && all(prediction$fitness_sd > 0))
  # two steps away is predicted below the average one-step neighbour
  expect_lt(max(prediction$fitness_mean[3:4]), mean(fit$summary$fitness_mean[!fit$summary$major]))
  # labels and matrices are interchangeable
  expect_equal(predict(fit, parse_karyotypes(two_steps))$fitness_mean, prediction$fitness_mean[3:4])
})

test_that("an observed karyotype outside the domain is placed at the level of observed neighbours", {
  sim <- simulate_sparse_joint_counts()
  counts <- rbind(sim$counts, "4.4.4.4" = c(1L, 0L), "4.4.4.3" = c(0L, 1L))
  fit <- fit_sparse_joint(counts, dt = sim$dt, beta = sim$beta, birth_rate = 1, major_frequency = 0.05, max_cn = 4)
  expect_equal(fit$dropped, c("4.4.4.4", "4.4.4.3"))
  expect_equal(fit$diagnostics$dropped_cells, 2)
  prediction <- predict(fit, c("4.4.4.4", "4.4.3.3"))
  expect_equal(prediction$in_domain, c(FALSE, FALSE))
  observed_minor <- !fit$summary$major & (fit$summary$count_t0 + fit$summary$count_t1) > 0
  self <- 2 * sparse_joint_self_probability(parse_karyotypes("4.4.4.4"), sim$beta) - 1
  expect_equal(prediction$fitness_mean[1L],
               self * exp(mean(fit$summary$log_birth_mean[observed_minor])), tolerance = 1e-10)
  # an unobserved karyotype equally far away is predicted lower
  expect_lt(prediction$fitness_mean[2L], prediction$fitness_mean[1L])
})

test_that("invalid inputs are rejected", {
  sim <- simulate_sparse_joint_counts()
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 2), "copy-number bounds")
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, birth_rate = 0), "birth_rate")
  expect_error(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, initial_frequency = c(0, -1)),
               "initial_frequency")
})
