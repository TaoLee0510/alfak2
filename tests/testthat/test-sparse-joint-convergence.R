fake_hessian_env <- function(H) {
  env <- new.env()
  env$last.par <- 0
  env$spHess <- function(par, random = TRUE) H
  env
}

test_that("inverse iteration finds the smallest eigenvalue of a sparse Hessian", {
  set.seed(5)
  n <- 60
  A <- Matrix::rsparsematrix(n, n, 0.05)
  H <- Matrix::forceSymmetric(Matrix::crossprod(A) + Matrix::Diagonal(n, 0.5))
  H <- methods::as(H, "CsparseMatrix")
  expected <- min(eigen(as.matrix(H), symmetric = TRUE, only.values = TRUE)$values)
  # Here the smallest eigenvalues are crowded, the worst case for inverse
  # iteration; the estimate is still well within what the check needs.
  estimate <- sparse_joint_inner_min_eigen(fake_hessian_env(H))
  expect_gte(estimate, expected - 1e-10)
  expect_equal(estimate, expected, tolerance = 1e-3)
  # A nearly flat direction, as on the ridges of a degenerate fit.
  v <- sin(seq_len(n) * 3)
  v <- v / sqrt(sum(v^2))
  values <- c(1e-7, seq(0.05, 3, length.out = n - 1L))
  Q <- qr.Q(qr(cbind(v, matrix(stats::rnorm(n * (n - 1L)), n))))
  flat <- methods::as(Matrix::forceSymmetric(Matrix::Matrix(Q %*% diag(values) %*% t(Q))), "CsparseMatrix")
  expect_equal(sparse_joint_inner_min_eigen(fake_hessian_env(flat)), 1e-7, tolerance = 1e-3)
  # Not positive definite: zero.
  indefinite <- methods::as(Matrix::forceSymmetric(Matrix::Matrix(Q %*% diag(c(-1, values[-1])) %*% t(Q))),
                            "CsparseMatrix")
  expect_equal(sparse_joint_inner_min_eigen(fake_hessian_env(indefinite)), 0)
  # It leaves the random number stream alone.
  set.seed(9)
  before <- stats::runif(1)
  set.seed(9)
  sparse_joint_inner_min_eigen(fake_hessian_env(H))
  expect_equal(stats::runif(1), before)
})

test_that("inverse iteration also gives the flat direction", {
  set.seed(5)
  n <- 40
  v <- sin(seq_len(n) * 3)
  v <- v / sqrt(sum(v^2))
  Q <- qr.Q(qr(cbind(v, matrix(stats::rnorm(n * (n - 1L)), n))))
  values <- c(1e-7, seq(0.05, 3, length.out = n - 1L))
  as_sparse <- function(M) methods::as(Matrix::forceSymmetric(Matrix::Matrix(M)), "CsparseMatrix")
  flat <- sparse_joint_inner_flattest(fake_hessian_env(as_sparse(Q %*% diag(values) %*% t(Q))))
  expect_equal(flat$value, 1e-7, tolerance = 1e-3)
  expect_equal(abs(sum(flat$vector * v)), 1, tolerance = 1e-8)
  # Not positive definite: the value stays zero, the direction is the most
  # negative eigenvalue's.
  flat <- sparse_joint_inner_flattest(fake_hessian_env(as_sparse(Q %*% diag(c(-0.5, values[-1])) %*% t(Q))))
  expect_equal(flat$value, 0)
  expect_equal(abs(sum(flat$vector * v)), 1, tolerance = 1e-6)
  expect_identical(sparse_joint_inner_min_eigen(fake_hessian_env(as_sparse(Q %*% diag(values) %*% t(Q)))),
                   sparse_joint_inner_flattest(fake_hessian_env(as_sparse(Q %*% diag(values) %*% t(Q))))$value)
})

test_that("the flat direction is named by karyotype, cells and share", {
  labels <- sparse_joint_latent_labels(
    c("zeta", "zeta", "zeta", "g", "g", "g", "coef", "coef", "theta1"),
    states = data.frame(karyotype = c("2.2", "2.3", "3.2"), cells = c("40 -> 50", "0 -> 2", "1 -> 0"),
                        major = c(TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
  )
  expect_equal(labels$variable, rep(c("initial abundance", "birth rate", "trend coefficient", "depth"),
                                    c(3, 3, 2, 1)))
  expect_equal(labels$karyotype, c("2.2", "2.3", "3.2", "2.2", "2.3", "3.2", NA, NA, NA))
  expect_equal(labels$cells[c(2, 6)], c("0 -> 2", "1 -> 0"))
  expect_equal(labels$major[1:3], c(TRUE, FALSE, FALSE))
  vector <- c(0.05, 0.95, 0, 0, 0, 0.3, 0.01, 0, 0)
  flat <- sparse_joint_flat_direction(vector, labels)
  expect_equal(flat$karyotype, c("2.3", "3.2"))
  expect_equal(flat$variable, c("initial abundance", "birth rate"))
  expect_equal(sum(flat$share), sum(vector[c(2, 6)]^2) / sum(vector^2))
  expect_equal(sparse_joint_describe_flat(flat),
               "along the initial abundance of karyotype 2.3 (not major, cells 0 -> 2, 91%)")
  composition <- sparse_joint_flat_composition(vector, labels)
  expect_equal(composition$variable, c("initial abundance", "birth rate"))
  expect_equal(composition$major, c(FALSE, FALSE))
  expect_equal(composition$share, c(sum(vector[1:2]^2) - vector[1]^2, vector[6]^2) / sum(vector^2))
  expect_equal(sum(sparse_joint_flat_composition(vector, labels, min_share = 0)$share), 1)
  expect_equal(sparse_joint_describe_flat(flat, composition),
               paste("in a direction made up of initial abundances of karyotypes that are not major (91%),",
                     "led by the initial abundance of karyotype 2.3 (not major, cells 0 -> 2, 91%)"))
  expect_null(sparse_joint_flat_direction(vector[-1], labels))
  expect_null(sparse_joint_flat_direction(NULL, labels))
  # A cohort: states carry their sample; shared birth rates their node.
  cohort <- sparse_joint_latent_labels(
    c("zeta", "zeta", "m", "m"),
    states = data.frame(karyotype = c("2.2", "2.2"), sample = c("A", "B"), cells = c("3 -> 4", "0 -> 1"),
                        major = c(TRUE, FALSE), stringsAsFactors = FALSE),
    nodes = data.frame(karyotype = c("2.2", "2.3"), major = c(TRUE, FALSE), stringsAsFactors = FALSE)
  )
  expect_equal(cohort$sample, c("A", "B", NA, NA))
  expect_equal(cohort$karyotype, c("2.2", "2.2", "2.2", "2.3"))
  expect_equal(sparse_joint_describe_flat(sparse_joint_flat_direction(c(0, 1, 0, 0), cohort)),
               "along the initial abundance of karyotype 2.2 in sample B (not major, cells 0 -> 1, 100%)")
})

test_that("finite-difference and projected gradients respect the bounds", {
  fn <- function(p) sum((p - c(1, -2))^2)
  g <- sparse_joint_fd_gradient(fn, c(0, 0), -5, 5, 1e-3)
  expect_equal(g, c(-2, 4), tolerance = 1e-6)
  # One-sided at a bound.
  g <- sparse_joint_fd_gradient(fn, c(-5, 0), -5, 5, 1e-3)
  expect_equal(g[[1L]], -12, tolerance = 1e-2)
  expect_equal(sparse_joint_projected_gradient(c(3, -1), c(-5, 0), -5, 5), c(0, -1))
  expect_equal(sparse_joint_projected_gradient(c(-3, -1), c(-5, 5), -5, 5), c(-3, 0))
})

test_that("compass search finds a bounded minimum and reports improvements", {
  fn <- function(p) sum((p - c(0.3, -7))^2)
  hits <- 0L
  run <- sparse_joint_compass_search(fn, c(-1, -1), lower = c(-2, -3), upper = c(2, 3),
                                     improved_hook = function() hits <<- hits + 1L)
  expect_lte(max(abs(run$par - c(0.3, -3))), 2e-3)
  expect_equal(run$convergence, 0L)
  expect_gt(hits, 1L)
  # Points where the objective is not finite are skipped.
  walled <- function(p) if (p[[1L]] > 0) NaN else fn(p)
  run <- sparse_joint_compass_search(walled, c(-1, -1), lower = -5, upper = 5)
  expect_lte(run$par[[1L]], 0)
  expect_true(is.finite(run$objective))
})

test_that("the candidate kept is the lowest finite objective, or the best accepted one on request", {
  expect_equal(sparse_joint_choose_candidate(c(10, 5, 7), c(0, 1, 0)), 2L)
  expect_equal(sparse_joint_choose_candidate(c(Inf, 5, NA), c(0, 1, 0)), 2L)
  expect_identical(sparse_joint_choose_candidate(c(Inf, NA), c(0, 0)), NA_integer_)
  old <- options(alfak2.choose_accepted = TRUE)
  on.exit(options(old), add = TRUE)
  expect_equal(sparse_joint_choose_candidate(c(10, 5, 7), c(0, 1, 0)), 3L)
  expect_equal(sparse_joint_choose_candidate(c(10, 5, 7), c(1, 1, 1)), 2L)
  expect_equal(sparse_joint_choose_candidate(c(Inf, 5, NA), c(0, 1, 0)), 2L)
})

test_that("checking the scales leaves TMB where the fit left it", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_alfak2_data(sim$counts, dt = sim$dt, beta = sim$beta)
  settings <- list(birth_rate = 1, birth_rate_log_sd = 0.05, major_frequency = 0.05, major_count = 3,
                   max_states = 5000, min_cn = 1L, max_cn = 4L, shell_power = 0.85, shell_sd = 0.5,
                   initial_frequency = c(mean = log(1e-5), sd = 3), initial_state = "balance", history = sim$dt,
                   balance_sd = 1, calibration = sparse_joint_default_calibration,
                   inner_control = list(maxit = 200, tol = 1e-10))
  problem <- sparse_joint_problem(data, settings)
  make <- function() TMB::MakeADFun(problem$tmb_data, problem$parameters, map = problem$map,
                                    random = problem$random, DLL = "alfak2", silent = TRUE,
                                    inner.control = settings$inner_control)
  plain <- make()
  optimum <- suppressWarnings(stats::nlminb(plain$par, plain$fn, plain$gr, lower = log(1e-4), upper = log(0.5)))
  plain$fn(optimum$par)
  checked <- make()
  suppressWarnings(stats::nlminb(checked$par, checked$fn, checked$gr, lower = log(1e-4), upper = log(0.5)))
  check <- sparse_joint_check_scales(checked, optimum$par, optimum$objective, log(1e-4), log(0.5),
                                     sparse_joint_convergence_tolerance)
  expect_named(check, c("pass", "reasons", "objective_difference", "joint_difference", "min_eigen",
                        "projected_gradient", "flat_direction", "flat_composition"))
  expect_true(identical(checked$env$random.start, plain$env$random.start))
  checked$fn(optimum$par)
  expect_identical(checked$env$last.par, plain$env$last.par)
})

test_that("a fit reports how its scale optimisation was checked", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                          history = sim$dt)
  check <- fit$diagnostics$optimisation
  expect_true(check$status %in% c("converged", "not converged"))
  expect_equal(fit$diagnostics$convergence, if (check$status == "converged") 0L else 1L)
  expect_true(all(c("stage", "reasons", "restarts", "evaluations", "nlminb_convergence", "min_eigen",
                    "projected_gradient") %in% names(check)))
  if (check$status == "converged") {
    expect_gte(check$min_eigen, sparse_joint_convergence_tolerance[["eigen"]])
    expect_lte(check$projected_gradient, sparse_joint_convergence_tolerance[["gradient"]])
    expect_null(check$flat_direction)
  } else {
    expect_s3_class(check$flat_direction, "data.frame")
  }
})

test_that("a fit that is not accepted warns once, with its reasons", {
  sim <- simulate_sparse_joint_counts()
  # An impossible gradient tolerance: nothing can be accepted.
  strict <- sparse_joint_convergence_tolerance
  strict[["gradient"]] <- -1
  local_mocked_bindings(sparse_joint_convergence_tolerance = strict)
  expect_warning(
    fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                            history = c(sim$dt, 300)),
    class = "alfak2_convergence_warning"
  )
  expect_equal(fit$diagnostics$convergence, 1L)
  expect_equal(fit$diagnostics$optimisation$status, "not converged")
  expect_true("gradient not small" %in% fit$diagnostics$optimisation$reasons)
  expect_true(all(fit$diagnostics$history_profile$convergence == 1))
  # No search again by default: the fit is nlminb's.
  expect_equal(fit$diagnostics$optimisation$stage, "first")
  expect_equal(fit$diagnostics$optimisation$restarts, 0L)
})

test_that("a fit with a near-singular inner Hessian names the flat direction", {
  sim <- simulate_sparse_joint_counts()
  # An impossible eigenvalue tolerance: every inner Hessian counts as
  # near-singular.
  strict <- sparse_joint_convergence_tolerance
  strict[["eigen"]] <- Inf
  local_mocked_bindings(sparse_joint_convergence_tolerance = strict)
  warning <- NULL
  fit <- withCallingHandlers(
    fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05, max_cn = 4,
                     history = sim$dt),
    alfak2_convergence_warning = function(w) {
      warning <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  flat <- fit$diagnostics$optimisation$flat_direction
  expect_s3_class(flat, "data.frame")
  expect_named(flat, c("variable", "karyotype", "sample", "cells", "major", "share"))
  composition <- fit$diagnostics$optimisation$flat_composition
  expect_named(composition, c("variable", "major", "share"))
  expect_lte(sum(composition$share), 1 + 1e-12)
  expect_gte(composition$share[[1L]], flat$share[[1L]])
  expect_true(all(diff(flat$share) <= 0))
  expect_lte(sum(flat$share), 1 + 1e-12)
  # Each karyotype named is a domain state, with its own counts.
  named <- flat[!is.na(flat$karyotype), , drop = FALSE]
  expect_gt(nrow(named), 0L)
  row <- match(named$karyotype, fit$summary$karyotype)
  expect_false(anyNA(row))
  expect_equal(named$cells, paste(fit$summary$count_t0[row], fit$summary$count_t1[row], sep = " -> "))
  expect_equal(named$major, fit$summary$major[row])
  expect_match(warning, "nearly flat in a direction made up of ", fixed = TRUE)
  expect_match(warning, if (is.na(flat$karyotype[[1L]])) flat$variable[[1L]] else flat$karyotype[[1L]],
               fixed = TRUE)
})

test_that("when nlminb fails, the prior-centred scales are kept, as before the check", {
  sim <- simulate_sparse_joint_counts()
  data <- prepare_alfak2_data(sim$counts, dt = sim$dt, beta = sim$beta)
  settings <- list(birth_rate = 1, birth_rate_log_sd = 0.05, major_frequency = 0.05, major_count = 3,
                   max_states = 5000, min_cn = 1L, max_cn = 4L, shell_power = 0.85, shell_sd = 0.5,
                   initial_frequency = c(mean = log(1e-5), sd = 3), initial_state = "balance", history = sim$dt,
                   balance_sd = 1, calibration = sparse_joint_default_calibration,
                   inner_control = list(maxit = 200, tol = 1e-10))
  problem <- sparse_joint_problem(data, settings)
  objective <- TMB::MakeADFun(problem$tmb_data, problem$parameters, map = problem$map, random = problem$random,
                              DLL = "alfak2", silent = TRUE, inner.control = settings$inner_control)
  # A gradient that fails makes nlminb fail.
  objective$gr <- function(...) stop("no gradient")
  optimum <- sparse_joint_optimise_scales(objective, list(eval.max = 200, iter.max = 100))
  expect_identical(optimum$par, objective$par)
  expect_equal(optimum$objective, as.numeric(objective$fn(objective$par)))
  expect_equal(optimum$convergence, 1L)
  expect_equal(optimum$check$stage, "none")
  expect_equal(optimum$check$restarts, 0L)
  expect_identical(optimum$check$reasons[[1L]], "nlminb failed")
  expect_identical(optimum$check$nlminb_message, "nlminb failed")
})

test_that("searching again is available on request", {
  sim <- simulate_sparse_joint_counts()
  strict <- sparse_joint_convergence_tolerance
  strict[["gradient"]] <- -1
  local_mocked_bindings(sparse_joint_convergence_tolerance = strict)
  old <- options(alfak2.convergence_restarts = TRUE)
  on.exit(options(old), add = TRUE)
  fit <- suppressWarnings(fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, major_frequency = 0.05,
                                           max_cn = 4, history = sim$dt))
  expect_gt(fit$diagnostics$optimisation$restarts, 0L)
  expect_false(identical(fit$diagnostics$optimisation$stage, "first"))
})
