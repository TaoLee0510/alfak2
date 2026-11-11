small_domain <- function() {
  grid <- as.matrix(expand.grid(1:4, 1:4, 1:4))
  grid[rowSums(abs(sweep(grid, 2L, c(2, 2, 2)))) <= 2, , drop = FALSE]
}

dense_reference <- function(log_n0, log_birth, kernel, dt) {
  birth <- exp(log_birth)
  net <- (2 * exp(kernel$log_self) - 1) * birth
  generator <- diag(net - max(net))
  for (e in seq_along(kernel$from)) {
    generator[kernel$to[e], kernel$from[e]] <- generator[kernel$to[e], kernel$from[e]] +
      2 * exp(kernel$log_prob[e]) * birth[kernel$from[e]]
  }
  log(as.numeric(Matrix::expm(Matrix::Matrix(generator * dt)) %*% exp(log_n0))) + max(net) * dt
}

test_that("sparse kernel reproduces the exact single-daughter kernel on one-step edges", {
  karyotypes <- small_domain()
  beta <- 0.003
  kernel <- sparse_division_kernel(karyotypes, beta)
  full <- alfak2_transition_operator_cpp(karyotypes, beta, "exact", division_process = TRUE)
  full_matrix <- matrix(0, nrow(karyotypes), nrow(karyotypes))
  full_matrix[cbind(full$from0 + 1L, full$to0 + 1L)] <- full$weight
  expect_equal(exp(kernel$log_self), diag(full_matrix), tolerance = 1e-12)
  expect_equal(exp(kernel$log_prob), full_matrix[cbind(kernel$from, kernel$to)], tolerance = 1e-12)
  one_step <- rowSums(abs(karyotypes[kernel$from, , drop = FALSE] - karyotypes[kernel$to, , drop = FALSE]))
  expect_true(all(one_step == 1))
  # every one-step pair inside the domain is an edge, in both directions
  distance <- as.matrix(stats::dist(karyotypes, method = "manhattan"))
  expect_equal(length(kernel$from), sum(distance == 1))
})

test_that("phi functions equal the time-ordered integrals they stand for", {
  phi1 <- function(a, b, dt) {
    stats::integrate(function(s) exp(a * (dt - s) + b * s), 0, dt, rel.tol = 1e-12)$value
  }
  phi2 <- function(a, b, c, dt) {
    inner <- function(s2) {
      exp(a * (dt - s2)) * stats::integrate(function(s1) exp(b * (s2 - s1) + c * s1), 0, s2, rel.tol = 1e-12)$value
    }
    stats::integrate(Vectorize(inner), 0, dt, rel.tol = 1e-11)$value
  }
  cases <- rbind(c(0.3, -0.2, 0.1), c(0.1, 0.1, 0.1), c(0.1, 0.1 + 1e-9, -0.4),
                 c(-0.5, 0.2, 0.2), c(0.05, -0.01, 0.049999), c(-1, -3, -2))
  for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    expect_equal(log_phi1(x[1], x[2], 7), log(phi1(x[1], x[2], 7)), tolerance = 1e-9)
    expect_equal(log_phi2(x[1], x[2], x[3], 7), log(phi2(x[1], x[2], x[3], 7)), tolerance = 1e-8)
  }
  # symmetric in its arguments, and finite where exp(a * dt) would overflow
  expect_equal(log_phi1(0.2, -0.7, 5), log_phi1(-0.7, 0.2, 5))
  expect_true(is.finite(log_phi1(3, 1, 1000)))
  expect_true(is.finite(log_phi2(3, 1, 2, 1000)))
})

test_that("uniformisation reference agrees with a dense matrix exponential", {
  set.seed(11)
  karyotypes <- small_domain()
  n <- nrow(karyotypes)
  for (dt in c(30, 360)) {
    kernel <- sparse_division_kernel(karyotypes, 5e-4)
    log_birth <- log(stats::runif(n, 0.4, 1))
    log_n0 <- log(stats::rexp(n))
    expect_equal(exact_forward_reference(log_n0, log_birth, kernel, dt),
                 dense_reference(log_n0, log_birth, kernel, dt), tolerance = 1e-10)
  }
})

test_that("closed form converges to the exact solution as the flux vanishes", {
  set.seed(12)
  karyotypes <- small_domain()
  n <- nrow(karyotypes)
  log_birth <- log(stats::runif(n, 0.85, 0.95))
  log_n0 <- log(stats::rexp(n))
  dt <- 60
  error <- function(beta, order) {
    kernel <- sparse_division_kernel(karyotypes, beta)
    max(abs(closed_form_forward_r(log_n0, log_birth, kernel, dt, order) -
              exact_forward_reference(log_n0, log_birth, kernel, dt)))
  }
  first <- c(error(2e-5, 1L), error(1e-5, 1L))
  second <- c(error(2e-5, 2L), error(1e-5, 2L))
  # the neglected term is one order higher in the flux
  expect_lt(first[1], 1e-3)
  expect_equal(first[1] / first[2], 4, tolerance = 0.2)
  expect_lt(second[1], first[1] / 50)
  expect_equal(second[1] / second[2], 8, tolerance = 0.3)
})

test_that("a deleterious neighbour settles at mutation-selection balance with its parent", {
  karyotypes <- rbind(c(2L, 2L), c(2L, 3L))
  beta <- 5e-5
  kernel <- sparse_division_kernel(karyotypes, beta)
  log_birth <- log(c(0.9, 0.8))
  # the neighbour starts absent
  log_n1 <- closed_form_forward_r(c(0, -Inf), log_birth, kernel, dt = 720, order = 1L)
  net <- (2 * exp(kernel$log_self) - 1) * exp(log_birth)
  edge <- which(kernel$from == 1L & kernel$to == 2L)
  flux <- 2 * exp(kernel$log_prob[edge]) * exp(log_birth[1L])
  expect_equal(exp(log_n1[2L] - log_n1[1L]), flux / (net[1L] - net[2L]), tolerance = 1e-6)
  # and the first-order result is then already exact to the order of the flux
  exact <- exact_forward_reference(c(0, -Inf), log_birth, kernel, dt = 720)
  expect_equal(log_n1[2L] - log_n1[1L], exact[2L] - exact[1L], tolerance = 1e-3)
})
