# Simulated counts shared by the sparse joint tests.

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

# A population whose neighbours have been fed by the three clones for a long
# time before the first sample, sampled deeply first and thinly afterwards:
# what the neighbours' fitness is can only be read from the first sample.
simulate_balanced_counts <- function(seed = 31, depth = c(30000, 300), history = 200, dt = 60, beta = 5e-4) {
  set.seed(seed)
  major <- rbind(c(2L, 2L, 2L, 2L), c(2L, 3L, 2L, 2L), c(2L, 2L, 1L, 2L))
  labels <- union(format_karyotypes(major), format_karyotypes(one_step_neighbours(major, 1L, 4L)))
  karyotypes <- parse_karyotypes(labels)
  n <- length(labels)
  birth <- c(1.00, 1.01, 0.99, stats::runif(n - 3L, 0.80, 0.97))
  kernel <- sparse_division_kernel(karyotypes, beta)
  log_n0 <- exact_forward_reference(log(c(0.60, 0.10, 0.30, rep(0, n - 3L))), log(birth), kernel, history)
  log_n1 <- exact_forward_reference(log_n0, log(birth), kernel, dt)
  draw <- function(log_n, size) stats::rmultinom(1L, size, exp(log_n - max(log_n)))[, 1L]
  counts <- cbind(draw(log_n0, depth[1L]), draw(log_n1, depth[2L]))
  dimnames(counts) <- list(labels, c("t0", "t1"))
  list(counts = counts[rowSums(counts) > 0, , drop = FALSE], labels = labels, birth = birth,
       net = (2 * exp(kernel$log_self) - 1) * birth, dt = dt, beta = beta, history = history)
}

# Two samples from one landscape with different abundant clones. Clone
# "2.3.2.2" is the fittest. It dominates sample B; sample A is thin and short,
# and has it as a neighbour of its dominant clone seen in a few cells at most.
# Neighbours start absent, as in `simulate_sparse_joint_counts()`.
simulate_cohort_counts <- function(seed = 41, depth = c(A = 300, B = 4000), dt = c(A = 10, B = 60),
                                   beta = 5e-4) {
  set.seed(seed)
  major <- rbind(c(2L, 2L, 2L, 2L), c(2L, 3L, 2L, 2L), c(2L, 2L, 1L, 2L), c(2L, 3L, 3L, 2L))
  labels <- union(format_karyotypes(major), format_karyotypes(one_step_neighbours(major, 1L, 4L)))
  karyotypes <- parse_karyotypes(labels)
  n <- length(labels)
  birth <- c(1.00, 1.03, 0.97, 1.01, stats::runif(n - 4L, 0.80, 0.94))
  kernel <- sparse_division_kernel(karyotypes, beta)
  sample_counts <- function(frequency0, depth, dt) {
    frequency0 <- c(frequency0, rep(0, n - 4L))
    log_n1 <- exact_forward_reference(log(frequency0), log(birth), kernel, dt)
    counts <- cbind(stats::rmultinom(1L, depth, frequency0)[, 1L],
                    stats::rmultinom(1L, depth, exp(log_n1 - max(log_n1)))[, 1L])
    dimnames(counts) <- list(labels, c("t0", "t1"))
    counts[rowSums(counts) > 0, , drop = FALSE]
  }
  list(samples = list(A = sample_counts(c(0.70, 0, 0.30, 0), depth[["A"]], dt[["A"]]),
                      B = sample_counts(c(0, 0.55, 0, 0.45), depth[["B"]], dt[["B"]])),
       labels = labels, net = (2 * exp(kernel$log_self) - 1) * birth, dt = unname(dt), beta = beta)
}
