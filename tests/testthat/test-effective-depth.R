test_that("effective-depth preprocessing preserves controlled column totals", {
  counts <- matrix(
    c(90, 9,
      10, 1,
      0, 90),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3", "3.3"), c("t0", "t1"))
  )
  out <- apply_effective_depth_counts(counts, effective_depth = 20, effective_depth_mode = "fixed")
  info <- attr(out, "effective_depth_info")

  expect_equal(colSums(out), c(t0 = 20L, t1 = 20L))
  expect_equal(info$effective_depth, c(20, 20))
  expect_equal(info$raw_depth, c(100, 100))
  expect_true(all(rowSums(out) > 0))
})

test_that("hash effective-depth rounding is stable under row reordering", {
  counts <- matrix(
    c(1, 1,
      1, 1,
      1, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3", "3.3"), c("t0", "t1"))
  )
  out1 <- apply_effective_depth_counts(
    counts,
    effective_depth = 2,
    effective_depth_mode = "fixed",
    effective_depth_rounding = "hash"
  )
  out2 <- apply_effective_depth_counts(
    counts[c("3.3", "2.2", "2.3"), , drop = FALSE],
    effective_depth = 2,
    effective_depth_mode = "fixed",
    effective_depth_rounding = "hash"
  )

  common <- intersect(rownames(out1), rownames(out2))
  expect_equal(out1[common, , drop = FALSE], out2[common, , drop = FALSE])
  expect_equal(sort(rownames(out1)), sort(rownames(out2)))
})

test_that("stochastic effective-depth rounding is seed-reproducible", {
  counts <- matrix(
    c(1, 1,
      1, 1,
      1, 1,
      1, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3", "3.2", "3.3"), c("t0", "t1"))
  )
  out1 <- apply_effective_depth_counts(
    counts,
    effective_depth = 2,
    effective_depth_mode = "fixed",
    effective_depth_rounding = "stochastic",
    effective_depth_seed = 11
  )
  out2 <- apply_effective_depth_counts(
    counts,
    effective_depth = 2,
    effective_depth_mode = "fixed",
    effective_depth_rounding = "stochastic",
    effective_depth_seed = 11
  )

  expect_equal(out1, out2)
})

test_that("effective-depth preprocessing preserves observation weights", {
  counts <- matrix(
    c(90, 9,
      10, 1,
      0, 90),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3", "3.3"), c("t0", "t1"))
  )
  weights <- matrix(c(1, 1, 0.2, 0.2, 1, 1), nrow = 3, byrow = TRUE,
                    dimnames = list(rownames(counts), c("t0", "t1")))
  attr(counts, "observation_weights") <- weights
  out <- apply_effective_depth_counts(counts, effective_depth = 20, effective_depth_mode = "fixed")
  out_weights <- attr(out, "observation_weights")

  expect_equal(out_weights[rownames(out), , drop = FALSE], weights[rownames(out), , drop = FALSE])
})

test_that("effective-depth preprocessing preserves zero-weight holdout targets", {
  counts <- matrix(
    c(99, 99,
      1, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("1.1", "3.3"), c("t0", "t1"))
  )
  weights <- matrix(1, nrow = 2, ncol = 2, dimnames = dimnames(counts))
  weights["3.3", ] <- 0
  attr(counts, "observation_weights") <- weights
  attr(counts, "holdout_mode") <- list(mode = "zero_observation_weight", labels = "3.3")

  dat <- prepare_counts_for_input_depth(
    counts,
    dt = 1,
    beta = 0.01,
    input_depth = "effective",
    effective_depth = 1,
    effective_depth_mode = "fixed"
  )
  graph <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3)

  expect_false("3.3" %in% dat$labels)
  expect_equal(dat$metadata$holdout_mode$labels, "3.3")
  expect_equal(graph$support_tier[match("3.3", graph$labels)], "weakly_supported")
  expect_false("2.3" %in% graph$labels)
})

test_that("effective-depth fit defaults to dirichlet-multinomial controls", {
  counts <- matrix(
    c(30, 20,
      10, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  fit <- fit_alfak2(
    counts,
    dt = 1,
    beta = 0.01,
    local_shell_depth = 1,
    global_extra_shell = 0,
    min_cn = 1,
    max_cn = 3,
    max_nodes = 100,
    lambda_l_grid = 1,
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05,
    input_depth = "effective",
    effective_depth = 20,
    effective_depth_mode = "fixed",
    control = list(eval.max = 120, iter.max = 120)
  )

  expect_s3_class(fit, "alfak2_fit")
  expect_equal(fit$local$diagnostics$observation_model, "dirichlet_multinomial")
  expect_equal(fit$diagnostics$dm_concentration, 50)
  expect_equal(colSums(fit$data$counts), c(t0 = 20L, t1 = 20L))
  expect_equal(fit$data$metadata$input_depth$effective_depth, c(20, 20))
  expect_true(all(is.finite(summarize_alfak2(fit)$fitness_mean)))
})

test_that("raw input depth keeps the existing multinomial default", {
  counts <- matrix(
    c(30, 20,
      10, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  fit <- fit_alfak2(
    counts,
    dt = 1,
    beta = 0.01,
    local_shell_depth = 1,
    global_extra_shell = 0,
    min_cn = 1,
    max_cn = 3,
    max_nodes = 100,
    lambda_l_grid = 1,
    lambda_e_grid = 0.1,
    sigma_obs_grid = 0.05,
    control = list(eval.max = 120, iter.max = 120)
  )

  expect_equal(fit$local$diagnostics$observation_model, "multinomial")
  expect_equal(fit$diagnostics$dm_concentration, 200)
  expect_equal(fit$diagnostics$input_depth$input_depth, "raw")
})
