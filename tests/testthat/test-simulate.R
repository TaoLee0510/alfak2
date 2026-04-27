test_that("toy landscape generator is reproducible", {
  a <- simulate_toy_landscape(n_chr = 3, min_cn = 1, max_cn = 3, seed = 42)
  b <- simulate_toy_landscape(n_chr = 3, min_cn = 1, max_cn = 3, seed = 42)
  expect_equal(a$karyotypes, b$karyotypes)
  expect_equal(a$fitness, b$fitness)
})

test_that("sparse simulator returns valid two-timepoint counts", {
  land <- simulate_toy_landscape(n_chr = 3, min_cn = 1, max_cn = 3, seed = 1)
  sim <- simulate_sparse_counts(land, n0 = 50, n1 = 60, seed = 2)
  expect_equal(ncol(sim$counts), 2)
  expect_true(all(sim$counts >= 0))
  expect_true(sum(sim$counts[, 1]) <= 50)
  expect_true(sum(sim$counts[, 2]) <= 60)
  expect_gt(sim$sparsity$observed_nodes, 0)
})
