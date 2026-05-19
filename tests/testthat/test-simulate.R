test_that("GRF oracle is reproducible by seed", {
  a <- simulate_grf_landscape(n_chr = 3, n_centroids = 6, seed = 42)
  b <- simulate_grf_landscape(n_chr = 3, n_centroids = 6, seed = 42)
  k <- rbind(c(2, 2, 2), c(2, 3, 2), c(1, 2, 3))
  expect_equal(predict_landscape_fitness(a, k), predict_landscape_fitness(b, k))
})

test_that("sparse simulator returns valid two-timepoint counts", {
  land <- simulate_grf_landscape(n_chr = 3, n_centroids = 6, seed = 1, max_cn = 3)
  sim <- simulate_sparse_counts(
    land,
    beta = 0.01,
    dt = 0.5,
    n0 = 50,
    n1 = 60,
    seed = 2,
    initial_population = 200,
    time_step = 0.25
  )
  expect_equal(ncol(sim$counts), 2)
  expect_true(all(sim$counts >= 0))
  expect_equal(sum(sim$counts[, 1]), 50)
  expect_equal(sum(sim$counts[, 2]), 60)
  expect_gt(sim$sparsity$observed_nodes, 0)
})
