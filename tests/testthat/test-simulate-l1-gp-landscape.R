test_that("paper-style GRF landscape is lightweight and reproducible", {
  x <- simulate_grf_landscape(n_chr = 22, n_centroids = 8, seed = 42)
  y <- simulate_grf_landscape(n_chr = 22, n_centroids = 8, seed = 42)
  z <- simulate_grf_landscape(n_chr = 22, n_centroids = 8, seed = 43)

  expect_s3_class(x, "alfak2_grf_landscape")
  expect_null(x$karyotypes)
  expect_null(x$fitness)
  expect_equal(x$centroids, y$centroids)
  expect_false(isTRUE(all.equal(x$centroids, z$centroids)))
})

test_that("GRF oracle anchors founder and is deterministic", {
  x <- simulate_grf_landscape(
    n_chr = 4,
    n_centroids = 10,
    founder = c(2, 2, 2, 2),
    founder_fitness = 1.25,
    scale = 0.8,
    seed = 1
  )
  k <- rbind(
    c(2, 2, 2, 2),
    c(2, 2, 3, 2),
    c(1, 2, 3, 2)
  )
  f1 <- predict_landscape_fitness(x, k)
  f2 <- predict_landscape_fitness(x, format_karyotypes(k))

  expect_equal(f1, f2)
  expect_equal(f1[1], 1.25, tolerance = 1e-12)
  expect_equal(f1, predict_landscape_fitness(x, k), tolerance = 1e-12)
})

test_that("legacy landscape constructor now returns lazy GRF", {
  x <- simulate_l1_gp_landscape(
    n_chr = 3,
    min_cn = 1,
    max_cn = 4,
    diploid_fitness = 1.1,
    seed = 1
  )
  expect_s3_class(x, "alfak2_grf_landscape")
  expect_equal(
    predict_landscape_fitness(x, matrix(c(2, 2, 2), nrow = 1)),
    1.1,
    tolerance = 1e-12
  )
})

test_that("lazy GRF sparse simulator returns valid observed truth", {
  land <- simulate_grf_landscape(
    n_chr = 22,
    n_centroids = 8,
    lambda = 0.8,
    seed = 1,
    max_cn = 4
  )
  sim <- simulate_sparse_counts(
    land,
    beta = 0.001,
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

  truth <- predict_landscape_fitness(land, rownames(sim$counts))
  expect_equal(sim$truth_observed, truth, tolerance = 1e-12)
})
