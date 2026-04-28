test_that("L1 GP landscape generator smoke test", {
  x <- simulate_l1_gp_landscape(
    n_chr = 3,
    min_cn = 1,
    max_cn = 4,
    seed = 1
  )
  expect_equal(length(x$fitness), 4^3)
  expect_equal(nrow(x$karyotypes), 4^3)
  expect_equal(length(x$labels), 4^3)
  expect_s3_class(x, "alfak2_landscape")
})

test_that("L1 GP landscape anchors diploid fitness", {
  x <- simulate_l1_gp_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = 1)
  dip <- which(rowSums(abs(x$karyotypes - x$diploid_cn)) == 0)
  expect_length(dip, 1)
  expect_equal(x$fitness[dip], x$diploid_fitness, tolerance = 1e-12)
})

test_that("L1 GP landscape stays inside requested range", {
  x <- simulate_l1_gp_landscape(
    n_chr = 3,
    min_cn = 1,
    max_cn = 4,
    lower = -3,
    upper = 4,
    seed = 1
  )
  expect_gte(min(x$fitness), x$lower - 1e-12)
  expect_lte(max(x$fitness), x$upper + 1e-12)
})

test_that("L1 GP landscape is reproducible by seed", {
  x1 <- simulate_l1_gp_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = 123)
  x2 <- simulate_l1_gp_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = 123)
  x3 <- simulate_l1_gp_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = 124)
  expect_equal(x1$fitness, x2$fitness)
  expect_false(isTRUE(all.equal(x1$fitness, x3$fitness)))
})

test_that("L1 GP landscape table is consistent", {
  x <- simulate_l1_gp_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = 1)
  expect_equal(x$table$fitness, x$fitness)
  expect_equal(x$table$label, x$labels)
  chr_cols <- paste0("chr", seq_len(x$n_chr))
  expect_equal(as.matrix(x$table[, chr_cols]), x$karyotypes)
})

test_that("L1 GP landscape can skip long table", {
  x <- simulate_l1_gp_landscape(
    n_chr = 3,
    min_cn = 1,
    max_cn = 4,
    seed = 1,
    include_table = FALSE
  )
  expect_null(x$table)
  expect_equal(length(x$fitness), 4^3)
  expect_equal(nrow(x$karyotypes), 4^3)
  expect_equal(length(x$labels), 4^3)
})

test_that("L1 GP landscape is compatible with sparse simulator", {
  land <- simulate_l1_gp_landscape(
    n_chr = 3,
    min_cn = 1,
    max_cn = 3,
    seed = 1,
    include_table = FALSE
  )
  sim <- simulate_sparse_counts(
    land,
    n0 = 50,
    n1 = 60,
    dt = 0.2,
    seed = 2
  )
  expect_equal(ncol(sim$counts), 2)
  expect_true(all(sim$counts >= 0))
  expect_gt(sim$sparsity$observed_nodes, 0)
})

test_that("default L1 GP benchmark size is correct", {
  skip_on_cran()
  x <- simulate_l1_gp_landscape(seed = 1, include_table = FALSE)
  expect_equal(length(x$fitness), 262144)
  expect_equal(nrow(x$karyotypes), 262144)
})
