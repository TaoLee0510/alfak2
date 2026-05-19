test_that("karyotype parsing round-trips", {
  labels <- c("2.2.1", "3.2.1", "2.3.1")
  mat <- alfak2:::parse_karyotypes(labels)
  expect_equal(as.character(alfak2:::format_karyotypes(mat)), labels)
})

test_that("shell neighbor graph contains expected nodes and edges", {
  counts <- matrix(c(10, 5), nrow = 1, dimnames = list("2.2", c("t0", "t1")))
  dat <- prepare_alfak2_data(counts, beta = 0.01)
  g <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3)
  expect_true(all(c("1.2", "2.1", "2.2", "2.3", "3.2") %in% g$labels))
  expect_gt(length(g$edge_from), 0)
  expect_true(any(g$support_tier == "local_borrowed"))
})

test_that("exact transition kernel uses chromosome-level probabilities", {
  g_exact <- build_karyotype_graph(c("2.2", "3.2"), beta = 0.1,
                                   transition_kernel = "exact",
                                   shell_depth = 0, min_cn = 1, max_cn = 3)
  g_linear <- build_karyotype_graph(c("2.2", "3.2"), beta = 0.1,
                                    transition_kernel = "linear",
                                    shell_depth = 0, min_cn = 1, max_cn = 3)

  expect_equal(g_exact$transition_kernel, "exact")
  expect_equal(g_linear$transition_kernel, "linear")
  expect_equal(g_exact$edge_weight[1], alfak2:::alfak2_pij_cpp(2, 3, 0.1) * alfak2:::alfak2_pij_cpp(2, 2, 0.1))
  expect_lt(g_exact$edge_weight[1], g_linear$edge_weight[1])
  expect_equal(as.numeric(tapply(g_exact$transition_weight, g_exact$transition_from0, sum)), c(1, 1))
})

test_that("zero-weight rows are retained as targets but not graph seeds", {
  counts <- matrix(
    c(10, 5,
      4, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("2.2", "2.3"), c("t0", "t1"))
  )
  weights <- matrix(1, nrow = 2, ncol = 2, dimnames = dimnames(counts))
  weights["2.3", ] <- 0
  attr(counts, "observation_weights") <- weights
  dat <- prepare_alfak2_data(counts, beta = 0.01)
  g <- build_karyotype_graph(dat, shell_depth = 1, min_cn = 1, max_cn = 3)

  expect_equal(g$support_tier[match("2.3", g$labels)], "local_borrowed")
  expect_false("3.3" %in% g$labels)
})
