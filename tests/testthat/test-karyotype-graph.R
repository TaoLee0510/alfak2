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
