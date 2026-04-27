test_that("legacy exports are absent", {
  ns <- getNamespaceExports("alfak2")
  expect_false(paste0("run_", "alfak2_", "bench", "mark") %in% ns)
  retired_exports <- c(
    paste0("al", "fak"),
    paste0("predict", "_", "evo"),
    paste0("al", "fak_", "cohort_", "transition")
  )
  retired_patterns <- paste(
    paste0("nn", "_", "prior"),
    paste0("cohort", "_", "transition"),
    paste0("K", "rig"),
    sep = "|"
  )
  expect_false(any(retired_exports %in% ns))
  expect_false(any(grepl(retired_patterns, ns)))
})
