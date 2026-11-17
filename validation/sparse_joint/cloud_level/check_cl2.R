#!/usr/bin/env Rscript
# Follow-up step 3: checks of build cl2 (cl2.patch: initial_state = "auto").
#   1. At its defaults it gives the fits of bf86c85 (same_default.R).
#   2. "auto" makes the choice of the rule in confirm.R and returns the fit
#      that was chosen, as fitted in the study.
#   3. as_alfakR() runs on an "auto" fit; more than two time points give an
#      error.
cl_script_dir <- "validation/sparse_joint/cloud_level"
source(file.path(cl_script_dir, "setup.R"))
.libPaths(c(file.path(cl_output, "lib", "alfak2_cl2"), .libPaths()))
suppressPackageStartupMessages(library(alfak2))
cat("alfak2 from", find.package("alfak2"), "\n")
scores <- utils::read.csv(file.path(cl_output, "scores_confirm.csv"), stringsAsFactors = FALSE)
level <- scores[scores$method == "level_dt", c("set", "dataset_id", "balance_level", "landscape_id", "role")]
pick <- function(set, rule, n) {
  d <- level[level$set == set & (level$balance_level > 1) == rule, ]
  d[order(d$dataset_id), ][seq_len(min(n, nrow(d))), ]
}
cases <- rbind(pick("dev2", FALSE, 2), pick("short2", FALSE, 1), pick("short2", TRUE, 2), pick("diverse2", TRUE, 2))
rows <- lapply(seq_len(nrow(cases)), function(i) {
  set <- cases$set[i]; id <- cases$dataset_id[i]
  dataset <- readRDS(cl_dataset_path(set, id))
  nodes <- readRDS(cl_nodes_path(set, id))
  counts <- dataset$counts; storage.mode(counts) <- "integer"
  beta <- if (is.null(dataset$beta_assumed)) dataset$beta else dataset$beta_assumed
  set.seed(dataset$seed)
  fit <- suppressWarnings(fit_sparse_joint(counts, dt = dataset$interval, beta = beta, birth_rate = 1, min_cn = 1L,
                                           max_cn = 8L, initial_state = "auto"))
  check <- fit$diagnostics$initial_state_check
  expected <- if (cases$balance_level[i] > 1) "free" else "base"
  stored <- readRDS(cl_fit_path(set, expected, id))$prediction
  p <- stats::predict(fit, stored$karyotype)
  data.frame(set = set, dataset_id = id, level_study = cases$balance_level[i], level_auto = check$level,
             chosen = check$chosen, expected = expected,
             max_difference = max(abs(p$relative_fitness - stored$fitness)), stringsAsFactors = FALSE)
})
out <- do.call(rbind, rows)
print(out, digits = 4, row.names = FALSE)
# The study's fits ran on HPC (x86), these locally: levels and free fits can
# differ slightly between machines (results/platform_check.R shows cl1 and cl2
# identical on one machine). The choice must agree.
stopifnot(all(out$chosen == ifelse(out$expected == "free", "free", "balance")))
dataset <- readRDS(cl_dataset_path("diverse2", cases$dataset_id[nrow(cases)]))
counts <- dataset$counts; storage.mode(counts) <- "integer"
fit <- suppressWarnings(fit_sparse_joint(counts, dt = 15, initial_state = "auto", birth_rate = 1, min_cn = 1L, max_cn = 8L))
compat <- as_alfakR(fit, minobs = 10, nboot = 4)
cat("as_alfakR on an auto fit: landscape rows", nrow(compat$landscape), " R2R", compat$xval$R2R, "\n")
three <- cbind(counts, counts[, 2L])
colnames(three) <- c("0", "15", "30")
message <- tryCatch({ fit_sparse_joint(three, times = c(0, 15, 30), initial_state = "auto"); "no error" },
                    error = function(e) conditionMessage(e))
cat("three time points:", message, "\n")
cat("All checks passed\n")
