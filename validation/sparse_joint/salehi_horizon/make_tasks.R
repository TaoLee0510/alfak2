#!/usr/bin/env Rscript
# Task table of forecast_score.R: every Salehi trajectory of the paper
# benchmark ("ends" design), for alfak2 (round 05a8855, the default
# initial_state = "balance"), alfakR (the same round's links) and alfak2 with
# initial_state = "free" (the alfak2-salehi-variants project), at minobs 5,
# 10 and 20. A task without a fit or without descendants records that and
# stops.
#   Rscript make_tasks.R OUT.tsv [all]   (from the paper benchmark's folder)
# With "all": the fits of every passage of a trajectory (trajectories of three
# or more passages): alfak2 anchored (the default), pairwise and state_space,
# alfakR, and alfak2 "free".
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
design <- if (length(args) > 1L) args[[2]] else "ends"
ids <- sort(salehi_trajectories(paths)$traj_id)
round <- "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits/salehi"
variants <- "/share/lab_crd/taoli/Project/alfak2-salehi-variants/results/fits/salehi"
roots <- c(alfak2 = round, alfakR = round, alfak2_free = variants)
if (design == "all") {
  roots <- c(alfak2 = round, alfak2_pairwise = round, alfak2_state_space = round, alfakR = round,
             alfak2_free = variants)
  ids <- ids[as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", ids)) >= 3L]
}
tasks <- expand.grid(id = ids, method = names(roots), minobs = c(5, 10, 20), stringsAsFactors = FALSE)
tasks$design <- design
tasks$root <- roots[tasks$method]
tasks$pm <- bench_settings$salehi$pm
tasks$nb <- bench_settings$salehi$nb
utils::write.table(tasks, args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("tasks", nrow(tasks), "\n")
