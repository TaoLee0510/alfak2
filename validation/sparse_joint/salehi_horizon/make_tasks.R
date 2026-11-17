#!/usr/bin/env Rscript
# Task table of forecast_score.R: every Salehi trajectory of the paper
# benchmark ("ends" design), for alfak2 (round 05a8855, the default
# initial_state = "balance"), alfakR (the same round's links) and alfak2 with
# initial_state = "free" (the alfak2-salehi-variants project), at minobs 5,
# 10 and 20. A task without a fit or without descendants records that and
# stops.
#   Rscript make_tasks.R OUT.tsv   (from the paper benchmark's folder)
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
ids <- sort(salehi_trajectories(paths)$traj_id)
roots <- c(alfak2 = "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits/salehi",
           alfakR = "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits/salehi",
           alfak2_free = "/share/lab_crd/taoli/Project/alfak2-salehi-variants/results/fits/salehi")
tasks <- expand.grid(id = ids, method = names(roots), minobs = c(5, 10, 20), stringsAsFactors = FALSE)
tasks$root <- roots[tasks$method]
tasks$pm <- bench_settings$salehi$pm
tasks$nb <- bench_settings$salehi$nb
utils::write.table(tasks, commandArgs(TRUE)[1], sep = "\t", quote = FALSE, row.names = FALSE)
cat("tasks", nrow(tasks), "\n")
