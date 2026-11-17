#!/usr/bin/env Rscript
# Corrected-rule study, step 0 (1): features of every two-time-point ABM and
# missegregation input of the ALFA-K paper benchmark, read only. From round
# 05a8855 (default, balance prior): interval, cells, majors, coverage, chosen
# history and its profile. From round 827c4b5: the level of the "auto" check
# and its convergence. The effect of switching against the truth, where it
# happened, is the benchmark's vs_truth_fits.tsv.
#   Rscript paper_features.R OUT.csv
bench <- "/share/lab_crd/taoli/Project/alfak2-paper-benchmark"
rows <- list()
for (src in c("abm", "misseg")) {
  # One fit per input: the fit does not depend on minobs (missegregation
  # inputs have minobs 10 only).
  dirs <- Sys.glob(file.path(bench, "results_05a8855", "fits", src, "*", "w*_ends", "alfak2", "minobs_*"))
  dirs <- dirs[!duplicated(dirname(dirs))]
  for (dir in dirs) {
    parts <- strsplit(dir, "/", fixed = TRUE)[[1L]]
    n <- length(parts)
    id <- parts[n - 3L]; design <- parts[n - 2L]
    fit_path <- file.path(dir, "alfak2_fit.Rds")
    new_run <- file.path(bench, "results_827c4b5", "fits", src, id, design, "alfak2", basename(dir), "run.Rds")
    if (!file.exists(fit_path)) next
    fit <- readRDS(fit_path)
    d <- fit$diagnostics
    counts <- fit$data$counts
    prof <- d$history_profile
    gap <- if (is.null(prof) || nrow(prof) < 2L) NA_real_ else diff(sort(prof$objective)[1:2])
    nr <- if (file.exists(new_run)) readRDS(new_run) else NULL
    chk <- if (is.null(nr)) NULL else nr$initial_state_check
    s <- fit$summary
    total <- rowSums(as.matrix(s[grep("^count_t", names(s))]))
    rows[[length(rows) + 1L]] <- data.frame(
      source = src, dataset = id, design = design, dt = fit$data$dt,
      cells_t0 = sum(counts[, 1L]), cells_t1 = sum(counts[, 2L]), karyotypes = nrow(counts),
      singleton_t0 = sum(counts[counts[, 1L] == 1, 1L]) / sum(counts[, 1L]),
      top_t0 = max(counts[, 1L]) / sum(counts[, 1L]),
      coverage_t0 = d$domain_coverage$fraction[[1L]], coverage_t1 = d$domain_coverage$fraction[[2L]],
      n_major = d$n_major, history = d$history, history_gap = gap, objective = d$objective,
      scale_residual = fit$trend$scales[["residual"]],
      shell_coef = fit$trend$coefficients[[length(fit$trend$coefficients)]],
      cloud_share_t0 = sum(s$count_t0[!s$major]) / sum(s$count_t0),
      level = if (is.null(chk)) NA_real_ else chk$level, level_sd = if (is.null(chk)) NA_real_ else chk$level_sd,
      check_convergence = if (is.null(chk)) NA_real_ else chk$check_convergence,
      chosen = if (is.null(chk)) NA_character_ else chk$chosen, stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
utils::write.csv(out, commandArgs(TRUE)[1], row.names = FALSE)
print(table(out$source, out$design))
