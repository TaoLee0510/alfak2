#!/usr/bin/env Rscript
# Follow-up step 2, read only: what the fits of the "all" design (every
# passage of a Salehi trajectory) recorded, for the round 05a8855 (anchored,
# pairwise, state_space) and the alfak2-salehi-variants project (free):
# samples, samples fitted, weight of the time points in between and the
# majors' dispersion it came from, history, objective, convergence.
#   Rscript all_design_runs.R OUT.csv
roots <- c(round = "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits/salehi",
           variants = "/share/lab_crd/taoli/Project/alfak2-salehi-variants/results/fits/salehi")
rows <- list()
for (r in names(roots)) {
  for (f in Sys.glob(file.path(roots[[r]], "*", "all", "alfak2*", "minobs_5", "run.Rds"))) {
    x <- readRDS(f)
    parts <- strsplit(f, "/", fixed = TRUE)[[1L]]
    n <- length(parts)
    v <- function(name) { y <- x[[name]]; if (is.null(y) || !length(y)) NA else as.numeric(y)[[1L]] }
    rows[[length(rows) + 1L]] <- data.frame(
      source = r, id = parts[n - 4L], method = parts[n - 2L], status = x$status,
      samples = v("samples"), samples_fitted = v("samples_fitted"), history = v("history"),
      intermediate_weight = v("intermediate_weight"), intermediate_dispersion = v("intermediate_dispersion"),
      drift = v("drift"), objective = v("objective"), convergence = v("convergence"),
      coverage_first = if (is.null(x$domain_coverage)) NA else x$domain_coverage$fraction[[1L]],
      coverage_last = if (is.null(x$domain_coverage)) NA else x$domain_coverage$fraction[[nrow(x$domain_coverage)]],
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
utils::write.csv(out, commandArgs(TRUE)[1], row.names = FALSE)
print(table(out$source, out$method, out$status))
