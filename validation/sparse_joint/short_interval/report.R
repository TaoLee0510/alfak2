#!/usr/bin/env Rscript
# Tables of the short-interval diagnosis from metrics_fit.csv and
# metrics_distance.csv (run.R --stage=score).
#
#   Rscript validation/sparse_joint/short_interval/report.R [--methods=a,b] [--output-dir=DIR]
#
# Cells are first-sample day (t0) x interval x wavelength class, five
# landscapes each; every number is a median over the datasets of a cell.

si_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(si_script_dir, "setup.R"))

si_report <- function(output_dir = si_default_output, methods = NULL, baseline = "balance") {
  fits <- utils::read.csv(file.path(output_dir, "metrics_fit.csv"), stringsAsFactors = FALSE)
  distance <- utils::read.csv(file.path(output_dir, "metrics_distance.csv"), stringsAsFactors = FALSE)
  if (is.null(methods)) methods <- setdiff(unique(fits$method), "zero")
  fits <- fits[fits$method %in% methods, ]
  options(width = 250)
  med <- function(x) round(stats::median(x, na.rm = TRUE), 3)
  cell <- c("t0", "interval", "wavelength_class")
  by_cell <- function(d, columns, keys = c("method", cell)) {
    groups <- split(d, d[keys], drop = TRUE)
    out <- do.call(rbind, lapply(groups, function(g) {
      cbind(g[1L, keys, drop = FALSE], n = nrow(g), as.data.frame(lapply(g[columns], med)))
    }))
    out[do.call(order, unname(as.list(out[keys]))), , drop = FALSE]
  }
  report_dir <- file.path(output_dir, "report")
  dir.create(report_dir, showWarnings = FALSE)
  show <- function(title, table, file) {
    cat("\n==", title, "\n")
    print(table, row.names = FALSE)
    utils::write.csv(table, file.path(report_dir, file), row.names = FALSE)
  }

  failed <- fits[fits$status != "ok", c("dataset_id", "method", "status")]
  cat("Fits:", nrow(fits), " failed:", nrow(failed), "\n")
  if (nrow(failed)) print(failed, row.names = FALSE)

  show("A. Held-out prediction of frequent karyotypes and the landscape against the truth",
       by_cell(fits, c("n_folds", "heldout_slope", "heldout_r", "full_slope", "full_r", "bias_fq", "between",
                       "landscape_r", "landscape_R2")), "cv_by_cell.csv")
  show("B. Fit: trend scales, bias of observed karyotypes (truth-scale, per day), convergence",
       by_cell(fits, c("scale_dosage", "scale_curvature", "scale_residual", "bias_major",
                       "bias_observed_nonmajor", "history", "n_major", "fit_seconds")), "fit_by_cell.csv")
  cat("Non-zero convergence codes:\n")
  print(table(method = fits$method, convergence = fits$convergence))

  keep <- distance$stratum %in% c("d0", "d1", "d2") & distance$method %in% c("zero", methods)
  d <- distance[keep, ]
  d$rmse_ratio <- d$rmse / d$rmse_zero
  show("C. Per-distance scores of the development benchmark (slope = prediction on truth)",
       by_cell(d, c("rmse_ratio", "slope", "spearman", "coverage95"), c("method", "stratum", "t0", "interval")),
       "distance_by_interval.csv")
  show("C2. The same by wavelength class, pooled over t0 and interval",
       by_cell(d, c("rmse_ratio", "slope", "spearman", "coverage95"), c("method", "stratum", "wavelength_class")),
       "distance_by_wavelength.csv")

  # Paired with the baseline on the same dataset.
  if (baseline %in% fits$method) {
    base <- fits[fits$method == baseline, ]
    others <- setdiff(methods, baseline)
    paired <- do.call(rbind, lapply(others, function(m) {
      x <- merge(base, fits[fits$method == m, ], by = c("dataset_id", cell), suffixes = c(".base", ".arm"))
      if (!nrow(x)) return(NULL)
      dd <- merge(d[d$method == baseline, c("dataset_id", "stratum", "rmse")],
                  d[d$method == m, c("dataset_id", "stratum", "rmse")], by = c("dataset_id", "stratum"))
      rmse_change <- function(s) {
        y <- dd[dd$stratum == s, ]
        stats::setNames(y$rmse.y / y$rmse.x - 1, y$dataset_id)
      }
      r1 <- rmse_change("d1")
      r2 <- rmse_change("d2")
      x$slope_gap_change <- abs(1 - x$heldout_slope.arm) - abs(1 - x$heldout_slope.base)
      x$heldout_r_change <- x$heldout_r.arm - x$heldout_r.base
      x$landscape_r_change <- x$landscape_r.arm - x$landscape_r.base
      x$d1_rmse_change <- r1[x$dataset_id]
      x$d2_rmse_change <- r2[x$dataset_id]
      x$method <- m
      x
    }))
    if (!is.null(paired)) {
      columns <- c("slope_gap_change", "heldout_r_change", "landscape_r_change", "d1_rmse_change", "d2_rmse_change")
      show(sprintf("D. Paired change against %s by t0 and interval (slope_gap = |1 - held-out slope|; rmse as ratio - 1)",
                   baseline),
           by_cell(paired, columns, c("method", "t0", "interval")), "paired_by_interval.csv")
      show("D2. The same by wavelength class", by_cell(paired, columns, c("method", "wavelength_class")),
           "paired_by_wavelength.csv")
    }
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else NULL
  }
  methods <- value("methods")
  si_report(if (is.null(value("output-dir"))) si_default_output else value("output-dir"),
            if (is.null(methods)) NULL else strsplit(methods, ",", fixed = TRUE)[[1L]])
}
