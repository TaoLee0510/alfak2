#!/usr/bin/env Rscript
# Summary of laplace_mean.R: posterior mode, mean and median against the
# truth, by kind of karyotype, and the coverage of the sampled 95% intervals.
#
#   Rscript validation/sparse_joint/unobserved_bias/laplace_mean_report.R [--source=rugged|diagnosis|dev]

if (!exists("ub_script_dir")) ub_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})

ub_laplace_report <- function(dir) {
  files <- list.files(dir, "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    d <- r$domain
    o <- r$outside
    kind <- c(ifelse(d$major, "major", ifelse(d$observed, "observed", "unobserved_in_domain")),
              paste0("outside_shell", pmin(o$shell, 3)))
    pick <- function(name) c(if (is.null(d[[name]])) rep(NA, nrow(d)) else d[[name]],
                             if (is.null(o[[name]])) rep(NA, nrow(o)) else o[[name]])
    data.frame(dataset_id = r$dataset_id, kind = kind, truth = pick("truth"), mode = pick("mode"),
               mean = pick("mean"), median = pick("median"), low = pick("low"), high = pick("high"),
               stringsAsFactors = FALSE)
  })
  x <- do.call(rbind, rows)
  x$covered <- x$truth >= x$low & x$truth <= x$high
  per <- do.call(rbind, lapply(split(x, list(x$dataset_id, x$kind), drop = TRUE), function(g) {
    rmse <- function(p) sqrt(mean((p - g$truth)^2))
    data.frame(dataset_id = g$dataset_id[1L], kind = g$kind[1L], truth = mean(g$truth), mode = mean(g$mode),
               mean = mean(g$mean), median = mean(g$median), rmse_mode = rmse(g$mode), rmse_mean = rmse(g$mean),
               rmse_median = rmse(g$median), covered = mean(g$covered), stringsAsFactors = FALSE)
  }))
  diag <- do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(dataset_id = r$dataset_id, n_latent = r$n_latent, acceptance = r$acceptance,
               ess = r$ess_unobserved_mean, seconds = r$seconds, mode_matches_fit = r$mode_matches_fit)
  }))
  list(per = per, diag = diag)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  src <- sub("^--source=", "", c(grep("^--source=", args, value = TRUE), "--source=rugged")[[1L]])
  dir <- file.path(ub_script_dir, "results", src, "laplace_mean")
  index <- utils::read.csv(file.path(ub_script_dir, "results", src, "dataset_index.csv"), stringsAsFactors = FALSE)
  out <- ub_laplace_report(dir)
  per <- merge(index[c("dataset_id", "wavelength_class", "interval")], out$per)
  options(width = 220)
  cat("Chains:", nrow(out$diag), " acceptance", paste(round(range(out$diag$acceptance), 2), collapse = "-"),
      " ESS of the unobserved mean, median", round(stats::median(out$diag$ess)), "min", round(min(out$diag$ess)),
      " seconds, median", round(stats::median(out$diag$seconds)), " mode = fit:", all(out$diag$mode_matches_fit), "\n\n")
  a <- stats::aggregate(cbind(truth, mode, mean, median, rmse_mode, rmse_mean, rmse_median, covered) ~
                          wavelength_class + kind, per, mean, na.action = stats::na.pass)
  a$ratio_mode <- a$mode / a$truth
  a$ratio_mean <- a$mean / a$truth
  a$ratio_median <- a$median / a$truth
  a$rmse_mean_change <- a$rmse_mean / a$rmse_mode - 1
  a$rmse_median_change <- a$rmse_median / a$rmse_mode - 1
  keep <- c("wavelength_class", "kind", "truth", "ratio_mode", "ratio_mean", "ratio_median", "rmse_mode",
            "rmse_mean_change", "rmse_median_change", "covered")
  a <- a[keep]
  a[-(1:2)] <- lapply(a[-(1:2)], round, 3)
  print(a[order(a$wavelength_class, a$kind), ], row.names = FALSE)
}
