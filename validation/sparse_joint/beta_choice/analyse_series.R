#!/usr/bin/env Rscript
# Choosing the missegregation rate, step 2: does the chosen rate keep the
# simulated sets of validation/sparse_joint/series_methods (build cv10)
# where they were? Their true rate is 5e-5, the default, so the fixed rate is
# the truth there.
#
# `ends_auto` and `anchored_auto` are composed by series_methods/run.R
# --stage=compose from one fit per candidate rate (the grid of
# beta = "auto"), and scored by --stage=score. Criterion, fixed before any
# score was seen (plan of 2026-11-15, memory alfak2-beta-auto-plan), for
# ends_auto against ends and anchored_auto against anchored, with the RMSE
# ratio auto / fixed paired by dataset:
#   e1: in at least 18 of the 24 cells (bottleneck x spacing x depth x time
#       points) the median ratio is at most 1.03 at every distance 0-2, and
#       the median over cells at most 1.02 at every distance;
#   e2, e3: the median ratio at most 1.03 at every distance.
#
#   Rscript validation/sparse_joint/beta_choice/analyse_series.R
#
# Writes validation/sparse_joint/beta_choice/results/series_report.md.

root <- "validation/sparse_joint/series_methods/results/cv10"
out <- "validation/sparse_joint/beta_choice/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pairs <- list(ends = "ends_auto", anchored = "anchored_auto")
cells <- c("bottleneck", "spacing", "depth", "samples")

read_set <- function(set) {
  x <- utils::read.csv(file.path(root, set, "metrics_distance.csv"), stringsAsFactors = FALSE)
  x$set <- set
  x[x$stratum %in% c("d0", "d1", "d2") & x$status == "ok", ]
}
ratios <- do.call(rbind, lapply(c("e1", "e2", "e3"), function(set) {
  x <- read_set(set)
  do.call(rbind, lapply(names(pairs), function(fixed) {
    a <- x[x$method == pairs[[fixed]], c("dataset_id", "stratum", "rmse", intersect(cells, names(x)))]
    f <- x[x$method == fixed, c("dataset_id", "stratum", "rmse")]
    m <- merge(a, f, by = c("dataset_id", "stratum"), suffixes = c("_auto", "_fixed"))
    if (!nrow(m)) return(NULL)
    m$ratio <- m$rmse_auto / m$rmse_fixed
    m$set <- set
    m$pair <- paste(pairs[[fixed]], "vs", fixed)
    m
  }))
}))

verdict <- do.call(rbind, lapply(unique(ratios$pair), function(p) {
  r <- ratios[ratios$pair == p, ]
  e1 <- r[r$set == "e1", ]
  per_cell <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples + stratum, e1, stats::median)
  worst <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples, per_cell, max)
  over_cells <- stats::aggregate(ratio ~ stratum, per_cell, stats::median)
  e1_pass <- sum(worst$ratio <= 1.03) >= 18 && all(over_cells$ratio <= 1.02)
  others <- stats::aggregate(ratio ~ set + stratum, r[r$set != "e1", ], stats::median)
  rest_pass <- all(others$ratio <= 1.03)
  data.frame(pair = p, e1_cells = sprintf("%d / %d", sum(worst$ratio <= 1.03), nrow(worst)),
             e1_over_cells = paste(sprintf("%.3f", over_cells$ratio[order(over_cells$stratum)]), collapse = "/"),
             e1_pass = e1_pass,
             e2 = paste(sprintf("%.3f", others$ratio[others$set == "e2"][order(others$stratum[others$set == "e2"])]),
                        collapse = "/"),
             e3 = paste(sprintf("%.3f", others$ratio[others$set == "e3"][order(others$stratum[others$set == "e3"])]),
                        collapse = "/"),
             e2_e3_pass = rest_pass, pass = e1_pass && rest_pass, stringsAsFactors = FALSE)
}))

# The rates chosen.
picks <- do.call(rbind, lapply(c("e1", "e2", "e3"), function(set) {
  do.call(rbind, lapply(unlist(pairs), function(m) {
    files <- list.files(file.path(root, set, "fits", m), pattern = "[.]rds$", full.names = TRUE)
    b <- vapply(files, function(f) {
      v <- readRDS(f)$beta_multiple
      if (is.null(v)) NA_real_ else v
    }, numeric(1))
    tab <- table(factor(b, levels = c(1 / 16, 1 / 4, 1, 4, 16, 64)))
    data.frame(set = set, method = m, t(as.integer(tab)), check.names = FALSE)
  }))
}))
names(picks)[-(1:2)] <- c("1/16", "1/4", "1", "4", "16", "64")

md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
by_set <- stats::aggregate(ratio ~ pair + set + stratum, ratios, stats::median)
by_set <- reshape(by_set, idvar = c("pair", "set"), timevar = "stratum", direction = "wide")
writeLines(c("# Choosing the missegregation rate: step 2 (series_methods cv10, true rate = default)", "",
             "## Verdict (criterion fixed before scoring); ratios are RMSE auto / fixed at distances 0/1/2", "",
             md(verdict), "", "## Median ratio by set", "", md(by_set), "",
             "## Rates chosen (multiples of 5e-5)", "", md(picks)),
           file.path(out, "series_report.md"))
cat(readLines(file.path(out, "series_report.md")), sep = "\n")
