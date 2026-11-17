#!/usr/bin/env Rscript
# Tables of stage 2 (README.md), from results/<set>/metrics_*.csv written by
# `run.R --stage=score`. A cell is one bottleneck, spacing, number of time
# points and depth; landscapes are the unit, and paired ratios are taken
# landscape by landscape.
#
#   Rscript validation/sparse_joint/multi_timepoint/report.R [--set=series]
#
# Gates of the plan: in every cell the median ratio of the RMSE on all time
# points to the first-and-last RMSE is at most 1.03 at distances 0, 1 and 2;
# the 95% intervals of the fit on all time points cover at least 90% of
# karyotypes, otherwise the cell is marked for recalibration or the drift
# model. Every multi-time-point arm present (`all`, `majors_any`, `majors1`,
# `majors`) is compared.

mt_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
args <- commandArgs(TRUE)
set <- if (any(grepl("^--set=", args))) sub("^--set=", "", grep("^--set=", args, value = TRUE)[[1L]]) else "series"
output_dir <- file.path(mt_script_dir, "results", set)
distance <- utils::read.csv(file.path(output_dir, "metrics_distance.csv"), stringsAsFactors = FALSE)
fits <- utils::read.csv(file.path(output_dir, "metrics_fit.csv"), stringsAsFactors = FALSE)
distance <- distance[distance$status == "ok" & distance$stratum %in% c("d0", "d1", "d2"), ]
cell_keys <- c("bottleneck", "spacing", "samples", "depth")
cells <- unique(distance[cell_keys])
cells <- cells[do.call(order, cells), ]

lines <- character()
add <- function(...) lines <<- c(lines, ...)
table_md <- function(x, digits = 2L) {
  # The `samples` column counts time points; the header says so.
  names(x)[names(x) == "samples"] <- "time points"
  x[] <- lapply(x, function(v) if (is.numeric(v) && any(v != round(v), na.rm = TRUE)) formatC(v, digits = digits, format = "f") else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
in_cell <- function(frame, cell) {
  keep <- rep(TRUE, nrow(frame))
  for (k in cell_keys) keep <- keep & frame[[k]] == cell[[k]]
  frame[keep, , drop = FALSE]
}
# Median over landscapes of method a's metric over method b's, and the number
# of landscapes where a is lower.
paired_ratio <- function(cell, stratum, a, b, metric = "rmse") {
  rows <- in_cell(distance, cell)
  rows <- rows[rows$stratum == stratum, ]
  m <- merge(rows[rows$method == a, c("landscape_id", metric)], rows[rows$method == b, c("landscape_id", metric)],
             by = "landscape_id")
  if (!nrow(m)) return(c(ratio = NA, better = NA, n = 0))
  r <- m[[2L]] / m[[3L]]
  c(ratio = stats::median(r, na.rm = TRUE), better = sum(r < 1, na.rm = TRUE), n = sum(is.finite(r)))
}
median_of <- function(cell, stratum, method, metric) {
  rows <- in_cell(distance, cell)
  v <- rows[[metric]][rows$method == method & rows$stratum == stratum]
  if (!length(v)) NA_real_ else stats::median(v, na.rm = TRUE)
}

add(sprintf("# Multi-timepoint stage 2 (%s): tables", set), "",
    sprintf("Generated %s. Medians over landscapes; ratios paired by landscape.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
status <- stats::aggregate(dataset_id ~ method + status, fits, length)
names(status)[3L] <- "fits"
add("## Fits", "", table_md(status, 0L), "")
failed <- fits[fits$status != "ok", intersect(c("dataset_id", "method", "detail"), names(fits))]
if (nrow(failed)) add("Failures:", "", table_md(utils::head(failed, 30)), "")

# Every multi-time-point arm against first and last (`ends`).
arms <- intersect(c("all", "majors_any", "majors1", "majors"), unique(fits$method))
summary_rows <- list()
for (arm in arms) {
  gate <- cells
  for (s in c("d0", "d1", "d2")) {
    r <- t(vapply(seq_len(nrow(cells)), function(i) paired_ratio(cells[i, ], s, arm, "ends"), numeric(3)))
    gate[[paste0(s, "_ratio")]] <- r[, "ratio"]
    gate[[paste0(s, "_better")]] <- sprintf("%d/%d", r[, "better"], r[, "n"])
  }
  complete <- is.finite(gate$d0_ratio) & is.finite(gate$d1_ratio) & is.finite(gate$d2_ratio)
  gate$rmse_gate <- ifelse(!complete, "NA",
                           ifelse(pmax(gate$d0_ratio, gate$d1_ratio, gate$d2_ratio) <= 1.03, "pass", "FAIL"))
  cover <- vapply(seq_len(nrow(cells)), function(i) {
    min(vapply(c("d0", "d1", "d2"), function(s) median_of(cells[i, ], s, arm, "coverage95"), numeric(1)))
  }, numeric(1))
  gate$coverage_min <- cover
  add(sprintf("## `%s` against first and last (`ends`): RMSE ratio", arm), "",
      "Ratio below 1: the multi-time-point fit is better. `better`: landscapes where it is better.",
      "`coverage_min`: the lowest median 95% coverage of the three distances.", "",
      table_md(gate, 3L), "",
      sprintf("Cells passing the RMSE gate: %d of %d; coverage at least 0.90 in %d.",
              sum(gate$rmse_gate == "pass"), sum(gate$rmse_gate != "NA"), sum(cover >= 0.90, na.rm = TRUE)), "")
  for (b in sort(unique(cells$bottleneck))) {
    g <- gate[gate$bottleneck == b, ]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      arm = arm, bottleneck = b, d0 = stats::median(g$d0_ratio, na.rm = TRUE),
      d1 = stats::median(g$d1_ratio, na.rm = TRUE), d2 = stats::median(g$d2_ratio, na.rm = TRUE),
      worst = max(pmax(g$d0_ratio, g$d1_ratio, g$d2_ratio), na.rm = TRUE),
      cells_pass = sprintf("%d/%d", sum(g$rmse_gate == "pass"), sum(g$rmse_gate != "NA")),
      coverage_ok = sprintf("%d/%d", sum(g$coverage_min >= 0.90, na.rm = TRUE), nrow(g)), stringsAsFactors = FALSE)
  }
}
summary_table <- do.call(rbind, summary_rows)
lines <- append(lines, c("## Summary: multi-time-point arms against `ends` by bottleneck", "",
                         "Median over cells of the median RMSE ratio; `worst`: the largest cell ratio.", "",
                         table_md(summary_table, 3L), ""), after = which(lines == "## Fits") - 1L)

coverage <- cells
for (method in c(arms, "ends")) for (s in c("d0", "d1", "d2")) {
  coverage[[paste(method, s)]] <- vapply(seq_len(nrow(cells)), function(i) median_of(cells[i, ], s, method, "coverage95"),
                                         numeric(1))
}
add("## 95% coverage", "", table_md(coverage, 3L), "")

skill <- cells
for (method in c(arms, "ends", "alfakR")) for (s in c("d0", "d1", "d2")) {
  skill[[paste(method, s)]] <- vapply(seq_len(nrow(cells)), function(i) median_of(cells[i, ], s, method, "skill"),
                                      numeric(1))
}
add("## Skill (1 - MSE / MSE of no difference)", "", table_md(skill, 2L), "")

spearman <- cells
for (method in c(arms, "ends", "alfakR")) for (s in c("d0", "d1", "d2")) {
  spearman[[paste(method, s)]] <- vapply(seq_len(nrow(cells)), function(i) median_of(cells[i, ], s, method, "spearman"),
                                         numeric(1))
}
add("## Rank correlation", "", table_md(spearman, 2L), "")

if (any(distance$method == "alfakR")) {
  versus <- cells[cells$bottleneck == max(cells$bottleneck), ]
  for (s in c("d0", "d1", "d2")) {
    r <- t(vapply(seq_len(nrow(versus)), function(i) paired_ratio(versus[i, ], s, "all", "alfakR"), numeric(3)))
    versus[[paste0(s, "_ratio")]] <- r[, "ratio"]
    versus[[paste0(s, "_better")]] <- sprintf("%d/%d", r[, "better"], r[, "n"])
  }
  add("## alfak2 all time points against alfakR all time points: RMSE ratio", "", table_md(versus, 3L), "")
}

# Dispersion of the majors' later counts about the fit on all time points, against
# the drift index K^2 N_s / (5 n0).
all_fits <- fits[fits$method == "all" & fits$status == "ok", ]
dispersion <- stats::aggregate(cbind(drift_index, dispersion) ~ bottleneck + spacing + samples + depth, all_fits,
                               stats::median)
add("## Dispersion of the majors (all time points) and the drift index", "", table_md(dispersion, 2L), "",
    sprintf("Spearman correlation of dispersion and drift index over fits: %.2f.",
            stats::cor(all_fits$dispersion, all_fits$drift_index, method = "spearman", use = "complete.obs")), "")

time <- stats::aggregate(seconds ~ method + samples + depth, fits[fits$status == "ok", ], stats::median)
add("## Seconds per fit (median)", "", table_md(time, 0L), "")

writeLines(lines, file.path(output_dir, "report.md"))
cat("Wrote", file.path(output_dir, "report.md"), "\n")
