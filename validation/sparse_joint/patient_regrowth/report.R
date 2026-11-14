#!/usr/bin/env Rscript
# Tables of the patient regrowth study (README.md), from results/metrics_*.csv
# written by `run.R --stage=score`. Landscapes are the unit: each cell is the
# median over the 12 development landscapes (patient 1), and paired
# differences are taken landscape by landscape.
#
#   Rscript validation/sparse_joint/patient_regrowth/report.R [--build=BUILD]
#
# Writes results/report.md, or results/BUILD/report.md for a build other
# than ff13cf4.

pr_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
build_arg <- grep("^--build=", commandArgs(TRUE), value = TRUE)
build <- if (length(build_arg)) sub("^--build=", "", build_arg[[1L]]) else "ff13cf4"
output_dir <- if (build == "ff13cf4") file.path(pr_script_dir, "results") else file.path(pr_script_dir, "results", build)
distance <- utils::read.csv(file.path(output_dir, "metrics_distance.csv"), stringsAsFactors = FALSE)
fits <- utils::read.csv(file.path(output_dir, "metrics_fit.csv"), stringsAsFactors = FALSE)
scenarios <- c("none", "s2000", "s0200", "s0020", "s0200_d50", "s0200_d75", "s0200_d875", "s0020_d875")
strata <- c(d0_ge5 = "observed, 5+ cells", d0_lt5 = "observed, <5 cells", d1 = "one step", d2 = "two steps")
distance <- distance[distance$status == "ok" & distance$stratum %in% names(strata), ]

lines <- character()
add <- function(...) lines <<- c(lines, ...)
table_md <- function(x, digits = 2L) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "f") else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
cell <- function(method, metric, stratum, scenario) {
  v <- distance[[metric]][distance$method == method & distance$stratum == stratum & distance$scenario == scenario]
  if (!length(v)) NA_real_ else stats::median(v, na.rm = TRUE)
}
grid <- function(method, metric, rows = scenarios) {
  out <- data.frame(scenario = rows, stringsAsFactors = FALSE)
  for (s in names(strata)) out[[s]] <- vapply(rows, function(r) cell(method, metric, s, r), numeric(1))
  out
}

add("# Patient regrowth: tables", "",
    sprintf("Generated %s from `metrics_distance.csv` (%d rows) and `metrics_fit.csv`.", format(Sys.time(), "%Y-%m-%d %H:%M"),
            nrow(distance)), "",
    "Truth is per calendar day (rate while growing x growing share). Medians over 12 landscapes.", "")

status <- stats::aggregate(dataset_id ~ method + status, fits, length)
names(status)[3] <- "fits"
add("## Fits", "", table_md(status, 0L), "")
failed <- fits[fits$status != "ok", intersect(c("dataset_id", "method", "detail"), names(fits)), drop = FALSE]
if (nrow(failed)) add("Failures:", "", table_md(failed), "")

for (metric in c("skill", "slope", "spearman", "coverage95")) {
  add(sprintf("## %s default (`base`): %s", build, metric), "", table_md(grid("base", metric)), "")
}

# Founder effect on the majors: RMSE of the fit's major karyotypes relative to
# the spread of their true fitness.
major <- fits[fits$method == "base" & fits$status == "ok", ]
major$major_ratio <- major$rmse_major / major$rmse_major_zero
add("## Majors (`base`): RMSE / RMSE of predicting no difference", "",
    table_md(data.frame(scenario = scenarios,
                        ratio = vapply(scenarios, function(s) stats::median(major$major_ratio[major$scenario == s],
                                                                            na.rm = TRUE), numeric(1)),
                        majors = vapply(scenarios, function(s) stats::median(major$n_major_scored[major$scenario == s]),
                                        numeric(1)))), "")

# Paired differences against `none`, landscape by landscape.
paired <- function(method, metric, baseline_method = "base", baseline_scenario = "none") {
  out <- data.frame(scenario = setdiff(scenarios, baseline_scenario), stringsAsFactors = FALSE)
  for (s in names(strata)) {
    out[[s]] <- vapply(out$scenario, function(r) {
      a <- distance[distance$method == method & distance$stratum == s & distance$scenario == r, ]
      b <- distance[distance$method == baseline_method & distance$stratum == s & distance$scenario == baseline_scenario, ]
      m <- merge(a[c("landscape_id", metric)], b[c("landscape_id", metric)], by = "landscape_id")
      d <- m[[2L]] - m[[3L]]
      sprintf("%+.2f (%d/%d worse)", stats::median(d, na.rm = TRUE), sum(d < 0, na.rm = TRUE), sum(is.finite(d)))
    }, character(1))
  }
  out
}
add("## Skill against `none`, paired by landscape (`base`)", "",
    "Median difference; in brackets, landscapes where the scenario is worse.", "", table_md(paired("base", "skill")), "")
add("## Rank correlation against `none`, paired by landscape (`base`)", "", table_md(paired("base", "spearman")), "")

# Dormant scenarios: what not knowing the dormant share costs.
dormant <- c("s0200", "s0200_d50", "s0200_d75", "s0200_d875", "s0020", "s0020_d875")
for (metric in c("skill", "slope", "spearman", "coverage95")) {
  rows <- list()
  for (method in c("base", "equivalent", "rates")) {
    g <- grid(method, metric, dormant)
    g$scenario <- paste(g$scenario, method)
    rows[[method]] <- g
  }
  add(sprintf("## Dormancy: %s by arm", metric), "",
      "`base` assumes one division per day, `equivalent` the growing share of it, `rates` the range 0.1-1 (point at 0.1).", "",
      table_md(do.call(rbind, rows)), "")
}

# Cohort: does fitting the 8 patients of a landscape reduce the cost of a
# bottleneck? Difference in differences against `none`.
cohort_rows <- list()
for (method in c("base", "cohort")) {
  g <- grid(method, "skill", c("none", "s0200", "s0020"))
  g$scenario <- paste(g$scenario, method)
  cohort_rows[[method]] <- g
}
add("## Cohort of 8 patients: skill", "", table_md(do.call(rbind, cohort_rows)), "")
add("## Cohort of 8 patients: skill against `none`, paired", "",
    table_md(rbind(paired("cohort", "skill", "cohort")[c(2, 3), ])), "")

writeLines(lines, file.path(output_dir, "report.md"))
cat("Wrote", file.path(output_dir, "report.md"), "\n")
