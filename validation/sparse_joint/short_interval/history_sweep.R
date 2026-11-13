#!/usr/bin/env Rscript
# How the error depends on the history, and what choosing it by marginal
# likelihood does, on the development benchmark copy (devbench.R): fits with
# history = interval (`hdt`) and 30, 100, 300, 1000, 3000 days (`h30` ...).
# A rule picks, per dataset, one of these fits by the lowest objective among
# its allowed candidates; its RMSE is compared with history = interval.
#
#   Rscript validation/sparse_joint/short_interval/history_sweep.R

si_script_dir <- "validation/sparse_joint/short_interval"
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_repo, "benchmark", "sparse_two_timepoint", "report.R"))
output_dir <- file.path(si_default_output, "devbench")
metrics <- utils::read.csv(file.path(output_dir, "metrics_dev.csv"), stringsAsFactors = FALSE)
index <- utils::read.csv(file.path(output_dir, "dataset_index_dev.csv"), stringsAsFactors = FALSE)
arms <- c("hdt", "h30", "h100", "h300", "h1000", "h3000")

fits <- do.call(rbind, lapply(index$dataset_id, function(id) {
  interval <- readRDS(file.path(output_dir, "datasets", paste0(id, ".rds")))$interval
  do.call(rbind, lapply(arms, function(arm) {
    path <- file.path(output_dir, "fits", sprintf("%s__%s.rds", id, arm))
    if (!file.exists(path)) return(NULL)
    d <- readRDS(path)$diagnostics
    if (is.null(d)) return(NULL)
    data.frame(dataset_id = id, arm = arm, interval = interval, history = d$history, objective = d$objective,
               stringsAsFactors = FALSE)
  }))
}))
scores <- metrics[metrics$method %in% arms & metrics$stratum %in% c("d0", "d1", "d2"),
                  c("dataset_id", "landscape_id", "scenario", "role", "method", "stratum", "rmse", "spearman")]
names(scores)[names(scores) == "method"] <- "arm"

# A rule: allowed candidates as a function of the interval; NULL = all.
rules <- list(
  interval = function(h, dt) h == dt,
  all = function(h, dt) rep(TRUE, length(h)),
  at_least_interval = function(h, dt) h >= dt,
  at_most_300 = function(h, dt) h == dt | h <= 300,
  at_least_interval_at_most_1000 = function(h, dt) h >= dt & h <= max(dt, 1000)
)
choose <- function(rule) {
  do.call(rbind, lapply(split(fits, fits$dataset_id), function(f) {
    ok <- rules[[rule]](f$history, f$interval)
    f[ok, ][which.min(f$objective[ok]), c("dataset_id", "arm", "history")]
  }))
}
options(width = 200)
cat("== Mean RMSE by fixed history, landscape-weighted (rows: scenario x distance)\n")
by_arm <- stats::aggregate(rmse ~ scenario + role + stratum + arm, scores, mean)
wide <- stats::reshape(by_arm, idvar = c("scenario", "role", "stratum"), timevar = "arm", direction = "wide")
names(wide) <- sub("rmse.", "", names(wide), fixed = TRUE)
wide <- wide[order(wide$role, wide$scenario, wide$stratum), c("role", "scenario", "stratum", arms)]
for (arm in arms[-1L]) wide[[arm]] <- round(wide[[arm]] / wide$hdt, 3)
wide$hdt <- round(wide$hdt, 4)
print(wide, row.names = FALSE)

cat("\n== Rules: chosen history and RMSE relative to history = interval\n")
for (rule in names(rules)[-1L]) {
  chosen <- choose(rule)
  picked <- merge(chosen[c("dataset_id", "arm")], scores, by = c("dataset_id", "arm"))
  base <- scores[scores$arm == "hdt", ]
  both <- merge(picked, base, by = c("dataset_id", "landscape_id", "scenario", "role", "stratum"),
                suffixes = c("", ".base"))
  per <- stats::aggregate(cbind(rmse, rmse.base, spearman, spearman.base) ~ scenario + role + stratum, both, mean)
  per$ratio <- per$rmse / per$rmse.base
  per$spearman_change <- per$spearman - per$spearman.base
  cat(sprintf("\n-- rule %s; chosen histories:\n", rule))
  chosen <- merge(chosen, index[c("dataset_id", "scenario")])
  print(table(chosen$scenario, chosen$history))
  core <- per[per$role == "core", ]
  cat(sprintf("core cells with ratio > 1.10: %d of %d; largest %.3f; mean core ratio %.3f\n",
              sum(core$ratio > 1.10), nrow(core), max(core$ratio), mean(core$ratio)))
  print(per[order(per$role, per$scenario, per$stratum), c("role", "scenario", "stratum", "ratio", "spearman_change")],
        row.names = FALSE, digits = 3)
}
