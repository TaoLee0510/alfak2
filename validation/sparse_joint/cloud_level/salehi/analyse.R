#!/usr/bin/env Rscript
# The arms of the cloud-level study on the Salehi passages, with the paper
# benchmark's Fig. 2d as the test: predicted against observed change
# (Wasserstein distances from the parent passage to the forecast and to the
# daughter passage, on the landscape's karyotypes). A confirmation: nothing is
# tuned on these results.
#
#   Rscript validation/sparse_joint/cloud_level/salehi/analyse.R [BENCHMARK_REPO]
#
# Reads the rows that the benchmark's scripts/fig2d_check.R wrote, one per
# fit, minobs and daughter passage:
#   <benchmark>/results_05a8855/analysis/fig2d_check/salehi   alfak2, alfakR
#   results/salehi/check            this study's arms (alfak2 build cl1)
#   results/salehi/variants_check   the alfak2-salehi-variants project's fits
# All comparisons are on the rows that every compared method has, without
# the paper's selection by cross-validation. Intervals: 95% percentile
# intervals of a bootstrap over lineage types (the paper's unit) and, beside
# it, over samples; 4000 draws, seed 1.
# Writes results/salehi/report/*.csv and prints them.

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
args <- commandArgs(TRUE)
benchmark <- if (length(args)) args[[1L]] else "~/Documents/GitHub/alfak2-paper-benchmark"
results <- file.path(script_dir, "..", "results", "salehi")
out_dir <- file.path(results, "report")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_chunks <- function(dir) {
  files <- list.files(dir, pattern = "^chunk_.*\\.Rds$", full.names = TRUE)
  rows <- lapply(files, readRDS)
  rows <- Filter(function(r) !is.null(r) && nrow(r) > 0L, rows)
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(columns, names(r))] <- NA
    r[columns]
  }))
}
sources <- list(file.path(benchmark, "results_05a8855", "analysis", "fig2d_check", "salehi"),
                file.path(results, "check"), file.path(results, "variants_check"))
tables <- Filter(Negate(is.null), lapply(sources[dir.exists(unlist(sources))], read_chunks))
columns <- unique(unlist(lapply(tables, names)))
x <- do.call(rbind, lapply(tables, function(r) {
  r[setdiff(columns, names(r))] <- NA
  r[columns]
}))
x$method <- sub("^alfak2_", "", x$method)
x$key <- paste(x$fi, x$min_obs, x$dj)
x$passages <- as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", x$fi))
x$sample <- sub("_X[0-9]+_l_.*$", "", x$fi)
cat("Rows by method:\n")
print(table(x$method))

sp <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 4L) NA_real_ else stats::cor(a[keep], b[keep], method = "spearman")
}
show <- function(table, title) {
  cat("\n==", title, "==\n")
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], function(v) signif(v, 3))
  print(table, row.names = FALSE)
}

# 1. The build at its defaults against the round's alfak2 fits.
if ("base" %in% x$method) {
  a <- x[x$method == "base", ]
  b <- x[x$method == "alfak2", ]
  m <- merge(a[c("key", "pred0", "var_f", "history")], b[c("key", "pred0", "var_f", "history")], by = "key",
             suffixes = c("_build", "_round"))
  same <- data.frame(rows = nrow(m), inputs = length(unique(a$fi)),
                     same_history = mean(m$history_build == m$history_round),
                     max_difference_variance = max(abs(m$var_f_build - m$var_f_round)),
                     max_difference_predicted_change = max(abs(m$pred0_build - m$pred0_round)))
  utils::write.csv(same, file.path(out_dir, "default_reproduces_round.csv"), row.names = FALSE)
  show(same, "Build cl1 at its defaults against the round's alfak2 fits")
}

# Derived arm `select`, as in the simulations (report.R): per input, the fit
# with the higher marginal likelihood of the round's alfak2 fit (the default)
# and level_dt. Objectives: the round's from step 0 (results/step0), level_dt's
# from its run.Rds (results/salehi/runs).
step0 <- utils::read.csv(file.path(script_dir, "..", "results", "step0", "salehi_fits.csv"), stringsAsFactors = FALSE)
step0 <- step0[step0$design == "ends" & step0$method == "alfak2", c("id", "objective")]
run_files <- Sys.glob(file.path(results, "runs", "salehi", "*", "ends", "alfak2_*", "minobs_5", "run.Rds"))
runs <- do.call(rbind, lapply(run_files, function(f) {
  r <- readRDS(f)
  parts <- strsplit(f, "/", fixed = TRUE)[[1L]]
  n <- length(parts)
  data.frame(id = parts[n - 4L], method = sub("^alfak2_", "", parts[n - 2L]), status = r$status,
             objective = if (is.null(r$objective)) NA_real_ else r$objective, stringsAsFactors = FALSE)
}))
if ("base" %in% runs$method) {
  same <- merge(runs[runs$method == "base", ], step0, by = "id", suffixes = c("_build", "_round"))
  cat("\nObjective of build cl1 at its defaults against the round's: max |difference|",
      signif(max(abs(same$objective_build - same$objective_round)), 3), "over", nrow(same), "inputs\n")
}
level_dt_objective <- runs[runs$method == "level_dt" & runs$status == "ok", c("id", "objective")]
choice <- merge(step0, level_dt_objective, by = "id", suffixes = c("_default", "_level_dt"))
choice$pick <- ifelse(choice$objective_level_dt < choice$objective_default, "level_dt", "alfak2")
choice$passages <- as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", choice$id))
picked <- as.data.frame(table(one_interval = choice$passages == 2L, pick = choice$pick))
utils::write.csv(choice, file.path(out_dir, "select_choice.csv"), row.names = FALSE)
show(picked, "Fits that select picks (Salehi, two time points)")
select_rows <- do.call(rbind, lapply(seq_len(nrow(choice)), function(i) {
  rows <- x[x$fi == choice$id[[i]] & x$method == choice$pick[[i]], ]
  rows$method <- rep("select", nrow(rows))
  rows
}))
x <- rbind(x, select_rows)

# 2. Rank agreement on common rows.
main <- intersect(c("alfakR", "alfak2", "level", "level_dt", "select", "free"), unique(x$method))
extra <- intersect(c("wide", "observed", "observed_free"), unique(x$method))
common_keys <- function(methods) Reduce(intersect, lapply(methods, function(m) x$key[x$method == m & is.finite(x$pred0) & is.finite(x$base)]))
subsets <- list(all = function(d) rep(TRUE, nrow(d)), `one interval` = function(d) d$passages == 2L,
                `more intervals` = function(d) d$passages > 2L)
boot_values <- function(d, cluster, statistic, draws = 4000L) {
  groups <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  vapply(seq_len(draws), function(i) statistic(d[unlist(groups[sample.int(length(groups), replace = TRUE)],
                                                        use.names = FALSE), , drop = FALSE]), numeric(1))
}
interval <- function(values) stats::quantile(values, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

rank_tables <- function(methods, label) {
  keys <- common_keys(methods)
  wide <- x[x$method == methods[[1L]] & x$key %in% keys, c("key", "fi", "type", "sample", "passages", "base")]
  for (m in methods) {
    rows <- x[x$method == m & x$key %in% keys, ]
    rows <- rows[match(wide$key, rows$key), ]
    wide[[paste0("pred0_", m)]] <- rows$pred0
    wide[[paste0("pred_", m)]] <- rows$pred
    wide[[paste0("var_", m)]] <- rows$var_f
  }
  per_method <- do.call(rbind, lapply(names(subsets), function(s) {
    d <- wide[subsets[[s]](wide), ]
    do.call(rbind, lapply(methods, function(m) {
      column <- paste0("pred0_", m)
      ci <- interval(boot_values(d, "type", function(b) sp(b$base, b[[column]])))
      data.frame(rows = label, subset = s, method = m, n = nrow(d), inputs = length(unique(d$fi)),
                 spearman = sp(d$base, d[[column]]), low = ci[1L], high = ci[2L],
                 predicted_change = stats::median(d[[column]]), observed_change = stats::median(d$base),
                 fitness_variance = stats::median(d[[paste0("var_", m)]]),
                 # The forecast is closer to the daughter passage than the parent is.
                 beats_parent = mean(d[[paste0("pred_", m)]] < d$base), stringsAsFactors = FALSE)
    }))
  }))
  paired <- do.call(rbind, lapply(names(subsets), function(s) {
    d <- wide[subsets[[s]](wide), ]
    do.call(rbind, lapply(setdiff(methods, c("alfakR", "alfak2")), function(m) {
      do.call(rbind, lapply(intersect(c("alfakR", "alfak2"), methods), function(reference) {
        do.call(rbind, lapply(c("type", "sample"), function(cluster) {
          statistic <- function(b) sp(b$base, b[[paste0("pred0_", m)]]) - sp(b$base, b[[paste0("pred0_", reference)]])
          values <- boot_values(d, cluster, statistic)
          ci <- interval(values)
          data.frame(rows = label, subset = s, method = m, against = reference,
                     cluster = if (cluster == "type") "lineage" else "sample", clusters = length(unique(d[[cluster]])),
                     n = nrow(d), difference = statistic(d), low = ci[1L], high = ci[2L],
                     share_above_0 = mean(values > 0, na.rm = TRUE), stringsAsFactors = FALSE)
        }))
      }))
    }))
  }))
  list(per_method = per_method, paired = paired)
}
main_tables <- rank_tables(main, "main arms")
utils::write.csv(main_tables$per_method, file.path(out_dir, "rank_by_method.csv"), row.names = FALSE)
utils::write.csv(main_tables$paired, file.path(out_dir, "rank_paired.csv"), row.names = FALSE)
show(main_tables$per_method, "Predicted against observed change, common rows of the main arms")
show(main_tables$paired, "Paired difference in Spearman")
if (length(extra)) {
  all_tables <- rank_tables(c(main, extra), "all arms")
  utils::write.csv(all_tables$per_method, file.path(out_dir, "rank_by_method_all_arms.csv"), row.names = FALSE)
  utils::write.csv(all_tables$paired, file.path(out_dir, "rank_paired_all_arms.csv"), row.names = FALSE)
  show(all_tables$per_method, "The same with the variants of alfak2-salehi-variants")
}

# 3. Does the scale follow the history? Trajectories with one interval of 15
# days, one row per fit and daughter passage (minobs 5).
one <- x[x$passages == 2L & x$min_obs == 5 & is.finite(x$fit_dt) & x$fit_dt == 15, ]
history <- do.call(rbind, lapply(intersect(c("alfak2", "level", "level_dt", "free"), unique(one$method)), function(m) {
  d <- one[one$method == m, ]
  data.frame(method = m, n = nrow(d), fits = length(unique(d$fi)),
             history_15 = sum(d$history == 15), history_30 = sum(d$history == 30), history_100 = sum(d$history >= 100),
             spearman_history_variance = sp(d$history, d$var_f), spearman_history_predicted = sp(d$history, d$pred0),
             spearman_history_observed = sp(d$history, d$base),
             fitness_variance = stats::median(d$var_f), predicted_change = stats::median(d$pred0),
             observed_change = stats::median(d$base),
             balance_level = stats::median(d$balance_level, na.rm = TRUE),
             shell_coef = stats::median(d$shell_coef, na.rm = TRUE),
             scale_residual = stats::median(d$scale_residual, na.rm = TRUE), stringsAsFactors = FALSE)
}))
utils::write.csv(history, file.path(out_dir, "history_one_interval_15d.csv"), row.names = FALSE)
show(history, "One interval of 15 days, minobs 5: history and scale")
by_history <- do.call(rbind, lapply(intersect(c("alfak2", "level"), unique(one$method)), function(m) {
  d <- one[one$method == m, ]
  d$group <- ifelse(d$history >= 100, "100+", as.character(d$history))
  rows <- stats::aggregate(cbind(var_f, pred0, base) ~ group, data = d, FUN = stats::median)
  cbind(method = m, rows, n = as.integer(table(d$group)[rows$group]))
}))
utils::write.csv(by_history, file.path(out_dir, "by_history_one_interval_15d.csv"), row.names = FALSE)
show(by_history, "The same by chosen history (medians)")

# 4. The level itself.
level <- x[x$method %in% c("level", "level_dt") & x$min_obs == 5, ]
level <- level[!duplicated(paste(level$method, level$fi)), ]
if (nrow(level)) {
  levels <- do.call(rbind, lapply(split(level, level$method), function(d) {
    q <- stats::quantile(d$balance_level, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE, names = FALSE)
    data.frame(method = d$method[[1L]], fits = nrow(d), q10 = q[1L], q25 = q[2L], median = q[3L], q75 = q[4L],
               q90 = q[5L], share_above_0 = mean(d$balance_level > 0, na.rm = TRUE))
  }))
  utils::write.csv(levels, file.path(out_dir, "balance_level.csv"), row.names = FALSE)
  show(levels, "Fitted level of the balance prior (log factor over the balance)")
}
