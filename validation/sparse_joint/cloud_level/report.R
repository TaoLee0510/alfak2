#!/usr/bin/env Rscript
# Tables of the cloud-level study from results/scores.csv (README.md).
#
#   Rscript validation/sparse_joint/cloud_level/report.R
#
# Writes results/report/*.csv and prints them. Every comparison is paired on
# the datasets where the arm and `base` both have a fit; intervals are 95%
# percentile intervals of a bootstrap over landscapes (2000 draws, seed 1).

cl_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
out_dir <- file.path(cl_script_dir, "results", "report")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scores <- utils::read.csv(file.path(cl_script_dir, "results", "scores.csv"), stringsAsFactors = FALSE)
arms <- c("base", "free", "level", "level_dt")
candidates <- setdiff(arms, "base")
draws <- 2000L

# Groups of datasets the tables are split by.
scores$group <- with(scores, ifelse(set == "dev", paste("dev", role), ifelse(set == "short", "short", "diverse")))
scores$subgroup <- with(scores, ifelse(set == "dev", scenario,
                                       ifelse(set == "short", sprintf("%02d days", interval), wavelength_class)))
status <- as.data.frame(table(set = scores$set, method = scores$method, status = scores$status), stringsAsFactors = FALSE)
status <- status[status$Freq > 0, ]
utils::write.csv(status, file.path(out_dir, "status.csv"), row.names = FALSE)
cat("\n== Fits ==\n")
print(status, row.names = FALSE)
ok <- scores[scores$status == "ok", ]

# Derived arms, no new fits: per dataset, the fit with the highest marginal
# likelihood (lowest objective) among
#   select   base and level_dt
#   select3  base, level and level_dt
# Added on 2026-11-16 after the interim scores of dev and short showed that a
# free level costs about a quarter of the accuracy there, and before the
# diverse and Salehi results were looked at.
derive <- function(name, from) {
  rows <- ok[ok$method %in% from, ]
  rows <- rows[is.finite(rows$objective), ]
  complete <- names(which(table(rows$dataset_id) == length(from)))
  rows <- rows[rows$dataset_id %in% complete, ]
  rows <- rows[order(rows$dataset_id, rows$objective), ]
  rows <- rows[!duplicated(rows$dataset_id), ]
  rows$chosen <- rows$method
  rows$method <- name
  rows
}
ok$chosen <- ok$method
ok <- rbind(ok, derive("select", c("base", "level_dt")), derive("select3", c("base", "level", "level_dt")))
arms <- c(arms, "select", "select3")
candidates <- setdiff(arms, "base")
chosen <- as.data.frame(table(set = ok$set[ok$method %in% c("select", "select3")],
                              arm = ok$method[ok$method %in% c("select", "select3")],
                              chosen = ok$chosen[ok$method %in% c("select", "select3")]), stringsAsFactors = FALSE)
chosen <- chosen[chosen$Freq > 0, ]
utils::write.csv(chosen, file.path(out_dir, "selected.csv"), row.names = FALSE)
cat("\n== Fits the derived arms chose ==\n")
print(chosen, row.names = FALSE)

wide <- function(data, column) {
  w <- stats::reshape(data[c("dataset_id", "landscape_id", "method", column)], idvar = c("dataset_id", "landscape_id"),
                      timevar = "method", direction = "wide")
  names(w) <- sub(paste0("^", column, "\\."), "", names(w))
  w
}

# Bootstrap over landscapes of a statistic of the rows of `w`.
boot <- function(w, statistic) {
  groups <- split(seq_len(nrow(w)), w$landscape_id)
  set.seed(1)
  values <- vapply(seq_len(draws), function(i) {
    statistic(w[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), , drop = FALSE])
  }, numeric(1))
  stats::quantile(values, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

# Ratio of mean RMSE to `base`, and difference in mean direction accuracy.
paired <- function(data, column, kind) {
  w <- wide(data, column)
  do.call(rbind, lapply(intersect(candidates, names(w)), function(arm) {
    use <- w[is.finite(w[[arm]]) & is.finite(w$base), , drop = FALSE]
    if (nrow(use) < 3L) return(NULL)
    statistic <- if (kind == "ratio") function(d) mean(d[[arm]]) / mean(d$base) else function(d) mean(d[[arm]]) - mean(d$base)
    interval <- boot(use, statistic)
    data.frame(score = column, arm = arm, n = nrow(use), landscapes = length(unique(use$landscape_id)),
               base = mean(use$base), value = mean(use[[arm]]), effect = statistic(use), low = interval[1L],
               high = interval[2L], stringsAsFactors = FALSE)
  }))
}

by_group <- function(columns, kind, split = "group") {
  rows <- lapply(split(ok, ok[[split]]), function(data) {
    rows <- do.call(rbind, lapply(columns, function(column) paired(data, column, kind)))
    if (is.null(rows)) return(NULL)
    cbind(group = data[[split]][[1L]], rows, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

show <- function(table, title, digits = 3) {
  cat("\n==", title, "==\n")
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], function(v) signif(v, digits))
  print(table, row.names = FALSE)
}

rmse <- by_group(paste0("rmse_d", 0:2), "ratio")
utils::write.csv(rmse, file.path(out_dir, "rmse_ratio.csv"), row.names = FALSE)
show(rmse, "RMSE of relative fitness: mean, and ratio to base (effect, low, high)")

ok$key <- paste(ok$group, ok$subgroup, sep = " | ")
rmse_sub <- by_group(paste0("rmse_d", 0:2), "ratio", "key")
utils::write.csv(rmse_sub, file.path(out_dir, "rmse_ratio_subgroups.csv"), row.names = FALSE)
show(rmse_sub[rmse_sub$score == "rmse_d1" & !grepl("^dev", rmse_sub$group), ],
     "RMSE at distance 1 by interval (short) and level (diverse)")

direction <- by_group(paste0("direction_d", 0:2), "difference")
utils::write.csv(direction, file.path(out_dir, "direction_difference.csv"), row.names = FALSE)
show(direction, "Direction accuracy: mean, and difference from base")
direction_sub <- by_group(paste0("direction_d", 0:2), "difference", "key")
utils::write.csv(direction_sub, file.path(out_dir, "direction_difference_subgroups.csv"), row.names = FALSE)

# Zero predictor beside the arms.
reference <- do.call(rbind, lapply(split(ok[ok$method == "base", ], ok$group[ok$method == "base"]), function(d) {
  data.frame(group = d$group[[1L]], n = nrow(d), rmse_zero_d0 = mean(d$rmse_zero_d0),
             rmse_zero_d1 = mean(d$rmse_zero_d1), rmse_zero_d2 = mean(d$rmse_zero_d2))
}))
utils::write.csv(reference, file.path(out_dir, "rmse_zero.csv"), row.names = FALSE)
show(reference, "RMSE of predicting no differences")

# Scale of the fitness differences and the Fig. 2d counterpart, per arm.
spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4L) NA_real_ else stats::cor(x[keep], y[keep], method = "spearman")
}
scale_rows <- function(data) {
  do.call(rbind, lapply(intersect(arms, unique(data$method)), function(arm) {
    d <- data[data$method == arm, ]
    ratio <- d$var_pred / d$var_true
    chosen <- d$history / d$interval
    data.frame(
      arm = arm, n = nrow(d),
      variance_ratio = stats::median(ratio, na.rm = TRUE),
      variance_ratio_q25 = stats::quantile(ratio, 0.25, na.rm = TRUE, names = FALSE),
      variance_ratio_q75 = stats::quantile(ratio, 0.75, na.rm = TRUE, names = FALSE),
      spearman_variance = spearman(d$var_pred, d$var_true),
      # Fig. 2d counterpart: predicted against observed change, and the ceiling
      # that the true fitness reaches on the same datasets.
      spearman_change = spearman(d$change_pred, d$change_observed),
      spearman_change_truth = spearman(d$change_true, d$change_observed),
      change_pred = stats::median(d$change_pred, na.rm = TRUE),
      change_true = stats::median(d$change_true, na.rm = TRUE),
      change_observed = stats::median(d$change_observed, na.rm = TRUE),
      beats_start = mean(d$error_pred < d$change_observed, na.rm = TRUE),
      beats_start_truth = mean(d$error_true < d$change_observed, na.rm = TRUE),
      # Does the scale follow the history that was chosen?
      history_over_interval = stats::median(chosen, na.rm = TRUE),
      share_history_longer = mean(chosen > 1.01, na.rm = TRUE),
      spearman_history_variance = spearman(chosen, ratio),
      spearman_history_shell = spearman(chosen, d$shell_coef),
      balance_level = stats::median(d$balance_level, na.rm = TRUE),
      shell_coef = stats::median(d$shell_coef, na.rm = TRUE),
      scale_residual = stats::median(d$scale_residual, na.rm = TRUE),
      seconds = stats::median(d$seconds, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
}
scale_table <- function(split) {
  do.call(rbind, lapply(split(ok, ok[[split]]), function(data) cbind(group = data[[split]][[1L]], scale_rows(data))))
}
scales <- scale_table("group")
utils::write.csv(scales, file.path(out_dir, "scale_and_change.csv"), row.names = FALSE)
show(scales[c("group", "arm", "n", "variance_ratio", "variance_ratio_q25", "variance_ratio_q75", "spearman_variance")],
     "Fitness variance among the first time point's cells: predicted over true")
show(scales[c("group", "arm", "spearman_change", "spearman_change_truth", "change_pred", "change_true",
              "change_observed", "beats_start", "beats_start_truth")],
     "Fig. 2d counterpart: predicted against observed change")
show(scales[c("group", "arm", "history_over_interval", "share_history_longer", "spearman_history_variance",
              "spearman_history_shell", "balance_level", "shell_coef", "scale_residual", "seconds")],
     "History, level and scales")
scales_sub <- scale_table("key")
utils::write.csv(scales_sub, file.path(out_dir, "scale_and_change_subgroups.csv"), row.names = FALSE)
show(scales_sub[!grepl("^dev", scales_sub$group),
                c("group", "arm", "n", "variance_ratio", "spearman_change", "spearman_change_truth", "beats_start",
                  "history_over_interval", "spearman_history_variance", "balance_level", "shell_coef")],
     "By interval (short) and level (diverse)")

# Criteria 1 and 2 of the README.
verdict <- function(groups, label) {
  do.call(rbind, lapply(candidates, function(arm) {
    if (is.null(rmse) || is.null(direction)) return(NULL)
    r <- rmse[rmse$group %in% groups & rmse$arm == arm, ]
    d <- direction[direction$group %in% groups & direction$arm == arm, ]
    if (!nrow(r)) return(NULL)
    data.frame(rule = label, arm = arm, worst_rmse_ratio = max(r$effect), worst_direction_difference = min(d$effect),
               pass = max(r$effect) <= 1.03 && min(d$effect) >= -0.01, stringsAsFactors = FALSE)
  }))
}
criteria <- rbind(verdict("dev core", "1 dev core"), verdict("short", "1 short"), verdict("diverse", "2 diverse"),
                  verdict("dev stress", "(dev stress, reported)"))
if (!is.null(criteria)) {
  utils::write.csv(criteria, file.path(out_dir, "criteria.csv"), row.names = FALSE)
  show(criteria, "Criteria 1 and 2")
}
