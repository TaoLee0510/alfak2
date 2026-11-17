#!/usr/bin/env Rscript
# Corrected-rule study, step 0: what tells the inputs where switching to
# initial_state = "free" helps from those where it hurts? (README.md)
#
# Inputs and the effect of switching against the truth:
#   (1) the ALFA-K paper benchmark: every two-time-point ABM and
#       missegregation input (paper_features.R), and the benchmark's paired
#       comparison of the switched fits with their balance fits
#       (results_827c4b5/analysis/initial_state/vs_truth_fits.tsv, copied to
#       results/step0/). The benchmark fits at the true rate, so these inputs
#       show only where switching must not happen.
#   (2) this study's simulations (scores.csv, scores_confirm.csv and the
#       paper-sized set paper3 in scores_paper3.csv): base, free and level_dt
#       on every input, scored against the truth.
# For each input the features a rule could use before or after the first fit:
#   level      the level of the level_dt fit (the "auto" check)
#   dt         the interval in days (one division per day is assumed)
#   history    the history the default chose, over the interval
#   coverage   the share of the first time point's cells in the domain
#   singleton  the share of its cells on karyotypes seen once
#   top        the share of its largest karyotype
#   n_major    the number of major karyotypes
# and the effect of switching: the change in landscape correlation with the
# truth (paper) or in RMSE at distance 1 and direction accuracy (study).
#   Rscript validation/sparse_joint/cloud_level/rule/step0.R
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- file.path(script_dir, "..", "results")
out_dir <- file.path(res, "rule_step0")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
q <- function(v) round(stats::quantile(v, c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE), 2)
fmt <- function(v) paste(q(v), collapse = " / ")

# (1) Paper benchmark.
pf <- utils::read.csv(file.path(res, "step0", "paper_features.csv"), stringsAsFactors = FALSE)
vt <- utils::read.delim(file.path(res, "step0", "vs_truth_fits.tsv"), stringsAsFactors = FALSE)
vt <- stats::aggregate(cbind(dr = r.new - r.old, dcv = r_heldout_true.new - r_heldout_true.old) ~ source + dataset + design,
                       data = vt[vt$chosen == "free", ], FUN = mean)
pf <- merge(pf, vt, by = c("source", "dataset", "design"), all.x = TRUE)
pf$hr <- pf$history / pf$dt
pf$group <- paste("paper", pf$source, pf$design)
pf$switched <- pf$chosen == "free"
pf$effect <- ifelse(pf$switched, ifelse(pf$dr < 0, "hurt", "helped"), "not switched")

# (2) This study.
read_scores <- function(file) {
  path <- file.path(res, file)
  if (!file.exists(path)) return(NULL)
  s <- utils::read.csv(path, stringsAsFactors = FALSE)
  s[s$status == "ok", ]
}
st <- do.call(rbind, lapply(c("scores_rule.csv", "scores_paper3.csv"), function(f) {
  s <- read_scores(f)
  if (is.null(s)) return(NULL)
  columns <- c("set", "dataset_id", "landscape_id", "wavelength_class", "role", "interval", "coverage_t0",
               "singleton_t0", "top_t0", "n_major", "method", "rmse_d1", "direction_d1", "landscape_r", "history",
               "balance_level")
  s[intersect(columns, names(s))]
}))
lev <- st[st$method == "level_dt", c("set", "dataset_id", "balance_level")]
names(lev)[3] <- "level"
# The level against the balance over the population's age (arm level_age).
lev_age <- st[st$method == "level_age", c("set", "dataset_id", "balance_level")]
names(lev_age)[3] <- "level_age"
lev <- merge(lev, lev_age, by = c("set", "dataset_id"), all.x = TRUE)
base <- st[st$method == "base", ]
free <- st[st$method == "free", c("set", "dataset_id", "rmse_d1", "direction_d1", "landscape_r")]
names(free)[3:5] <- paste0(names(free)[3:5], "_free")
sd <- merge(merge(base, lev, by = c("set", "dataset_id")), free, by = c("set", "dataset_id"))
sd$hr <- sd$history / sd$interval
sd$drmse <- sd$rmse_d1_free / sd$rmse_d1 - 1
sd$ddir <- sd$direction_d1_free - sd$direction_d1
sd$dr <- sd$landscape_r_free - sd$landscape_r
sd$group <- ifelse(sd$set %in% c("dev", "dev2"), paste(sd$set, ifelse(sd$role == "core", "core", "stress")),
                   ifelse(sd$set %in% c("diverse", "diverse2", "paper3"), paste(sd$set, sd$role), sd$set))
sd$group <- ifelse(sd$set %in% c("short", "short2"), paste(sd$set, sprintf("%02dd", sd$interval)), sd$group)
sd$group <- ifelse(sd$set == "paper3", paste(sd$set, sd$role, sprintf("%02dd", sd$interval)), sd$group)
# Switching helps where free is clearly better against the truth.
sd$effect <- ifelse(sd$drmse < -0.05 | sd$dr > 0.05, "helped", ifelse(sd$drmse > 0.05 | sd$dr < -0.05, "hurt", "neutral"))

feat_table <- function(d, features) {
  do.call(rbind, lapply(split(d, d$group), function(g) {
    data.frame(group = g$group[1], inputs = nrow(g),
               t(vapply(features, function(f) fmt(g[[f]]), character(1))), stringsAsFactors = FALSE)
  }))
}
features <- c("level", "coverage_t0", "singleton_t0", "top_t0", "n_major", "hr")
pf_t <- feat_table(pf, features)
cat("== Paper benchmark: features (10% / 50% / 90%) ==\n"); print(pf_t, row.names = FALSE)
eff_p <- stats::aggregate(cbind(n = 1, dr, dcv) ~ group + effect, data = transform(pf, dr = ifelse(is.na(dr), 0, dr),
                                                                                      dcv = ifelse(is.na(dcv), 0, dcv)),
                          FUN = function(v) if (all(v == 1)) length(v) else round(stats::median(v), 3))
cat("\n== Paper benchmark: effect of the switch (median change in landscape r and held-out r) ==\n"); print(eff_p, row.names = FALSE)
sd_t <- feat_table(sd, c("level_age", features))
sd_e <- do.call(rbind, lapply(split(sd, sd$group), function(g) data.frame(
  group = g$group[1], inputs = nrow(g), helped = mean(g$effect == "helped"), hurt = mean(g$effect == "hurt"),
  rmse_change = round(stats::median(g$drmse), 3), direction_change = round(stats::median(g$ddir, na.rm = TRUE), 3),
  landscape_r_change = round(stats::median(g$dr, na.rm = TRUE), 3))))
cat("\n== Study simulations: features ==\n"); print(sd_t, row.names = FALSE)
cat("\n== Study simulations: effect of free against base ==\n"); print(sd_e, row.names = FALSE)
utils::write.csv(pf_t, file.path(out_dir, "paper_features.csv"), row.names = FALSE)
utils::write.csv(sd_t, file.path(out_dir, "study_features.csv"), row.names = FALSE)
utils::write.csv(sd_e, file.path(out_dir, "study_effects.csv"), row.names = FALSE)
saveRDS(list(paper = pf, study = sd), file.path(out_dir, "inputs.Rds"))

# How well does each level separate the inputs where switching helped from
# those where it hurt? Area under the curve (helped against hurt), and the
# share of each kind above a threshold.
auc <- function(score, positive) {
  ok <- is.finite(score) & !is.na(positive)
  score <- score[ok]; positive <- positive[ok]
  n1 <- sum(positive); n0 <- sum(!positive)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score)[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
lab <- sd[sd$effect != "neutral", ]
lab$helped <- lab$effect == "helped"
sep <- data.frame(feature = c("level", "level_age", "coverage_t0", "singleton_t0"),
                  auc = c(auc(lab$level, lab$helped), auc(lab$level_age, lab$helped),
                          auc(-lab$coverage_t0, lab$helped), auc(lab$singleton_t0, lab$helped)))
cat("\n== Separation of helped from hurt inputs (study simulations; AUC, 0.5 = none) ==\n"); print(sep, digits = 3)
thr <- do.call(rbind, lapply(c(0, 0.5, 1, 1.5), function(t) do.call(rbind, lapply(split(sd, sd$group), function(g) data.frame(
  threshold = t, group = g$group[1], inputs = nrow(g), switched = mean(g$level_age > t, na.rm = TRUE),
  helped_among_switched = if (any(g$level_age > t, na.rm = TRUE)) mean(g$effect[which(g$level_age > t)] == "helped") else NA,
  hurt_among_switched = if (any(g$level_age > t, na.rm = TRUE)) mean(g$effect[which(g$level_age > t)] == "hurt") else NA)))))
cat("\n== Switching when level_age exceeds a threshold ==\n"); print(thr, digits = 2, row.names = FALSE)
utils::write.csv(sep, file.path(out_dir, "separation.csv"), row.names = FALSE)
utils::write.csv(thr, file.path(out_dir, "threshold_level_age.csv"), row.names = FALSE)
saveRDS(list(paper = pf, study = sd), file.path(out_dir, "inputs.Rds"))
