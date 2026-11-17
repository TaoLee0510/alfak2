#!/usr/bin/env Rscript
# Choosing the missegregation rate by marginal likelihood, step 0 verdict.
#
# Reads the fits of sim_fit.R (results/fits/, one file per input and rate).
# For each input the rate with the lowest objective among the six of the grid
# (5e-5 x 1/16, 1/4, 1, 4, 16, 64) is the pick. Criteria, fixed before any
# fit was scored (plan of 2026-11-15, memory alfak2-beta-auto-plan):
#   (1) at least 70% of picks within 4 times of the true rate, and the median
#       |log2(pick / true)| at most 1;
#   (2) the paired median of r(pick) - r(true rate) at least -0.02, and on
#       inputs whose true rate is at least 4 times the default the paired
#       median of r(pick) - r(5e-5) above 0.
# r is the landscape Pearson r against the true fitness (the benchmark's
# landscape_scores on its minobs-10 karyotype set).
#
#   Rscript validation/sparse_joint/beta_choice/analyse_sim.R
#
# Writes results/sim_choice.rds and results/sim_report.md.

base <- "validation/sparse_joint/beta_choice/results"
files <- list.files(file.path(base, "fits"), pattern = "\\.rds$", full.names = TRUE)
rows <- lapply(files, function(f) readRDS(f)$summary)
cols <- unique(unlist(lapply(rows, names)))
s <- do.call(rbind, lapply(rows, function(r) {
  for (c in setdiff(cols, names(r))) r[[c]] <- NA
  r[cols]
}))
grid <- c(1 / 16, 1 / 4, 1, 4, 16, 64)
s$key <- paste(s$source, s$sim_id, s$window)
s$on_grid <- vapply(s$beta_multiple, function(m) any(abs(grid - m) < 1e-9), logical(1))
s$is_true <- abs(s$beta - s$true_beta) < 1e-12

per_input <- do.call(rbind, lapply(split(s, s$key), function(d) {
  g <- d[d$on_grid & is.finite(d$objective), ]
  tr <- d[d$is_true & is.finite(d$objective), ]
  fx <- d[abs(d$beta_multiple - 1) < 1e-9 & is.finite(d$objective), ]
  if (nrow(g) < length(grid) || !nrow(tr) || !nrow(fx)) {
    return(data.frame(key = d$key[1], source = d$source[1], window = d$window[1], w = d$w[1],
                      true_beta = d$true_beta[1], complete = FALSE, stringsAsFactors = FALSE))
  }
  pick <- g[which.min(g$objective), ]
  data.frame(key = d$key[1], source = d$source[1], window = d$window[1], w = d$w[1], true_beta = d$true_beta[1],
             complete = TRUE, pick_beta = pick$beta, log2_ratio = log2(pick$beta / d$true_beta[1]),
             r_pick = pick$r, r_true = tr$r[1], r_fixed = fx$r[1], R_pick = pick$R, R_true = tr$R[1], R_fixed = fx$R[1],
             objective_gain_vs_fixed = fx$objective[1] - pick$objective,
             objective_gain_true_vs_fixed = fx$objective[1] - tr$objective[1], stringsAsFactors = FALSE)
}))
saveRDS(list(fits = s, inputs = per_input), file.path(base, "sim_choice.rds"))
ok <- per_input[per_input$complete, ]

within4 <- mean(abs(ok$log2_ratio) <= 2 + 1e-9)
med_abs <- stats::median(abs(ok$log2_ratio))
d_true <- stats::median(ok$r_pick - ok$r_true, na.rm = TRUE)
high <- ok[ok$true_beta >= 4 * 5e-5, ]
d_fixed_high <- stats::median(high$r_pick - high$r_fixed, na.rm = TRUE)
c1 <- within4 >= 0.70 && med_abs <= 1
c2 <- d_true >= -0.02 && d_fixed_high > 0
verdict <- data.frame(
  check = c("picks within 4x of the truth", "median |log2(pick/true)|", "criterion 1",
            "median r(pick) - r(true)", "median r(pick) - r(5e-5), truth >= 4x", "criterion 2", "PASS"),
  value = c(sprintf("%.2f (>= 0.70)", within4), sprintf("%.2f (<= 1)", med_abs), c1,
            sprintf("%.3f (>= -0.02)", d_true), sprintf("%.3f (> 0, %d inputs)", d_fixed_high, nrow(high)), c2,
            c1 && c2))

md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
ok$group <- ifelse(ok$source == "abm", paste0("ABM w=", ok$w), paste0("misseg ", format(ok$true_beta)))
picks <- as.data.frame.matrix(table(ok$group, factor(ok$pick_beta / 5e-5, levels = grid)))
picks <- cbind(group = rownames(picks), picks)
by_group <- do.call(rbind, lapply(split(ok, list(ok$group, ok$window), drop = TRUE), function(d) data.frame(
  group = d$group[1], window = d$window[1], n = nrow(d), median_log2 = stats::median(d$log2_ratio),
  within4 = mean(abs(d$log2_ratio) <= 2 + 1e-9), r_pick = stats::median(d$r_pick, na.rm = TRUE),
  r_true = stats::median(d$r_true, na.rm = TRUE), r_fixed = stats::median(d$r_fixed, na.rm = TRUE),
  d_true = stats::median(d$r_pick - d$r_true, na.rm = TRUE), d_fixed = stats::median(d$r_pick - d$r_fixed, na.rm = TRUE))))
by_group <- by_group[order(by_group$group, by_group$window), ]
failed <- sum(!is.na(s$error))

writeLines(c("# Choosing the missegregation rate by marginal likelihood: step 0", "",
             sprintf("Inputs complete: %d of %d; failed fits: %d.", nrow(ok), nrow(per_input), failed), "",
             "## Verdict (criteria fixed before scoring)", "", md(verdict), "",
             "## Picks (multiples of 5e-5), all windows", "", md(picks), "",
             "## By group and window", "", md(by_group)),
           file.path(base, "sim_report.md"))
cat(readLines(file.path(base, "sim_report.md")), sep = "\n")
