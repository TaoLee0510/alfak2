#!/usr/bin/env Rscript
# Corrected-rule study, step 0: a threshold on the level of the "auto" check
# (history at the interval) that falls with the interval,
#   switch to "free" if level > c0 - c1 * log(dt / 6),
# on the study's simulations (effect of free against base known) and on the
# paper benchmark's inputs (where every switch is unwarranted). Exploratory:
# the grid is searched on the same data it is shown on.
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
inp <- readRDS(file.path(script_dir, "..", "results", "rule_step0", "inputs.Rds"))
sd <- inp$study; pf <- inp$paper
sd$dt <- sd$interval
rule <- function(level, dt, c0, c1) is.finite(level) & level > c0 - c1 * log(dt / 6)
grid <- expand.grid(c0 = seq(1.5, 3.5, by = 0.25), c1 = seq(0, 1.5, by = 0.25))
score <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  s <- rule(sd$level, sd$dt, g$c0, g$c1)
  p <- rule(pf$level, pf$dt, g$c0, g$c1)
  data.frame(g, helped_switched = sum(s & sd$effect == "helped") / sum(sd$effect == "helped"),
             hurt_switched = sum(s & sd$effect == "hurt") / sum(sd$effect == "hurt"),
             paper_switched = mean(p, na.rm = TRUE))
}))
score$net <- score$helped_switched - score$hurt_switched - score$paper_switched
best <- score[order(-score$net), ][1:8, ]
cat("== Best thresholds (share of helped and hurt study inputs switched; share of paper inputs switched) ==\n")
print(best, digits = 2, row.names = FALSE)
b <- best[1, ]
sd$switch <- rule(sd$level, sd$dt, b$c0, b$c1)
pf$switch <- rule(pf$level, pf$dt, b$c0, b$c1)
cat(sprintf("\n== By group at c0 = %.2f, c1 = %.2f ==\n", b$c0, b$c1))
grp <- do.call(rbind, lapply(split(sd, sd$group), function(g) data.frame(
  group = g$group[1], inputs = nrow(g), switched = mean(g$switch),
  helped_if_switched = if (any(g$switch)) mean(g$effect[g$switch] == "helped") else NA,
  hurt_if_switched = if (any(g$switch)) mean(g$effect[g$switch] == "hurt") else NA,
  rmse_change_rule = stats::median(ifelse(g$switch, g$drmse, 0)))))
print(grp, digits = 2, row.names = FALSE)
pg <- do.call(rbind, lapply(split(pf, paste(pf$source, pf$design)), function(g) data.frame(
  group = paste("paper", g$source[1], g$design[1]), inputs = nrow(g), switched = mean(g$switch, na.rm = TRUE))))
print(pg, digits = 2, row.names = FALSE)
