#!/usr/bin/env Rscript
# Salehi forecasts against later passages: who forecasts better on real data?
# (README.md). Reads the rows of forecast_score.R (results/*.Rds, copied from
# HPC) and compares alfak2 (default, balance prior), alfakR and alfak2 with
# initial_state = "free" on the rows that the compared methods share (same
# trajectory, minobs and later passage).
#
# Main score: the mean log-likelihood per cell of the later passage's counts
# under the forecast, minus that under the parent's frequencies
# (persistence), each mixed with a uniform distribution over the landscape
# with weight 0.01. Sampling noise of the later passage enters both alike, so
# a forecast gains only by putting probability where the cells go. Also: how
# often the forecast is closer than the parent in Wasserstein distance, and the
# sampling noise of the later passage (two random halves).
#
# Intervals: 95% percentile intervals of a bootstrap over the six lineage
# types (the paper's unit), and over samples (data sets); 4000 draws, seed 1.
#   Rscript validation/sparse_joint/salehi_horizon/analyse.R
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- file.path(script_dir, "results")
out_dir <- file.path(res, "report")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
files <- list.files(file.path(res, "rows"), pattern = "[.]Rds$", full.names = TRUE)
parts <- lapply(files, readRDS)
status <- as.data.frame(table(method = vapply(parts, function(p) p$task$method, character(1)),
                              status = vapply(parts, `[[`, character(1), "status")))
cat("Tasks:\n"); print(status[status$Freq > 0, ], row.names = FALSE)
x <- do.call(rbind, lapply(parts, `[[`, "rows"))
x$method <- c(alfak2 = "alfak2", alfakR = "alfakR", alfak2_free = "alfak2 free")[x$method]
# Lineage types as the paper benchmark assigns them (salehi_assign_lineage()).
lineage_map <- c(SA1035T = "SA1035", SA1035U = "SA1035", SA532 = "SA532", SA609 = "SA609", SA609UnBU = "SA609",
                 SA906a = "p53 k.o", SA906b = "p53 k.o", SA039U = "p53 w.t", SA535_CISPLATIN_CombinedT = "SA535",
                 SA535_CISPLATIN_CombinedU = "SA535", SA000 = "SA609", SA609R2 = "SA609",
                 SA535_CISPLATIN_CombinedH = "SA535", SA001 = "SA609", SA609H2 = "SA609", SA609H1 = "SA609",
                 SA609R2b = "SA609", SA609TIIDT = "SA609")
x$sample <- sub("_X[0-9]+_l_.*$", "", x$id)
x$type <- unname(lineage_map[x$sample])
x$passages <- as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", x$id))
x$treatment <- ifelse(x$treat_parent == x$treat_desc, "unchanged", "changed")
x$gain <- x$ll_fc_e2 - x$ll_p0_e2
x$gain3 <- x$ll_fc_e3 - x$ll_p0_e3
x$closer <- x$w_fc < x$w_p0
x$key <- paste(x$id, x$minobs, x$uid)
x$horizon <- ifelse(x$depth >= 3L, "3-4", as.character(x$depth))
cat("\nRows by method:\n"); print(table(x$method, x$depth))

boot <- function(d, statistic, cluster, draws = 4000L) {
  g <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  v <- vapply(seq_len(draws), function(i) statistic(d[unlist(g[sample.int(length(g), replace = TRUE)],
                                                        use.names = FALSE), , drop = FALSE]), numeric(1))
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
wide <- function(methods) {
  keys <- Reduce(intersect, lapply(methods, function(m) x$key[x$method == m & is.finite(x$gain)]))
  w <- x[x$method == methods[[1L]] & x$key %in% keys,
         c("key", "id", "minobs", "uid", "depth", "horizon", "days", "type", "sample", "treatment", "passages",
           "coverage", "cells_desc", "noise", "w_p0", "ll_p0_e2")]
  for (m in methods) {
    r <- x[x$method == m & x$key %in% keys, ]
    r <- r[match(w$key, r$key), ]
    w[[paste0("gain|", m)]] <- r$gain
    w[[paste0("gain3|", m)]] <- r$gain3
    w[[paste0("ll|", m)]] <- r$ll_fc_e2
    w[[paste0("closer|", m)]] <- r$closer
    w[[paste0("w|", m)]] <- r$w_fc
  }
  w
}

summarise <- function(w, methods, by) {
  groups <- split(w, w[[by]])
  do.call(rbind, lapply(names(groups), function(gname) {
    d <- groups[[gname]]
    rows <- lapply(methods, function(m) {
      ci <- boot(d, function(b) mean(b[[paste0("gain|", m)]]), "type")
      data.frame(subset = gname, method = m, rows = nrow(d), trajectories = length(unique(d$id)),
                 lineage_types = length(unique(d$type)), gain = mean(d[[paste0("gain|", m)]]), low = ci[1L],
                 high = ci[2L], gain_eps_0.001 = mean(d[[paste0("gain3|", m)]]),
                 closer_than_parent = mean(d[[paste0("closer|", m)]]),
                 excess_over_noise = stats::median(d[[paste0("w|", m)]] - d$noise, na.rm = TRUE),
                 parent_excess_over_noise = stats::median(d$w_p0 - d$noise, na.rm = TRUE), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }))
}
paired <- function(w, a, b, by) {
  groups <- c(list(all = w), split(w, w[[by]]))
  do.call(rbind, lapply(names(groups), function(gname) {
    d <- groups[[gname]]
    stat <- function(z) mean(z[[paste0("ll|", a)]] - z[[paste0("ll|", b)]])
    rbind(
      data.frame(subset = gname, comparison = paste(a, "-", b), cluster = "lineage type", rows = nrow(d),
                 clusters = length(unique(d$type)), difference = stat(d), low = boot(d, stat, "type")[1L],
                 high = boot(d, stat, "type")[2L], share_a_better = mean(d[[paste0("ll|", a)]] > d[[paste0("ll|", b)]]),
                 stringsAsFactors = FALSE),
      data.frame(subset = gname, comparison = paste(a, "-", b), cluster = "sample", rows = nrow(d),
                 clusters = length(unique(d$sample)), difference = stat(d), low = boot(d, stat, "sample")[1L],
                 high = boot(d, stat, "sample")[2L], share_a_better = mean(d[[paste0("ll|", a)]] > d[[paste0("ll|", b)]]),
                 stringsAsFactors = FALSE))
  }))
}
show <- function(t, title) {
  cat("\n==", title, "==\n")
  num <- vapply(t, is.numeric, logical(1))
  t[num] <- lapply(t[num], function(v) signif(v, 3))
  print(t, row.names = FALSE)
}

# 1. alfak2 (default) against alfakR.
two <- c("alfak2", "alfakR")
w2 <- wide(two)
s_h <- summarise(w2, two, "horizon")
p_h <- paired(w2, "alfak2", "alfakR", "horizon")
show(s_h, "Log-score gain over persistence per cell, by horizon (passages ahead)")
show(p_h, "alfak2 minus alfakR, log score per cell")
utils::write.csv(s_h, file.path(out_dir, "gain_by_horizon.csv"), row.names = FALSE)
utils::write.csv(p_h, file.path(out_dir, "paired_by_horizon.csv"), row.names = FALSE)
p_m <- paired(w2, "alfak2", "alfakR", "minobs")
show(p_m, "alfak2 minus alfakR, by minobs")
utils::write.csv(p_m, file.path(out_dir, "paired_by_minobs.csv"), row.names = FALSE)
p_t <- paired(w2, "alfak2", "alfakR", "treatment")
show(p_t, "alfak2 minus alfakR, by treatment between parent and later passage")
utils::write.csv(p_t, file.path(out_dir, "paired_by_treatment.csv"), row.names = FALSE)
w2$input_passages <- ifelse(w2$passages == 2L, "one interval", "more intervals")
p_p <- paired(w2, "alfak2", "alfakR", "input_passages")
show(p_p, "alfak2 minus alfakR, by the input's intervals")
utils::write.csv(p_p, file.path(out_dir, "paired_by_input.csv"), row.names = FALSE)

# 2. The three arms on the rows all three share.
three <- c("alfak2", "alfakR", "alfak2 free")
w3 <- wide(three)
s3 <- summarise(w3, three, "horizon")
show(s3, "Three arms, common rows")
p3 <- rbind(paired(w3, "alfak2 free", "alfak2", "horizon"), paired(w3, "alfak2 free", "alfakR", "horizon"))
show(p3, "alfak2 free against the others")
utils::write.csv(s3, file.path(out_dir, "gain_three_arms.csv"), row.names = FALSE)
utils::write.csv(p3, file.path(out_dir, "paired_free.csv"), row.names = FALSE)
saveRDS(list(two = w2, three = w3), file.path(out_dir, "rows.Rds"))

# 3. Figure: gain over persistence by horizon, all methods in one panel; and
# the paired difference alfak2 - alfakR.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  theme_set(theme_bw(base_size = 11))
  s3$horizon <- factor(s3$subset, levels = c("1", "2", "3-4"))
  a <- ggplot(s3, aes(horizon, gain, colour = method, group = method)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4)) +
    labs(x = "passages ahead of the parent", y = "log score over persistence (per cell)", colour = NULL,
         title = "a  Forecasts against later passages")
  d <- p_h[p_h$cluster == "lineage type" & p_h$subset != "all", ]
  d$horizon <- factor(d$subset, levels = c("1", "2", "3-4"))
  b <- ggplot(d, aes(horizon, difference)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = low, ymax = high)) +
    labs(x = "passages ahead of the parent", y = "alfak2 minus alfakR (per cell)", title = "b  Paired difference")
  for (ext in c("png", "pdf")) {
    f <- file.path(out_dir, "figures", paste0("horizon.", ext))
    if (ext == "png") grDevices::png(f, width = 10, height = 4, units = "in", res = 150) else grDevices::pdf(f, width = 10, height = 4)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(0.6, 0.4), "npc"))))
    print(a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grDevices::dev.off()
  }
}
