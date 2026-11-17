#!/usr/bin/env Rscript
# Salehi forecasts against later passages, for fits of every passage of a
# trajectory ("all" design, trajectories of three or more passages): alfak2
# anchored (the default), pairwise and state_space, alfakR, and alfak2 "free",
# and beside them the two-time-point ("ends") fits of the same trajectories.
# Score and intervals as in analyse.R: log-likelihood per cell of the later
# passage under the forecast minus under persistence (mixture weight 0.01),
# paired on the rows (trajectory, minobs, later passage) the compared fits
# share; bootstrap over lineage types and over samples.
#   Rscript validation/sparse_joint/salehi_horizon/analyse_all.R
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- file.path(script_dir, "results")
out_dir <- file.path(res, "report_all")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
read_rows <- function(dir, design) {
  parts <- lapply(list.files(dir, pattern = "[.]Rds$", full.names = TRUE), readRDS)
  status <- data.frame(design = design, method = vapply(parts, function(p) p$task$method, character(1)),
                       status = vapply(parts, `[[`, character(1), "status"))
  rows <- do.call(rbind, lapply(parts, function(p) {
    r <- p$rows
    if (is.null(r)) return(NULL)
    r$design <- design
    r[c("id", "design", "method", "minobs", "uid", "depth", "days", "treat_parent", "treat_desc", "cells_desc",
        "coverage", "ll_fc_e2", "ll_p0_e2", "ll_fc_e3", "ll_p0_e3", "w_fc", "w_p0", "noise")]
  }))
  list(rows = rows, status = status)
}
a <- read_rows(file.path(res, "rows_all"), "all")
e <- read_rows(file.path(res, "rows"), "ends")
st <- rbind(a$status, e$status)
cat("Tasks:\n"); print(stats::aggregate(n ~ design + method + status, data = transform(st, n = 1), FUN = sum), row.names = FALSE)
x <- rbind(a$rows, e$rows)
label <- c(alfak2 = "anchored", alfak2_pairwise = "pairwise", alfak2_state_space = "state_space", alfakR = "alfakR",
           alfak2_free = "free")
x$arm <- paste(x$design, label[x$method], sep = ":")
lineage_map <- c(SA1035T = "SA1035", SA1035U = "SA1035", SA532 = "SA532", SA609 = "SA609", SA609UnBU = "SA609",
                 SA906a = "p53 k.o", SA906b = "p53 k.o", SA039U = "p53 w.t", SA535_CISPLATIN_CombinedT = "SA535",
                 SA535_CISPLATIN_CombinedU = "SA535", SA000 = "SA609", SA609R2 = "SA609",
                 SA535_CISPLATIN_CombinedH = "SA535", SA001 = "SA609", SA609H2 = "SA609", SA609H1 = "SA609",
                 SA609R2b = "SA609", SA609TIIDT = "SA609")
x$sample <- sub("_X[0-9]+_l_.*$", "", x$id)
x$type <- unname(lineage_map[x$sample])
x$gain <- x$ll_fc_e2 - x$ll_p0_e2
x$gain3 <- x$ll_fc_e3 - x$ll_p0_e3
x$closer <- x$w_fc < x$w_p0
x$horizon <- ifelse(x$depth >= 3L, "3-4", as.character(x$depth))
x$key <- paste(x$id, x$minobs, x$uid)
x <- x[as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", x$id)) >= 3L, ]
cat("\nRows by arm and horizon:\n"); print(table(x$arm, x$horizon))

boot <- function(d, statistic, cluster, draws = 4000L) {
  g <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  v <- vapply(seq_len(draws), function(i) statistic(d[unlist(g[sample.int(length(g), replace = TRUE)],
                                                        use.names = FALSE), , drop = FALSE]), numeric(1))
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
wide <- function(arms) {
  keys <- Reduce(intersect, lapply(arms, function(m) x$key[x$arm == m & is.finite(x$gain)]))
  w <- x[x$arm == arms[[1L]] & x$key %in% keys, c("key", "id", "minobs", "uid", "horizon", "type", "sample")]
  for (m in arms) {
    r <- x[x$arm == m & x$key %in% keys, ]
    r <- r[match(w$key, r$key), ]
    w[[paste0("gain|", m)]] <- r$gain
    w[[paste0("ll|", m)]] <- r$ll_fc_e2
    w[[paste0("closer|", m)]] <- r$closer
  }
  w
}
summarise <- function(w, arms) {
  do.call(rbind, lapply(c("all", sort(unique(w$horizon))), function(h) {
    d <- if (h == "all") w else w[w$horizon == h, ]
    do.call(rbind, lapply(arms, function(m) {
      ci <- boot(d, function(b) mean(b[[paste0("gain|", m)]]), "type")
      data.frame(horizon = h, arm = m, rows = nrow(d), trajectories = length(unique(d$id)),
                 gain = mean(d[[paste0("gain|", m)]]), low = ci[1L], high = ci[2L],
                 closer_than_parent = mean(d[[paste0("closer|", m)]]), stringsAsFactors = FALSE)
    }))
  }))
}
paired <- function(w, pairs) {
  do.call(rbind, lapply(c("all", sort(unique(w$horizon))), function(h) {
    d <- if (h == "all") w else w[w$horizon == h, ]
    do.call(rbind, lapply(pairs, function(p) {
      # Gains over persistence, not log-likelihoods: the two designs choose
      # their landscapes' karyotypes from different time points, so their
      # likelihoods are on different supports; within one design the two are
      # the same difference.
      stat <- function(z) mean(z[[paste0("gain|", p[1])]] - z[[paste0("gain|", p[2])]])
      lt <- boot(d, stat, "type"); sa <- boot(d, stat, "sample")
      data.frame(horizon = h, comparison = paste(p[1], "-", p[2]), rows = nrow(d), difference = stat(d),
                 low_lineage = lt[1L], high_lineage = lt[2L], low_sample = sa[1L], high_sample = sa[2L],
                 share_first_better = mean(d[[paste0("gain|", p[1])]] > d[[paste0("gain|", p[2])]]),
                 stringsAsFactors = FALSE)
    }))
  }))
}
show <- function(t, title) {
  cat("\n==", title, "==\n")
  num <- vapply(t, is.numeric, logical(1))
  t[num] <- lapply(t[num], function(v) signif(v, 3))
  print(t, row.names = FALSE)
}

# 1. The fits of every passage: alfak2's three methods against alfakR.
arms <- c("all:anchored", "all:pairwise", "all:state_space", "all:alfakR")
w <- wide(arms)
s1 <- summarise(w, arms)
p1 <- paired(w, list(c("all:anchored", "all:alfakR"), c("all:pairwise", "all:alfakR"),
                     c("all:state_space", "all:alfakR"), c("all:state_space", "all:anchored"),
                     c("all:pairwise", "all:anchored")))
show(s1, "All time points: log score over persistence per cell")
show(p1, "All time points: paired differences")
utils::write.csv(s1, file.path(out_dir, "gain_all.csv"), row.names = FALSE)
utils::write.csv(p1, file.path(out_dir, "paired_all.csv"), row.names = FALSE)

# 2. Does fitting every passage help? Each method against its two-time-point
# fit of the same trajectory, and the best of each design.
arms2 <- c("all:anchored", "all:state_space", "all:alfakR", "ends:anchored", "ends:alfakR")
w2 <- wide(arms2)
s2 <- summarise(w2, arms2)
p2 <- paired(w2, list(c("all:anchored", "ends:anchored"), c("all:state_space", "ends:anchored"),
                      c("all:alfakR", "ends:alfakR"), c("ends:anchored", "ends:alfakR")))
show(s2, "All time points against two time points (same trajectories)")
show(p2, "Paired: every passage against two time points")
utils::write.csv(s2, file.path(out_dir, "gain_all_vs_ends.csv"), row.names = FALSE)
utils::write.csv(p2, file.path(out_dir, "paired_all_vs_ends.csv"), row.names = FALSE)

# 3. alfak2 "free" (every passage), where its fits exist.
arms3 <- c("all:anchored", "all:free", "all:alfakR")
w3 <- wide(arms3)
p3 <- paired(w3, list(c("all:free", "all:anchored"), c("all:free", "all:alfakR")))
show(p3, "alfak2 free, every passage")
utils::write.csv(p3, file.path(out_dir, "paired_all_free.csv"), row.names = FALSE)
saveRDS(list(all = w, all_vs_ends = w2, free = w3), file.path(out_dir, "rows.Rds"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  theme_set(theme_bw(base_size = 11))
  s <- s2[s2$horizon != "all", ]
  s$horizon <- factor(s$horizon, levels = c("1", "2", "3-4"))
  s$arm <- factor(s$arm, levels = arms2)
  f1 <- ggplot(s, aes(horizon, gain, colour = arm, group = arm)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    labs(x = "passages ahead of the parent", y = "log score over persistence (per cell)", colour = NULL,
         title = "Trajectories of 3+ passages: every passage (all) and two (ends)")
  ggsave(file.path(out_dir, "figures", "horizon_all.png"), f1, width = 8, height = 4.5, dpi = 150)
  ggsave(file.path(out_dir, "figures", "horizon_all.pdf"), f1, width = 8, height = 4.5)
}
