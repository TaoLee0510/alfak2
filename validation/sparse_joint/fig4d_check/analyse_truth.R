#!/usr/bin/env Rscript
# Fig. 4d check, step 2 summary: estimated against true delta f between
# frequent karyotypes and their one-chromosome neighbours in the paper's
# simulations (rows from abm_truth.R, copied to results/truth/).
#
#   Rscript validation/sparse_joint/fig4d_check/analyse_truth.R
#
# Writes results/truth.md.

base <- "validation/sparse_joint/fig4d_check/results"
files <- list.files(file.path(base, "truth"), pattern = "^truth_rows_.*\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(files, readRDS))
x$window <- sub("^w([0-9]+)_.*", "\\1", x$design)
x$kind <- ifelse(grepl("_ends$", x$design), "2 time points", "all time points")

md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
summ <- function(d) data.frame(n = nrow(d), est = stats::median(d$est), true = stats::median(d$true),
                               abs_est = stats::median(abs(d$est)), abs_true = stats::median(abs(d$true)),
                               bias = stats::median(d$est - d$true),
                               slope = unname(stats::coef(stats::lm(est ~ true, d))[2]),
                               r = stats::cor(d$est, d$true))
group <- function(d, by) {
  keys <- unique(d[by])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    m <- Reduce(`&`, lapply(by, function(b) d[[b]] == keys[[b]][i]))
    if (sum(m) < 20) return(NULL)
    cbind(keys[i, , drop = FALSE], summ(d[m, ]))
  }))
  out[do.call(order, unname(as.list(out[by]))), ]
}

# The alfak2 source of each neighbour, for the alfakR rows too (same fit design).
src <- unique(x[x$method == "alfak2", c("source_set", "sim_id", "design", "k_id", "nn", "source")])
names(src)[names(src) == "source"] <- "alfak2_source"
x <- merge(x, src, by = c("source_set", "sim_id", "design", "k_id", "nn"), all.x = TRUE)
both <- x[!is.na(x$alfak2_source), ]

t1 <- group(both, c("source_set", "kind", "method"))
t2 <- group(both[both$kind == "2 time points", ], c("source_set", "alfak2_source", "method"))
t3 <- group(both[both$kind == "2 time points" & both$source_set == "abm", ], c("window", "alfak2_source", "method"))
fits <- unique(x[x$method == "alfak2", c("source_set", "design", "sim_id", "history", "dt_total", "cells_mean")])
t4 <- aggregate(cbind(history, dt_total, cells_mean) ~ source_set + design, fits, stats::median)

writeLines(c("# Fig. 4d check, step 2: estimated against true delta f (minobs 10)", "",
             "delta f = f(frequent) - f(neighbour). est and true are medians; bias the median of est - true;",
             "slope and r from est ~ true over all pairs of a group.", "",
             "## By method and design", "", md(t1), "",
             "## Two time points, by where alfak2's value of the neighbour comes from", "", md(t2), "",
             "## ABM, two time points, by window", "", md(t3), "",
             "## The alfak2 fits", "", md(t4)),
           file.path(base, "truth.md"))
cat(readLines(file.path(base, "truth.md")), sep = "\n")
