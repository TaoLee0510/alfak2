#!/usr/bin/env Rscript
# Follow-up step 2: sampling noise against the observed change of Fig. 2d
# ("all" design, the 324 rows every method has). Reads report_all/rows.Rds
# (all_design_diagnosis.R) and report_all/noise_check.csv (noise_check.R).
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
d <- file.path(script_dir, "..", "results", "salehi", "report_all")
w <- readRDS(file.path(d, "rows.Rds"))
n <- utils::read.csv(file.path(d, "noise_check.csv"), stringsAsFactors = FALSE)
n$key <- paste(n$fi, n$min_obs, n$dj)
w <- merge(w, n[c("key", "cells_parent", "cells_daughter", "noise", "noise_all", "noise_daughter")], by = "key")
sp <- function(a, b) { k <- is.finite(a) & is.finite(b); stats::cor(a[k], b[k], method = "spearman") }
partial <- function(y, a, z) {
  ok <- is.finite(y) & is.finite(a) & stats::complete.cases(z)
  z <- as.data.frame(lapply(as.data.frame(z)[ok, , drop = FALSE], rank))
  ry <- stats::residuals(stats::lm(rank(y[ok]) ~ ., data = z)); ra <- stats::residuals(stats::lm(rank(a[ok]) ~ ., data = z))
  stats::cor(ry, ra)
}
# Lineage-type bootstrap of a Spearman.
boot <- function(a, b, draws = 4000L) {
  g <- split(seq_len(nrow(w)), w$type); set.seed(1)
  v <- vapply(seq_len(draws), function(i) { s <- unlist(g[sample.int(length(g), replace = TRUE)], use.names = FALSE); sp(a[s], b[s]) }, numeric(1))
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
noise_mean <- (w$noise + w$noise_daughter) / 2
pred <- c("all:alfakR", "all:alfak2", "all:state_space", "all:free")
tab <- data.frame(
  predictor = c("noise: halves of the parent passage", "noise: halves of the daughter passage", "noise: mean of the two",
                "recent change (passage before -> parent)", paste("forecast", pred)),
  spearman = c(sp(w$base, w$noise), sp(w$base, w$noise_daughter), sp(w$base, noise_mean), sp(w$base, w$base_prev),
               sapply(pred, function(m) sp(w$base, w[[paste0("pred0|", m)]]))))
ci <- t(sapply(list(w$noise, w$noise_daughter, noise_mean, w$base_prev, w[["pred0|all:alfakR"]], w[["pred0|all:alfak2"]],
                    w[["pred0|all:state_space"]], w[["pred0|all:free"]]), function(p) boot(w$base, p)))
tab$low <- ci[, 1]; tab$high <- ci[, 2]
cat("Rows:", nrow(w), "\n\n== Spearman with the observed change ==\n"); print(tab, digits = 3, row.names = FALSE)
cat("\nRecent change vs parent noise:", round(sp(w$base_prev, w$noise), 3), "\n")
cat("alfakR forecast size vs parent noise:", round(sp(w[["pred0|all:alfakR"]], w$noise), 3),
    "  alfak2:", round(sp(w[["pred0|all:alfak2"]], w$noise), 3), "\n")
pt <- data.frame(method = pred,
  given_noise = sapply(pred, function(m) partial(w$base, w[[paste0("pred0|", m)]], data.frame(noise_mean))),
  given_noise_and_recent = sapply(pred, function(m) partial(w$base, w[[paste0("pred0|", m)]], data.frame(noise_mean, w$base_prev))))
cat("\n== Partial Spearman with the observed change, sampling noise held fixed ==\n"); print(pt, digits = 3, row.names = FALSE)
cat("\nObserved change over noise (median):", round(stats::median(w$base / noise_mean, na.rm = TRUE), 2),
    "  share of rows where the observed change is within 1.5x the noise:", round(mean(w$base < 1.5 * noise_mean, na.rm = TRUE), 2), "\n")
utils::write.csv(tab, file.path(d, "noise_spearman.csv"), row.names = FALSE)
utils::write.csv(pt, file.path(d, "noise_partial.csv"), row.names = FALSE)
