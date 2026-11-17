#!/usr/bin/env Rscript
# Follow-up step 2: the sampling-noise check on the trajectories with one
# interval (two passages, "ends" design), for the arms of the first round.
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- file.path(script_dir, "..", "results", "salehi")
rd <- function(dir) {
  f <- list.files(dir, pattern = "^chunk_.*[.]Rds$", full.names = TRUE)
  r <- Filter(function(z) !is.null(z) && nrow(z) > 0, lapply(f, readRDS))
  do.call(rbind, lapply(r, function(z) z[c("fi", "method", "min_obs", "dj", "type", "base", "pred0")]))
}
x <- rbind(rd("~/Documents/GitHub/alfak2-paper-benchmark/results_05a8855/analysis/fig2d_check/salehi"),
           rd(file.path(res, "check")), rd(file.path(res, "variants_check")))
x$method <- sub("^alfak2_", "", x$method)
x$key <- paste(x$fi, x$min_obs, x$dj)
n <- utils::read.csv(file.path(res, "report_all", "noise_one_interval.csv"), stringsAsFactors = FALSE)
n$key <- paste(n$fi, n$min_obs, n$dj)
methods <- c("alfakR", "alfak2", "free", "level_dt")
keys <- Reduce(intersect, c(list(n$key), lapply(methods, function(m) x$key[x$method == m & is.finite(x$pred0)])))
w <- n[match(keys, n$key), ]
w$base <- x$base[x$method == "alfakR"][match(keys, x$key[x$method == "alfakR"])]
w$type <- x$type[x$method == "alfakR"][match(keys, x$key[x$method == "alfakR"])]
for (m in methods) w[[m]] <- x$pred0[x$method == m][match(keys, x$key[x$method == m])]
sp <- function(a, b) { k <- is.finite(a) & is.finite(b); stats::cor(a[k], b[k], method = "spearman") }
partial <- function(y, a, z) {
  ok <- is.finite(y) & is.finite(a) & is.finite(z)
  ry <- stats::residuals(stats::lm(rank(y[ok]) ~ rank(z[ok]))); ra <- stats::residuals(stats::lm(rank(a[ok]) ~ rank(z[ok])))
  stats::cor(ry, ra)
}
noise <- (w$noise + w$noise_daughter) / 2
cat("Rows:", nrow(w), " trajectories:", length(unique(w$fi)), "\n")
cat("Spearman with the observed change: noise parent", round(sp(w$base, w$noise), 3), " noise daughter",
    round(sp(w$base, w$noise_daughter), 3), " noise mean", round(sp(w$base, noise), 3), "\n")
tab <- data.frame(method = methods, spearman = sapply(methods, function(m) sp(w$base, w[[m]])),
                  with_noise = sapply(methods, function(m) sp(w[[m]], noise)),
                  partial_given_noise = sapply(methods, function(m) partial(w$base, w[[m]], noise)))
print(tab, digits = 3, row.names = FALSE)
cat("Observed change over noise (median):", round(stats::median(w$base / noise, na.rm = TRUE), 2), "\n")
utils::write.csv(tab, file.path(res, "report_all", "noise_one_interval_summary.csv"), row.names = FALSE)

# Lineage-type bootstrap of the difference in partial Spearman (noise held
# fixed) between each arm and alfakR, here and on the all-design rows.
boot_diff <- function(d, y, a, b, z, cluster = "type", draws = 4000L) {
  g <- split(seq_len(nrow(d)), d[[cluster]]); set.seed(1)
  f <- function(s) partial(d[[y]][s], d[[a]][s], z[s]) - partial(d[[y]][s], d[[b]][s], z[s])
  v <- vapply(seq_len(draws), function(i) f(unlist(g[sample.int(length(g), replace = TRUE)], use.names = FALSE)), numeric(1))
  c(difference = f(seq_len(nrow(d))), stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}
cat("\nOne interval, partial given noise, difference from alfakR (lineage bootstrap):\n")
print(round(t(sapply(c("alfak2", "free", "level_dt"), function(m) boot_diff(w, "base", m, "alfakR", noise))), 3))
a <- readRDS(file.path(res, "report_all", "rows.Rds"))
na <- utils::read.csv(file.path(res, "report_all", "noise_check.csv"), stringsAsFactors = FALSE)
na$key <- paste(na$fi, na$min_obs, na$dj)
a <- merge(a, na[c("key", "noise", "noise_daughter")], by = "key")
za <- (a$noise + a$noise_daughter) / 2
names(a) <- sub("^pred0[|]all:", "all_", names(a))
cat("\nAll time points, partial given noise, difference from alfakR (lineage bootstrap):\n")
print(round(t(sapply(c("all_alfak2", "all_free", "all_pairwise", "all_state_space"), function(m) boot_diff(a, "base", m, "all_alfakR", za))), 3))
