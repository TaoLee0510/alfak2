#!/usr/bin/env Rscript
# Follow-up step 2, read only: why alfak2 ranks the size of change worse than
# alfakR when every passage of a Salehi trajectory is fitted ("all" design).
# Rows from the benchmark's scripts/fig2d_check.R (FIG2D_DESIGN=all), one per
# fit, minobs and daughter passage, on the rows all methods share.
#   Rscript validation/sparse_joint/cloud_level/salehi/all_design_diagnosis.R
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- file.path(script_dir, "..", "results", "salehi")
out_dir <- file.path(res, "report_all")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(dir) {
  f <- list.files(dir, pattern = "^chunk_.*\\.Rds$", full.names = TRUE)
  r <- Filter(function(z) !is.null(z) && nrow(z) > 0, lapply(f, readRDS))
  cc <- unique(unlist(lapply(r, names)))
  do.call(rbind, lapply(r, function(z) { z[setdiff(cc, names(z))] <- NA; z[cc] }))
}
bind <- function(tabs) {
  cc <- unique(unlist(lapply(tabs, names)))
  do.call(rbind, lapply(tabs, function(z) { z[setdiff(cc, names(z))] <- NA; z[cc] }))
}
all <- bind(list(rd(file.path(res, "check_all_round")), rd(file.path(res, "check_all_variants"))))
ends <- rd(file.path(res, "check_ends_round"))
all$design <- "all"; ends$design <- "ends"
x <- bind(list(all, ends))
x$method <- paste(x$design, sub("^alfak2_", "", x$method), sep = ":")
x$key <- paste(x$fi, x$min_obs, x$dj)
sp <- function(a, b) { k <- is.finite(a) & is.finite(b); if (sum(k) < 5) NA_real_ else stats::cor(a[k], b[k], method = "spearman") }
methods <- c("all:alfakR", "all:alfak2", "all:pairwise", "all:state_space", "all:free", "ends:alfakR", "ends:alfak2")
keys <- Reduce(intersect, lapply(methods, function(m) x$key[x$method == m & is.finite(x$pred0) & is.finite(x$base)]))
cat("Common rows:", length(keys), " trajectories:", length(unique(sub(" .*$", "", keys))), "\n")
w <- x[x$method == "all:alfakR" & x$key %in% keys, c("key", "fi", "type", "passages", "base", "base_full", "base_prev",
                                                      "base_prev_all", "span_all", "entropy_parent", "cov_parent")]
for (m in methods) {
  r <- x[x$method == m & x$key %in% keys, ]
  r <- r[match(w$key, r$key), ]
  w[[paste0("pred0|", m)]] <- r$pred0
  w[[paste0("base|", m)]] <- r$base
  w[[paste0("var|", m)]] <- r$var_f
  if (m == "all:alfak2") { w$weight <- r$intermediate_weight; w$fitted <- r$samples_fitted; w$history <- r$history }
}
# 1. Is the observed-change axis the same for all methods?
cat("\nObserved-change axis identical to alfakR's (all design): ")
print(sapply(methods, function(m) mean(abs(w[[paste0("base|", m)]] - w$base) < 1e-9)))
# 2. Spearman with the observed change, by method, and simple predictors.
tab <- data.frame(predictor = c(paste("forecast", methods), "recent change (landscape)", "recent change (all cells)",
                                "total change of the trajectory", "parent entropy", "parent coverage"),
                  spearman = c(sapply(methods, function(m) sp(w$base, w[[paste0("pred0|", m)]])),
                               sp(w$base, w$base_prev), sp(w$base, w$base_prev_all), sp(w$base, w$span_all),
                               sp(w$base, w$entropy_parent), sp(w$base, w$cov_parent)))
utils::write.csv(tab, file.path(out_dir, "spearman.csv"), row.names = FALSE)
cat("\n== Spearman with the observed change (common rows) ==\n"); print(tab, digits = 3, row.names = FALSE)
# 3. What does each forecast follow?
feat <- c("base_prev", "base_prev_all", "span_all", "entropy_parent", "passages")
follow <- do.call(rbind, lapply(methods, function(m) {
  p <- w[[paste0("pred0|", m)]]
  data.frame(method = m, t(sapply(feat, function(f) round(sp(p, w[[f]]), 3))),
             with_alfakR_all = round(sp(p, w[["pred0|all:alfakR"]]), 3), with_alfak2_all = round(sp(p, w[["pred0|all:alfak2"]]), 3))
}))
utils::write.csv(follow, file.path(out_dir, "forecast_follows.csv"), row.names = FALSE)
cat("\n== Spearman of each forecast's size with features ==\n"); print(follow, row.names = FALSE)
# 4. Partial: alfakR's advantage after the recent change.
partial <- function(y, a, z) {
  ok <- is.finite(y) & is.finite(a) & is.finite(z)
  ry <- stats::residuals(stats::lm(rank(y[ok]) ~ rank(z[ok]))); ra <- stats::residuals(stats::lm(rank(a[ok]) ~ rank(z[ok])))
  stats::cor(ry, ra)
}
pt <- data.frame(method = methods, partial_given_recent = sapply(methods, function(m) partial(w$base, w[[paste0("pred0|", m)]], w$base_prev_all)),
                 partial_given_total = sapply(methods, function(m) partial(w$base, w[[paste0("pred0|", m)]], w$span_all)))
cat("\n== Partial Spearman with the observed change ==\n"); print(pt, digits = 3, row.names = FALSE)
# 5. alfak2 anchored by intermediate weight and samples dropped.
w$wgroup <- cut(w$weight, c(-Inf, 0.5, 0.99, Inf), labels = c("< 0.5", "0.5-0.99", "1"))
w$dropped <- w$fitted < w$passages
g <- do.call(rbind, lapply(split(w, list(w$wgroup, w$dropped), drop = TRUE), function(d) data.frame(
  weight = d$wgroup[1], dropped_samples = d$dropped[1], n = nrow(d), trajectories = length(unique(d$fi)),
  alfak2 = sp(d$base, d[["pred0|all:alfak2"]]), free = sp(d$base, d[["pred0|all:free"]]),
  state_space = sp(d$base, d[["pred0|all:state_space"]]), alfakR = sp(d$base, d[["pred0|all:alfakR"]]))))
cat("\n== By alfak2's intermediate weight ==\n"); print(g, digits = 3, row.names = FALSE)
# 6. By number of passages.
w$plen <- cut(w$passages, c(2, 3, 4, 20), labels = c("3", "4", "5+"))
h <- do.call(rbind, lapply(split(w, w$plen), function(d) data.frame(passages = d$plen[1], n = nrow(d),
  t(sapply(methods, function(m) round(sp(d$base, d[[paste0("pred0|", m)]]), 3))), check.names = FALSE)))
cat("\n== By passages ==\n"); print(h, row.names = FALSE)
saveRDS(w, file.path(out_dir, "rows.Rds"))
