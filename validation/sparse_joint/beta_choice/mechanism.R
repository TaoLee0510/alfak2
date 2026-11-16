#!/usr/bin/env Rscript
# Short-interval plan, step 1: what makes a dataset prefer a higher
# missegregation rate when the true rate is the default?
#
# Reads the candidate fits of the two-time-point method (`ends`, `ends_b4`,
# `ends_b16`, ...) on series_methods sets e1-e3 (build cv10), whose true rate
# is the default. For each dataset: gain4 = objective(x1) - objective(x4) (> 0:
# the likelihood prefers four times the true rate), and features of the
# counts of its first and last time point:
#   cells           mean cells per time point
#   n_kary          karyotypes seen
#   singleton_share share of karyotypes seen in exactly one cell over both
#   new_share       share of the last time point's karyotypes not seen at the
#                   first
#   top_share       share of cells on the most frequent karyotype
#   turnover        total-variation distance between the two compositions
# plus the design (bottleneck, spacing, interval) and the stage-2 drift index.
#
#   Rscript validation/sparse_joint/beta_choice/mechanism.R
#
# Writes validation/sparse_joint/beta_choice/results/mechanism.md.

root <- "validation/sparse_joint/series_methods/results/cv10"
out <- "validation/sparse_joint/beta_choice/results"
mult <- c(ends_b1o16 = 1 / 16, ends_b1o4 = 1 / 4, ends = 1, ends_b4 = 4, ends_b16 = 16, ends_b64 = 64)

rows <- do.call(rbind, lapply(c("e1", "e2", "e3"), function(set) {
  idx <- utils::read.csv(file.path(root, set, "dataset_index.csv"), stringsAsFactors = FALSE)
  do.call(rbind, lapply(idx$dataset_id, function(id) {
    obj <- vapply(names(mult), function(m) {
      f <- readRDS(file.path(root, set, "fits", m, paste0(id, ".rds")))
      o <- f$diagnostics$objective
      if (identical(f$status, "ok") && length(o)) as.numeric(o) else NA_real_
    }, numeric(1))
    d <- readRDS(file.path(root, set, "datasets", paste0(id, ".rds")))
    x <- d$counts[, c(1L, ncol(d$counts)), drop = FALSE]
    x <- x[rowSums(x) > 0, , drop = FALSE]
    p0 <- x[, 1] / sum(x[, 1])
    p1 <- x[, 2] / sum(x[, 2])
    i <- idx[idx$dataset_id == id, ]
    data.frame(set = set, dataset_id = id, landscape = d$landscape_id, bottleneck = d$bottleneck,
               spacing = d$spacing, interval = max(d$times), depth = d$depth,
               drift_index = if ("drift_index" %in% names(i)) i$drift_index else NA_real_,
               cells = mean(colSums(x)), n_kary = nrow(x), singleton_share = mean(rowSums(x) == 1),
               new_share = mean(x[, 1] == 0 & x[, 2] > 0) / mean(x[, 2] > 0),
               top_share = max(rowSums(x)) / sum(x), turnover = sum(abs(p0 - p1)) / 2,
               gain4 = obj[["ends"]] - obj[["ends_b4"]], gain16 = obj[["ends"]] - obj[["ends_b16"]],
               pick = mult[[which.min(obj)]], stringsAsFactors = FALSE)
  }))
}))
saveRDS(rows, file.path(out, "mechanism_rows.rds"))

md <- function(d) {
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
rows$high <- rows$pick >= 4
by_design <- stats::aggregate(cbind(high, gain4, singleton_share, new_share, turnover) ~ set + depth + spacing + bottleneck,
                              rows, stats::median)
by_design$high <- stats::aggregate(high ~ set + depth + spacing + bottleneck, rows, mean)$high
# Standardised effects on gain4 (all sets together).
z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)
fit_all <- stats::lm(gain4 ~ z(log(cells)) + z(log(interval)) + z(log(bottleneck)) + z(singleton_share) +
                       z(new_share) + z(top_share) + z(turnover) + z(log(n_kary)), rows)
coef_all <- as.data.frame(summary(fit_all)$coefficients)
coef_all <- cbind(term = rownames(coef_all), coef_all)
# The same within the shallow, short-spacing cells, where most wrong picks are.
sub <- rows[rows$depth == 200 & rows$spacing == 6, ]
fit_sub <- stats::lm(gain4 ~ z(log(bottleneck)) + z(singleton_share) + z(new_share) + z(top_share) + z(turnover) +
                       z(log(n_kary)) + z(log(interval)), sub)
coef_sub <- as.data.frame(summary(fit_sub)$coefficients)
coef_sub <- cbind(term = rownames(coef_sub), coef_sub)

writeLines(c("# Short-interval plan, step 1: what drives a preference for a higher rate (true rate = default)", "",
             "gain4 = objective at the true rate minus at 4x (> 0 prefers 4x); high = share picking >= 4x.", "",
             "## By design (medians; high = share)", "", md(by_design), "",
             sprintf("## gain4 on standardised features, all sets (R2 %.2f, n %d)", summary(fit_all)$r.squared, nrow(rows)),
             "", md(coef_all), "",
             sprintf("## The same within 200 cells and 6-day spacing (R2 %.2f, n %d)", summary(fit_sub)$r.squared, nrow(sub)),
             "", md(coef_sub)),
           file.path(out, "mechanism.md"))
cat(readLines(file.path(out, "mechanism.md")), sep = "\n")
