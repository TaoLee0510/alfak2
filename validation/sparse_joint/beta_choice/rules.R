#!/usr/bin/env Rscript
# Short-interval plan, step 3: rules that decide when to leave the default
# missegregation rate, evaluated offline on the candidate fits.
#
# Each dataset has a fit at 5e-5 x {1/16, 1/4, 1, 4, 16, 64} (objective and
# accuracy). A rule maps the candidate objectives and two features of the
# counts (interval between the first and last time point, mean cells per
# time point) to one multiple; abstaining is x1, the default.
#   R0  lowest objective (beta = "auto" of ac36fd1)
#   A   R0 only if interval >= I days and cells >= C, else x1
#   B   R0 only if its objective is at least T below x1's, else x1
#   C   lowest objective + (log(multiple) / s)^2 / 2 (lognormal prior on the
#       rate around the default, sd s on the log scale)
# Data:
#   null         series_methods e1-e3 (true rate = default), methods ends and
#                anchored with their candidates (build cv10)
#   alternatives the paper's simulations of step 0 (sim_fit.R; landscape r),
#                and series_methods e4 (true rate 4x or 16x; RMSE), ends only
#                (the two-time-point fit; anchored was not fitted on e4)
# Tuning half: landscapes mt_L01-L06, step-0 replicates 01-12; test half:
# mt_L07-L12, replicates 13-25.
# Criteria, fixed before scoring (memory alfak2-beta-shortinterval-plan):
#   (1) null: e1 >= 18 of 24 cells with median RMSE ratio (rule / x1) <= 1.03
#       at every distance 0-2 and the median over cells <= 1.02; e2 and e3
#       medians <= 1.03 at every distance; for ends and for anchored;
#   (2) alternatives: >= 70% of picks within 4x of the truth; where the truth
#       is >= 4x the default, accuracy better than x1 (step 0: paired median
#       r(rule) - r(x1) > 0; e4: median RMSE ratio < 1 at every distance).
# Among the rules passing on the tuning half, the one with the best
# alternatives (lowest e4 median ratio, then highest step-0 gain) is chosen,
# simpler first on a tie (A, B, C before combinations; coarser parameters
# first), and judged on the test half.
#
#   Rscript validation/sparse_joint/beta_choice/rules.R
#
# Writes validation/sparse_joint/beta_choice/results/rules.md and rules.rds.

root <- "validation/sparse_joint/series_methods/results/cv10"
out <- "validation/sparse_joint/beta_choice/results"
grid <- c(1 / 16, 1 / 4, 1, 4, 16, 64)
suffix <- c("_b1o16", "_b1o4", "", "_b4", "_b16", "_b64")

# ---- candidate tables ------------------------------------------------------------
series_candidates <- function(set, base) {
  if (!dir.exists(file.path(root, set))) return(NULL)
  idx <- utils::read.csv(file.path(root, set, "dataset_index.csv"), stringsAsFactors = FALSE)
  dist <- utils::read.csv(file.path(root, set, "metrics_distance.csv"), stringsAsFactors = FALSE)
  dist <- dist[dist$stratum %in% c("d0", "d1", "d2") & dist$status == "ok", c("dataset_id", "method", "stratum", "rmse")]
  rows <- do.call(rbind, lapply(idx$dataset_id, function(id) {
    d <- readRDS(file.path(root, set, "datasets", paste0(id, ".rds")))
    obj <- vapply(paste0(base, suffix), function(m) {
      p <- file.path(root, set, "fits", m, paste0(id, ".rds"))
      if (!file.exists(p)) return(NA_real_)
      f <- readRDS(p)
      o <- f$diagnostics$objective
      if (identical(f$status, "ok") && length(o)) as.numeric(o) else NA_real_
    }, numeric(1))
    x <- d$counts[, c(1L, ncol(d$counts)), drop = FALSE]
    data.frame(source = set, base = base, id = id, landscape = d$landscape_id, interval = max(d$times),
               cells = mean(colSums(x)), true_multiple = d$beta / 5e-5, multiple = grid, objective = obj,
               method = paste0(base, suffix), stringsAsFactors = FALSE)
  }))
  rmse <- reshape(dist, idvar = c("dataset_id", "method"), timevar = "stratum", direction = "wide")
  names(rmse) <- sub("^rmse[.]", "rmse_", names(rmse))
  merge(rows, rmse, by.x = c("id", "method"), by.y = c("dataset_id", "method"), all.x = TRUE)
}
step0_candidates <- function() {
  files <- list.files(file.path(out, "fits"), pattern = "[.]rds$", full.names = TRUE)
  s <- do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)$summary
    r[intersect(c("source", "sim_id", "window", "true_beta", "beta_multiple", "objective", "r", "cells_mean", "dt"),
                names(r))]
  }))
  s <- s[s$beta_multiple %in% grid, ]
  data.frame(source = paste0("step0_", s$source), base = "ends", id = paste(s$source, s$sim_id, s$window),
             landscape = sub(".*_rep_", "rep", s$sim_id), interval = s$dt, cells = s$cells_mean,
             true_multiple = s$true_beta / 5e-5, multiple = s$beta_multiple, objective = s$objective, r = s$r,
             stringsAsFactors = FALSE)
}

null <- do.call(rbind, lapply(c("e1", "e2", "e3"), function(set) rbind(series_candidates(set, "ends"),
                                                                        series_candidates(set, "anchored"))))
e4 <- series_candidates("e4", "ends")
step0 <- step0_candidates()
half <- function(landscape, source) {
  if (grepl("^step0", source[1])) ifelse(as.integer(sub("rep", "", landscape)) <= 12L, "tune", "test") else
    ifelse(as.integer(sub("mt_L", "", landscape)) <= 6L, "tune", "test")
}

# ---- rules -----------------------------------------------------------------------
rules <- c(list(R0 = list(type = "R0")),
           unlist(lapply(c(20, 50, 100, 200), function(i) lapply(c(0, 300, 500, 1000), function(cc)
             list(type = "A", interval = i, cells = cc))), recursive = FALSE),
           lapply(c(2, 5, 10, 20, 40, 80), function(t) list(type = "B", threshold = t)),
           lapply(log(c(2, 4, 8, 16)), function(s) list(type = "C", sd = s)),
           with(expand.grid(t = c(5, 20), cc = c(300, 500), i = c(20, 50, 100)),
                Map(function(i, cc, t) list(type = "AB", interval = i, cells = cc, threshold = t), i, cc, t)))
rule_name <- function(r) switch(r$type, R0 = "R0", A = sprintf("A(I>=%g,C>=%g)", r$interval, r$cells),
                                B = sprintf("B(T=%g)", r$threshold), C = sprintf("C(sd=log %g)", exp(r$sd)),
                                AB = sprintf("AB(I>=%g,C>=%g,T=%g)", r$interval, r$cells, r$threshold))
names(rules) <- vapply(rules, rule_name, character(1))

pick_one <- function(d, r) {
  ok <- is.finite(d$objective)
  if (!any(ok)) return(1)
  base <- d$objective[d$multiple == 1]
  if (r$type == "C") {
    pen <- d$objective + (log(d$multiple) / r$sd)^2 / 2
    return(d$multiple[which.min(ifelse(ok, pen, Inf))])
  }
  best <- d$multiple[which.min(ifelse(ok, d$objective, Inf))]
  gate <- if (r$type %in% c("A", "AB")) d$interval[1] >= r$interval && d$cells[1] >= r$cells else TRUE
  thr <- if (r$type %in% c("B", "AB")) is.finite(base) && (base - min(d$objective[ok])) >= r$threshold else TRUE
  if (gate && thr) best else 1
}
apply_rule <- function(cand, r) {
  if (is.null(cand) || !nrow(cand)) return(NULL)
  key <- paste(cand$source, cand$base, cand$id)
  do.call(rbind, lapply(split(cand, key), function(d) {
    m <- pick_one(d, r)
    chosen <- d[d$multiple == m, ][1, ]
    fixed <- d[d$multiple == 1, ][1, ]
    out <- data.frame(source = d$source[1], base = d$base[1], id = d$id[1], landscape = d$landscape[1],
                      half = half(d$landscape[1], d$source[1]), true_multiple = d$true_multiple[1], pick = m,
                      stringsAsFactors = FALSE)
    for (s in c("d0", "d1", "d2")) {
      col <- paste0("rmse_", s)
      out[[paste0("ratio_", s)]] <- if (col %in% names(d)) chosen[[col]] / fixed[[col]] else NA_real_
    }
    out$dr <- if ("r" %in% names(d)) chosen$r - fixed$r else NA_real_
    out
  }))
}

# ---- criteria ----------------------------------------------------------------------
judge <- function(null_p, alt_e4, alt_s0, which_half) {
  n <- null_p[null_p$half == which_half, ]
  res <- list()
  for (b in c("ends", "anchored")) {
    e1 <- n[n$source == "e1" & n$base == b, ]
    if (!nrow(e1)) next
    idx <- utils::read.csv(file.path(root, "e1", "dataset_index.csv"), stringsAsFactors = FALSE)
    e1 <- merge(e1, idx[c("dataset_id", "bottleneck", "spacing", "depth", "samples")], by.x = "id", by.y = "dataset_id")
    cell <- stats::aggregate(cbind(ratio_d0, ratio_d1, ratio_d2) ~ bottleneck + spacing + depth + samples, e1,
                             stats::median)
    worst <- pmax(cell$ratio_d0, cell$ratio_d1, cell$ratio_d2)
    over <- c(stats::median(cell$ratio_d0), stats::median(cell$ratio_d1), stats::median(cell$ratio_d2))
    other <- vapply(c("e2", "e3"), function(s) {
      x <- n[n$source == s & n$base == b, ]
      max(stats::median(x$ratio_d0, na.rm = TRUE), stats::median(x$ratio_d1, na.rm = TRUE),
          stats::median(x$ratio_d2, na.rm = TRUE))
    }, numeric(1))
    res[[b]] <- c(e1_cells = sum(worst <= 1.03), e1_over = max(over), e2 = other[["e2"]], e3 = other[["e3"]])
  }
  null_pass <- all(vapply(res, function(v) v[["e1_cells"]] >= 18 && v[["e1_over"]] <= 1.02 && v[["e2"]] <= 1.03 &&
                            v[["e3"]] <= 1.03, logical(1)))
  a4 <- alt_e4[alt_e4$half == which_half, ]
  a0 <- alt_s0[alt_s0$half == which_half, ]
  alt <- rbind(a4[c("true_multiple", "pick")], a0[c("true_multiple", "pick")])
  within4 <- mean(abs(log2(alt$pick / alt$true_multiple)) <= 2 + 1e-9)
  # Datasets whose fits failed at every candidate (no cell in the model domain)
  # have no ratio; a rule cannot change them, so they are left out.
  e4_high <- if (nrow(a4)) max(stats::median(a4$ratio_d0, na.rm = TRUE), stats::median(a4$ratio_d1, na.rm = TRUE),
                               stats::median(a4$ratio_d2, na.rm = TRUE)) else NA
  s0_high <- stats::median(a0$dr[a0$true_multiple >= 4], na.rm = TRUE)
  alt_pass <- within4 >= 0.70 && (!nrow(a4) || isTRUE(e4_high < 1)) && isTRUE(s0_high > 0)
  list(null = res, null_pass = null_pass, within4 = within4, e4_ratio = e4_high, s0_gain = s0_high,
       alt_pass = alt_pass, pass = null_pass && alt_pass)
}

results <- lapply(names(rules), function(nm) {
  r <- rules[[nm]]
  np <- apply_rule(null, r)
  ap <- apply_rule(e4, r)
  sp <- apply_rule(step0, r)
  list(name = nm, rule = r, null = np, e4 = ap, step0 = sp, tune = judge(np, ap, sp, "tune"),
       test = judge(np, ap, sp, "test"))
})
names(results) <- names(rules)
saveRDS(results, file.path(out, "rules.rds"))

row_of <- function(x, h) {
  j <- x[[h]]
  data.frame(rule = x$name, null_pass = j$null_pass,
             e1_cells_ends = j$null$ends[["e1_cells"]], e1_cells_anch = j$null$anchored[["e1_cells"]],
             e2_max = max(j$null$ends[["e2"]], j$null$anchored[["e2"]]),
             within4 = j$within4, e4_ratio = j$e4_ratio, s0_gain = j$s0_gain, alt_pass = j$alt_pass, pass = j$pass,
             stringsAsFactors = FALSE)
}
tune <- do.call(rbind, lapply(results, row_of, h = "tune"))
passing <- tune[tune$pass, ]
chosen <- if (nrow(passing)) {
  passing <- passing[order(passing$e4_ratio, -passing$s0_gain), ]
  passing$rule[[1L]]
} else NA_character_
test <- if (!is.na(chosen)) row_of(results[[chosen]], "test") else NULL
md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
# Why the alternatives fail (descriptive, tuning half): how often a rule leaves
# the default where the truth is >= 4x, and what that does to the error
# (ratio at distance 1; better < 0.97, worse > 1.03).
why <- do.call(rbind, lapply(results, function(x) {
  a4 <- x$e4[x$e4$half == "tune" & is.finite(x$e4$ratio_d1), ]
  a0 <- x$step0[x$step0$half == "tune" & x$step0$true_multiple >= 4, ]
  n0 <- x$null[x$null$half == "tune" & x$null$source == "e1" & x$null$base == "ends", ]
  data.frame(rule = x$name, e1_moved = mean(n0$pick != 1), e4_moved = mean(a4$pick != 1),
             e4_better = mean(a4$ratio_d1 < 0.97), e4_worse = mean(a4$ratio_d1 > 1.03),
             e4_true_pick = mean(a4$pick == a4$true_multiple), step0_high_moved = mean(a0$pick != 1),
             stringsAsFactors = FALSE)
}))

# Picks on the Salehi trajectories of the Fig. 4d check (fits at the six rates;
# interval and cells as there): in full under R0 and the chosen rule, and as
# shares under every rule (criterion 3).
salehi_dir <- file.path("validation/sparse_joint/fig4d_check/results", c("beta", "beta_ext"))
salehi_files <- list.files(salehi_dir, pattern = "[.]rds$", full.names = TRUE)
salehi <- if (length(salehi_files)) {
  x <- do.call(rbind, lapply(salehi_files, function(f) readRDS(f)[c("fi", "context", "beta_multiple", "objective",
                                                                     "cells_mean", "dt")]))
  x <- unique(x[x$beta_multiple %in% grid, ])
  data.frame(source = "salehi", base = "ends", id = x$fi, landscape = x$context, interval = x$dt, cells = x$cells_mean,
             true_multiple = NA_real_, multiple = x$beta_multiple, objective = x$objective, stringsAsFactors = FALSE)
} else NULL
salehi_picks <- if (!is.null(salehi)) do.call(rbind, lapply(unique(c("R0", stats::na.omit(chosen))), function(nm) {
  p <- apply_rule(salehi, rules[[nm]])
  tab <- as.data.frame.matrix(table(p$landscape, factor(p$pick, levels = grid)))
  names(tab) <- c("1/16", "1/4", "1", "4", "16", "64")
  cbind(rule = nm, context = rownames(tab), tab, stringsAsFactors = FALSE)
})) else NULL
# Every rule, compactly: share of trajectories moved off the default, up (>= 4x)
# and down (<= 1/4x), by context.
salehi_all <- if (!is.null(salehi)) do.call(rbind, lapply(names(rules), function(nm) {
  p <- apply_rule(salehi, rules[[nm]])
  ctx <- sort(unique(p$landscape))
  row <- data.frame(rule = nm, stringsAsFactors = FALSE)
  for (cx in ctx) {
    q <- p[p$landscape == cx, ]
    row[[paste0(cx, "_up")]] <- mean(q$pick >= 4)
    row[[paste0(cx, "_down")]] <- mean(q$pick <= 1 / 4)
  }
  row
})) else NULL
salehi_range <- if (!is.null(salehi)) {
  r <- unique(salehi[c("id", "landscape", "interval", "cells")])
  stats::aggregate(cbind(interval, cells) ~ landscape, r, function(v) stats::median(v))
} else NULL

writeLines(c("# Short-interval plan, step 3: rules for leaving the default rate", "",
             sprintf("Datasets: null %d, e4 %d, step 0 %d (each with its candidates).",
                     length(unique(paste(null$source, null$base, null$id))),
                     length(unique(paste(e4$source, e4$base, e4$id))), length(unique(step0$id))), "",
             "## Tuning half", "", md(tune), "",
             sprintf("Chosen on the tuning half: %s", chosen), "",
             "## Test half (the chosen rule)", "", md(test), "",
             "## Why the alternatives fail (tuning half; e1 = true rate x1, ends)", "",
             "Share of datasets moved off the default; e4 better / worse = RMSE ratio at distance 1 below 0.97 / above 1.03.",
             "", md(why), "",
             "## Salehi trajectories (Fig. 4d check): rates picked, count", "", md(salehi_picks), "",
             "Median interval (days) and cells per time point by context:", "", md(salehi_range), "",
             "## Salehi trajectories under every rule: share moved up (>= 4x) and down (<= 1/4x)", "",
             md(salehi_all)),
           file.path(out, "rules.md"))
cat(readLines(file.path(out, "rules.md")), sep = "\n")
