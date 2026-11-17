#!/usr/bin/env Rscript
# Tables, gates and one overview figure for a build's fits on e1-e3
# (README.md). Reads results/<build>/<set>/metrics_*.csv (run.R
# --stage=score); writes results/<build>/report.md, gates.csv and
# overview.pdf.
#
#   Rscript validation/sparse_joint/series_methods/report.R --build=BUILD [--methods=ends,anchored,...]
#
# The gates were fixed before any fit of the three methods was scored (plan
# of 2026-11-14, handoff):
#   e1  not worse than the first and last time point alone: in at least 18 of
#       the 24 cells (bottleneck x spacing x depth x number of time points) the
#       median RMSE ratio to `ends` is at most 1.03 at every distance, and the median
#       over cells is at most 1.02 at every distance;
#   e3  false positives controlled: median share of karyotypes called
#       different at most 0.05, and at most 0.10 at every bottleneck;
#   all no more failed fits than `ends` in any set.
# Among the methods that pass, the default is the one with the lower RMSE
# ratio on e2 (median over distances 0-2), then the higher power at effect
# 0.05 there; within 2% of each other, the simpler (anchored, then pairwise,
# then state_space). The paper benchmark (E4) is weighed before the choice.

sm_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(sm_script_dir, "setup.R"))

rp_read <- function(build, set, part) {
  path <- file.path(sm_default_output, build, set, sprintf("metrics_%s.csv", part))
  if (!file.exists(path)) return(NULL)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  x$set <- set
  x
}

# RMSE ratio of each method to `ends`, paired by dataset, then medians.
rp_ratios <- function(distance, methods) {
  d <- distance[distance$stratum %in% c("d0", "d1", "d2") & distance$status == "ok", ]
  ends <- d[d$method == "ends", c("set", "dataset_id", "stratum", "rmse")]
  names(ends)[names(ends) == "rmse"] <- "rmse_ends"
  x <- merge(d[d$method %in% methods, ], ends, by = c("set", "dataset_id", "stratum"))
  x$ratio <- x$rmse / x$rmse_ends
  x
}

rp_cell <- function(x) {
  cell <- c("set", "method", "bottleneck", "spacing", "depth", "samples")
  if ("start" %in% names(x)) cell <- c(cell, "start")
  cell <- intersect(cell, names(x))
  cell <- cell[vapply(cell, function(column) !all(is.na(x[[column]])), logical(1))]
  stats::aggregate(ratio ~ ., x[c(cell, "stratum", "ratio")], stats::median)
}

rp_gates <- function(ratios, power, fits, methods) {
  rows <- lapply(setdiff(methods, "ends"), function(m) {
    e1 <- ratios[ratios$set == "e1" & ratios$method == m, ]
    cells <- if (nrow(e1)) rp_cell(e1) else NULL
    e1_pass <- if (is.null(cells)) NA else {
      per_cell <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples, cells, max)
      overall <- stats::aggregate(ratio ~ stratum, cells, stats::median)
      sum(per_cell$ratio <= 1.03) >= 18 && all(overall$ratio <= 1.02)
    }
    e3 <- power[power$set == "e3" & power$method == m & power$stratum == "all", ]
    e3_pass <- if (!nrow(e3)) NA else {
      by_bottleneck <- stats::aggregate(false_positive ~ bottleneck, e3, stats::median)
      stats::median(e3$false_positive, na.rm = TRUE) <= 0.05 && all(by_bottleneck$false_positive <= 0.10)
    }
    failed <- vapply(sm_sets, function(s) {
      f <- fits[fits$set == s, ]
      if (!nrow(f)) return(NA)
      sum(f$method == m & f$status != "ok") <= sum(f$method == "ends" & f$status != "ok")
    }, logical(1))
    e2 <- ratios[ratios$set == "e2" & ratios$method == m, ]
    e2_power <- power[power$set == "e2" & power$method == m & power$stratum == "all", ]
    data.frame(method = m,
               e1_cells_pass = if (is.null(cells)) NA else {
                 per_cell <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples, cells, max)
                 sprintf("%d / %d", sum(per_cell$ratio <= 1.03), nrow(per_cell))
               },
               e1_pass = e1_pass, e3_false_positive = if (nrow(e3)) stats::median(e3$false_positive, na.rm = TRUE) else NA,
               e3_pass = e3_pass, failures_pass = all(failed, na.rm = TRUE),
               e2_ratio = if (nrow(e2)) stats::median(e2$ratio) else NA,
               e2_power_0.05 = if (nrow(e2_power)) stats::median(e2_power$power_0.05, na.rm = TRUE) else NA,
               # Means beside the pre-set medians, which are zero when most
               # datasets have no call at all.
               e3_false_positive_mean = if (nrow(e3)) mean(e3$false_positive, na.rm = TRUE) else NA,
               e2_power_0.05_mean = if (nrow(e2_power)) mean(e2_power$power_0.05, na.rm = TRUE) else NA,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$passes <- out$e1_pass %in% TRUE & out$e3_pass %in% TRUE & out$failures_pass
  out
}

# Power at a 5% false-positive rate: for each method and distance, the
# smallest threshold of |estimate| / model sd that at most 5% of e3's
# karyotypes exceed (counts pooled over datasets); then, on e1 and e2, the
# share of karyotypes with |truth| >= e beyond it on the right side, and the
# direction accuracy (right side at threshold 0).
rp_controlled_power <- function(counts, fpr = 0.05) {
  if (is.null(counts) || !any(counts$set == "e3")) return(NULL)
  null <- stats::aggregate(cbind(exceed, n) ~ method + stratum + z,
                           counts[counts$set == "e3" & counts$kind == "null", ], sum)
  null$rate <- null$exceed / null$n
  threshold <- do.call(rbind, lapply(split(null, list(null$method, null$stratum), drop = TRUE), function(d) {
    d <- d[order(d$z), ]
    data.frame(method = d$method[[1L]], stratum = d$stratum[[1L]], z_star = d$z[which(d$rate <= fpr)[1L]],
               stringsAsFactors = FALSE)
  }))
  effect <- stats::aggregate(cbind(exceed, n) ~ set + method + stratum + kind + z,
                             counts[counts$set %in% c("e1", "e2") & counts$kind != "null", ], sum)
  effect <- merge(effect, threshold, by = c("method", "stratum"))
  at <- effect[abs(effect$z - effect$z_star) < 1e-9, ]
  at$power <- at$exceed / at$n
  direction <- effect[effect$z == 0, ]
  direction$direction <- direction$exceed / direction$n
  wide <- function(x, value) {
    w <- stats::reshape(x[c("set", "method", "stratum", "kind", value)], idvar = c("set", "method", "stratum"),
                        timevar = "kind", direction = "wide")
    names(w) <- sub(paste0("^", value, "\\.effect_"), paste0(value, "_"), names(w))
    w
  }
  out <- merge(merge(threshold, wide(at, "power"), by = c("method", "stratum")),
               wide(direction, "direction")[c("set", "method", "stratum", "direction_0.05")],
               by = c("set", "method", "stratum"))
  out[order(out$set, out$stratum, out$method), ]
}

rp_markdown <- function(x, digits = 3) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "g") else as.character(v))
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

rp_figure <- function(path, ratios, power, fits, methods) {
  grDevices::pdf(path, width = 11, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(5, 4.5, 3, 1))
  colours <- c(ends = "grey40", anchored = "#1b6ca8", pairwise = "#d1495b", pairwise_same = "#edae49",
               state_space = "#00798c")
  shown <- intersect(methods, names(colours))
  for (set in sm_sets) {
    x <- ratios[ratios$set == set & ratios$method %in% setdiff(shown, "ends"), ]
    if (!nrow(x)) { graphics::plot.new(); next }
    m <- stats::aggregate(ratio ~ method + stratum, x, stats::median)
    tab <- stats::xtabs(ratio ~ method + stratum, m)
    tab <- tab[intersect(setdiff(shown, "ends"), rownames(tab)), , drop = FALSE]
    graphics::barplot(tab, beside = TRUE, col = colours[rownames(tab)], ylim = c(0, max(1.5, max(tab) * 1.1)),
                      main = sprintf("%s: RMSE / ends", set), ylab = "median ratio")
    graphics::abline(h = 1, lty = 2)
  }
  legend_methods <- setdiff(shown, "ends")
  graphics::legend("topright", legend = legend_methods, fill = colours[legend_methods], bty = "n", cex = 0.9)
  # Power at effect 0.05 on e1 and e2; false positives on e3.
  p <- power[power$stratum == "all" & power$method %in% shown, ]
  for (set in c("e1", "e2")) {
    q <- p[p$set == set, ]
    if (!nrow(q)) { graphics::plot.new(); next }
    m <- tapply(q$power_0.05, q$method, stats::median, na.rm = TRUE)[shown]
    graphics::barplot(m, col = colours[shown], ylim = c(0, 1), main = sprintf("%s: power at 0.05", set), las = 2)
  }
  q <- p[p$set == "e3", ]
  if (nrow(q)) {
    m <- tapply(q$false_positive, q$method, stats::median, na.rm = TRUE)[shown]
    graphics::barplot(m, col = colours[shown], ylim = c(0, max(0.2, max(m, na.rm = TRUE) * 1.2)),
                      main = "e3: false positives", las = 2)
    graphics::abline(h = 0.05, lty = 2)
  } else graphics::plot.new()
}

rp_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    args[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  build <- if (is.null(args$build)) stop("--build is required.", call. = FALSE) else args$build
  methods <- if (is.null(args$methods)) c("ends", "anchored", "pairwise", "pairwise_same", "state_space") else
    strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  # The sets' indices differ in their columns: bind on the union.
  read_all <- function(part) {
    parts <- Filter(Negate(is.null), lapply(sm_sets, function(s) rp_read(build, s, part)))
    columns <- unique(unlist(lapply(parts, names)))
    do.call(rbind, lapply(parts, function(x) {
      for (column in setdiff(columns, names(x))) x[[column]] <- NA
      x[columns]
    }))
  }
  distance <- read_all("distance")
  power <- read_all("power")
  fits <- read_all("fit")
  counts <- do.call(rbind, lapply(sm_sets, function(s) {
    path <- file.path(sm_default_output, build, s, "power_counts.rds")
    if (!file.exists(path)) return(NULL)
    x <- readRDS(path)
    x$set <- s
    x[c("set", "dataset_id", "method", "stratum", "kind", "z", "n", "exceed")]
  }))
  controlled <- rp_controlled_power(counts[counts$method %in% methods, ])
  if (!is.null(controlled)) {
    utils::write.csv(controlled, file.path(sm_default_output, build, "power_controlled.csv"), row.names = FALSE)
  }
  ratios <- rp_ratios(distance, methods)
  gates <- rp_gates(ratios, power, fits, methods)
  out_dir <- file.path(sm_default_output, build)
  utils::write.csv(gates, file.path(out_dir, "gates.csv"), row.names = FALSE)
  rp_figure(file.path(out_dir, "overview.pdf"), ratios, power, fits, methods)
  accuracy <- stats::aggregate(ratio ~ set + method + stratum, ratios, stats::median)
  accuracy <- stats::reshape(accuracy, idvar = c("set", "method"), timevar = "stratum", direction = "wide")
  coverage <- stats::aggregate(coverage95 ~ set + method + stratum,
                               distance[distance$stratum %in% c("d0", "d1", "d2") & distance$method %in% methods, ],
                               stats::median)
  coverage <- stats::reshape(coverage, idvar = c("set", "method"), timevar = "stratum", direction = "wide")
  pw <- power[power$stratum == "all" & power$method %in% methods, ]
  # Mean over datasets: most datasets call few karyotypes, so medians are 0.
  power_table <- stats::aggregate(cbind(power_0.02, power_0.05, power_0.1, false_positive) ~ set + method, pw,
                                  function(v) mean(v, na.rm = TRUE), na.action = stats::na.pass)
  power_distance <- stats::aggregate(cbind(power_0.05, false_positive) ~ set + method + stratum,
                                     power[power$method %in% methods & power$stratum != "all", ],
                                     function(v) mean(v, na.rm = TRUE), na.action = stats::na.pass)
  f <- fits[fits$method %in% methods, ]
  fit_table <- stats::aggregate(cbind(failed = f$status != "ok", not_accepted = f$convergence %in% 1,
                                      seconds = f$seconds, peak_mb = f$peak_mb) ~ set + method,
                                data.frame(f[c("set", "method")], failed = f$status != "ok",
                                           not_accepted = f$convergence %in% 1, seconds = f$seconds,
                                           peak_mb = f$peak_mb),
                                function(v) if (all(is.na(v))) NA else stats::median(v, na.rm = TRUE),
                                na.action = stats::na.pass)
  counts <- stats::aggregate(cbind(fits = 1, failed = f$status != "ok", not_accepted = f$convergence %in% 1) ~
                               set + method, data.frame(f[c("set", "method")], fits = 1, failed = f$status != "ok",
                                                        not_accepted = f$convergence %in% 1), sum)
  lines <- c(sprintf("# Series methods, build %s", build), "",
             "Median RMSE ratio to `ends` (paired by dataset), by distance:", "", rp_markdown(accuracy), "",
             "Median 95% coverage by distance:", "", rp_markdown(coverage), "",
             "Mean power (effect 0.02, 0.05, 0.1) and false positives over datasets:", "", rp_markdown(power_table), "",
             "Same by distance (power at 0.05, false positives):", "", rp_markdown(power_distance), "",
             paste("Power at a 5% false-positive rate (threshold z_star on |estimate| / model sd from e3;",
                   "direction: right sign at |truth| >= 0.05):"), "",
             if (is.null(controlled)) "(no counts)" else rp_markdown(controlled), "",
             "Fits (counts; median seconds and peak MB):", "", rp_markdown(counts), "",
             rp_markdown(fit_table[c("set", "method", "seconds", "peak_mb")]), "",
             "Gates (header of report.R):", "", rp_markdown(gates), "")
  writeLines(lines, file.path(out_dir, "report.md"))
  cat(lines, sep = "\n")
}

if (sys.nframe() == 0L) rp_main()
