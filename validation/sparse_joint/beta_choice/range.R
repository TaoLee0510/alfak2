#!/usr/bin/env Rscript
# Short-interval plan, fallback: when the missegregation rate cannot be chosen
# from the data, fit over a range of rates and report how much the answer
# moves (as fit_sparse_joint_rates() does for the birth rate).
#
# Offline check on the two-time-point candidate fits of series_methods (build
# cv10): e1 (true rate x1) and e4 (true rate x4 or x16). For a range of three
# rates the prediction is
#   point     the fit at the lowest rate of the range (or at its centre);
#   interval  the union of the 95% intervals of the three fits, written as
#             the standard deviation whose symmetric interval about the point
#             contains it (as predict.alfak2_sparse_joint_rate_fit()).
# Each is scored against the truth exactly as a fit (stt_score_predictions),
# and compared with the fit at the default rate alone.
#
#   Rscript validation/sparse_joint/beta_choice/range.R [--workers=3] [--report-only]
#
# --report-only rebuilds the report from results/range.rds.
#
# Writes validation/sparse_joint/beta_choice/results/range.md and range.rds.

local({
  here <- "validation/sparse_joint/series_methods"
  sm_script_dir <<- normalizePath(here)
  source(file.path(here, "setup.R"))
})
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(grep("^--workers=", args))) as.integer(sub("^--workers=", "", grep("^--workers=", args, value = TRUE))) else 3L
out <- "validation/sparse_joint/beta_choice/results"
z <- stats::qnorm(0.975)
suffix <- c(`1/16` = "_b1o16", `1/4` = "_b1o4", `1` = "", `4` = "_b4", `16` = "_b16", `64` = "_b64")
ranges <- list(
  `default only` = list(rates = "1", point = "1"),
  `1/4-4, lowest` = list(rates = c("1/4", "1", "4"), point = "1/4"),
  `1/4-4, centre` = list(rates = c("1/4", "1", "4"), point = "1"),
  `1-16, lowest` = list(rates = c("1", "4", "16"), point = "1"),
  `1-16, centre` = list(rates = c("1", "4", "16"), point = "4"),
  `1/16-16, lowest` = list(rates = c("1/16", "1/4", "1", "4", "16"), point = "1/16"),
  `1/16-16, centre` = list(rates = c("1/16", "1/4", "1", "4", "16"), point = "1")
)

score_set <- function(set) {
  config <- sm_config(set, "cv10")
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  rows <- si_parallel(index$dataset_id, function(dataset_id) {
    fits <- lapply(suffix, function(s) {
      p <- si_fit_path(config, dataset_id, paste0("ends", s))
      if (!file.exists(p)) return(NULL)
      f <- readRDS(p)
      if (identical(f$status, "ok")) f$prediction else NULL
    })
    if (any(vapply(fits, is.null, logical(1)))) return(NULL)
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- mt_dataset_nodes(config, dataset)
    truth <- sm_true_growth_rate(nodes$karyotype, dataset$truth)
    do.call(rbind, lapply(names(ranges), function(nm) {
      r <- ranges[[nm]]
      at <- fits[r$rates]
      mean_at <- vapply(at, function(p) p$fitness, numeric(nrow(at[[1L]])))
      sd_at <- vapply(at, function(p) p$fitness_sd, numeric(nrow(at[[1L]])))
      mean_at <- matrix(mean_at, ncol = length(at))
      sd_at <- matrix(sd_at, ncol = length(at))
      point <- fits[[r$point]]$fitness
      low <- apply(mean_at - z * sd_at, 1L, min)
      high <- apply(mean_at + z * sd_at, 1L, max)
      prediction <- fits[[r$point]]
      prediction$fitness <- point
      prediction$fitness_sd <- pmax(high - point, point - low) / z
      score <- stt_score_predictions(nodes, truth, prediction, config$scoring$reference_min_count)
      m <- score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse", "coverage95",
                                                                           "interval_score95", "spearman")]
      cbind(data.frame(set = set, dataset_id = dataset_id, range = nm, stringsAsFactors = FALSE), m)
    }))
  }, workers)
  x <- do.call(rbind, rows)
  merge(x, index[c("dataset_id", "landscape_id", "bottleneck", "spacing", "depth", "samples")], by = "dataset_id")
}

if ("--report-only" %in% args) {
  x <- readRDS(file.path(out, "range.rds"))
} else {
  x <- rbind(score_set("e1"), score_set("e4"))
  x$truth <- ifelse(x$set == "e1", "x1", ifelse(as.integer(sub("mt_L", "", x$landscape_id)) %% 2L == 1L, "x4", "x16"))
  base <- x[x$range == "default only", c("dataset_id", "set", "stratum", "rmse", "interval_score95")]
  x <- merge(x, base, by = c("dataset_id", "set", "stratum"), suffixes = c("", "_default"))
  x$rmse_ratio <- x$rmse / x$rmse_default
  x$is_ratio <- x$interval_score95 / x$interval_score95_default
  saveRDS(x, file.path(out, "range.rds"))
}

md <- function(d) {
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
order_range <- function(d) {
  d[order(match(d$truth, c("x1", "x4", "x16")), match(d$range, names(ranges)),
          if ("spacing" %in% names(d)) d$spacing else rep(0, nrow(d))), ]
}
agg <- stats::aggregate(cbind(rmse_ratio, coverage95, is_ratio) ~ truth + range + stratum, x, stats::median)
wide <- reshape(agg, idvar = c("truth", "range"), timevar = "stratum", direction = "wide")
by_spacing <- stats::aggregate(cbind(rmse_ratio, coverage95) ~ truth + range + spacing, x[x$stratum == "d1", ],
                               stats::median)
writeLines(c("# Fallback: fitting over a range of missegregation rates (two time points, cv10)", "",
             "Point = the fit at the named end of the range; interval = union of the fits' 95% intervals.",
             "rmse_ratio and is_ratio (interval score) are relative to the fit at the default rate alone;",
             "medians over datasets. Truth x1 = e1, x4 and x16 = e4.", "",
             "## By truth and distance", "", md(order_range(wide)), "",
             "## Distance 1, by spacing (days between time points)", "", md(order_range(by_spacing))),
           file.path(out, "range.md"))
cat(readLines(file.path(out, "range.md")), sep = "\n")
