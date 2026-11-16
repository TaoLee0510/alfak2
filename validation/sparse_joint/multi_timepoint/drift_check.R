#!/usr/bin/env Rscript
# Stage 2 diagnostic: is the fit on all time points worse than first+last
# because of passaging drift? The pilot (landscape mt_L05) found the fit on all
# time points worse even at bottleneck 2e4, mostly on karyotypes seen only at
# intermediate time points. This refits that landscape's windows at bottleneck 2e5, near the
# ALFA-K paper's ABM (n0 = 2e5), with the same seeds, spacings, depths and
# numbers of time points, and compares all time points with first+last there.
#
#   Rscript validation/sparse_joint/multi_timepoint/drift_check.R --stage=datasets|fit|score [--workers=4]
#     [--methods=all,ends]   (one build per run: `majors` on its own)
#
# Results in results/drift_check/ (git-ignored); truth populations linked from
# results/series/truth.

mt_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(mt_script_dir, "setup.R"))

dc_config <- function() {
  config <- mt_config(file.path(mt_default_output, "drift_check"), "series")
  config$landscapes <- config$landscapes[config$landscapes$landscape_id == "mt_L05", ]
  config$bottlenecks <- 2e5
  config
}

dc_link_truth <- function(config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  link <- file.path(config$output_dir, "truth")
  if (!file.exists(link)) file.symlink("../series/truth", link)
}

# Errors of all and ends by group of scored karyotypes: observed at the first
# or last time point, observed only in between, one and two steps away. Pilot
# datasets of the same landscape at bottleneck 2e4 for comparison.
dc_score <- function(config) {
  series <- mt_config(NULL, "series")
  sets <- list(b200000 = config, b20000 = series)
  rows <- list()
  for (set in names(sets)) {
    cfg <- sets[[set]]
    index <- utils::read.csv(si_index_path(cfg), stringsAsFactors = FALSE)
    ids <- grep(if (set == "b20000") "^mt_L05_b20000_" else "^mt_L05_", index$dataset_id, value = TRUE)
    for (id in ids) {
      dataset <- readRDS(stt_dataset_path(cfg, id))
      nodes <- mt_dataset_nodes(cfg, dataset)
      truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
      counts <- dataset$counts
      ends_seen <- rownames(counts)[rowSums(counts[, c(1L, ncol(counts)), drop = FALSE]) > 0]
      group <- ifelse(nodes$distance == 0, ifelse(nodes$karyotype %in% ends_seen, "d0_ends", "d0_middle"),
                      paste0("d", nodes$distance))
      for (method in c("all", "ends", "majors", "majors_ends")) {
        path <- si_fit_path(cfg, id, method)
        if (!file.exists(path)) next
        fit <- readRDS(path)
        if (!identical(fit$status, "ok")) next
        s <- stt_score_predictions(nodes, truth, fit$prediction, cfg$scoring$reference_min_count)
        err <- s$nodes$pred_relative - s$nodes$truth_relative
        for (g in unique(group)) {
          rows[[length(rows) + 1L]] <- data.frame(set = set, id = id, method = method, group = g,
                                                  n = sum(group == g), rmse = sqrt(mean(err[group == g]^2)),
                                                  bias = mean(err[group == g]), stringsAsFactors = FALSE)
        }
      }
    }
  }
  x <- do.call(rbind, rows)
  w <- stats::reshape(x, idvar = c("set", "id", "group", "n"), timevar = "method", direction = "wide")
  w$ratio <- w$rmse.all / w$rmse.ends
  w$ratio_majors <- if ("rmse.majors" %in% names(w)) w$rmse.majors / w$rmse.ends else NA_real_
  w$ratio_majors_ends <- if ("rmse.majors_ends" %in% names(w)) w$rmse.majors_ends / w$rmse.ends else NA_real_
  utils::write.csv(w, file.path(config$output_dir, "groups.csv"), row.names = FALSE)
  med <- function(v) stats::median(v, na.rm = TRUE)
  out <- stats::aggregate(w[c("ratio", "ratio_majors", "ratio_majors_ends", "n")], w[c("set", "group")], med)
  print(out, digits = 3)
  utils::write.csv(out, file.path(config$output_dir, "summary.csv"), row.names = FALSE)
}

dc_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    args[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- dc_config()
  workers <- if (is.null(args$workers)) 4L else as.integer(args$workers)
  switch(if (is.null(args$stage)) "" else args$stage,
    datasets = { dc_link_truth(config); mt_stage_datasets(config, workers) },
    fit = mt_stage_fit(config, workers, if (is.null(args$methods)) c("all", "ends") else
      strsplit(args$methods, ",", fixed = TRUE)[[1L]], NULL),
    score = dc_score(config),
    stop("--stage must be datasets, fit or score.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) dc_main()
