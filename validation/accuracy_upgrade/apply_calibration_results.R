#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "interval_calibration.R"))
source(file.path(script_dir, "metrics.R"))
args <- commandArgs(TRUE)
if (length(args) < 3L) stop("Usage: apply_calibration_results.R DATA RAW_RESULTS FROZEN_OBJECTS")
data_dir <- args[[1L]]; results_dir <- args[[2L]]; calibration_dir <- args[[3L]]
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
frozen <- utils::read.csv(file.path(calibration_dir, "frozen_calibration_manifest.csv"), stringsAsFactors = FALSE)
for (mi in seq_len(nrow(frozen))) {
  model <- frozen$model[mi]
  object_file <- file.path(calibration_dir, basename(frozen$object_file[mi]))
  if (!identical(unname(tools::md5sum(object_file)), frozen$md5[mi])) stop("Frozen calibration checksum changed")
  calibration <- readRDS(object_file)
  status_files <- list.files(results_dir, "_status\\.csv$", full.names = TRUE)
  status_tables <- lapply(status_files, utils::read.csv, stringsAsFactors = FALSE)
  files <- status_files[vapply(status_tables, function(x) identical(unique(x$model), model), logical(1))]
  if (!length(files)) next
  raw_nodes <- do.call(rbind, lapply(sub("_status\\.csv$", "_nodes.csv", files), utils::read.csv, stringsAsFactors = FALSE))
  raw_nodes <- raw_nodes[!duplicated(raw_nodes[, c("dataset_id", "model", "karyotype")]), , drop = FALSE]
  raw_status <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
  raw_status <- raw_status[!duplicated(raw_status$dataset_id), , drop = FALSE]
  raw_metrics <- do.call(rbind, lapply(sub("_status\\.csv$", "_metrics.csv", files), utils::read.csv, stringsAsFactors = FALSE))
  raw_metrics <- raw_metrics[!duplicated(raw_metrics[, c("dataset_id", "model", "distance_group", "metric")]), , drop = FALSE]
  provenance <- lapply(sub("_status\\.csv$", "_provenance.rds", files), readRDS)
  keys <- lapply(provenance, alfak2_calibration_model_key)
  if (!all(vapply(keys, identical, logical(1), calibration$model_key))) stop("Test/calibration model fingerprint mismatch")
  new_model <- paste0(model, "_calibrated")
  metrics <- nodes <- statuses <- list()
  for (di in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[di])
    x <- raw_nodes[raw_nodes$dataset_id == dataset$dataset_id, , drop = FALSE]
    if (!nrow(x)) stop(paste("Missing raw prediction rows for", dataset$dataset_id))
    applied <- apply_alfak2_interval_calibration(x, calibration, calibration$model_key, alfak2_calibration_scenario(dataset))
    stopifnot(identical(applied$pred_relative, x$pred_relative), identical(applied$sd_relative, x$sd_relative))
    applied$model <- new_model
    scoring <- applied
    for (level in calibration$levels) {
      suffix <- sprintf("%02d", round(100 * level))
      scoring[[paste0("relative_low_", suffix)]] <- applied[[paste0("calibrated_low_", suffix)]]
      scoring[[paste0("relative_high_", suffix)]] <- applied[[paste0("calibrated_high_", suffix)]]
    }
    # Missing calibrated intervals are failures of interval availability, not
    # successes. Only this scoring view masks SDs; exported predictions/SDs stay intact.
    scoring$sd_relative[scoring$calibration_status != "empirically_calibrated"] <- NA_real_
    status <- raw_status[match(dataset$dataset_id, raw_status$dataset_id), , drop = FALSE]
    status$model <- new_model
    status$warnings <- paste(status$warnings, "Empirical calibration only; point/SD unchanged; unvalidated interval rows:",
                              sum(applied$calibration_status != "empirically_calibrated"))
    scored <- score_dataset(scoring, dataset, status$runtime_seconds, status$fit_status, status$covariance_status)
    unchanged_metrics <- raw_metrics[raw_metrics$dataset_id == dataset$dataset_id &
                                      raw_metrics$distance_group %in% c("heldout_counts", "local_target"), , drop = FALSE]
    unchanged_metrics$model <- new_model
    if (nrow(unchanged_metrics)) scored <- rbind(scored, unchanged_metrics)
    metrics[[di]] <- scored; nodes[[di]] <- applied; statuses[[di]] <- status
  }
  utils::write.csv(do.call(rbind, metrics), file.path(results_dir, paste0(new_model, "_metrics.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, nodes), file.path(results_dir, paste0(new_model, "_nodes.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, statuses), file.path(results_dir, paste0(new_model, "_status.csv")), row.names = FALSE)
}
metric_files <- list.files(results_dir, "_status\\.csv$", full.names = TRUE)
all_metrics <- do.call(rbind, lapply(sub("_status\\.csv$", "_metrics.csv", metric_files), utils::read.csv, stringsAsFactors = FALSE))
all_metrics <- all_metrics[!duplicated(all_metrics[, c("dataset_id", "model", "distance_group", "metric")]), , drop = FALSE]
paired <- paired_landscape_summary(all_metrics, baseline_model = "baseline_calibrated")
utils::write.csv(paired, file.path(results_dir, "paired_against_calibrated_baseline.csv"), row.names = FALSE)
cat("Applied frozen empirical intervals; all reported fitness point estimates and SDs unchanged.\n")
