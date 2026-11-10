#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "interval_calibration.R"))
args <- commandArgs(TRUE)
if (length(args) < 4L) stop("Usage: freeze_interval_calibration.R CAL_DATA CAL_RESULTS OUTPUT SEALED_TEST_DATA")
data_dir <- args[[1L]]; results_dir <- args[[2L]]; output <- args[[3L]]; test_data <- args[[4L]]
if (file.exists(file.path(test_data, "dataset_manifest.csv"))) stop("Test data already generated: cannot certify pre-test calibration freeze")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output, "_calibration\\.rds$"))) stop("Calibration objects already frozen; refusing to overwrite")
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) > 0, all(manifest$split == "calibration"))
datasets <- lapply(manifest$dataset_file, readRDS)
scenario <- alfak2_calibration_scenario(datasets[[1L]])
stopifnot(all(vapply(datasets, function(x) identical(alfak2_calibration_scenario(x), scenario), logical(1))))
node_files <- list.files(results_dir, "_nodes\\.csv$", full.names = TRUE)
tables <- lapply(node_files, utils::read.csv, stringsAsFactors = FALSE)
nodes <- do.call(rbind, tables)
nodes <- nodes[!duplicated(nodes[, c("dataset_id", "model", "karyotype")]), , drop = FALSE]
# Check every model before writing any immutable object, so a still-running
# worker cannot leave a partially frozen calibration collection.
for (model in unique(nodes$model)) {
  x <- nodes[nodes$model == model, , drop = FALSE]
  if (!setequal(unique(x$dataset_id), manifest$dataset_id)) stop(paste("Incomplete calibration datasets for", model))
  expected <- stats::setNames(vapply(datasets, function(x) length(x$target_labels), integer(1)), manifest$dataset_id)
  actual <- table(x$dataset_id)
  if (any(as.integer(actual[names(expected)]) != expected)) stop(paste("Incomplete calibration target rows for", model))
}
objects <- list()
for (model in unique(nodes$model)) {
  x <- nodes[nodes$model == model, , drop = FALSE]
  if (!setequal(unique(x$dataset_id), manifest$dataset_id)) stop(paste("Incomplete calibration datasets for", model))
  files <- node_files[vapply(tables, function(z) identical(unique(z$model), model), logical(1))]
  provenance <- lapply(sub("_nodes\\.csv$", "_provenance.rds", files), readRDS)
  keys <- lapply(provenance, alfak2_calibration_model_key)
  if (!all(vapply(keys, identical, logical(1), keys[[1L]]))) stop("Different model implementations/arguments within calibration")
  calibration <- fit_alfak2_interval_calibration(x, keys[[1L]], scenario)
  calibration$freeze_provenance <- list(
    frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), test_manifest_present = FALSE,
    sealed_test_directory = test_data, calibration_manifest_md5 = tools::md5sum(file.path(data_dir, "dataset_manifest.csv")),
    calibration_dataset_md5 = tools::md5sum(manifest$dataset_file), calibration_nodes_md5 = tools::md5sum(files),
    calibration_source_md5 = tools::md5sum(file.path(script_dir, c("interval_calibration.R", "freeze_interval_calibration.R"))))
  file <- file.path(output, paste0(model, "_calibration.rds"))
  saveRDS(calibration, file)
  utils::write.csv(calibration$table, file.path(output, paste0(model, "_quantiles.csv")), row.names = FALSE)
  objects[[length(objects) + 1L]] <- data.frame(model = model, object_file = normalizePath(file), md5 = unname(tools::md5sum(file)),
    n_landscapes = length(calibration$calibration_landscapes), n_datasets = length(calibration$calibration_datasets),
    test_generated_at_freeze = FALSE, stringsAsFactors = FALSE)
}
object_manifest <- do.call(rbind, objects)
utils::write.csv(object_manifest, file.path(output, "frozen_calibration_manifest.csv"), row.names = FALSE)
saveRDS(list(scenario = scenario, object_manifest = object_manifest,
              calibration_landscape_ids = sort(unique(manifest$landscape_id)),
              frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), test_generated = FALSE),
         file.path(output, "freeze_record.rds"))
print(object_manifest, row.names = FALSE)
