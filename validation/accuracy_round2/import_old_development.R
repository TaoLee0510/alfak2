#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
root <- round2_arg("old-root", "/private/tmp/alfak2_accuracy_20261104_0415")
output <- round2_arg("output", "/private/tmp/alfak2_iteration2_20261104_0518/development81-data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (phase in c("cal", "test")) {
  data_dir <- file.path(root, paste0(phase, "-exact-data"))
  results_dir <- file.path(root, paste0(phase, "-exact-results"))
  manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
  manifest <- manifest[grepl("_A1_S1$", manifest$dataset_id), ]
  for (i in seq_len(nrow(manifest))) {
    original <- readRDS(manifest$dataset_file[i]); dataset <- original
    dataset$original_round1_split <- original$split
    dataset$split <- "development"
    dataset$round2_role <- "old_exact_revealed_development"
    dataset$domain <- "cn4_1to3"
    file <- file.path(output, basename(manifest$dataset_file[i]))
    if (file.exists(file) && !identical(readRDS(file), dataset)) stop("Immutable imported development dataset differs")
    if (!file.exists(file)) saveRDS(dataset, file)
    stopifnot(identical(original$counts, dataset$counts), identical(original$gauge_weight, dataset$gauge_weight),
              identical(original$truth_relative, dataset$truth_relative), identical(original$target_distance, dataset$target_distance))
    cache <- file.path(results_dir, paste0("division_v1_", dataset$dataset_id, ".rds"))
    frozen <- readRDS(cache)
    ii <- match(dataset$target_labels, frozen$nodes$karyotype)
    stopifnot(!anyNA(ii), max(abs(frozen$nodes$truth_relative[ii] - dataset$truth_relative)) < 1e-12,
              identical(as.numeric(frozen$nodes$distance[ii]), as.numeric(dataset$target_distance)))
    row <- manifest[i, ]
    row$original_split <- row$split; row$split <- "development"; row$mode <- "old_exact_revealed_development"
    row$domain <- "cn4_1to3"; row$original_dataset_file <- row$dataset_file
    row$original_dataset_sha256 <- unname(round2_sha256(row$original_dataset_file))
    row$dataset_file <- normalizePath(file); row$dataset_sha256 <- unname(round2_sha256(file))
    row$v1_reference_cache <- normalizePath(cache); row$v1_reference_cache_sha256 <- unname(round2_sha256(cache))
    rows[[length(rows) + 1L]] <- row
  }
}
manifest <- do.call(rbind, rows)
stopifnot(nrow(manifest) == 18L, length(unique(manifest$landscape_id)) == 18L)
round2_write_csv(manifest, file.path(output, "dataset_manifest.csv"))
cat("Imported 18 pre-fixed A1/S1 old exact landscapes as development; original files unchanged.\n")
