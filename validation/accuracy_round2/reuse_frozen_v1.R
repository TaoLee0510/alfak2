#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "score_models.R"))
data_dir <- round2_arg("data", "/private/tmp/alfak2_iteration2_20261104_0518/development81-data")
output <- round2_arg("output", "/private/tmp/alfak2_iteration2_20261104_0518/development81-v1-reuse")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
all_nodes <- metrics <- statuses <- audit <- list()
expected_scales <- c(baseline = .5, main = .07, curvature = .04, interaction = .04, residual = .125, edge = .25)
for (i in seq_len(nrow(manifest))) {
  dataset <- readRDS(manifest$dataset_file[i]); original <- readRDS(manifest$original_dataset_file[i])
  frozen <- readRDS(manifest$v1_reference_cache[i]); nodes <- frozen$nodes
  stopifnot(identical(dataset$counts, original$counts), identical(dataset$gauge_weight, original$gauge_weight),
    identical(dataset$truth_relative, original$truth_relative), dataset$split == "development",
    identical(unname(round2_sha256(manifest$v1_reference_cache[i])), manifest$v1_reference_cache_sha256[i]),
    identical(frozen$local_diagnostics$forward_model, "division_process"),
    identical(frozen$local_diagnostics$inference_method, "laplace"),
    max(abs(frozen$global_hyperparameters[names(expected_scales)] - expected_scales)) < 1e-12)
  index <- match(dataset$target_labels, nodes$karyotype)
  stopifnot(!anyNA(index), max(abs(nodes$truth_relative[index] - dataset$truth_relative)) < 1e-12)
  nodes <- nodes[index, ]; nodes$model <- "v1_joint"; nodes$split <- "development"
  nodes$distance <- as.numeric(dataset$target_distance)
  score <- round2_score_nodes(nodes, dataset)
  status <- data.frame(dataset_id = dataset$dataset_id, model = "v1_joint", family = dataset$family,
    landscape_id = dataset$landscape_id, fit_status = "reused_frozen_v1",
    covariance_status = frozen$status$covariance_status, elapsed_seconds = 0,
    warnings = "No new fit; identical original counts/reference; source cache SHA and fixed prior scales checked")
  all_nodes[[i]] <- nodes; metrics[[i]] <- score; statuses[[i]] <- status
  audit[[i]] <- data.frame(dataset_id = dataset$dataset_id, original_split = manifest$original_split[i],
    original_dataset_sha256 = manifest$original_dataset_sha256[i],
    reference_cache_sha256 = manifest$v1_reference_cache_sha256[i], reused_without_refit = TRUE,
    unchanged_counts = TRUE, unchanged_reference = TRUE, original_implementation_md5 = frozen$cache_key$implementation_hash)
  saveRDS(list(nodes = nodes, metrics = score, status = status, original_cache_key = frozen$cache_key,
    source_cache_sha256 = manifest$v1_reference_cache_sha256[i], original_local_diagnostics = frozen$local_diagnostics),
    file.path(output, paste0("v1_joint_", dataset$dataset_id, ".rds")))
}
round2_write_csv(do.call(rbind, metrics), file.path(output, "case_metrics.csv"))
round2_write_csv(do.call(rbind, statuses), file.path(output, "fit_status.csv"))
round2_write_csv(do.call(rbind, audit), file.path(output, "reuse_audit.csv"))
saveRDS(all_nodes, file.path(output, "node_tables.rds"))
cat("Reused 18 frozen V1 posterior summaries with exact count/reference checks; no new fit.\n")
