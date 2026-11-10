#!/usr/bin/env Rscript
# Conditional sampling-depth sensitivity on frozen biological trajectories.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: prepare_depth_sensitivity.R PILOT_DATA OUTPUT")
input <- args[[1L]]; output <- args[[2L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
manifest <- utils::read.csv(file.path(input, "dataset_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(all(manifest$split == "train"))
rows <- list()
for (i in seq_len(nrow(manifest))) {
  original <- readRDS(manifest$dataset_file[i])
  k <- original$target_karyotypes
  distance_matrix <- as.matrix(stats::dist(k, method = "manhattan"))
  for (depth in c(150L, 1000L, 5000L)) {
    dataset <- original
    set.seed(original$sampling_seed)
    y0 <- as.numeric(stats::rmultinom(1, depth, original$truth_pi0))
    y1 <- as.numeric(stats::rmultinom(1, depth, original$truth_pi1))
    counts <- cbind(t0 = y0, t1 = y1)
    rownames(counts) <- original$target_labels
    observed <- which(rowSums(counts) > 0)
    counts <- counts[observed, , drop = FALSE]
    # Preserve the exact frozen pilot draw. Re-normalizing stored probabilities
    # can change RNG consumption at probability-one binomial boundaries.
    if (depth == 150L) {
      counts <- original$counts
      observed <- match(rownames(counts), original$target_labels)
    }
    dataset$dataset_id <- paste0("depth", depth, "_", original$dataset_id)
    dataset$counts <- counts
    dataset$actual_target_distance <- apply(distance_matrix[, observed, drop = FALSE], 1L, min)
    # Keep the estimand and strata fixed across depths. Do not recenter each run.
    dataset$target_distance <- original$target_distance
    dataset$gauge_weight <- original$gauge_weight
    dataset$config$n0 <- dataset$config$n1 <- depth
    dataset$config$depth_sensitivity <- list(reference_depth = 150L, sampling_depth = depth,
      fixed = c("latent populations", "all seeds", "target nodes", "reference-distance strata", "reference contrast"),
      interpretation = "independent multinomial observations conditional on a fixed biological trajectory")
    file <- file.path(output, paste0(dataset$dataset_id, ".rds"))
    saveRDS(dataset, file)
    row <- manifest[i, ]
    row$dataset_id <- dataset$dataset_id
    row$observed_nodes <- nrow(counts)
    row$max_distance <- max(dataset$target_distance)
    row$actual_max_distance <- max(dataset$actual_target_distance)
    row$sampling_depth <- depth
    row$dataset_file <- normalizePath(file)
    rows[[length(rows) + 1L]] <- row
  }
}
utils::write.csv(do.call(rbind, rows), file.path(output, "dataset_manifest.csv"), row.names = FALSE)
file.copy(file.path(input, "landscape_seed_manifest.csv"), file.path(output, "landscape_seed_manifest.csv"), overwrite = TRUE)
cat(sprintf("Prepared %d development depth datasets; calibration/test untouched.\n", length(rows)))
