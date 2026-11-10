#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "simulation_backends.R")); source(file.path(script_dir, "truth_and_data.R"))
config <- accuracy_round2_config()
domain <- round2_arg("domain", "cn3_1to4"); split <- round2_arg("split", "development")
profile <- round2_arg("profile", "standard")
output <- round2_arg("output", file.path(script_dir, "data", paste(domain, split, profile, sep = "_")))
if (!domain %in% names(config$domains)) stop("Unknown domain")
if (!profile %in% c("standard", "support_limit")) stop("Unknown profile")
replicates <- as.integer(strsplit(round2_arg("landscapes", "1,2"), ",", fixed = TRUE)[[1L]])
manifest <- accuracy_round2_seed_manifest(config)
split_name <- if (split == "test") "test_reserved" else split
if (!split_name %in% unique(manifest$split)) stop("Unknown split")
manifest <- manifest[manifest$domain == domain & manifest$split == split_name & manifest$landscape_replicate %in% replicates, ]
if (!nrow(manifest)) stop("No predeclared landscapes selected")
if (split != "development" && identical(config$status, "development_only_new_splits_sealed"))
  stop("SEALED: round-2 development stage is not cleared; provisional freeze fields do not authorize new calibration/test generation")
if (split != "development") {
  freeze_file <- round2_arg("generation-freeze-rds")
  if (is.null(freeze_file)) stop("SEALED: candidate/protocol freeze required before calibration; model/calibration/analysis/power freeze required before test")
  freeze <- readRDS(freeze_file)
  required <- c("candidate_source_sha256", "candidate_arguments_sha256", "protocol_sha256")
  if (!all(required %in% names(freeze))) stop("Incomplete candidate freeze")
  if (split == "test" && (!identical(freeze$power_status, "FINAL_N_FROZEN_BEFORE_NEW_TEST") ||
      is.null(freeze$calibration_objects_sha256) || is.null(freeze$analysis_sha256))) stop("SEALED: test power/calibration/analysis freeze incomplete")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(config = config, selected_landscapes = manifest,
              generator_sha256 = round2_sha256(file.path(script_dir, c("truth_and_data.R", "simulation_backends.R", "prepare_data.R"))),
              phase = split, model_family_input_forbidden = TRUE), file.path(output, "generation_protocol.rds"))
rows <- list()
for (li in seq_len(nrow(manifest))) {
  landscape <- manifest[li, ]
  landscape$split <- split
  lattice <- round2_lattice(config$domains[[domain]])
  truth <- round2_truth(domain, config$domains[[domain]], landscape$family, landscape$landscape_seed, lattice)
  for (a in config$abm_replicates) {
    first <- round2_prepare_dataset(config, domain, landscape, a, 1L, profile, truth = truth)
    latent <- list(pop0 = first$latent_population_t0, pop1 = first$latent_population_t1)
    for (s in config$sampling_replicates) {
      dataset <- if (s == 1L) first else round2_prepare_dataset(config, domain, landscape, a, s, profile, latent, truth)
      file <- file.path(output, paste0(dataset$dataset_id, ".rds"))
      if (file.exists(file) && !identical(readRDS(file), dataset)) stop("Immutable dataset differs; refusing overwrite")
      if (!file.exists(file)) saveRDS(dataset, file)
      rows[[length(rows) + 1L]] <- data.frame(dataset_id = dataset$dataset_id, domain = domain, profile = profile,
        family = landscape$family, landscape_id = landscape$landscape_id, split = split,
        landscape_seed = landscape$landscape_seed, abm_seed = dataset$abm_seed, sampling_seed = dataset$sampling_seed,
        target_nodes = length(dataset$target_labels), observed_nodes = nrow(dataset$counts),
        actual_max_distance = max(dataset$actual_target_distance), auxiliary_max_distance = max(dataset$auxiliary_stencil_distance),
        actual_has_d4 = any(dataset$actual_target_distance == 4L), actual_has_d5 = any(dataset$actual_target_distance == 5L),
        dataset_file = normalizePath(file), dataset_sha256 = unname(round2_sha256(file)))
    }
  }
}
round2_write_csv(do.call(rbind, rows), file.path(output, "dataset_manifest.csv"))
cat(sprintf("Prepared %d %s datasets; actual distances retained, no seed replacement.\n", length(rows), split))
