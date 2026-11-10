#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "config.R"))
output <- round2_arg("output", file.path(script_dir, "results", "protocol_preparation"))
v1_source <- round2_arg("v1-source", "/private/tmp/alfak2_iteration2_20261104_0518/v1-source")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
config <- accuracy_round2_config(); manifest <- accuracy_round2_seed_manifest(config)
stopifnot(!anyDuplicated(manifest$landscape_seed), all(manifest$landscape_replicate[manifest$split != "development"] > 8L))
round2_write_csv(manifest, file.path(output, "reserved_landscape_seed_manifest.csv"))
saveRDS(config, file.path(output, "protocol_config.rds"))
for (model in c("original_baseline", "current_default_two_layer", "v1_joint", "A", "B", "zero_point_predictor"))
  saveRDS(accuracy_round2_model_arguments(model), file.path(output, paste0(model, "_arguments.rds")))
source_fingerprint <- round2_source_fingerprint(v1_source)
round2_write_csv(source_fingerprint$manifest, file.path(output, "v1_source_sha256.csv"))
round2_write_csv(data.frame(item = c("v1_source", "protocol_config", "reserved_seed_manifest"),
                           sha256 = c(source_fingerprint$sha256,
                             unname(round2_sha256(file.path(output, "protocol_config.rds"))),
                             unname(round2_sha256(file.path(output, "reserved_landscape_seed_manifest.csv"))))),
                 file.path(output, "freeze_preparation_sha256.csv"))
round2_write_csv(data.frame(domain = names(config$domains), calibration_landscapes = 18L,
  provisional_test_landscapes = 36L, maximum_test_landscapes = 60L, datasets_per_landscape = 4L,
  calibration_fits_per_model = 72L, provisional_test_fits_per_model = 144L,
  max_test_fits_per_model = 240L, maximum_concurrent_large_fits = 2L,
  power_status = c("deferred_candidate_v2_pair_variance", "unvalidated_new_domain_farfield_variance")),
  file.path(output, "resource_budget.csv"))
cat("Reserved manifests and source SHA256 written. No calibration/test landscapes generated.\n")
