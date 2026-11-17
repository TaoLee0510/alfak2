#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "config.R")); source(file.path(script_dir, "truth_and_data.R"))
config <- accuracy_round2_config(); manifest <- accuracy_round2_seed_manifest(config)
stopifnot(!anyDuplicated(manifest$landscape_seed),
  all(manifest$landscape_replicate[manifest$split %in% c("calibration", "test_reserved")] > 8L),
  all(vapply(config$domains, function(x) identical(x$primary_distance, "actual_observed_support"), logical(1))))
old <- readRDS(file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "node_matrices.rds"))
old_manifest <- utils::read.csv(file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "test_dataset_manifest.csv"), stringsAsFactors = FALSE)
domain <- config$domains$cn4_1to3; lattice <- round2_lattice(domain)
for (landscape in unique(old_manifest$landscape_id)) {
  column <- which(old_manifest$landscape_id == landscape)[1L]
  family <- old_manifest$family[column]
  replicate <- as.integer(sub(".*_L([0-9]+)$", "\\1", landscape))
  seed <- 110000L + 1000L * match(family, config$families) + replicate
  truth <- round2_truth("cn4_1to3", domain, family, seed, lattice)
  stopifnot(max(abs(truth$birth - old$truth_birth_rate[, column])) < 1e-12)
}
domain <- config$domains$cn3_1to4; lattice <- round2_lattice(domain)
founder <- which(rowSums(abs(lattice$karyotypes - 2L)) == 0)
stencil <- which(lattice$distance_matrix[founder, ] <= 1L)
auxiliary <- apply(lattice$distance_matrix[, stencil, drop = FALSE], 1, min)
stopifnot(length(lattice$labels) == 64L, max(auxiliary) == 5L)
for (fi in seq_along(config$families)) {
  truth <- round2_truth("cn3_1to4", domain, config$families[fi], domain$seed_namespace + 1000L * fi + 1L, lattice)
  stopifnot(all(truth$birth >= 0.5 & truth$birth <= 2), truth$birth[founder] == 1)
}
for (split in c("calibration", "test")) {
  output <- tempfile("sealed_generation_must_fail_")
  result <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(script_dir, "prepare_data.R")), "--domain", "cn3_1to4", "--split", split,
      "--landscapes", if (split == "calibration") "9" else "15", "--output", shQuote(output)), stdout = TRUE, stderr = TRUE))
  stopifnot(!is.null(attr(result, "status")), attr(result, "status") != 0L,
    any(grepl("SEALED", result)), !file.exists(file.path(output, "dataset_manifest.csv")))
}
fake_freeze <- tempfile(fileext = ".rds")
saveRDS(list(candidate_source_sha256 = "provisional", candidate_arguments_sha256 = "provisional",
  protocol_sha256 = "provisional", calibration_objects_sha256 = "provisional", analysis_sha256 = "provisional",
  power_status = "FINAL_N_FROZEN_BEFORE_NEW_TEST"), fake_freeze)
output <- tempfile("sealed_provisional_freeze_must_fail_")
result <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
  c(shQuote(file.path(script_dir, "prepare_data.R")), "--domain", "cn3_1to4", "--split", "test", "--landscapes", "15",
    "--generation-freeze-rds", shQuote(fake_freeze), "--output", shQuote(output)), stdout = TRUE, stderr = TRUE))
stopifnot(!is.null(attr(result, "status")), attr(result, "status") != 0L, any(grepl("SEALED", result)),
          !file.exists(file.path(output, "dataset_manifest.csv")))
unlink(fake_freeze)
cat("PASS: 4chr truth continuity; 3chr positive bounded link; disjoint seed namespaces; actual-distance primary; sealed calibration/test generation blocked\n")
