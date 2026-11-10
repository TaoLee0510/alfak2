#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
root <- normalizePath(file.path(script_dir, "..", ".."))
work <- round2_arg("work", "/private/tmp/alfak2_iteration2_20261104_0518")
output <- file.path(script_dir, "results", "development81_comparison", "archive")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source_snapshot <- file.path(work, "development81-exact-source")
reference <- readRDS(file.path(work, "development81-B", "B_provenance.rds"))
manifest <- reference$implementation$manifest
source_choices <- list()
for (i in seq_len(nrow(manifest))) {
  relative <- manifest$file[i]
  choices <- if (relative == "binary/alfak2.so") file.path(work, "current-lib", "alfak2", "libs", "alfak2.so") else
    c(file.path(root, relative), file.path(work, "v1-source", relative))
  choices <- choices[file.exists(choices)]
  matches <- choices[unname(round2_sha256(choices)) == manifest$sha256[i]]
  if (!length(matches)) stop(paste("Cannot recover exact fitted source:", relative))
  destination <- file.path(source_snapshot, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination)) {
    if (unname(round2_sha256(destination)) != manifest$sha256[i]) stop("Refusing to overwrite a different source snapshot")
  } else if (!file.copy(matches[1L], destination)) stop("Source archive copy failed")
  source_choices[[i]] <- data.frame(file = relative, sha256 = manifest$sha256[i], recovered_from = matches[1L])
}
fingerprint <- round2_source_fingerprint(source_snapshot, file.path(source_snapshot, "binary", "alfak2.so"))
stopifnot(identical(fingerprint$sha256, reference$implementation$sha256))
round2_write_csv(do.call(rbind, source_choices), file.path(output, "fitted_source_manifest.csv"))

for (model in c("A", "B", "current_default_two_layer")) {
  from <- file.path(work, paste0("development81-", if (model == "current_default_two_layer") "currentdefault" else model), paste0(model, "_provenance.rds"))
  p <- readRDS(from)
  if (model %in% c("A", "B")) stopifnot(identical(p$implementation$sha256, reference$implementation$sha256))
  file.copy(from, file.path(output, basename(from)), overwrite = TRUE)
}
scoring_files <- c(file.path(script_dir, c("run_models.R", "score_models.R", "common.R", "config.R")),
  file.path(script_dir, "..", "accuracy_upgrade", c("metrics.R", "posterior_sampling.R")))
scoring_hashes <- round2_sha256(scoring_files)
stopifnot(identical(round2_hash_object(scoring_hashes), reference$scoring_sha256))
round2_write_csv(data.frame(original_file = names(scoring_hashes), sha256 = unname(scoring_hashes)),
  file.path(output, "fitted_scoring_manifest.csv"))
for (file in scoring_files) {
  relative <- if (grepl("accuracy_upgrade", file, fixed = TRUE)) file.path("accuracy_upgrade", basename(file)) else file.path("accuracy_round2", basename(file))
  destination <- file.path(output, "scoring_source", relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  file.copy(file, destination, overwrite = TRUE)
}
file.copy(file.path(work, "development81-v1-reuse", "reuse_audit.csv"), file.path(output, "frozen_v1_reuse_audit.csv"), overwrite = TRUE)
file.copy(file.path(work, "development81-data", "dataset_manifest.csv"), file.path(output, "imported_dataset_manifest.csv"), overwrite = TRUE)

source_tar <- file.path(output, "fitted_R_CPP_source.tar.gz")
if (!file.exists(source_tar)) {
  previous <- getwd(); setwd(source_snapshot)
  utils::tar(source_tar, files = manifest$file[manifest$file != "binary/alfak2.so"], compression = "gzip", tar = "internal")
  setwd(previous)
}
data_tar <- file.path(output, "development81_datasets.tar.gz")
if (!file.exists(data_tar)) {
  previous <- getwd(); setwd(file.path(work, "development81-data"))
  utils::tar(data_tar, files = list.files(), compression = "gzip", tar = "internal")
  setwd(previous)
}
nodes <- list()
for (entry in c("v1-reuse", "A", "B", "currentdefault", "zero"))
  nodes <- c(nodes, readRDS(file.path(work, paste0("development81-", entry), "node_tables.rds")))
saveRDS(nodes, file.path(output, "node_tables.rds"), version = 2)
record <- list(role = "development_only_no_new_calibration_or_test",
  source_sha256 = reference$implementation$sha256,
  fitted_scoring_sha256 = reference$scoring_sha256,
  exact_source_snapshot = source_snapshot,
  exact_dll = file.path(source_snapshot, "binary", "alfak2.so"),
  dll_sha256 = manifest$sha256[manifest$file == "binary/alfak2.so"],
  dll_archived_in_repository = FALSE,
  source_note = "All source files match fitted SHA256; fit-alfak2.R recovered from frozen v1-source because current API docs changed after fits",
  stage = "no_single_winner_frozen; nested calibration and joint power deferred",
  generated_at = format(Sys.time(), tz = "America/New_York", usetz = TRUE))
saveRDS(record, file.path(output, "archive_record.rds"), version = 2)
artifacts <- list.files(output, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "artifact_manifest.csv"]
round2_write_csv(data.frame(file = substring(artifacts, nchar(output) + 2L), bytes = file.info(artifacts)$size,
  sha256 = unname(round2_sha256(artifacts))), file.path(output, "artifact_manifest.csv"))
cat("PASS: exact fitted source/DLL fingerprint recovered; scoring SHA matched; 90 node tables and 18 datasets archived; no new splits opened\n")
