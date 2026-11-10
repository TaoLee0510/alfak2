#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "config.R")); source(file.path(script_dir, "score_models.R"))
source(file.path(script_dir, "..", "accuracy_upgrade", "posterior_sampling.R"))
data_dir <- round2_arg("data"); output <- round2_arg("output")
if (is.null(data_dir) || is.null(output)) stop("Supply --data and --output")
models <- strsplit(round2_arg("models", "v1_joint"), ",", fixed = TRUE)[[1L]]
source_dir <- round2_arg("source"); dll <- round2_arg("dll"); library_dir <- round2_arg("library")
limit <- as.integer(round2_arg("limit", "2147483647")); mode <- round2_arg("select-mode")
config <- accuracy_round2_config(); slot_dir <- round2_arg("fit-slot-dir", config$slot_directory)
dir.create(output, recursive = TRUE, showWarnings = FALSE); dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
if (!is.null(mode)) manifest <- manifest[manifest$mode == mode, ]
manifest <- head(manifest, limit)
for (i in seq_len(nrow(manifest))) if (!file.exists(manifest$dataset_file[i])) manifest$dataset_file[i] <- file.path(data_dir, basename(manifest$dataset_file[i]))
runtime <- new.env(parent = globalenv())
if (!all(models == "zero_point_predictor")) {
  requireNamespace("Rcpp"); requireNamespace("TMB")
  if (!is.null(library_dir)) {
    library(alfak2, lib.loc = library_dir)
    runtime$fit_alfak2 <- alfak2::fit_alfak2
    implementation <- list(sha256 = round2_hash_object(round2_sha256(list.files(find.package("alfak2"), "[.](rdb|rdx|so)$", recursive = TRUE, full.names = TRUE))))
  } else {
    if (is.null(source_dir) || is.null(dll)) stop("Supply --source + --dll, or --library")
    implementation <- round2_source_fingerprint(source_dir, dll)
    loaded <- dyn.load(normalizePath(dll))
    routines <- getDLLRegisteredRoutines(loaded)$`.Call`
    for (name in names(routines)) assign(name, routines[[name]], runtime)
    for (file in sort(list.files(file.path(source_dir, "R"), "[.]R$", full.names = TRUE))) sys.source(file, runtime)
    stopifnot(identical(round2_source_fingerprint(source_dir, dll)$sha256, implementation$sha256))
  }
} else implementation <- list(sha256 = "analytic_zero_no_fit")
legacy_scoring <- new.env(parent = globalenv())
sys.source(file.path(script_dir, "..", "accuracy_upgrade", "metrics.R"), legacy_scoring)
scoring_sha <- round2_hash_object(round2_sha256(c(file.path(script_dir, c("run_models.R", "score_models.R", "common.R", "config.R")),
  file.path(script_dir, "..", "accuracy_upgrade", c("metrics.R", "posterior_sampling.R")))))
all_nodes <- all_metrics <- statuses <- list()
for (model in models) {
  arguments <- accuracy_round2_model_arguments(model)
  if (!is.null(round2_arg("fit-arguments-rds"))) arguments <- readRDS(round2_arg("fit-arguments-rds"))
  if (any(c("family", "truth", "truth_relative", "oracle_initial_probability", "population_growth") %in% names(arguments))) stop("Oracle/generating-family fit inputs forbidden")
  saveRDS(list(model = model, arguments = arguments, implementation = implementation, arguments_sha256 = round2_hash_object(arguments),
    scoring_sha256 = scoring_sha, manifest_sha256 = round2_sha256(file.path(data_dir, "dataset_manifest.csv")),
    session_info = sessionInfo(), scoring_distance = "native_target_distance when available; otherwise actual target_distance",
    oracle_information_passed = FALSE), file.path(output, paste0(model, "_provenance.rds")))
  for (i in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[i])
    if (!is.null(dataset$native_target_distance)) {
      dataset$auxiliary_common_distance <- dataset$target_distance
      dataset$target_distance <- dataset$native_target_distance
    }
    key <- list(model = model, implementation_sha256 = implementation$sha256, scoring_sha256 = scoring_sha,
      arguments_sha256 = round2_hash_object(arguments), dataset_sha256 = unname(round2_sha256(manifest$dataset_file[i])))
    cache <- file.path(output, paste0(model, "_", dataset$dataset_id, ".rds"))
    previous <- if (file.exists(cache)) readRDS(cache) else NULL
    if (!is.null(previous) && identical(previous$cache_key, key)) result <- previous else {
      fit <- NULL; started <- proc.time()[["elapsed"]]; warnings <- character()
      if (model == "zero_point_predictor") {
        nodes <- data.frame(dataset_id = dataset$dataset_id, family = dataset$family, landscape_id = dataset$landscape_id,
          split = dataset$split, model = model, karyotype = dataset$target_labels, distance = dataset$target_distance,
          truth_raw = dataset$truth, truth_relative = dataset$truth_relative, pred_raw = 1, pred_relative = 0,
          sd_raw = NA_real_, sd_relative = NA_real_, covariance_source = "point_only_no_intervals", interval_source = "unavailable")
        covariance_status <- "not_applicable_point_only"; fit_status <- "analytic_zero"
      } else {
        slot <- NULL
        for (number in 1:2) if (dir.create(file.path(slot_dir, paste0("slot", number)), showWarnings = FALSE)) { slot <- file.path(slot_dir, paste0("slot", number)); break }
        if (is.null(slot)) stop("Two fit slots already occupied; no third fit started")
        saveRDS(list(pid = Sys.getpid(), model = model, dataset = dataset$dataset_id), file.path(slot, "owner.rds"))
        fit_arguments <- list(counts = dataset$counts, dt = dataset$config$dt, beta = dataset$config$beta,
          min_cn = dataset$config$min_cn, max_cn = dataset$config$max_cn,
          max_nodes = length(dataset$target_labels), local_shell_depth = 1L,
          global_extra_shell = ncol(dataset$target_karyotypes) * (dataset$config$max_cn - dataset$config$min_cn))
        if ("target_karyotypes" %in% names(formals(runtime$fit_alfak2))) fit_arguments$target_karyotypes <- dataset$target_karyotypes
        for (name in names(arguments)) fit_arguments[[name]] <- arguments[[name]]
        fit <- tryCatch(withCallingHandlers(do.call(runtime$fit_alfak2, fit_arguments), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        }), error = identity)
        unlink(slot, recursive = TRUE)
        fit_status <- if (inherits(fit, "error")) conditionMessage(fit) else "ok"
        covariance_status <- if (inherits(fit, "error")) "fit_failed" else fit$local$diagnostics$covariance_status
        if (inherits(fit, "error")) {
          nodes <- data.frame(dataset_id = dataset$dataset_id, family = dataset$family, landscape_id = dataset$landscape_id,
            split = dataset$split, model = model, karyotype = dataset$target_labels, distance = dataset$target_distance,
            truth_raw = dataset$truth, truth_relative = dataset$truth_relative, pred_raw = NA_real_, pred_relative = NA_real_,
            sd_raw = NA_real_, sd_relative = NA_real_, covariance_source = "fit_failed", interval_source = "fit_failed")
          fit <- NULL
        } else nodes <- round2_extract_nodes(fit, dataset, model, legacy_scoring)
      }
      elapsed <- proc.time()[["elapsed"]] - started
      result <- list(cache_key = key, fit = fit, nodes = nodes, metrics = round2_score_nodes(nodes, dataset),
        status = data.frame(dataset_id = dataset$dataset_id, model = model, family = dataset$family, landscape_id = dataset$landscape_id,
          fit_status = fit_status, covariance_status = covariance_status, elapsed_seconds = elapsed,
          warnings = paste(unique(warnings), collapse = " | ")))
      saveRDS(result, cache)
    }
    id <- paste(model, dataset$dataset_id)
    all_nodes[[id]] <- result$nodes; all_metrics[[id]] <- result$metrics; statuses[[id]] <- result$status
    cat(sprintf("%s %d/%d %s %.3fs %s covariance=%s\n", model, i, nrow(manifest), dataset$dataset_id,
      result$status$elapsed_seconds, result$status$fit_status, result$status$covariance_status))
    round2_write_csv(do.call(rbind, all_metrics), file.path(output, "case_metrics.csv"))
    round2_write_csv(do.call(rbind, statuses), file.path(output, "fit_status.csv"))
  }
}
saveRDS(all_nodes, file.path(output, "node_tables.rds"))
