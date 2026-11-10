#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "metrics.R"))
args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
data_dir <- get_arg("data", file.path(script_dir, "data"))
output_dir <- get_arg("output", file.path(script_dir, "results"))
model <- get_arg("model", "candidate")
library_dir <- get_arg("library")
source_dir <- get_arg("source")
dll <- get_arg("dll")
extra_args <- get_arg("fit-args-rds")
limit <- as.integer(get_arg("limit", "2147483647"))
ids <- get_arg("ids")
workers <- as.integer(get_arg("workers", "1"))
worker_index <- as.integer(get_arg("worker-index", "1"))
stopifnot(workers >= 1L, worker_index >= 1L, worker_index <= workers)
file_prefix <- if (workers > 1L) sprintf("%s_worker%02d", model, worker_index) else model
if (!is.null(library_dir)) {
  library(alfak2, lib.loc = normalizePath(library_dir), character.only = FALSE)
  fit_function <- alfak2::fit_alfak2
  implementation <- normalizePath(find.package("alfak2"))
  implementation_files <- list.files(implementation, recursive = TRUE, full.names = TRUE)
  implementation_files <- implementation_files[grepl("\\.(so|dll|rdb|rdx)$|/DESCRIPTION$", implementation_files)]
} else if (!is.null(source_dir) && !is.null(dll)) {
  requireNamespace("Rcpp"); requireNamespace("TMB")
  dll_info <- dyn.load(normalizePath(dll))
  env <- new.env(parent = globalenv())
  routines <- getDLLRegisteredRoutines(dll_info)$`.Call`
  for (name in names(routines)) assign(name, routines[[name]], envir = env)
  for (file in sort(list.files(file.path(source_dir, "R"), "\\.R$", full.names = TRUE))) sys.source(file, env)
  fit_function <- env$fit_alfak2
  implementation <- normalizePath(source_dir)
  implementation_files <- c(list.files(file.path(source_dir, "R"), "\\.R$", full.names = TRUE),
                            list.files(file.path(source_dir, "src"), "\\.(cpp|h)$", full.names = TRUE), normalizePath(dll))
} else stop("Supply --library DIR or --source DIR --dll FILE")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(alfak2.warning_log_dir = file.path(output_dir, "warnings"))
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
if (!is.null(ids)) manifest <- manifest[manifest$dataset_id %in% strsplit(ids, ",", fixed = TRUE)[[1L]], , drop = FALSE]
manifest <- head(manifest, limit)
if (workers > 1L) manifest <- manifest[(seq_len(nrow(manifest)) - 1L) %% workers == worker_index - 1L, , drop = FALSE]
user_args <- if (is.null(extra_args)) list() else readRDS(extra_args)
fingerprint <- function(x) {
  temporary <- tempfile("alfak2_benchmark_fingerprint_")
  on.exit(unlink(temporary))
  saveRDS(x, temporary, version = 2)
  unname(tools::md5sum(temporary))
}
implementation_hash <- fingerprint(tools::md5sum(sort(implementation_files)))
arguments_hash <- fingerprint(user_args)
scoring_hash <- fingerprint(tools::md5sum(file.path(script_dir, c("metrics.R", "run_pilot.R"))))
saveRDS(list(model = model, implementation = implementation,
             implementation_hash = implementation_hash, arguments_hash = arguments_hash, scoring_hash = scoring_hash,
             session_info = utils::sessionInfo(), extra_fit_arguments = user_args,
             data_manifest_md5 = unname(tools::md5sum(file.path(data_dir, "dataset_manifest.csv"))),
             data_md5 = tools::md5sum(manifest$dataset_file)), file.path(output_dir, paste0(file_prefix, "_provenance.rds")))
all_metrics <- all_nodes <- statuses <- list()
for (i in seq_len(nrow(manifest))) {
  dataset <- readRDS(manifest$dataset_file[i])
  started <- proc.time()[["elapsed"]]
  cache <- file.path(output_dir, paste0(model, "_", dataset$dataset_id, ".rds"))
  cache_key <- list(model = model, implementation_hash = implementation_hash, arguments_hash = arguments_hash,
                    scoring_hash = scoring_hash, dataset_hash = unname(tools::md5sum(manifest$dataset_file[i])))
  cached <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(cached) && identical(cached$scoring_schema_version, 3L) && identical(cached$cache_key, cache_key)) {
    result <- cached
    cat(sprintf("%s %d/%d %s cached\n", model, i, nrow(manifest), dataset$dataset_id))
  } else {
    fit_args <- c(list(counts = dataset$counts, dt = dataset$config$dt, beta = dataset$config$beta,
                       min_cn = dataset$config$min_cn, max_cn = dataset$config$max_cn,
                       max_nodes = dataset$config$max_nodes,
                       local_shell_depth = dataset$config$local_shell_depth,
                       global_extra_shell = dataset$config$n_chr * (dataset$config$max_cn - dataset$config$min_cn)),
                  dataset$config$fit_arguments)
    for (name in names(user_args)) fit_args[[name]] <- user_args[[name]]
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(do.call(fit_function, fit_args), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    runtime <- proc.time()[["elapsed"]] - started
    fit_status <- if (inherits(fit, "error")) conditionMessage(fit) else "ok"
    covariance_status <- if (inherits(fit, "error")) "fit_failed" else fit$local$diagnostics$covariance_status
    if (is.null(covariance_status)) covariance_status <- "not_reported"
    if (!inherits(fit, "error")) {
      nodes <- tryCatch(extract_target_predictions(fit, dataset, model, recover_baseline = model == "baseline"), error = function(e) e)
      if (inherits(nodes, "error")) {
        warnings <- c(warnings, paste("Covariance extraction failed:", conditionMessage(nodes)))
        fit$global$fitness_covariance <- NULL; fit$global$covariance <- NULL
        nodes <- extract_target_predictions(fit, dataset, model, recover_baseline = FALSE)
      }
      frequency <- fit$local$summary
      idx <- match(dataset$target_labels, frequency$karyotype)
      frequency_metrics <- data.frame(
        dataset_id = dataset$dataset_id, family = dataset$family, landscape_id = dataset$landscape_id,
        split = dataset$split, model = model, distance_group = "local_target",
        fit_status = fit_status, covariance_status = covariance_status,
        metric = c("pi0_rmse_local", "pi1_rmse_local", "pi0_target_mass", "pi1_target_mass"),
        value = c(rmse(frequency$pi0[idx] - dataset$truth_pi0), rmse(frequency$pi1[idx] - dataset$truth_pi1),
                  sum(dataset$truth_pi0[!is.na(idx)]), sum(dataset$truth_pi1[!is.na(idx)])), stringsAsFactors = FALSE)
      predictive_metrics <- heldout_count_metrics(fit, dataset, model, fit_status, covariance_status)
    } else {
      nodes <- data.frame(dataset_id = dataset$dataset_id, family = dataset$family,
                          landscape_id = dataset$landscape_id, split = dataset$split, model = model,
                          karyotype = dataset$target_labels, distance = dataset$target_distance,
                          truth_raw = dataset$truth, truth_relative = dataset$truth_relative,
                          pred_raw = NA_real_, pred_relative = NA_real_, sd_raw = NA_real_, sd_relative = NA_real_,
                          covariance_source = "fit_failed", stringsAsFactors = FALSE)
      for (level in dataset$config$interval_levels) {
        suffix <- sprintf("%02d", round(100 * level))
        nodes[[paste0("relative_low_", suffix)]] <- NA_real_
        nodes[[paste0("relative_high_", suffix)]] <- NA_real_
      }
      nodes$interval_source <- "fit_failed"
      frequency_metrics <- NULL
      predictive_metrics <- NULL
    }
    metrics <- score_dataset(nodes, dataset, runtime, fit_status, covariance_status)
    if (!is.null(frequency_metrics)) metrics <- rbind(metrics, frequency_metrics)
    if (!is.null(predictive_metrics)) metrics <- rbind(metrics, predictive_metrics)
    result <- list(scoring_schema_version = 3L, cache_key = cache_key, metrics = metrics, nodes = nodes,
                   status = data.frame(dataset_id = dataset$dataset_id, model = model,
                                        fit_status = fit_status, covariance_status = covariance_status,
                                        runtime_seconds = runtime, warnings = paste(unique(warnings), collapse = " | "),
                                        stringsAsFactors = FALSE),
                   diagnostics = if (!inherits(fit, "error")) fit$diagnostics else NULL,
                   local_diagnostics = if (!inherits(fit, "error")) fit$local$diagnostics else NULL,
                   global_hyperparameters = if (!inherits(fit, "error")) fit$global$hyperparameters else NULL)
    saveRDS(result, cache)
    cat(sprintf("%s %d/%d %s %.2fs %s cov=%s\n", model, i, nrow(manifest), dataset$dataset_id,
                runtime, fit_status, covariance_status))
  }
  all_metrics[[i]] <- result$metrics; all_nodes[[i]] <- result$nodes; statuses[[i]] <- result$status
  utils::write.csv(do.call(rbind, all_metrics), file.path(output_dir, paste0(file_prefix, "_metrics.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, all_nodes), file.path(output_dir, paste0(file_prefix, "_nodes.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, statuses), file.path(output_dir, paste0(file_prefix, "_status.csv")), row.names = FALSE)
}
