#!/usr/bin/env Rscript
# One active fit at a time. The matched policy fixes the eta prior across the
# paired count samples; it is an approximation diagnostic, not an SBC claim.
args <- commandArgs(TRUE)
argument <- function(name, default = NULL) {
  index <- match(paste0('--', name), args)
  if (is.na(index)) default else args[[index + 1L]]
}
repo <- argument('source', '/Users/4482173/Documents/GitHub/alfak2')
frozen <- argument('v1-source', '/private/tmp/alfak2_iteration2_20261104_0518/v1-source')
dll_path <- argument('dll', file.path(repo, 'src/alfak2.so'))
data_dir <- argument('data', '/private/tmp/alfak2_iteration2_20261104_0518/process-data')
output_dir <- argument('output', '/private/tmp/alfak2_iteration2_20261104_0518/process-results')
policy <- argument('prior', 'matched')
stopifnot(policy %in% c('matched', 'deployment'))
limit <- as.integer(argument('limit', '2147483647'))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(alfak2.warning_log_dir = file.path(output_dir, 'condition-logs'))
requireNamespace('Rcpp'); requireNamespace('TMB')
dll <- dyn.load(normalizePath(dll_path))
runtime <- new.env(parent = globalenv())
for (name in names(getDLLRegisteredRoutines(dll)$`.Call`))
  assign(name, getDLLRegisteredRoutines(dll)$`.Call`[[name]], runtime)
source_files <- sort(list.files(file.path(frozen, 'R'), '[.]R$', full.names = TRUE))
for (path in source_files) sys.source(path, runtime)
# Only the local wrapper gains the new diagnostic fields required by the DLL;
# oracle_initial_probability remains NULL, so no oracle information enters.
current_local_wrapper <- file.path(repo, 'R/local-fit.R')
sys.source(current_local_wrapper, runtime)
# Freeze the statistical prior at v1 while allowing the diagnostic-only pi0
# input fields (disabled here) required by the newly compiled objective.
frozen_prior <- new.env(parent = runtime)
sys.source(file.path(frozen, 'R/joint-landscape.R'), frozen_prior)
runtime$build_joint_landscape_prior <- frozen_prior$build_joint_landscape_prior
scoring <- new.env(parent = globalenv())
sys.source(file.path(repo, 'validation/accuracy_upgrade/metrics.R'), scoring)
manifest <- head(read.csv(file.path(data_dir, 'dataset_manifest.csv'), stringsAsFactors = FALSE), limit)
implementation_md5 <- tools::md5sum(c(source_files, current_local_wrapper, normalizePath(dll_path),
  file.path(frozen, 'R/joint-landscape.R'),
  file.path(repo, 'validation/accuracy_upgrade/metrics.R'),
  file.path(repo, 'validation/accuracy_round2/process_run_fits.R')))
saveRDS(list(implementation_md5 = implementation_md5, prior_policy = policy,
  purpose = 'development-only process attribution; not SBC, not candidate promotion',
  data_md5 = tools::md5sum(manifest$dataset_file), session_info = sessionInfo()),
  file.path(output_dir, paste0(policy, '_provenance.rds')))
all_metrics <- all_nodes <- all_status <- list()
for (i in seq_len(nrow(manifest))) {
  dataset <- readRDS(manifest$dataset_file[[i]])
  model <- paste0('process_', policy, '_', dataset$mode)
  cache <- file.path(output_dir, paste0(policy, '_', dataset$dataset_id, '.rds'))
  key <- list(implementation_md5 = implementation_md5, policy = policy,
              dataset_md5 = tools::md5sum(manifest$dataset_file[[i]]))
  cached <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(cached) && identical(cached$cache_key, key)) {
    result <- cached
  } else {
    fit_args <- list(counts = dataset$counts, dt = dataset$config$dt, beta = dataset$config$beta,
      min_cn = dataset$config$min_cn, max_cn = dataset$config$max_cn, max_nodes = 27L,
      target_karyotypes = dataset$target_karyotypes, pipeline = 'joint_count',
      forward_model = 'division_process', observation_model = 'multinomial',
      joint_control = dataset$config$joint_control)
    if (policy == 'matched') {
      fit_args$eta_prior_sd <- 5
      fit_args$eta_borrowed_prior_mean <- 0
      fit_args$eta_borrowed_prior_sd <- 5
      fit_args$eta_distance_penalty <- 0
    }
    warnings <- character()
    started <- proc.time()[['elapsed']]
    fit <- tryCatch(withCallingHandlers(do.call(runtime$fit_alfak2, fit_args), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    }), error = function(e) e)
    elapsed <- proc.time()[['elapsed']] - started
    if (inherits(fit, 'error')) {
      status <- data.frame(dataset_id = dataset$dataset_id, pair_id = dataset$pair_id,
        family = dataset$family, landscape_id = dataset$landscape_id, mode = dataset$mode,
        model = model, prior_policy = policy, fit_status = conditionMessage(fit),
        covariance_status = 'fit_failed', inner_gradient = NA_real_,
        interval_valid_fraction = 0, elapsed_seconds = elapsed, warnings = paste(warnings, collapse = ' | '))
      nodes <- metrics <- NULL
      fit <- NULL
    } else {
      nodes <- scoring$extract_target_predictions(fit, dataset, model)
      nodes$native_distance <- dataset$native_target_distance
      nodes$mode <- dataset$mode; nodes$pair_id <- dataset$pair_id
      metrics <- scoring$score_dataset(nodes, dataset, elapsed, 'ok', fit$local$diagnostics$covariance_status)
      metrics$mode <- dataset$mode; metrics$pair_id <- dataset$pair_id
      metrics$prior_policy <- policy
      status <- data.frame(dataset_id = dataset$dataset_id, pair_id = dataset$pair_id,
        family = dataset$family, landscape_id = dataset$landscape_id, mode = dataset$mode,
        model = model, prior_policy = policy, fit_status = 'ok',
        covariance_status = fit$local$diagnostics$covariance_status,
        inner_gradient = fit$local$diagnostics$conditional_gradient_norm,
        interval_valid_fraction = mean(fit$global$summary$interval_valid),
        elapsed_seconds = elapsed, warnings = paste(warnings, collapse = ' | '))
    }
    result <- list(cache_key = key, nodes = nodes, metrics = metrics, status = status, fit = fit,
                   dataset_id = dataset$dataset_id, fit_arguments = fit_args)
    saveRDS(result, cache)
  }
  all_metrics[[i]] <- result$metrics; all_nodes[[i]] <- result$nodes; all_status[[i]] <- result$status
  write.csv(do.call(rbind, all_status), file.path(output_dir, paste0(policy, '_status.csv')), row.names = FALSE)
  if (length(Filter(Negate(is.null), all_metrics)))
    write.csv(do.call(rbind, all_metrics), file.path(output_dir, paste0(policy, '_metrics.csv')), row.names = FALSE)
  if (length(Filter(Negate(is.null), all_nodes)))
    write.csv(do.call(rbind, all_nodes), file.path(output_dir, paste0(policy, '_nodes.csv')), row.names = FALSE)
  cat(sprintf('%s %d/%d %s %.3fs %s cov=%s\n', policy, i, nrow(manifest), dataset$dataset_id,
    result$status$elapsed_seconds, result$status$fit_status, result$status$covariance_status))
}
