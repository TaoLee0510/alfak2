# Stages of the short-interval diagnosis. Sourced by run.R.

si_parallel <- function(items, fun, workers) {
  results <- parallel::mclapply(items, function(item) {
    tryCatch(fun(item), error = function(e) structure(conditionMessage(e), class = "si_task_error"))
  }, mc.cores = max(1L, workers), mc.preschedule = FALSE)
  failed <- vapply(results, function(x) inherits(x, "si_task_error") || inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    stop(sum(failed), " task(s) failed; first error: ", as.character(results[[which(failed)[1L]]]), call. = FALSE)
  }
  results
}

si_stage_truth <- function(config, workers) {
  landscapes <- config$landscapes
  tasks <- lapply(seq_len(nrow(landscapes)), function(i) landscapes[i, ])
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  paths <- si_parallel(tasks, function(landscape) stt_generate_truth(config, landscape, 1L), workers)
  cat("Truth populations available:", length(paths), "\n")
}

si_stage_datasets <- function(config, workers) {
  plan <- si_dataset_plan(config)
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  paths <- si_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
    path <- stt_make_dataset(config, truth, item$scenario, item$scenario$t0, item$seed)
    stt_dataset_nodes_si(config, readRDS(path))
    path
  }, workers)
  index <- do.call(rbind, lapply(seq_along(plan), function(i) {
    dataset <- readRDS(paths[[i]])
    data.frame(
      dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
      wavelength_class = dataset$wavelength_class, scenario = dataset$scenario,
      t0 = dataset$t0, interval = dataset$interval,
      t(stt_count_profile(dataset$counts)),
      selection_signal = si_selection_signal(dataset),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(index, si_index_path(config), row.names = FALSE)
  cat("Datasets available:", nrow(index), " index:", si_index_path(config), "\n")
}

# Scored nodes, cached as in the development benchmark.
stt_dataset_nodes_si <- function(config, dataset) {
  path <- file.path(config$output_dir, "nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  scoring <- config$scoring
  nodes <- stt_scoring_nodes(dataset$counts, scoring$min_cn, scoring$max_cn, scoring$d1_cap,
                             scoring$d2_cap, dataset$seed)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(nodes, temporary)
  file.rename(temporary, path)
  nodes
}

si_stage_fit <- function(config, workers, method_names, datasets = NULL) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  if (!is.null(datasets)) index <- index[grepl(datasets, index$dataset_id), , drop = FALSE]
  builds <- unique(unlist(lapply(config$methods[method_names], function(m) m$build)))
  if (length(builds) > 1L) stop("Methods of one run must share an alfak2 build.", call. = FALSE)
  if (length(builds)) {
    library_dir <- si_library(builds, if (is.null(config$library_root)) config$output_dir else config$library_root)
    if (!dir.exists(file.path(library_dir, "alfak2"))) stop("No alfak2 build in ", library_dir, call. = FALSE)
    .libPaths(c(library_dir, .libPaths()))
    loadNamespace("alfak2")
    cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  }
  tasks <- list()
  for (dataset_id in index$dataset_id) for (method_name in method_names) {
    if (!file.exists(si_fit_path(config, dataset_id, method_name))) {
      tasks[[length(tasks) + 1L]] <- list(dataset_id = dataset_id, method_name = method_name)
    }
  }
  cat("Fits to run:", length(tasks), "\n")
  if (!length(tasks)) return(invisible())
  si_parallel(tasks, function(task) {
    dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
    nodes <- stt_dataset_nodes_si(config, dataset)
    set.seed(dataset$seed)
    result <- si_run_method(config, dataset, nodes, config$methods[[task$method_name]])
    result$dataset_id <- task$dataset_id
    result$method <- task$method_name
    path <- si_fit_path(config, task$dataset_id, task$method_name)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(result, temporary)
    file.rename(temporary, path)
    TRUE
  }, workers)
  invisible()
}

# Two tables over every fit present: per-distance scores (long) and one row of
# fit, cross-validation and landscape metrics per dataset and method.
si_stage_score <- function(config, workers) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  method_names <- intersect(names(config$methods), list.files(file.path(config$output_dir, "fits")))
  scored <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- stt_dataset_nodes_si(config, dataset)
    rows <- lapply(method_names, function(method_name) {
      path <- si_fit_path(config, dataset_id, method_name)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      score <- si_score_fit(config, dataset, nodes, fit)
      key <- data.frame(dataset_id = dataset_id, method = method_name, status = fit$status,
                        stringsAsFactors = FALSE)
      scales <- if (is.null(fit$scales)) c(dosage = NA, curvature = NA, residual = NA) else
        fit$scales[c("dosage", "curvature", "residual")]
      names(scales) <- paste0("scale_", names(scales))
      d <- fit$diagnostics
      fit_row <- cbind(key, t(c(score$extra, scales,
                                objective = if (is.null(d)) NA else d$objective,
                                convergence = if (is.null(d)) NA else d$convergence,
                                n_states = if (is.null(d)) NA else d$n_states,
                                n_major = if (is.null(d)) NA else d$n_major,
                                history = if (is.null(d)) NA else d$history,
                                fit_seconds = if (is.null(fit$fit_seconds)) NA else fit$fit_seconds,
                                cv_seconds = if (is.null(fit$cv_seconds)) NA else fit$cv_seconds)))
      if (!is.null(fit$cv)) fit_row <- cbind(fit_row, t(fit$cv))
      fit_row$trend_terms <- if (is.null(fit$trend_terms)) NA_character_ else fit$trend_terms
      list(distance = cbind(key, score$metrics), fit = fit_row)
    })
    rows <- Filter(Negate(is.null), rows)
    list(distance = si_bind(lapply(rows, `[[`, "distance")), fit = si_bind(lapply(rows, `[[`, "fit")))
  }, workers)
  distance <- merge(index, si_bind(lapply(scored, `[[`, "distance")), by = "dataset_id")
  fits <- merge(index, si_bind(lapply(scored, `[[`, "fit")), by = "dataset_id")
  utils::write.csv(distance, file.path(config$output_dir, "metrics_distance.csv"), row.names = FALSE)
  utils::write.csv(fits, file.path(config$output_dir, "metrics_fit.csv"), row.names = FALSE)
  cat("Scored:", nrow(fits), "fits in", file.path(config$output_dir, "metrics_fit.csv"), "\n")
}

si_bind <- function(rows) {
  rows <- Filter(function(r) !is.null(r) && nrow(r) > 0L, rows)
  if (!length(rows)) return(NULL)
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(columns, names(r))] <- NA
    r[columns]
  }))
}
