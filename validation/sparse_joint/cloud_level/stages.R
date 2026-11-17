# Stages of the cloud-level study. Sourced by run.R.

cl_parallel <- function(items, fun, workers) {
  results <- parallel::mclapply(items, function(item) {
    tryCatch(fun(item), error = function(e) structure(conditionMessage(e), class = "cl_task_error"))
  }, mc.cores = max(1L, workers), mc.preschedule = FALSE)
  failed <- vapply(results, function(x) inherits(x, "cl_task_error") || inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    stop(sum(failed), " task(s) failed; first error: ", as.character(results[[which(failed)[1L]]]), call. = FALSE)
  }
  results
}

cl_use_build <- function() {
  library_dir <- cl_library()
  if (!dir.exists(file.path(library_dir, "alfak2"))) stop("No alfak2 build in ", library_dir, call. = FALSE)
  .libPaths(c(library_dir, .libPaths()))
  loadNamespace("alfak2")
  cat("alfak2 from", dirname(find.package("alfak2")), "\n")
}

cl_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(object, temporary)
  file.rename(temporary, path)
  invisible(path)
}

# Copies of the development scoreboard and of the short-interval diagnosis
# set: datasets, scored nodes and one index with the columns the report uses.
cl_stage_inputs <- function() {
  cl_use_build()
  sources <- list(
    dev = list(dir = file.path(cl_repo, "benchmark", "results", "sparse_two_timepoint"), index = "dataset_index_dev.csv"),
    short = list(dir = file.path(cl_repo, "validation", "sparse_joint", "short_interval", "results"),
                 index = "dataset_index.csv")
  )
  for (set in names(sources)) {
    source_dir <- sources[[set]]$dir
    index <- utils::read.csv(file.path(source_dir, sources[[set]]$index), stringsAsFactors = FALSE)
    if (set == "dev") index <- index[index$scoreboard, , drop = FALSE]
    for (part in c("datasets", "nodes")) {
      dir.create(file.path(cl_input_dir(set), part), recursive = TRUE, showWarnings = FALSE)
      from <- file.path(source_dir, part, paste0(index$dataset_id, ".rds"))
      if (!all(file.exists(from))) stop("Missing ", part, " of set ", set, call. = FALSE)
      file.copy(from, file.path(cl_input_dir(set), part), overwrite = FALSE)
    }
    rows <- do.call(rbind, lapply(index$dataset_id, function(id) {
      dataset <- readRDS(cl_dataset_path(set, id))
      data.frame(dataset_id = id, landscape_id = dataset$landscape_id, wavelength_class = dataset$wavelength_class,
                 scenario = dataset$scenario, interval = dataset$interval,
                 role = if (set == "dev") index$role[index$dataset_id == id] else "short",
                 t(cl_diversity(dataset$counts)), stringsAsFactors = FALSE)
    }))
    utils::write.csv(rows, cl_index_path(set), row.names = FALSE)
    cat("Set", set, ":", nrow(rows), "datasets\n")
  }
}

# The populations of a simulated set (diverse, or a confirmation set), one per
# landscape.
cl_stage_truth <- function(workers, landscapes = NULL, set = "diverse") {
  config <- cl_set_config(set)
  rows <- config$landscapes
  if (!is.null(landscapes)) rows <- rows[grepl(landscapes, rows$landscape_id), , drop = FALSE]
  stt_load_abm(config, file.path(cl_repo, "benchmark", "sparse_two_timepoint"))
  cl_parallel(lapply(seq_len(nrow(rows)), function(i) rows[i, ]), function(landscape) {
    started <- proc.time()[["elapsed"]]
    level_config <- config
    if (!is.null(landscape$beta_multiple)) level_config$beta <- cl_assumed_beta * landscape$beta_multiple
    path <- stt_generate_truth(level_config, landscape, 1L)
    cat(landscape$landscape_id, "truth in", round(proc.time()[["elapsed"]] - started), "s\n")
    path
  }, workers)
  invisible()
}

cl_stage_datasets <- function(workers, landscapes = NULL, set = "diverse") {
  if (!is.null(cl_confirm_sets[[set]]) && cl_confirm_sets[[set]]$kind != "diverse") {
    return(cl_stage_scenario_datasets(workers, set))
  }
  cl_use_build()
  config <- cl_set_config(set)
  rows <- config$landscapes
  if (!is.null(landscapes)) rows <- rows[grepl(landscapes, rows$landscape_id), , drop = FALSE]
  rows <- rows[file.exists(stt_truth_path(config, rows$landscape_id, 1L)), , drop = FALSE]
  stt_load_abm(config, file.path(cl_repo, "benchmark", "sparse_two_timepoint"))
  scenario <- config$scenarios[1L, ]
  index <- cl_parallel(lapply(seq_len(nrow(rows)), function(i) rows[i, ]), function(landscape) {
    truth <- readRDS(stt_truth_path(config, landscape$landscape_id, 1L))
    seed <- as.integer(landscape$landscape_seed %% 100000L * 1000L + 21L)
    path <- stt_make_dataset(config, truth, scenario, scenario$t0, seed)
    dataset <- readRDS(path)
    if (is.null(dataset$level)) {
      dataset$level <- landscape$level
      dataset$mechanism <- landscape$mechanism
      dataset$beta_assumed <- cl_assumed_beta
      dataset$call_error <- landscape$call_error
      if (landscape$call_error > 0) {
        dataset$counts_before_errors <- dataset$counts
        dataset$counts <- cl_call_errors(dataset$counts, landscape$call_error, seed + 7L,
                                         config$abm$min_cn, config$abm$max_cn)
      }
      saveRDS(dataset, path)
    }
    nodes_path <- cl_nodes_path(set, dataset$dataset_id)
    if (!file.exists(nodes_path)) {
      scoring <- config$scoring
      cl_save(stt_scoring_nodes(dataset$counts, scoring$min_cn, scoring$max_cn, scoring$d1_cap, scoring$d2_cap,
                                dataset$seed), nodes_path)
    }
    data.frame(dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
               wavelength_class = dataset$level, scenario = dataset$scenario, interval = dataset$interval,
               role = dataset$mechanism, t(cl_diversity(dataset$counts)), stringsAsFactors = FALSE)
  }, workers)
  index <- do.call(rbind, index)
  utils::write.csv(index, cl_index_path(set), row.names = FALSE)
  cat("Datasets of", set, ":", nrow(index), "\n")
  print(stats::aggregate(cbind(singleton_t0, top_t0, coverage_t0, coverage_t1, n_major, n_karyotypes) ~ wavelength_class,
                         data = index, FUN = function(v) round(stats::median(v), 2)))
}

# Datasets of a confirmation set with the dev or short design: one per
# landscape and scenario, seeded as the development benchmark seeds its own.
cl_stage_scenario_datasets <- function(workers, set) {
  cl_use_build()
  config <- cl_set_config(set)
  stt_load_abm(config, file.path(cl_repo, "benchmark", "sparse_two_timepoint"))
  plan <- list()
  for (i in seq_len(nrow(config$landscapes))) for (s in seq_len(nrow(config$scenarios))) {
    plan[[length(plan) + 1L]] <- list(landscape = config$landscapes[i, ], scenario = config$scenarios[s, ],
                                      seed = as.integer(config$landscapes$landscape_seed[i] %% 100000L * 1000L + 20L + s))
  }
  index <- cl_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
    t0 <- if (is.null(item$scenario$t0)) config$t0 else item$scenario$t0
    path <- stt_make_dataset(config, truth, item$scenario, t0, item$seed)
    dataset <- readRDS(path)
    # Conditions of the paper-sized set: a true rate other than the assumed
    # one, or calling errors added to the sampled cells.
    landscape <- item$landscape
    if (!is.null(landscape$level) && is.null(dataset$level)) {
      dataset$level <- landscape$level
      dataset$mechanism <- landscape$mechanism
      dataset$beta_assumed <- cl_assumed_beta
      dataset$call_error <- landscape$call_error
      if (landscape$call_error > 0) {
        dataset$counts_before_errors <- dataset$counts
        dataset$counts <- cl_call_errors(dataset$counts, landscape$call_error, item$seed + 7L,
                                         config$abm$min_cn, config$abm$max_cn)
      }
      saveRDS(dataset, path)
    }
    nodes_path <- cl_nodes_path(set, dataset$dataset_id)
    if (!file.exists(nodes_path)) {
      scoring <- config$scoring
      cl_save(stt_scoring_nodes(dataset$counts, scoring$min_cn, scoring$max_cn, scoring$d1_cap, scoring$d2_cap,
                                dataset$seed), nodes_path)
    }
    data.frame(dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
               wavelength_class = dataset$wavelength_class, scenario = dataset$scenario, interval = dataset$interval,
               role = if (is.null(dataset$level)) item$scenario$role else dataset$level,
               t(cl_diversity(dataset$counts)), stringsAsFactors = FALSE)
  }, workers)
  index <- do.call(rbind, index)
  utils::write.csv(index, cl_index_path(set), row.names = FALSE)
  cat("Datasets of", set, ":", nrow(index), "\n")
}

# One line per fit: set, dataset and method.
cl_stage_tasks <- function(method_names, sets = cl_sets, tasks_file = "tasks.tsv") {
  tasks <- do.call(rbind, lapply(sets, function(set) {
    index <- utils::read.csv(cl_index_path(set), stringsAsFactors = FALSE)
    expand.grid(set = set, dataset_id = index$dataset_id, method = method_names, stringsAsFactors = FALSE)
  }))
  utils::write.table(tasks, cl_tasks_path(tasks_file), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Tasks:", nrow(tasks), "in", cl_tasks_path(tasks_file), "\n")
  print(table(tasks$set, tasks$method))
}

cl_fit_task <- function(set, dataset_id, method) {
  path <- cl_fit_path(set, method, dataset_id)
  if (file.exists(path)) return(invisible(path))
  dataset <- readRDS(cl_dataset_path(set, dataset_id))
  nodes <- readRDS(cl_nodes_path(set, dataset_id))
  set.seed(dataset$seed)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(cl_fit(dataset, nodes, cl_methods[[method]]),
                     error = function(e) list(status = "failed", detail = conditionMessage(e)))
  result$set <- set
  result$dataset_id <- dataset_id
  result$method <- method
  result$total_seconds <- proc.time()[["elapsed"]] - started
  cl_save(result, path)
  cat(set, dataset_id, method, result$status, round(result$total_seconds), "s\n")
  invisible(path)
}

# `task`: a line of tasks.tsv (one Slurm array task). Otherwise every missing
# fit of `sets` and `method_names`, optionally those whose dataset matches
# `datasets`.
cl_stage_fit <- function(workers, task = NULL, method_names = NULL, sets = cl_sets, datasets = NULL,
                         tasks_file = "tasks.tsv") {
  cl_use_build()
  if (!is.null(task)) {
    row <- utils::read.delim(cl_tasks_path(tasks_file), stringsAsFactors = FALSE)[as.integer(task), ]
    if (is.na(row$set)) stop("No task ", task, call. = FALSE)
    return(cl_fit_task(row$set, row$dataset_id, row$method))
  }
  tasks <- list()
  for (set in sets) {
    ids <- utils::read.csv(cl_index_path(set), stringsAsFactors = FALSE)$dataset_id
    if (!is.null(datasets)) ids <- ids[grepl(datasets, ids)]
    for (id in ids) for (method in method_names) {
      if (!file.exists(cl_fit_path(set, method, id))) tasks[[length(tasks) + 1L]] <- list(set, id, method)
    }
  }
  cat("Fits to run:", length(tasks), "\n")
  cl_parallel(tasks, function(t) cl_fit_task(t[[1L]], t[[2L]], t[[3L]]), workers)
  invisible()
}

# One row per fit present.
cl_stage_score <- function(workers, sets = cl_sets, out_file = "scores.csv") {
  cl_use_build()
  rows <- list()
  for (set in sets) {
    if (!file.exists(cl_index_path(set))) next
    index <- utils::read.csv(cl_index_path(set), stringsAsFactors = FALSE)
    methods <- intersect(names(cl_methods), list.files(file.path(cl_output, "fits", set)))
    scored <- cl_parallel(index$dataset_id, function(id) {
      dataset <- readRDS(cl_dataset_path(set, id))
      nodes <- readRDS(cl_nodes_path(set, id))
      do.call(rbind, lapply(methods, function(method) {
        path <- cl_fit_path(set, method, id)
        if (!file.exists(path)) return(NULL)
        fit <- readRDS(path)
        key <- data.frame(set = set, dataset_id = id, method = method, status = fit$status, stringsAsFactors = FALSE)
        if (fit$status != "ok") return(key)
        cbind(key, t(cl_score(dataset, nodes, fit)))
      }))
    }, workers)
    scored <- Filter(Negate(is.null), scored)
    columns <- unique(unlist(lapply(scored, names)))
    scored <- do.call(rbind, lapply(scored, function(r) {
      r[setdiff(columns, names(r))] <- NA
      r[columns]
    }))
    rows[[set]] <- merge(index, scored, by = "dataset_id")
  }
  columns <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(r) {
    r[setdiff(columns, names(r))] <- NA
    r[columns]
  }))
  utils::write.csv(out, file.path(cl_output, out_file), row.names = FALSE)
  cat("Scored:", nrow(out), "fits in", file.path(cl_output, out_file), "\n")
  print(table(out$set, out$method, out$status))
}
