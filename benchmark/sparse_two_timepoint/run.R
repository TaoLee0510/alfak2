#!/usr/bin/env Rscript
# Sparse two-timepoint benchmark runner.
#
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=validate
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=truth    [--split=dev]
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=datasets [--split=dev]
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=fit      [--split=dev] [--methods=a,b]
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=score    [--split=dev]
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=report   [--split=dev]
#   Rscript benchmark/sparse_two_timepoint/run.R --stage=profile  [--split=dev] [--real-dir=DIR]
#
# Options: --workers=N (default 2), --output-dir=DIR, --cohort (also build the cohort
# datasets of every patient), --unseal-heldout (required for --split=heldout).
# Set STT_R_LIB to a library directory holding the alfak2 and alfakR builds
# under test; otherwise the installed packages are used.

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R",
               "report.R", "validate_abm.R")) {
  source(file.path(script_dir, part))
}

stt_parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) stop("Unrecognised argument: ", arg, call. = FALSE)
    body <- substring(arg, 3L)
    key <- sub("=.*$", "", body)
    out[[gsub("-", "_", key)]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  out
}

stt_parallel <- function(items, fun, workers) {
  results <- parallel::mclapply(items, function(item) {
    tryCatch(fun(item), error = function(e) structure(conditionMessage(e), class = "stt_task_error"))
  }, mc.cores = max(1L, workers), mc.preschedule = FALSE)
  failed <- vapply(results, function(x) inherits(x, "stt_task_error") || inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    stop(sum(failed), " task(s) failed; first error: ", as.character(results[[which(failed)[1L]]]), call. = FALSE)
  }
  results
}

stt_selected_landscapes <- function(config, split) {
  config$landscapes[config$landscapes$split == split, , drop = FALSE]
}

stt_dataset_plan <- function(config, split, cohort) {
  landscapes <- stt_selected_landscapes(config, split)
  plan <- list()
  for (i in seq_len(nrow(landscapes))) {
    for (patient in seq_len(config$patients_per_landscape)) {
      on_scoreboard <- patient %in% config$scoreboard_patients
      for (s in seq_len(nrow(config$scenarios))) {
        scenario <- config$scenarios[s, ]
        in_cohort <- cohort && identical(scenario$scenario, config$cohort_scenario)
        if (!on_scoreboard && !in_cohort) next
        plan[[length(plan) + 1L]] <- list(
          landscape = landscapes[i, ], patient = patient, scenario = scenario,
          scoreboard = on_scoreboard,
          seed = as.integer(landscapes$landscape_seed[i] %% 100000L * 1000L + patient * 20L + s)
        )
      }
    }
  }
  plan
}

stt_stage_truth <- function(config, split, workers, cohort) {
  landscapes <- stt_selected_landscapes(config, split)
  patients <- if (cohort) seq_len(config$patients_per_landscape) else config$scoreboard_patients
  tasks <- list()
  for (i in seq_len(nrow(landscapes))) for (patient in patients) {
    tasks[[length(tasks) + 1L]] <- list(landscape = landscapes[i, ], patient = patient)
  }
  stt_load_abm(config, script_dir)
  paths <- stt_parallel(tasks, function(task) stt_generate_truth(config, task$landscape, task$patient), workers)
  cat("Truth populations available:", length(paths), "\n")
}

stt_stage_datasets <- function(config, split, workers, cohort) {
  plan <- stt_dataset_plan(config, split, cohort)
  stt_load_abm(config, script_dir)  # the bottleneck scenario restarts the simulator
  paths <- stt_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, item$patient))
    stt_make_dataset(config, truth, item$scenario, config$t0, item$seed)
  }, workers)
  index <- do.call(rbind, lapply(seq_along(plan), function(i) {
    dataset <- readRDS(paths[[i]])
    data.frame(
      dataset_id = dataset$dataset_id, split = dataset$split, landscape_id = dataset$landscape_id,
      wavelength_class = dataset$wavelength_class, patient = dataset$patient,
      scenario = dataset$scenario, scoreboard = plan[[i]]$scoreboard,
      t(stt_count_profile(dataset$counts)),
      stringsAsFactors = FALSE
    )
  }))
  stt_parallel(index$dataset_id[index$scoreboard], function(dataset_id) {
    stt_dataset_nodes(config, readRDS(stt_dataset_path(config, dataset_id)))
    TRUE
  }, workers)
  index_path <- file.path(config$output_dir, sprintf("dataset_index_%s.csv", split))
  utils::write.csv(index, index_path, row.names = FALSE)
  cat("Datasets available:", nrow(index), " index:", index_path, "\n")
}

stt_fit_path <- function(config, dataset_id, method_name) {
  file.path(config$output_dir, "fits", sprintf("%s__%s.rds", dataset_id, method_name))
}

# Scored nodes of a dataset. Enumerating the two-step shell of a 22-chromosome
# sample is slow and memory-hungry, so it is done once and cached on disk.
stt_dataset_nodes <- function(config, dataset) {
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

stt_stage_fit <- function(config, split, workers, method_names) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)),
                           stringsAsFactors = FALSE)
  index <- index[index$scoreboard, , drop = FALSE]
  tasks <- list()
  for (dataset_id in index$dataset_id) for (method_name in method_names) {
    if (!file.exists(stt_fit_path(config, dataset_id, method_name))) {
      tasks[[length(tasks) + 1L]] <- list(dataset_id = dataset_id, method_name = method_name)
    }
  }
  cat("Fits to run:", length(tasks), "\n")
  if (!length(tasks)) return(invisible())
  dir.create(file.path(config$output_dir, "fits"), recursive = TRUE, showWarnings = FALSE)
  stt_parallel(tasks, function(task) {
    dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
    nodes <- stt_dataset_nodes(config, dataset)
    work_dir <- file.path(tempdir(), paste0(task$dataset_id, "__", task$method_name))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
    set.seed(dataset$seed)
    result <- stt_run_method(config$methods[[task$method_name]], dataset, nodes, work_dir)
    result$dataset_id <- task$dataset_id
    result$method <- task$method_name
    saveRDS(result, stt_fit_path(config, task$dataset_id, task$method_name))
    TRUE
  }, workers)
  invisible()
}

stt_stage_score <- function(config, split, workers, method_names) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)),
                           stringsAsFactors = FALSE)
  index <- index[index$scoreboard, , drop = FALSE]
  scored <- stt_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- stt_dataset_nodes(config, dataset)
    truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
    rows <- lapply(method_names, function(method_name) {
      path <- stt_fit_path(config, dataset_id, method_name)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      prediction <- if (is.null(fit$prediction)) {
        data.frame(karyotype = character(), fitness = numeric(), stringsAsFactors = FALSE)
      } else fit$prediction
      score <- stt_score_predictions(nodes, truth_rate, prediction, config$scoring$reference_min_count)
      cbind(
        data.frame(dataset_id = dataset_id, method = method_name, status = fit$status,
                   seconds = fit$seconds, stringsAsFactors = FALSE),
        score$metrics
      )
    })
    do.call(rbind, rows)
  }, workers)
  metrics <- merge(index, do.call(rbind, scored), by = "dataset_id")
  path <- file.path(config$output_dir, sprintf("metrics_%s.csv", split))
  utils::write.csv(metrics, path, row.names = FALSE)
  cat("Scored rows:", nrow(metrics), " file:", path, "\n")
}

# Compare the observable profile of the simulated datasets with real inputs.
stt_stage_profile <- function(config, split, real_dir) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)),
                           stringsAsFactors = FALSE)
  cat("Simulated datasets, scenario", config$cohort_scenario, "\n")
  print(stt_profile_quantiles(index[index$scenario == config$cohort_scenario, ]), digits = 2, row.names = FALSE)
  if (!is.null(real_dir)) {
    real <- stt_profile_directory(real_dir)
    real <- real[pmin(real$depth_t0, real$depth_t1) >= 100, ]
    cat("\nReal inputs with at least 100 cells at both timepoints (n =", nrow(real), ")\n")
    print(stt_profile_quantiles(real), digits = 2, row.names = FALSE)
  }
}

main <- function() {
  args <- stt_parse_args(commandArgs(trailingOnly = TRUE))
  library_dir <- Sys.getenv("STT_R_LIB", unset = "")
  if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- if (is.null(args$output_dir)) stt_config() else stt_config(args$output_dir)
  split <- if (is.null(args$split)) "dev" else args$split
  if (!split %in% config$landscapes$split) stop("Unknown split: ", split, call. = FALSE)
  if (identical(split, "heldout") && !isTRUE(args$unseal_heldout)) {
    stop("Held-out landscapes are reserved for the final confirmation; pass --unseal-heldout to use them.",
         call. = FALSE)
  }
  # Fits at 22 chromosomes need 1-4.5 GB each, so the default is conservative.
  workers <- if (is.null(args$workers)) 2L else as.integer(args$workers)
  cohort <- isTRUE(args$cohort)
  method_names <- if (is.null(args$methods)) names(config$methods) else strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  stage <- if (is.null(args$stage)) stop("--stage is required.", call. = FALSE) else args$stage
  switch(stage,
    validate = {
      stt_load_abm(config, script_dir)
      checks <- stt_validate_abm(config)
      print(checks, row.names = FALSE, right = FALSE)
      if (!all(checks$passed)) stop("Simulator validation failed.", call. = FALSE)
    },
    profile = stt_stage_profile(config, split, args$real_dir),
    truth = stt_stage_truth(config, split, workers, cohort),
    datasets = stt_stage_datasets(config, split, workers, cohort),
    fit = stt_stage_fit(config, split, workers, method_names),
    score = stt_stage_score(config, split, workers, method_names),
    report = stt_stage_report(config, split),
    stop("Unknown stage: ", stage, call. = FALSE)
  )
}

if (sys.nframe() == 0L) main()
