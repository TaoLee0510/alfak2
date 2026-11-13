#!/usr/bin/env Rscript
# The development benchmark (benchmark/sparse_two_timepoint, development
# landscapes, 8 core and 6 stress scenarios) for the arms of the
# short-interval study, in a directory of its own: datasets and scored nodes
# are copied from the benchmark's results, which are not written to.
#
#   R=validation/sparse_joint/short_interval/devbench.R
#   Rscript $R --stage=setup
#   Rscript $R --stage=fit --methods=zero,balance [--workers=3] [--task=N]
#   Rscript $R --stage=count --methods=zero,balance   (number of tasks)
#   Rscript $R --stage=score
#   Rscript $R --stage=report --candidates=a,b
#
# --split=heldout --unseal-heldout runs any stage on the held-out landscapes
# (results/devbench_heldout), for the final confirmation only.
#
# The acceptance check of the study (README.md): against `balance`, the
# current default, no core scenario may have an RMSE more than 10% higher at
# distance 0, 1 or 2.

si_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_script_dir, "stages.R"))
for (part in c("methods.R", "report.R")) source(file.path(si_repo, "benchmark", "sparse_two_timepoint", part))

si_devbench_source <- file.path(si_repo, "benchmark", "results", "sparse_two_timepoint")

stt_fit_path <- function(config, dataset_id, method_name) {
  file.path(config$output_dir, "fits", sprintf("%s__%s.rds", dataset_id, method_name))
}

# `split = "heldout"` uses the benchmark's held-out landscapes, reserved for
# the final confirmation, in a directory of its own.
si_devbench_config <- function(split = "dev") {
  output_dir <- file.path(si_default_output, if (split == "dev") "devbench" else "devbench_heldout")
  config <- stt_config(output_dir)
  study <- si_config()
  sparse <- function(build, ...) {
    list(type = "alfak2_sparse_joint", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    balance = sparse("8327282"),
    history_grid = sparse("si1", history = "grid"),
    trend_centre_low = sparse("si1", trend_scale_centre = 0.001),
    trend_exponential = sparse("si1", trend_scale_prior = "exponential"),
    trend_select = sparse("si1", trend_terms = "select"),
    grid_exponential = sparse("si1", history = "grid", trend_scale_prior = "exponential"),
    grid_select = sparse("si1", history = "grid", trend_terms = "select")
  )
  # One fixed history each (days): with `balance` (history = interval) they
  # show how the error depends on the history and which one the marginal
  # likelihood picks (their objectives, saved by si_devbench_sparse_joint()).
  for (h in c(30, 100, 300, 1000, 3000)) config$methods[[paste0("h", h)]] <- sparse("si1", history = h)
  config$methods$hdt <- sparse("si1")
  config$methods$history_profile <- sparse("si1", history = "profile")
  config$methods$history_estimate <- sparse("si2", history = "estimate")
  # The package with the new default (build si4), at its defaults.
  config$methods$default_new <- sparse("si4")
  # Cohort fit of the eight patients of a landscape with the balance prior
  # (the benchmark's alfak2_sparse_joint_cohort_balance), scored on patient 1:
  # history = interval (8327282) and each sample's chosen history (si4).
  cohort <- function(build) {
    list(type = "alfak2_sparse_joint_cohort", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
         members = "landscape", dataset_dir = file.path(output_dir, "datasets"),
         landscapes = config$landscapes[c("landscape_id", "split")], cohort_size = 8L,
         arguments = list(trend = "select", initial_state = "balance"))
  }
  config$methods$cohort_balance <- cohort("8327282")
  config$methods$cohort_balance_new <- cohort("si4")
  config$library_root <- study$output_dir
  config$split <- split
  config
}

si_devbench_setup <- function(config) {
  full <- utils::read.csv(file.path(si_devbench_source, sprintf("dataset_index_%s.csv", config$split)),
                          stringsAsFactors = FALSE)
  index <- full[full$scoreboard, , drop = FALSE]
  for (part in c("datasets", "nodes")) dir.create(file.path(config$output_dir, part), recursive = TRUE,
                                                  showWarnings = FALSE)
  # Datasets of every patient (the cohort members as well); scored nodes of
  # the scoreboard datasets.
  for (part in c("datasets", "nodes")) {
    files <- paste0((if (part == "datasets") full else index)$dataset_id, ".rds")
    ok <- file.copy(file.path(si_devbench_source, part, files), file.path(config$output_dir, part, files),
                    overwrite = FALSE)
    present <- file.exists(file.path(config$output_dir, part, files))
    if (!all(present)) stop("Missing ", part, " for ", sum(!present), " datasets.", call. = FALSE)
  }
  utils::write.csv(index, file.path(config$output_dir, sprintf("dataset_index_%s.csv", config$split)), row.names = FALSE)
  cat("Copied", nrow(index), "datasets and their nodes to", config$output_dir, "\n")
}

# stt_method_alfak2_sparse_joint() of the benchmark, also keeping the fit's
# diagnostics: objective, history (and its profile), trend terms, scales.
si_devbench_sparse_joint <- function(dataset, nodes, method) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  started <- proc.time()[["elapsed"]]
  out <- tryCatch({
    fit <- do.call(alfak2::fit_sparse_joint, c(
      list(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
           min_cn = method$min_cn, max_cn = method$max_cn),
      method$arguments
    ))
    prediction <- stats::predict(fit, nodes$karyotype)
    d <- fit$diagnostics
    list(status = "ok", detail = NA_character_,
         prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                                 fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                                 shell = prediction$shell, stringsAsFactors = FALSE),
         diagnostics = list(objective = d$objective, convergence = d$convergence, history = d$history,
                            history_profile = d$history_profile, trend_terms = fit$trend$terms,
                            trend_profile = d$trend_profile, scales = fit$trend$scales,
                            n_major = d$n_major, n_states = d$n_states))
  }, error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e)))
  out$seconds <- proc.time()[["elapsed"]] - started
  out
}

# Every dataset x method, in a fixed order, for one Slurm array task each.
si_devbench_tasks <- function(config, method_names) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", config$split)), stringsAsFactors = FALSE)
  grid <- expand.grid(method_name = method_names, dataset_id = index$dataset_id, stringsAsFactors = FALSE)
  # Cohort methods exist for the cohort scenario only.
  cohort <- vapply(config$methods[grid$method_name], function(m) identical(m$type, "alfak2_sparse_joint_cohort"),
                   logical(1))
  grid <- grid[!cohort | grepl(paste0("_", config$cohort_scenario, "$"), grid$dataset_id), , drop = FALSE]
  lapply(seq_len(nrow(grid)), function(i) list(dataset_id = grid$dataset_id[i], method_name = grid$method_name[i]))
}

si_devbench_fit <- function(config, workers, method_names, task = NULL) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", config$split)), stringsAsFactors = FALSE)
  builds <- unique(unlist(lapply(config$methods[method_names], function(m) m$build)))
  if (length(builds) > 1L) stop("Methods of one run must share an alfak2 build.", call. = FALSE)
  if (length(builds)) {
    .libPaths(c(si_library(builds, config$library_root), .libPaths()))
    loadNamespace("alfak2")
    cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  }
  tasks <- si_devbench_tasks(config, method_names)
  if (!is.null(task)) {
    if (as.integer(task) > length(tasks)) stop("Task ", task, " beyond the ", length(tasks), " tasks.", call. = FALSE)
    tasks <- tasks[as.integer(task)]
  }
  tasks <- Filter(function(t) !file.exists(stt_fit_path(config, t$dataset_id, t$method_name)), tasks)
  cat("Fits to run:", length(tasks), "\n")
  if (!length(tasks)) return(invisible())
  dir.create(file.path(config$output_dir, "fits"), recursive = TRUE, showWarnings = FALSE)
  si_parallel(tasks, function(task) {
    dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
    nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(task$dataset_id, ".rds")))
    method <- config$methods[[task$method_name]]
    if (identical(method$arguments$history, "grid")) method$arguments$history <- si_history_grid(dataset$interval)
    if (identical(method$arguments$history, "profile")) {
      method$arguments$history <- si_history_profile(dataset$interval)
    }
    set.seed(dataset$seed)
    result <- if (identical(method$type, "alfak2_sparse_joint") && !identical(method$build, "8327282")) {
      si_devbench_sparse_joint(dataset, nodes, method)
    } else stt_run_method(method, dataset, nodes, tempdir())
    result$dataset_id <- task$dataset_id
    result$method <- task$method_name
    path <- stt_fit_path(config, task$dataset_id, task$method_name)
    temporary <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(result, temporary)
    file.rename(temporary, path)
    TRUE
  }, workers)
  invisible()
}

si_devbench_score <- function(config, workers) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", config$split)), stringsAsFactors = FALSE)
  method_names <- names(config$methods)
  scored <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(dataset_id, ".rds")))
    truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
    si_bind(lapply(method_names, function(method_name) {
      path <- stt_fit_path(config, dataset_id, method_name)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      prediction <- if (is.null(fit$prediction)) {
        data.frame(karyotype = character(), fitness = numeric(), stringsAsFactors = FALSE)
      } else fit$prediction
      score <- stt_score_predictions(nodes, truth_rate, prediction, config$scoring$reference_min_count)
      cbind(data.frame(dataset_id = dataset_id, method = method_name, status = fit$status,
                       seconds = fit$seconds, stringsAsFactors = FALSE), score$metrics)
    }))
  }, workers)
  metrics <- merge(index, si_bind(scored), by = "dataset_id")
  path <- file.path(config$output_dir, sprintf("metrics_%s.csv", config$split))
  utils::write.csv(metrics, path, row.names = FALSE)
  cat("Scored rows:", nrow(metrics), " file:", path, "\n")
}

# Per scenario and distance: RMSE and Spearman, landscape-weighted as in the
# benchmark, of each candidate against `balance`; pooled paired differences
# with landscape bootstrap intervals.
si_devbench_report <- function(config, candidates, baseline = "balance") {
  metrics <- utils::read.csv(file.path(config$output_dir, sprintf("metrics_%s.csv", config$split)), stringsAsFactors = FALSE)
  cells <- stt_landscape_scenario_means(metrics)
  options(width = 220)
  report_dir <- file.path(config$output_dir, "report")
  dir.create(report_dir, showWarnings = FALSE)
  role <- unique(metrics[c("scenario", "role")])
  per_fit <- metrics[!duplicated(metrics[c("dataset_id", "method")]), ]
  failures <- stats::aggregate(list(failed = per_fit$status != "ok"), per_fit["method"], sum)
  cat("Failed fits by method:\n")
  print(failures, row.names = FALSE)
  for (candidate in candidates) {
    safety <- merge(stt_safety_cells(cells, candidate, baseline, "rmse"), role, by = "scenario")
    spearman <- stt_safety_cells(cells, candidate, baseline, "spearman")
    safety$spearman_baseline <- spearman$baseline_value[match(paste(safety$scenario, safety$stratum),
                                                              paste(spearman$scenario, spearman$stratum))]
    safety$spearman_candidate <- spearman$candidate_value[match(paste(safety$scenario, safety$stratum),
                                                                paste(spearman$scenario, spearman$stratum))]
    safety <- safety[order(safety$role, safety$scenario, safety$stratum), ]
    core <- safety[safety$role == "core", ]
    cat(sprintf("\n== %s against %s: RMSE ratio by scenario (acceptance: core ratio <= 1.10)\n", candidate, baseline))
    print(safety[c("role", "scenario", "stratum", "baseline_value", "candidate_value", "ratio",
                   "spearman_baseline", "spearman_candidate")], row.names = FALSE, digits = 3)
    cat(sprintf("Core cells above 1.10: %d of %d; largest ratio %.3f\n", sum(core$ratio > 1.10), nrow(core),
                max(core$ratio)))
    utils::write.csv(safety, file.path(report_dir, sprintf("safety_%s.csv", candidate)), row.names = FALSE)
    for (scope in c("core", "stress")) {
      part <- cells[cells$scenario %in% role$scenario[role$role == scope], ]
      pooled <- rbind(stt_paired_comparison(part, candidate, baseline, "rmse"),
                      stt_paired_comparison(part, candidate, baseline, "spearman"),
                      stt_paired_comparison(part, candidate, baseline, "coverage95"))
      cat(sprintf("Pooled over %s scenarios (paired by landscape, 95%% bootstrap interval):\n", scope))
      print(pooled[c("metric", "stratum", "n_landscapes", "baseline_mean", "candidate_mean", "difference",
                     "ci_low", "ci_high")], row.names = FALSE, digits = 3)
    }
  }
}

si_devbench_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    args[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  split <- if (is.null(args$split)) "dev" else args$split
  if (!split %in% c("dev", "heldout")) stop("Unknown split: ", split, call. = FALSE)
  if (split == "heldout" && !isTRUE(args$unseal_heldout)) {
    stop("Held-out landscapes are reserved for the final confirmation; pass --unseal-heldout.", call. = FALSE)
  }
  config <- si_devbench_config(split)
  workers <- if (is.null(args$workers)) 3L else as.integer(args$workers)
  split_names <- function(x) if (is.null(x)) character() else strsplit(x, ",", fixed = TRUE)[[1L]]
  method_names <- split_names(args$methods)
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  switch(if (is.null(args$stage)) "" else args$stage,
    setup = si_devbench_setup(config),
    fit = si_devbench_fit(config, workers, method_names, args$task),
    count = cat(length(si_devbench_tasks(config, method_names)), "\n"),
    score = si_devbench_score(config, workers),
    report = si_devbench_report(config, split_names(args$candidates)),
    stop("--stage must be setup, fit, count, score or report.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) si_devbench_main()
