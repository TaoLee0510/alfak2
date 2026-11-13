#!/usr/bin/env Rscript
# Cohort model at a short interval: the eight patients of a diagnosis
# landscape (wavelengths 0.8 and 1.6, first sample at day 120, 6 days apart,
# 1000 cells each), fitted together with the balance prior and scored on
# patient 1 against the truth, with history = interval (8327282) and with each
# sample's history chosen as fit_sparse_joint() chooses it (si4).
#
#   R=validation/sparse_joint/short_interval/cohort_short.R
#   Rscript $R --stage=datasets
#   Rscript $R --stage=fit --method=old|new --task=N      (N: landscape 1-10)
#   Rscript $R --stage=score
#
# Patient 1 is the diagnosis dataset; patients 2-8 are new trajectories on the
# same landscapes. Output: results/cohort_short.

si_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(si_script_dir, "setup.R"))

si_cohort_config <- function() {
  config <- si_config()
  config$diagnosis_dir <- config$output_dir
  config$output_dir <- file.path(si_default_output, "cohort_short")
  config$landscapes <- config$landscapes[config$landscapes$wavelength_class %in% c("w0.8", "w1.6"), ]
  config$scenario <- config$scenarios[config$scenarios$scenario == "t0120_06d", ]
  config$scenario_index <- which(config$scenarios$scenario == "t0120_06d")
  config$cohort_size <- 8L
  config$builds <- c(old = "8327282", new = "si4")
  config
}

si_cohort_member_path <- function(config, landscape_id, patient) {
  file.path(config$output_dir, "datasets", sprintf("%s_P%02d_%s.rds", landscape_id, patient, config$scenario$scenario))
}

si_cohort_datasets <- function(config) {
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  dir.create(file.path(config$output_dir, "datasets"), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(config$landscapes))) {
    landscape <- config$landscapes[i, ]
    first <- sprintf("%s_P01_%s.rds", landscape$landscape_id, config$scenario$scenario)
    file.copy(file.path(config$diagnosis_dir, "datasets", first), si_cohort_member_path(config, landscape$landscape_id, 1L),
              overwrite = FALSE)
    for (patient in 2:config$cohort_size) {
      truth <- readRDS(stt_generate_truth(config, landscape, patient))
      seed <- as.integer(landscape$landscape_seed %% 100000L * 1000L + patient * 20L + config$scenario_index)
      stt_make_dataset(config, truth, config$scenario, config$scenario$t0, seed)
    }
  }
  cat("Cohort datasets:", length(list.files(file.path(config$output_dir, "datasets"))), "\n")
}

si_cohort_fit <- function(config, method, task) {
  landscape <- config$landscapes[as.integer(task), ]
  path <- file.path(config$output_dir, "fits", sprintf("%s__%s.rds", landscape$landscape_id, method))
  if (file.exists(path)) return(invisible())
  .libPaths(c(si_library(config$builds[[method]]), .libPaths()))
  members <- lapply(seq_len(config$cohort_size), function(p) readRDS(si_cohort_member_path(config, landscape$landscape_id, p)))
  samples <- lapply(members, function(m) {
    counts <- m$counts
    storage.mode(counts) <- "integer"
    counts
  })
  names(samples) <- vapply(members, function(m) m$dataset_id, character(1))
  nodes <- readRDS(file.path(config$diagnosis_dir, "nodes", paste0(members[[1L]]$dataset_id, ".rds")))
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    fit <- alfak2::fit_sparse_joint_cohort(samples, dt = vapply(members, function(m) m$interval, numeric(1)),
                                           beta = members[[1L]]$beta, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                                           trend = "select", initial_state = "balance")
    prediction <- stats::predict(fit, nodes$karyotype, sample = 1L)
    list(status = "ok", history = fit$diagnostics$history, trend = fit$diagnostics$trend,
         prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                                 fitness_sd = prediction$relative_sd, stringsAsFactors = FALSE))
  }, error = function(e) list(status = "failed", detail = conditionMessage(e)))
  result$seconds <- proc.time()[["elapsed"]] - started
  result$build <- config$builds[[method]]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, path)
  cat("Saved", path, result$status, "\n")
}

si_cohort_score <- function(config) {
  rows <- list()
  for (landscape_id in config$landscapes$landscape_id) {
    dataset <- readRDS(si_cohort_member_path(config, landscape_id, 1L))
    nodes <- readRDS(file.path(config$diagnosis_dir, "nodes", paste0(dataset$dataset_id, ".rds")))
    truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
    for (method in names(config$builds)) {
      path <- file.path(config$output_dir, "fits", sprintf("%s__%s.rds", landscape_id, method))
      if (!file.exists(path)) next
      fit <- readRDS(path)
      if (!identical(fit$status, "ok")) {
        rows[[length(rows) + 1L]] <- data.frame(landscape_id = landscape_id, method = method, status = fit$status)
        next
      }
      score <- stt_score_predictions(nodes, truth_rate, fit$prediction, config$scoring$reference_min_count)
      m <- score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), ]
      rows[[length(rows) + 1L]] <- data.frame(landscape_id = landscape_id, wavelength_class = dataset$wavelength_class,
                                              method = method, status = "ok", stratum = m$stratum,
                                              rmse_ratio = m$rmse / m$rmse_zero, spearman = m$spearman,
                                              coverage95 = m$coverage95, seconds = fit$seconds,
                                              median_history = if (is.null(fit$history)) NA_real_ else stats::median(fit$history))
    }
  }
  out <- si_bind(rows)
  utils::write.csv(out, file.path(config$output_dir, "metrics.csv"), row.names = FALSE)
  options(width = 200)
  ok <- out[out$status == "ok", ]
  print(stats::aggregate(cbind(rmse_ratio, spearman, coverage95, seconds) ~ method + stratum, ok,
                         function(x) round(stats::median(x), 3)), row.names = FALSE)
  if (any(out$status != "ok")) print(out[out$status != "ok", ], row.names = FALSE)
}

si_bind <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(columns, names(r))] <- NA
    r[columns]
  }))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else NULL
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- si_cohort_config()
  switch(value("stage"),
    datasets = si_cohort_datasets(config),
    fit = si_cohort_fit(config, value("method"), value("task")),
    score = si_cohort_score(config),
    stop("--stage must be datasets, fit or score.", call. = FALSE)
  )
}
