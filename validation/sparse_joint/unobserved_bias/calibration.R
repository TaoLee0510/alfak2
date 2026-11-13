#!/usr/bin/env Rscript
# Calibration of the reported uncertainty of fit_sparse_joint() across the
# regimes studied so far (README.md in this directory).
#
# The shipped constants (floor 0.009, deficit 0.6) were fitted to the core
# scenarios of the development benchmark. Here the default fits (af362af, or
# its identical predecessor) of three sets are pooled:
#
#   dev        the development benchmark's development landscapes, copied in
#              ../short_interval/results/devbench (method history_profile)
#   diagnosis  dev-sized populations at 6-42 day intervals, wavelengths
#              0.2-1.6 (../short_interval, method history_profile)
#   rugged     paper-sized young populations, wavelengths 0.2-0.8
#              (../lattice_residual, method base)
#
#   Rscript validation/sparse_joint/unobserved_bias/calibration.R collect [--workers=6]
#
# `collect` writes one row per scored karyotype to
# results/calibration_nodes.rds; the forms are fitted and compared by
# calibration_forms.R.

if (!exists("ub_script_dir")) ub_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(ub_script_dir, "setup.R"))

ub_calibration_sources <- function() {
  si <- si_default_output
  list(
    dev = list(dir = file.path(si, "devbench"), index = "dataset_index_dev.csv",
               fit = function(id) file.path(si, "devbench", "fits", sprintf("%s__history_profile.rds", id))),
    diagnosis = list(dir = si, index = "dataset_index.csv",
                     fit = function(id) file.path(si, "fits", "history_profile", paste0(id, ".rds"))),
    rugged = list(dir = file.path(lr_default_output, "rugged"), index = "dataset_index.csv",
                  fit = function(id) file.path(lr_default_output, "rugged", "fits", "base", paste0(id, ".rds")))
  )
}

# Scored karyotypes of one dataset with the fit's prediction, its model sd and
# shell, and the truth, centred as the development benchmark centres them.
ub_calibration_nodes <- function(source, name, id, reference_min_count = 5) {
  fit <- readRDS(source$fit(id))
  if (!identical(fit$status, "ok")) return(NULL)
  dataset <- readRDS(file.path(source$dir, "datasets", paste0(id, ".rds")))
  scored <- readRDS(file.path(source$dir, "nodes", paste0(id, ".rds")))
  truth <- stt_true_growth_rate(scored$karyotype, dataset$truth)
  table <- stt_score_predictions(scored, truth, fit$prediction, reference_min_count)$nodes
  at <- match(table$karyotype, fit$prediction$karyotype)
  data.frame(source = name, dataset_id = id, distance = table$distance,
             total = table$count_t0 + table$count_t1, shell = fit$prediction$shell[at],
             pred = table$pred_relative, truth = table$truth_relative,
             model_sd = fit$prediction$model_sd[at], fitness_sd = fit$prediction$fitness_sd[at],
             stringsAsFactors = FALSE)
}

ub_calibration_collect <- function(workers) {
  sources <- ub_calibration_sources()
  out <- lapply(names(sources), function(name) {
    source <- sources[[name]]
    index <- utils::read.csv(file.path(source$dir, source$index), stringsAsFactors = FALSE)
    index <- index[file.exists(vapply(index$dataset_id, source$fit, character(1))), ]
    rows <- si_parallel(index$dataset_id, function(id) ub_calibration_nodes(source, name, id), workers)
    nodes <- do.call(rbind, rows)
    info <- lapply(index$dataset_id, function(id) {
      d <- readRDS(file.path(source$dir, "datasets", paste0(id, ".rds")))
      data.frame(dataset_id = id, wavelength_class = d$wavelength_class, interval = d$interval,
                 t0 = if (is.null(d$t0)) NA else d$t0, scenario = d$scenario, stringsAsFactors = FALSE)
    })
    info <- do.call(rbind, info)
    info$role <- if (!is.null(index$role)) index$role[match(info$dataset_id, index$dataset_id)] else name
    merge(nodes, info)
  })
  nodes <- do.call(rbind, out)
  nodes <- nodes[is.finite(nodes$model_sd) & is.finite(nodes$pred), ]
  nodes$kind <- ifelse(nodes$distance == 0L, ifelse(nodes$total >= 5, "observed_5plus", "observed_lt5"),
                       paste0("unobserved_d", nodes$distance))
  path <- file.path(ub_default_output, "calibration_nodes.rds")
  saveRDS(nodes, path)
  cat("Collected", nrow(nodes), "karyotypes of", length(unique(nodes$dataset_id)), "datasets in", path, "\n")
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workers <- as.integer(sub("^--workers=", "", c(grep("^--workers=", args, value = TRUE), "--workers=6")[[1L]]))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  switch(if (length(args)) args[[1L]] else "",
    collect = ub_calibration_collect(workers),
    stop("Usage: calibration.R collect [--workers=6]", call. = FALSE)
  )
}
