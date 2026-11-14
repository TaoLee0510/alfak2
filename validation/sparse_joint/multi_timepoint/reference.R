#!/usr/bin/env Rscript
# Reference fits of an alfak2 build, and their comparison between builds
# (stage 0 and the regression check of stage 1 of the multi-timepoint work;
# README.md in this directory).
#
#   R=validation/sparse_joint/multi_timepoint/reference.R
#   Rscript $R --stage=fit --build=ff13cf4 [--workers=6] [--parts=REGEX] [--task=N]
#   Rscript $R --stage=count [--parts=REGEX]
#   Rscript $R --stage=compare --builds=ff13cf4,mt1
#
# `--parts` keeps the parts whose name matches; `--task=N` fits only the Nth
# item of that list (one task of a job array, see hpc/), `count` prints its
# length.
#
# A build is the library results/lib/alfak2_<build>. Fits go to
# results/reference/<build>/ for simulated data and to
# ~/Documents/alfak2_private/multi_timepoint/reference/<build>/ for the GBM
# inputs, which are patient data and stay outside the repository. The GBM
# inputs and the development benchmark's datasets are read in place.
#
# What is fitted, all two-timepoint:
#   dev      fit_sparse_joint() at its defaults (birth rate 1, copy numbers
#            1-8) on patient 1 of the 12 development landscapes in the
#            benchmark's 8 core and 2 short-interval scenarios;
#   entry    fit_alfak2() on two of them;
#   free     initial_state = "free" on twelve of them;
#   rates    fit_sparse_joint_rates() over 0.25-4 on two of them;
#   cohort   fit_sparse_joint_cohort() of the 8 patients of two landscapes
#            (typical scenario), free and balance;
#   compat   as_alfakR() and alfak_compat() on four of them;
#   gbm      fit_sparse_joint() on the 20 GBM inputs at birth rate
#            sqrt(0.02 * 1) per day, the centre of the range real_cohort.R used.
# Each item keeps the fit's numbers and its predictions at fixed karyotypes,
# computed by the build that made the fit.

if (!exists("mt_script_dir")) mt_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
mt_repo <- normalizePath(file.path(mt_script_dir, "..", "..", ".."))
for (part in c("landscape.R", "config.R", "score.R")) {
  source(file.path(mt_repo, "benchmark", "sparse_two_timepoint", part))
}
mt_results <- file.path(mt_script_dir, "results")
mt_benchmark <- file.path(mt_repo, "benchmark", "results", "sparse_two_timepoint")
mt_gbm_inputs <- "/Users/4482173/Documents/GitHub/PANcanKFLs/Results/data/processed/ALFA_K/GBM/ALFA-K_inputs"
mt_private <- file.path(path.expand("~"), "Documents", "alfak2_private", "multi_timepoint", "reference")

mt_library <- function(build) file.path(mt_results, "lib", paste0("alfak2_", build))

mt_reference_dir <- function(build, part) {
  if (part == "gbm") file.path(mt_private, build, "gbm") else file.path(mt_results, "reference", build, part)
}

mt_scenarios <- c("very_sparse_50_50_180d", "sparse_200_200_180d", "sparse_200_200_360d", "typical_500_700_360d",
                  "unbalanced_100_3000_360d", "unbalanced_3000_100_360d", "long_500_700_1080d",
                  "bottleneck_500_700_360d", "short_500_700_30d", "short_500_700_60d")

mt_dev_ids <- function() {
  config <- stt_config(mt_benchmark)
  landscapes <- config$landscapes$landscape_id[config$landscapes$split == "dev"]
  as.vector(outer(sprintf("%s_P01", landscapes), mt_scenarios, paste, sep = "_"))
}

mt_dataset <- function(id) readRDS(file.path(mt_benchmark, "datasets", paste0(id, ".rds")))
mt_nodes <- function(id) readRDS(file.path(mt_benchmark, "nodes", paste0(id, ".rds")))
mt_counts <- function(dataset) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  counts
}

# Everything numeric a sparse joint fit holds that a change could move.
mt_fit_record <- function(fit, karyotypes) {
  d <- fit$diagnostics
  list(summary = fit$summary, posterior = fit$posterior,
       trend = fit$trend[c("coefficients", "covariance", "centre", "scales", "observed_level", "reference")],
       diagnostics = d[intersect(names(d), c("convergence", "message", "optimisation", "objective", "n_states",
                                              "n_major", "n_edges",
                                              "dropped_karyotypes", "dropped_cells", "history", "history_profile",
                                              "domain_coverage", "fitted_birth_rate", "finite_uncertainty"))],
       prediction = if (length(karyotypes)) stats::predict(fit, karyotypes) else NULL)
}

mt_cohort_record <- function(fit, karyotypes) {
  d <- fit$diagnostics
  list(landscape = fit$landscape, samples = fit$samples,
       diagnostics = d[intersect(names(d), c("convergence", "message", "optimisation", "objective", "trend",
                                              "trend_objective", "trend_convergence", "n_states",
                                              "n_nodes", "shared_nodes", "domain_coverage", "history_profiles"))],
       prediction = stats::predict(fit, karyotypes, sample = 1L))
}

mt_gbm_inputs_list <- function() {
  files <- list.files(mt_gbm_inputs, pattern = "[.][Rr]ds$", full.names = TRUE)
  stats::setNames(files, sub("[.][Rr]ds$", "", basename(files)))
}

# The tasks of the fit stage: one per item, each writing one file.
mt_tasks <- function() {
  dev <- mt_dev_ids()
  typical <- grep("typical", dev, value = TRUE)
  tasks <- c(
    lapply(dev, function(id) list(part = "dev", id = id)),
    lapply(typical[c(1, 7)], function(id) list(part = "entry", id = id)),
    lapply(typical, function(id) list(part = "free", id = id)),
    lapply(typical[c(4, 10)], function(id) list(part = "rates", id = id)),
    lapply(c("dev_L01", "dev_L07"), function(l) list(part = "cohort", id = paste0(l, "_free"), landscape = l,
                                                     initial_state = "free")),
    lapply(c("dev_L01", "dev_L07"), function(l) list(part = "cohort", id = paste0(l, "_balance"), landscape = l,
                                                     initial_state = "balance")),
    lapply(typical[c(1, 4, 7, 10)], function(id) list(part = "compat", id = id)),
    lapply(names(mt_gbm_inputs_list()), function(id) list(part = "gbm", id = id))
  )
  tasks
}

mt_run_task <- function(task) {
  fit_dev <- function(dataset, ...) {
    alfak2::fit_sparse_joint(mt_counts(dataset), dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                             min_cn = 1L, max_cn = 8L, ...)
  }
  switch(task$part,
    dev = , free = {
      dataset <- mt_dataset(task$id)
      fit <- if (task$part == "dev") fit_dev(dataset) else fit_dev(dataset, initial_state = "free")
      mt_fit_record(fit, mt_nodes(task$id)$karyotype)
    },
    entry = {
      dataset <- mt_dataset(task$id)
      fit <- alfak2::fit_alfak2(mt_counts(dataset), dt = dataset$interval, beta = dataset$beta,
                                sparse_joint_control = list(birth_rate = 1, min_cn = 1L, max_cn = 8L))
      mt_fit_record(fit, mt_nodes(task$id)$karyotype)
    },
    rates = {
      dataset <- mt_dataset(task$id)
      fit <- alfak2::fit_sparse_joint_rates(mt_counts(dataset), dt = dataset$interval, beta = dataset$beta,
                                            min_cn = 1L, max_cn = 8L, birth_rate_range = c(0.25, 4))
      prediction <- stats::predict(fit, mt_nodes(task$id)$karyotype)
      list(rates = fit$rates, prediction = prediction, rank_agreement = attr(prediction, "rank_agreement"),
           fits = lapply(fit$fits, function(f) mt_fit_record(f, character())))
    },
    cohort = {
      ids <- sprintf("%s_P%02d_typical_500_700_360d", task$landscape, 1:8)
      members <- lapply(ids, mt_dataset)
      samples <- stats::setNames(lapply(members, mt_counts), ids)
      fit <- alfak2::fit_sparse_joint_cohort(samples, dt = vapply(members, `[[`, numeric(1), "interval"),
                                             beta = members[[1L]]$beta, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                                             initial_state = task$initial_state)
      mt_cohort_record(fit, mt_nodes(ids[[1L]])$karyotype)
    },
    compat = {
      dataset <- mt_dataset(task$id)
      fit <- fit_dev(dataset)
      set.seed(1)
      compat <- alfak2::as_alfakR(fit, minobs = 10, nboot = 20, n0 = 5e4, nb = 2e6, correct_efflux = TRUE)
      yi <- list(x = mt_counts(dataset), dt = 1)
      colnames(yi$x) <- c("0", as.character(dataset$interval))
      outdir <- tempfile("alfak_compat_")
      set.seed(1)
      r2r <- alfak2::alfak_compat(yi, outdir, minobs = 10, nboot = 20, n0 = 5e4, nb = 2e6, pm = dataset$beta,
                                  correct_efflux = TRUE, birth_rate = 1, min_cn = 1L, max_cn = 8L)
      files <- c("bootstrap_res.Rds", "landscape.Rds", "landscape_posterior_samples.Rds", "xval.Rds")
      written <- stats::setNames(lapply(files, function(f) readRDS(file.path(outdir, f))), files)
      unlink(outdir, recursive = TRUE)
      list(as_alfakR = compat[c("bootstrap_res", "landscape", "landscape_posterior_samples", "xval")],
           alfak_compat = list(r2r = r2r, files = written))
    },
    gbm = {
      input <- readRDS(mt_gbm_inputs_list()[[task$id]])
      counts <- as.matrix(input$x)
      storage.mode(counts) <- "integer"
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      interval <- diff(as.numeric(colnames(counts))) * (if (is.null(input$dt)) 1 else input$dt)
      cn <- alfak2:::parse_karyotypes(rownames(counts))
      min_cn <- min(1L, min(cn))
      max_cn <- max(cn) + 1L
      fit <- alfak2::fit_sparse_joint(counts, dt = interval, beta = 5e-5, birth_rate = sqrt(0.02), min_cn = min_cn,
                                      max_cn = max_cn)
      nodes <- stt_scoring_nodes(counts, min_cn, max_cn, d1_cap = 500L, d2_cap = 500L, seed = 1L)
      mt_fit_record(fit, nodes$karyotype)
    },
    stop("Unknown part ", task$part, call. = FALSE)
  )
}

# The items, in a fixed order, of the parts matching `parts`.
mt_selected_tasks <- function(parts = NULL) {
  tasks <- mt_tasks()
  if (!is.null(parts)) tasks <- Filter(function(task) grepl(parts, task$part), tasks)
  tasks
}

mt_stage_fit <- function(build, workers, parts = NULL, task = NULL) {
  library_dir <- mt_library(build)
  if (!dir.exists(file.path(library_dir, "alfak2"))) stop("No alfak2 build in ", library_dir, call. = FALSE)
  .libPaths(c(library_dir, .libPaths()))
  loadNamespace("alfak2")
  cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  tasks <- mt_selected_tasks(parts)
  if (!is.null(task)) {
    task <- as.integer(task)
    if (is.na(task) || task < 1L || task > length(tasks)) stop("--task must be 1 to ", length(tasks), call. = FALSE)
    tasks <- tasks[task]
    cat("Task", task, ":", tasks[[1L]]$part, tasks[[1L]]$id, "\n")
  }
  path_of <- function(task) file.path(mt_reference_dir(build, task$part), paste0(task$id, ".rds"))
  tasks <- Filter(function(task) !file.exists(path_of(task)), tasks)
  cat("Items to fit:", length(tasks), "\n")
  results <- parallel::mclapply(tasks, function(task) {
    started <- proc.time()[["elapsed"]]
    record <- tryCatch(mt_run_task(task), error = function(e) list(error = conditionMessage(e),
                                                                   error_class = class(e)[[1L]]))
    record$seconds <- proc.time()[["elapsed"]] - started
    path <- path_of(task)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(record, temporary)
    file.rename(temporary, path)
    if (is.null(record$error)) "ok" else paste(task$part, task$id, record$error)
  }, mc.cores = workers, mc.preschedule = FALSE)
  failed <- unlist(Filter(function(x) !identical(x, "ok"), results))
  cat("Done.", length(failed), "item(s) recorded an error.\n")
  if (length(failed)) cat(utils::head(failed, 10), sep = "\n")
}

# Largest difference between two records, relative to the scale of each
# numeric leaf (its largest absolute value), and whether every non-numeric
# leaf is identical. Leaves are matched by their path in the record.
mt_leaves <- function(x, path = "") {
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x) && !is.null(names(x)) || (is.list(x) && length(x))) {
    keys <- if (is.null(names(x))) as.character(seq_along(x)) else names(x)
    return(unlist(lapply(seq_along(x), function(i) mt_leaves(x[[i]], paste0(path, "/", keys[[i]]))),
                  recursive = FALSE))
  }
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  stats::setNames(list(x), path)
}

mt_compare_records <- function(a, b, ignore = "/seconds$") {
  la <- mt_leaves(a)
  lb <- mt_leaves(b)
  la <- la[!grepl(ignore, names(la))]
  lb <- lb[!grepl(ignore, names(lb))]
  missing <- union(setdiff(names(la), names(lb)), setdiff(names(lb), names(la)))
  common <- intersect(names(la), names(lb))
  worst <- 0
  worst_at <- NA_character_
  mismatched <- character()
  for (key in common) {
    x <- la[[key]]
    y <- lb[[key]]
    if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      same_na <- identical(is.na(x), is.na(y))
      if (!same_na) {
        mismatched <- c(mismatched, key)
        next
      }
      ok <- !is.na(x)
      if (!any(ok)) next
      infinite <- !is.finite(x[ok]) | !is.finite(y[ok])
      if (any(infinite) && !identical(x[ok][infinite], y[ok][infinite])) {
        mismatched <- c(mismatched, key)
        next
      }
      xf <- x[ok][!infinite]
      yf <- y[ok][!infinite]
      if (!length(xf)) next
      scale <- max(abs(c(xf, yf)), 1e-12)
      relative <- max(abs(xf - yf)) / scale
      if (relative > worst) {
        worst <- relative
        worst_at <- key
      }
    } else if (!identical(x, y)) {
      mismatched <- c(mismatched, key)
    }
  }
  list(worst = worst, worst_at = worst_at, mismatched = mismatched, missing = missing)
}

mt_stage_compare <- function(builds, tolerance = 1e-8) {
  stopifnot(length(builds) == 2L)
  rows <- list()
  for (part in c("dev", "entry", "free", "rates", "cohort", "compat", "gbm")) {
    dirs <- vapply(builds, mt_reference_dir, character(1), part = part)
    files <- intersect(list.files(dirs[[1L]]), list.files(dirs[[2L]]))
    for (file in files) {
      result <- mt_compare_records(readRDS(file.path(dirs[[1L]], file)), readRDS(file.path(dirs[[2L]], file)))
      rows[[length(rows) + 1L]] <- data.frame(
        part = part, item = sub("[.]rds$", "", file), worst_relative = result$worst, worst_at = result$worst_at,
        mismatched = paste(utils::head(result$mismatched, 3), collapse = " "),
        missing = paste(utils::head(result$missing, 3), collapse = " "), stringsAsFactors = FALSE
      )
    }
  }
  table <- do.call(rbind, rows)
  table$pass <- table$worst_relative <= tolerance & !nzchar(table$mismatched) & !nzchar(table$missing)
  out <- file.path(mt_results, "reference", sprintf("compare_%s_%s.csv", builds[[1L]], builds[[2L]]))
  utils::write.csv(table, out, row.names = FALSE)
  print(stats::aggregate(cbind(items = 1, pass = pass, worst_relative = worst_relative) ~ part, table,
                         function(x) if (all(x %in% c(0, 1)) && length(unique(x)) <= 2) sum(x) else max(x)))
  cat("Items:", nrow(table), " passing:", sum(table$pass), " largest relative difference:",
      format(max(table$worst_relative)), "\n", "Table:", out, "\n")
  invisible(table)
}

mt_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  workers <- if (is.null(args$workers)) 6L else as.integer(args$workers)
  switch(if (is.null(args$stage)) "" else args$stage,
    fit = mt_stage_fit(args$build, workers, args$parts, args$task),
    count = cat(length(mt_selected_tasks(args$parts)), "\n"),
    compare = mt_stage_compare(strsplit(args$builds, ",", fixed = TRUE)[[1L]]),
    stop("--stage must be fit, count or compare.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) mt_main()
