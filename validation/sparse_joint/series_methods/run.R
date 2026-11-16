#!/usr/bin/env Rscript
# The three multi-timepoint methods on simulated sets e1-e3 (README.md).
#
#   R=validation/sparse_joint/series_methods/run.R
#   Rscript $R --stage=truth    --set=e4 [--landscapes=REGEX] [--workers=3]
#   Rscript $R --stage=datasets --set=e2|e3|e4 [--workers=6]
#   Rscript $R --stage=sweeps   --set=e2 --build=BUILD     (e2's gate)
#   Rscript $R --stage=fit      --set=SET --build=BUILD --methods=ends,anchored [--datasets=REGEX] [--workers=6]
#   Rscript $R --stage=compose  --set=SET --build=BUILD   (ends_auto, anchored_auto)
#   Rscript $R --stage=score    --set=SET --build=BUILD [--workers=6]
#
# This file only dispatches, so it can stay unchanged while a run is in
# progress.

sm_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(sm_script_dir, "setup.R"))

sm_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  set <- if (is.null(args$set)) stop("--set is required.", call. = FALSE) else args$set
  stage <- if (is.null(args$stage)) "" else args$stage
  workers <- if (is.null(args$workers)) 6L else as.integer(args$workers)
  if (stage == "datasets") return(sm_stage_datasets(sm_config(set), workers))
  if (stage == "truth") return(sm_stage_truth(sm_config(set), workers, args$landscapes))
  if (is.null(args$build)) stop("--build is required for stage ", stage, ".", call. = FALSE)
  config <- sm_config(set, args$build)
  method_names <- if (is.null(args$methods)) character() else strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  switch(stage,
    sweeps = {
      sm_link(config)
      .libPaths(c(si_library(config$build, config$library_root), .libPaths()))
      loadNamespace("alfak2")
      table <- sm_sweep_table(config)
      utils::write.csv(table, file.path(config$output_dir, "sweeps.csv"), row.names = FALSE)
      multi <- table[table$samples > 2L, ]
      print(stats::aggregate(cbind(any_no_major = multi$middle_no_major > 0,
                                   any_no_domain = multi$middle_no_domain > 0,
                                   end_no_domain = multi$ends_no_domain > 0) ~ start + bottleneck + samples,
                             multi, mean), digits = 2)
    },
    fit = sm_stage_fit(config, workers, method_names, args$datasets),
    score = sm_stage_score(config, workers),
    compose = sm_stage_compose(config),
    stop("--stage must be datasets, sweeps, fit, compose or score.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) sm_main()
