#!/usr/bin/env Rscript
# Patient regrowth study of fit_sparse_joint() at ff13cf4 (README.md).
#
#   R=validation/sparse_joint/patient_regrowth/run.R
#   Rscript $R --stage=link
#   Rscript $R --stage=datasets [--workers=6]
#   Rscript $R --stage=fit --methods=zero,base [--datasets=REGEX] [--workers=6]
#   Rscript $R --stage=score [--workers=6]
#   Rscript $R --stage=symmetry --dataset=dev_L01_P01_s0200_d875
#
# `--build=BUILD` (default ff13cf4) fits with that alfak2 build of the
# multi-timepoint work, into results/BUILD/; run `--stage=link --build=BUILD`
# first, after the datasets.
#
# `equivalent` is meant for the dormant scenarios (--datasets=_d[0-9]+$) and
# `cohort` for the cohort scenarios (--datasets="_(none|s0200|s0020)$").
# This file only dispatches, so it can stay unchanged while a run is in
# progress.

pr_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(pr_script_dir, "setup.R"))

pr_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- pr_config(build = if (is.null(args$build)) "ff13cf4" else args$build)
  workers <- if (is.null(args$workers)) 6L else as.integer(args$workers)
  method_names <- if (is.null(args$methods)) character() else strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  switch(if (is.null(args$stage)) "" else args$stage,
    link = pr_stage_link(config),
    datasets = pr_stage_datasets(config, workers),
    fit = si_stage_fit(config, workers, method_names, args$datasets),
    score = pr_stage_score(config, workers),
    symmetry = pr_stage_symmetry(config, args$dataset),
    stop("--stage must be link, datasets, fit, score or symmetry.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) pr_main()
