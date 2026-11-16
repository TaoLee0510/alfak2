#!/usr/bin/env Rscript
# Level of the balance prior (README.md).
#
#   R=validation/sparse_joint/cloud_level/run.R
#   Rscript $R --stage=inputs                       copy the dev and short sets
#   Rscript $R --stage=truth    [--workers=3] [--landscapes=REGEX] [--set=diverse|dev2|short2|diverse2]
#   Rscript $R --stage=datasets [--workers=3] [--landscapes=REGEX] [--set=...]
#   Rscript $R --stage=tasks    --methods=base,free,level,level_dt [--sets=dev,short,diverse]
#   Rscript $R --stage=fit      --task=N [--tasks-file=tasks.tsv]   one line of results/<tasks file>
#   Rscript $R --stage=fit      --methods=... [--sets=...] [--datasets=REGEX] [--workers=3]
#   Rscript $R --stage=score    [--workers=3] [--sets=...] [--out=scores.csv]
#
# The stage code lives in setup.R and stages.R; this file only dispatches, so
# it can stay unchanged while a run is in progress.

cl_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(cl_script_dir, "setup.R"))
source(file.path(cl_script_dir, "stages.R"))

cl_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  workers <- if (is.null(args$workers)) 3L else as.integer(args$workers)
  split <- function(x) if (is.null(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1L]]
  methods <- split(args$methods)
  unknown <- setdiff(methods, names(cl_methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  sets <- if (is.null(args$sets)) cl_sets else split(args$sets)
  tasks_file <- if (is.null(args$tasks_file)) "tasks.tsv" else args$tasks_file
  switch(if (is.null(args$stage)) "" else args$stage,
    inputs = cl_stage_inputs(),
    truth = cl_stage_truth(workers, args$landscapes, if (is.null(args$set)) "diverse" else args$set),
    datasets = cl_stage_datasets(workers, args$landscapes, if (is.null(args$set)) "diverse" else args$set),
    tasks = cl_stage_tasks(methods, sets, tasks_file),
    fit = cl_stage_fit(workers, args$task, methods, sets, args$datasets, tasks_file),
    score = cl_stage_score(workers, sets, if (is.null(args$out)) "scores.csv" else args$out),
    stop("--stage must be inputs, truth, datasets, tasks, fit or score.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) cl_main()
