#!/usr/bin/env Rscript
# Unobserved-karyotype bias study (README.md).
#
#   R=validation/sparse_joint/unobserved_bias/run.R
#   Rscript $R --source=rugged|diagnosis --stage=link
#   Rscript $R --source=rugged|diagnosis --stage=fit --methods=zero,base [--datasets=REGEX] [--workers=6]
#   Rscript $R --source=rugged|diagnosis --stage=score [--workers=6]
#
# This file only dispatches, so it can stay unchanged while a run is in
# progress.

ub_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(ub_script_dir, "setup.R"))

ub_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- ub_config(if (is.null(args$source)) "rugged" else args$source)
  workers <- if (is.null(args$workers)) 3L else as.integer(args$workers)
  method_names <- if (is.null(args$methods)) character() else strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  switch(if (is.null(args$stage)) "" else args$stage,
    link = ub_link(config),
    fit = si_stage_fit(config, workers, method_names, args$datasets),
    score = ub_stage_score(config, workers),
    stop("--stage must be link, fit or score.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) ub_main()
