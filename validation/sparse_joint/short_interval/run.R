#!/usr/bin/env Rscript
# Short-interval diagnosis of fit_sparse_joint() (README.md).
#
#   R=validation/sparse_joint/short_interval/run.R
#   Rscript $R --stage=truth    [--workers=3]
#   Rscript $R --stage=datasets [--workers=3]
#   Rscript $R --stage=fit      --methods=zero,balance [--datasets=REGEX] [--workers=3]
#   Rscript $R --stage=score    [--workers=3]
#
# Options: --output-dir=DIR (default: results/ beside this script);
# --set=heldout --unseal-heldout for the held-out landscapes (seed base 94000,
# results/heldout_short), reserved for the final confirmation.
# The stage code lives in setup.R and stages.R; this file only dispatches, so
# it can stay unchanged while a run is in progress.

si_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_script_dir, "stages.R"))

si_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  set <- if (is.null(args$set)) "diagnosis" else args$set
  if (set == "heldout" && !isTRUE(args$unseal_heldout)) {
    stop("The held-out short-interval set is for the final confirmation; pass --unseal-heldout.", call. = FALSE)
  }
  config <- si_config(args$output_dir, set)
  workers <- if (is.null(args$workers)) 3L else as.integer(args$workers)
  method_names <- if (is.null(args$methods)) character() else strsplit(args$methods, ",", fixed = TRUE)[[1L]]
  unknown <- setdiff(method_names, names(config$methods))
  if (length(unknown)) stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  switch(if (is.null(args$stage)) "" else args$stage,
    truth = si_stage_truth(config, workers),
    datasets = si_stage_datasets(config, workers),
    fit = si_stage_fit(config, workers, method_names, args$datasets),
    score = si_stage_score(config, workers),
    stop("--stage must be truth, datasets, fit or score.", call. = FALSE)
  )
}

if (sys.nframe() == 0L) si_main()
