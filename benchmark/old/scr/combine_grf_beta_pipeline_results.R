#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Combine beta-calibrated alfak2/alfakR GRF benchmark outputs.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/combine_grf_beta_pipeline_results.R \\\n",
    "    --pipeline-dir=benchmark/results/grf_alfak2_beta_calibrated_pipeline\n\n",
    "Options:\n",
    "  --pipeline-dir=PATH\n",
    "  --output-dir=PATH       # default: <pipeline-dir>/combined\n",
    "  --include-node-accuracy=true\n",
    "  --path-map=/old/root=/new/root  # optional comma-separated runtime path remapping\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) stop("Unexpected positional argument: ", arg, call. = FALSE)
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      val <- "true"
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- val
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || !length(value) || !nzchar(as.character(value[[1L]]))) default else value[[1L]]
}

arg_logical <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  value <- tolower(trimws(as.character(value)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Expected a boolean for --", gsub("_", "-", name), call. = FALSE)
}

parse_path_map <- function(args) {
  value <- arg_value(args, "path_map", Sys.getenv("ALFAK_BENCH_PATH_MAP", ""))
  if (is.null(value) || !nzchar(as.character(value))) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  entries <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  entries <- entries[nzchar(entries)]
  if (!length(entries)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(entries, function(entry) {
    eq <- regexpr("=", entry, fixed = TRUE)[[1L]]
    if (eq < 2L) stop("Invalid --path-map entry: ", entry, call. = FALSE)
    data.frame(
      from = substr(entry, 1L, eq - 1L),
      to = substr(entry, eq + 1L, nchar(entry)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(nchar(out$from), decreasing = TRUE), , drop = FALSE]
}

apply_path_map_chr <- function(x, path_map) {
  if (is.null(path_map) || !nrow(path_map) || !length(x)) return(x)
  out <- as.character(x)
  na <- is.na(out)
  out[na] <- ""
  for (i in seq_len(nrow(path_map))) {
    from <- path_map$from[[i]]
    to <- path_map$to[[i]]
    hit <- startsWith(out, from)
    out[hit] <- paste0(to, substring(out[hit], nchar(from) + 1L))
  }
  out[na] <- NA_character_
  out
}

read_tsv <- function(path) {
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

normalize_path <- function(path, base = getwd()) {
  if (grepl("^/", path)) normalizePath(path, winslash = "/", mustWork = FALSE) else
    normalizePath(file.path(base, path), winslash = "/", mustWork = FALSE)
}

read_optional_tsv <- function(path, meta) {
  if (!file.exists(path)) return(NULL)
  x <- read_tsv(path)
  if (!nrow(x)) return(x)
  for (nm in names(meta)) {
    if (!nm %in% names(x)) x[[nm]] <- meta[[nm]]
  }
  x
}

rbind_fill <- function(xs) {
  xs <- Filter(Negate(is.null), xs)
  if (!length(xs)) return(data.frame())
  nms <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(nms, names(x))
    for (nm in missing) x[[nm]] <- rep(NA, nrow(x))
    x[, nms, drop = FALSE]
  })
  do.call(rbind, xs)
}

discover_manifest <- function(pipeline_dir) {
  bench_root <- file.path(pipeline_dir, "benchmark")
  beta_dirs <- list.dirs(bench_root, full.names = TRUE, recursive = FALSE)
  beta_dirs <- beta_dirs[grepl("^beta_", basename(beta_dirs))]
  data.frame(
    pm = NA_real_,
    fit_beta_label = sub("^beta_", "", basename(beta_dirs)),
    calibration_dir = file.path(pipeline_dir, "calibration", basename(beta_dirs)),
    benchmark_dir = beta_dirs,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  pipeline_dir <- normalize_path(arg_value(args, "pipeline_dir", "benchmark/results/grf_alfak2_beta_calibrated_pipeline"))
  output_dir <- normalize_path(arg_value(args, "output_dir", file.path(pipeline_dir, "combined")))
  tables_dir <- file.path(output_dir, "tables")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_path <- file.path(pipeline_dir, "tables", "beta_manifest.tsv")
  manifest <- if (file.exists(manifest_path)) read_tsv(manifest_path) else discover_manifest(pipeline_dir)
  if (!nrow(manifest)) stop("No beta manifest rows or benchmark beta directories found.", call. = FALSE)

  path_map <- parse_path_map(args)
  for (nm in c("calibration_dir", "benchmark_dir")) {
    manifest[[nm]] <- vapply(manifest[[nm]], normalize_path, character(1), base = pipeline_dir)
    manifest[[nm]] <- apply_path_map_chr(manifest[[nm]], path_map)
  }
  if (!"fit_beta_label" %in% names(manifest)) {
    manifest$fit_beta_label <- sub("^beta_", "", basename(manifest$benchmark_dir))
  }

  best_rows <- list()
  table_names <- c(
    "fit_results.tsv",
    "summary_by_lambda_time_minobs_method.tsv",
    "alfak2_delta_vs_alfakR.tsv",
    "missing_fit_tasks.tsv"
  )
  if (arg_logical(args, "include_node_accuracy", TRUE)) {
    table_names <- append(table_names, "node_accuracy.tsv", after = 1L)
  }
  bench_rows <- setNames(vector("list", length(table_names)), table_names)
  for (nm in table_names) bench_rows[[nm]] <- list()

  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    meta <- list(
      pipeline_pm = if ("pm" %in% names(row)) row$pm[[1L]] else NA_real_,
      pipeline_beta_label = row$fit_beta_label[[1L]],
      calibration_dir = row$calibration_dir[[1L]],
      benchmark_dir = row$benchmark_dir[[1L]]
    )
    best_rows[[i]] <- read_optional_tsv(
      file.path(row$calibration_dir[[1L]], "tables", "best_params.tsv"),
      meta
    )
    for (nm in table_names) {
      bench_rows[[nm]][[i]] <- read_optional_tsv(
        file.path(row$benchmark_dir[[1L]], "tables", nm),
        meta
      )
    }
  }

  best <- rbind_fill(best_rows)
  write_tsv(manifest, file.path(tables_dir, "beta_manifest.tsv"))
  write_tsv(best, file.path(tables_dir, "best_params_by_pm.tsv"))
  for (nm in table_names) {
    out <- rbind_fill(bench_rows[[nm]])
    write_tsv(out, file.path(tables_dir, nm))
  }
  message("Wrote combined beta-calibrated outputs under: ", tables_dir)
  invisible(list(manifest = manifest, best = best))
}

if (identical(environment(), globalenv())) main()
