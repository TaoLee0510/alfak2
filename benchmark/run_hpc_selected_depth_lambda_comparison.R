#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

cap_threads <- function(x, max_threads = 60L) {
  x <- suppressWarnings(as.integer(x))
  if (!is.finite(x) || x < 1L) x <- 1L
  min(as.integer(x), as.integer(max_threads))
}

hpc_worker_threads <- function(max_threads = 60L) {
  cap_threads(Sys.getenv("SLURM_CPUS_PER_TASK", "1"), max_threads = max_threads)
}

parallel_lapply <- function(x, fun, threads = hpc_worker_threads(), ...) {
  threads <- cap_threads(threads)
  if (threads <= 1L || length(x) <= 1L || .Platform$OS.type == "windows") {
    return(lapply(x, fun, ...))
  }
  parallel::mclapply(
    x,
    fun,
    ...,
    mc.cores = min(threads, length(x)),
    mc.preschedule = FALSE
  )
}

parse_hpc_args <- function(args) {
  out <- list(
    init = FALSE,
    combine = FALSE,
    resume = FALSE,
    submit = FALSE,
    output_dir = "benchmark/results/hpc_22chr_9method_depth_lambda_comparison",
    qos = "xxlarge",
    mem = "8G",
    cpus_per_task = 1L,
    time = "04:00:00",
    summary_mem = "256G",
    summary_cpus_per_task = 60L,
    summary_time = "08:00:00",
    summary_dependency = "afterany"
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--init")) out$init <- TRUE
    else if (identical(arg, "--combine")) out$combine <- TRUE
    else if (identical(arg, "--resume")) out$resume <- TRUE
    else if (identical(arg, "--submit")) out$submit <- TRUE
    else if (identical(arg, "--task-id")) {
      i <- i + 1L
      out$task_id <- as.integer(args[[i]])
    } else if (grepl("^--task-id=", arg)) out$task_id <- as.integer(sub("^--task-id=", "", arg))
    else if (identical(arg, "--output-dir")) {
      i <- i + 1L
      out$output_dir <- args[[i]]
    } else if (grepl("^--output-dir=", arg)) out$output_dir <- sub("^--output-dir=", "", arg)
    else if (identical(arg, "--alfakR-repo")) {
      i <- i + 1L
      out$alfakR_repo <- args[[i]]
    } else if (grepl("^--alfakR-repo=", arg)) out$alfakR_repo <- sub("^--alfakR-repo=", "", arg)
    else if (identical(arg, "--qos")) {
      i <- i + 1L
      out$qos <- args[[i]]
    } else if (grepl("^--qos=", arg)) out$qos <- sub("^--qos=", "", arg)
    else if (identical(arg, "--mem")) {
      i <- i + 1L
      out$mem <- args[[i]]
    } else if (grepl("^--mem=", arg)) out$mem <- sub("^--mem=", "", arg)
    else if (identical(arg, "--cpus-per-task")) {
      i <- i + 1L
      out$cpus_per_task <- as.integer(args[[i]])
    } else if (grepl("^--cpus-per-task=", arg)) out$cpus_per_task <- as.integer(sub("^--cpus-per-task=", "", arg))
    else if (identical(arg, "--time")) {
      i <- i + 1L
      out$time <- args[[i]]
    } else if (grepl("^--time=", arg)) out$time <- sub("^--time=", "", arg)
    else if (identical(arg, "--summary-mem")) {
      i <- i + 1L
      out$summary_mem <- args[[i]]
    } else if (grepl("^--summary-mem=", arg)) out$summary_mem <- sub("^--summary-mem=", "", arg)
    else if (identical(arg, "--summary-cpus-per-task")) {
      i <- i + 1L
      out$summary_cpus_per_task <- cap_threads(args[[i]])
    } else if (grepl("^--summary-cpus-per-task=", arg)) out$summary_cpus_per_task <- cap_threads(sub("^--summary-cpus-per-task=", "", arg))
    else if (identical(arg, "--summary-time")) {
      i <- i + 1L
      out$summary_time <- args[[i]]
    } else if (grepl("^--summary-time=", arg)) out$summary_time <- sub("^--summary-time=", "", arg)
    else if (identical(arg, "--summary-dependency")) {
      i <- i + 1L
      out$summary_dependency <- args[[i]]
    } else if (grepl("^--summary-dependency=", arg)) out$summary_dependency <- sub("^--summary-dependency=", "", arg)
    else if (identical(arg, "--help") || identical(arg, "-h")) {
      cat(
        "Usage:\n",
        "  Rscript benchmark/run_hpc_selected_depth_lambda_comparison.R --init [--output-dir DIR]\n",
        "  Rscript benchmark/run_hpc_selected_depth_lambda_comparison.R --task-id N --resume [--output-dir DIR]\n",
        "  Rscript benchmark/run_hpc_selected_depth_lambda_comparison.R --combine [--output-dir DIR]\n",
        "  Rscript benchmark/run_hpc_selected_depth_lambda_comparison.R --init --submit [--output-dir DIR]\n",
        "\n",
        "When --submit is used, the script submits both the task array and a dependent\n",
        "summary job. The summary job defaults to --summary-dependency=afterany so it\n",
        "runs after the array finishes even if some array tasks fail. Summary resources\n",
        "default to --summary-mem=256G and --summary-cpus-per-task=60.\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }
  out
}

repo_root <- function() {
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "benchmark/run_hpc_selected_depth_lambda_comparison.R"), ".."), mustWork = TRUE)
}

load_runner_env <- function(repo_dir) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo_dir, "benchmark/run_full_second_layer_comparison.R"), envir = env)
  env
}

hpc_selected_alfak2_grid <- function() {
  out <- expand.grid(
    input_mode = c("full", "minobs_matched", "soft_minobs"),
    extrapolation_method = alfak2:::second_layer_alfak2_methods(),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out$selection_source <- "all_9_extrapolation_methods"
  out
}

hpc_selected_build_run_index <- function() {
  n_chr <- 22L
  lambdas <- c(0.2, 0.6, 0.8, 1.0, 1.2)
  sample_depths <- c(200L, 2000L)
  landscape_reps <- 1:10
  fit_repeats <- 1:5
  landscape_grid <- expand.grid(
    sample_depth_index = seq_along(sample_depths),
    lambda_index = seq_along(lambdas),
    landscape_rep = landscape_reps,
    fit_repeat = fit_repeats,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  landscape_grid$sample_depth <- sample_depths[landscape_grid$sample_depth_index]
  landscape_grid$grf_lambda <- lambdas[landscape_grid$lambda_index]
  landscape_grid$n_chr <- n_chr
  landscape_grid$landscape_id <- sprintf(
    "chr%s_depth%s_lambda%s_rep%s",
    landscape_grid$n_chr,
    landscape_grid$sample_depth,
    gsub("\\.", "p", as.character(landscape_grid$grf_lambda)),
    landscape_grid$landscape_rep
  )
  alfak2 <- merge(
    landscape_grid,
    data.frame(package = "alfak2", hpc_selected_alfak2_grid(), stringsAsFactors = FALSE),
    all = TRUE
  )
  alfak2$minobs <- NA_integer_
  alfak2$NN_prior_slot <- NA_character_
  alfak2$NN_prior <- NA_character_
  alfak2$NN_prior_value <- NA_character_

  slots <- alfak2:::second_layer_alfakR_slots()
  alfakR_modes <- expand.grid(
    package = "alfakR",
    minobs = c(5L, 10L, 20L),
    slot_row = seq_len(nrow(slots)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  alfakR_modes$NN_prior_slot <- slots$NN_prior_slot[alfakR_modes$slot_row]
  alfakR_modes$NN_prior <- slots$NN_prior[alfakR_modes$slot_row]
  alfakR_modes$NN_prior_value <- slots$NN_prior_value[alfakR_modes$slot_row]
  alfakR_modes$slot_row <- NULL
  alfakR <- merge(landscape_grid, alfakR_modes, all = TRUE)
  alfakR$input_mode <- NA_character_
  alfakR$extrapolation_method <- NA_character_
  alfakR$selection_source <- "all_alfakR_settings"

  cols <- c(
    "package", "n_chr", "sample_depth", "sample_depth_index", "grf_lambda", "lambda_index",
    "landscape_id", "landscape_rep", "fit_repeat", "input_mode",
    "extrapolation_method", "selection_source", "minobs", "NN_prior_slot",
    "NN_prior", "NN_prior_value"
  )
  out <- rbind(alfak2[, cols], alfakR[, cols])
  out <- out[order(out$sample_depth, out$grf_lambda, out$landscape_rep, out$fit_repeat, out$package, out$input_mode, out$extrapolation_method, out$minobs, out$NN_prior_slot), ]
  out$task_id <- seq_len(nrow(out))
  out$run_id <- sprintf("hpc_chr22_9method_%05d", out$task_id)
  out$landscape_seed <- 100000L + 10000L * as.integer(out$lambda_index) + 100L * as.integer(out$landscape_rep)
  out$count_seed <- 200000L + 100000L * as.integer(out$sample_depth_index) +
    10000L * as.integer(out$lambda_index) + 100L * as.integer(out$landscape_rep) +
    as.integer(out$fit_repeat)
  out$fit_seed <- 300000L + 100000L * as.integer(out$sample_depth_index) +
    10000L * as.integer(out$lambda_index) + 100L * as.integer(out$landscape_rep) +
    as.integer(out$fit_repeat)
  out[, c("task_id", "run_id", setdiff(names(out), c("task_id", "run_id")))]
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

write_slurm_script <- function(repo_dir, output_dir, alfakR_repo, n_tasks, qos, mem, cpus_per_task, time) {
  slurm_dir <- file.path(output_dir, "slurm")
  log_dir <- file.path(output_dir, "logs")
  dir.create(slurm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(slurm_dir, "run_hpc_selected_depth_lambda_comparison.sbatch")
  lines <- c(
    "#!/bin/bash",
    "#SBATCH --job-name=alfak2_hpc_selected",
    sprintf("#SBATCH --array=1-%d", as.integer(n_tasks)),
    sprintf("#SBATCH --qos=%s", qos),
    sprintf("#SBATCH --cpus-per-task=%d", as.integer(cpus_per_task)),
    sprintf("#SBATCH --mem=%s", mem),
    sprintf("#SBATCH --time=%s", time),
    sprintf("#SBATCH --output=%s/%%A_%%a.out", log_dir),
    sprintf("#SBATCH --error=%s/%%A_%%a.err", log_dir),
    "",
    "source /etc/profile || true",
    "module load R/4.4.2-gfbf-2024a",
    "set -euo pipefail",
    sprintf("REPO_DIR=${REPO_DIR:-%s}", shQuote(repo_dir)),
    sprintf("OUTPUT_DIR=${OUTPUT_DIR:-%s}", shQuote(output_dir)),
    sprintf("ALFAKR_REPO=${ALFAKR_REPO:-%s}", shQuote(alfakR_repo)),
    "R_BIN=${R_BIN:-Rscript}",
    "cd \"${REPO_DIR}\"",
    "\"${R_BIN}\" benchmark/run_hpc_selected_depth_lambda_comparison.R \\",
    "  --task-id \"${SLURM_ARRAY_TASK_ID}\" \\",
    "  --output-dir \"${OUTPUT_DIR}\" \\",
    "  --alfakR-repo \"${ALFAKR_REPO}\" \\",
    "  --resume"
  )
  writeLines(lines, script)
  Sys.chmod(script, mode = "0755")
  script
}

write_summary_slurm_script <- function(repo_dir,
                                       output_dir,
                                       alfakR_repo,
                                       qos,
                                       mem,
                                       cpus_per_task,
                                       time) {
  slurm_dir <- file.path(output_dir, "slurm")
  log_dir <- file.path(output_dir, "logs")
  dir.create(slurm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(slurm_dir, "summarize_hpc_selected_depth_lambda_comparison.sbatch")
  lines <- c(
    "#!/bin/bash",
    "#SBATCH --job-name=alfak2_hpc_summary",
    sprintf("#SBATCH --qos=%s", qos),
    sprintf("#SBATCH --cpus-per-task=%d", as.integer(cpus_per_task)),
    sprintf("#SBATCH --mem=%s", mem),
    sprintf("#SBATCH --time=%s", time),
    sprintf("#SBATCH --output=%s/summary_%%j.out", log_dir),
    sprintf("#SBATCH --error=%s/summary_%%j.err", log_dir),
    "",
    "source /etc/profile || true",
    "module load R/4.4.2-gfbf-2024a",
    "set -euo pipefail",
    sprintf("REPO_DIR=${REPO_DIR:-%s}", shQuote(repo_dir)),
    sprintf("OUTPUT_DIR=${OUTPUT_DIR:-%s}", shQuote(output_dir)),
    sprintf("ALFAKR_REPO=${ALFAKR_REPO:-%s}", shQuote(alfakR_repo)),
    "R_BIN=${R_BIN:-Rscript}",
    "cd \"${REPO_DIR}\"",
    "\"${R_BIN}\" benchmark/run_hpc_selected_depth_lambda_comparison.R \\",
    "  --combine \\",
    "  --output-dir \"${OUTPUT_DIR}\" \\",
    "  --alfakR-repo \"${ALFAKR_REPO}\"",
    "\"${R_BIN}\" benchmark/scr/summarize_hpc_selected_depth_lambda_comparison.R \\",
    "  --input-dir=\"${OUTPUT_DIR}\" \\",
    "  --threads=\"${SLURM_CPUS_PER_TASK:-1}\""
  )
  writeLines(lines, script)
  Sys.chmod(script, mode = "0755")
  script
}

submit_sbatch <- function(script, extra_args = character()) {
  cmd_out <- tryCatch(
    system2("sbatch", c("--parsable", extra_args, script), stdout = TRUE, stderr = TRUE),
    error = function(e) {
      structure(conditionMessage(e), status = 1L)
    }
  )
  status <- attr(cmd_out, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    stop("`sbatch` failed with status ", status, ":\n", paste(cmd_out, collapse = "\n"), call. = FALSE)
  }
  first_line <- if (length(cmd_out)) cmd_out[[1L]] else ""
  job_id <- sub(";.*$", "", trimws(first_line))
  if (!nzchar(job_id)) {
    stop("Could not parse Slurm job id from sbatch output:\n", paste(cmd_out, collapse = "\n"), call. = FALSE)
  }
  list(job_id = job_id, output = cmd_out)
}

load_repos_for_task <- function(runner, repo_dir, alfakR_repo) {
  runner$load_repos(repo_dir, alfakR_repo)
}

run_task <- function(args, runner, repo_dir, output_dir, alfakR_repo) {
  run_index_path <- file.path(output_dir, "run_index.csv")
  if (!file.exists(run_index_path)) {
    stop("Missing run_index.csv. Run --init before submitting tasks.", call. = FALSE)
  }
  run_index <- read.csv(run_index_path, stringsAsFactors = FALSE)
  task_id <- as.integer(args$task_id)
  if (!is.finite(task_id) || task_id < 1L || task_id > nrow(run_index)) {
    stop("`--task-id` is outside the run index: ", task_id, call. = FALSE)
  }
  alfakR_status <- load_repos_for_task(runner, repo_dir, alfakR_repo)
  row <- run_index[task_id, , drop = FALSE]
  result <- runner$run_one(row, output_dir, resume = args$resume, alfakR_loaded = identical(alfakR_status, "loaded"))
  cat(
    "task_id=", task_id,
    " run_id=", row$run_id,
    " package=", row$package,
    " status=", result$result$status,
    " dependency_status=", result$result$dependency_status %||% NA_character_,
    "\n",
    sep = ""
  )
}

combine_results <- function(runner, output_dir, threads = hpc_worker_threads()) {
  threads <- cap_threads(threads)
  run_index_path <- file.path(output_dir, "run_index.csv")
  if (!file.exists(run_index_path)) stop("Missing run_index.csv.", call. = FALSE)
  run_index <- read.csv(run_index_path, stringsAsFactors = FALSE)
  cache_paths <- file.path(output_dir, "run_cache", paste0(run_index$run_id, ".rds"))
  message("Combining ", nrow(run_index), " run caches with ", threads, " worker(s).")
  read_one_result <- function(i) {
    if (file.exists(cache_paths[[i]])) {
      readRDS(cache_paths[[i]])
    } else {
      row <- run_index[i, , drop = FALSE]
      list(
        row = row,
        result = list(
          status = "missing",
          failure_status = "missing",
          error_message = "Missing run cache.",
          runtime_seconds = NA_real_,
          fit_path = NA_character_,
          dependency_status = "missing",
          diagnostics = list()
        ),
        metrics = data.frame()
      )
    }
  }
  results <- parallel_lapply(seq_len(nrow(run_index)), read_one_result, threads = threads)
  run_index_done <- runner$rbind_fill(lapply(results, function(x) {
    data.frame(
      as.data.frame(x$row, stringsAsFactors = FALSE),
      fit_status = x$result$status,
      failure_status = x$result$failure_status,
      runtime_seconds = x$result$diagnostics$runtime_seconds %||% x$result$runtime_seconds %||% NA_real_,
      dependency_status = x$result$dependency_status %||% NA_character_,
      error_message = x$result$error_message %||% NA_character_,
      output_path = x$result$fit_path %||% NA_character_,
      cache_key = x$result$fit_path %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(run_index_done, file.path(output_dir, "run_index.csv"))

  metrics <- runner$rbind_fill(lapply(results, `[[`, "metrics"))
  write_csv(metrics, file.path(output_dir, "metrics_by_run.csv"))
  failures <- run_index_done[
    run_index_done$fit_status != "ok" | runner$dependency_requires_attention(run_index_done$dependency_status),
    ,
    drop = FALSE
  ]
  write_csv(failures, file.path(output_dir, "failures.csv"))
  dependency_cols <- runner$present_cols(run_index_done, c(
    "task_id", "run_id", "package", "n_chr", "sample_depth", "grf_lambda", "landscape_id",
    "landscape_rep", "fit_repeat", "input_mode", "extrapolation_method",
    "minobs", "NN_prior_slot", "dependency_status", "fit_status",
    "failure_status", "runtime_seconds", "output_path"
  ))
  dependency <- run_index_done[, dependency_cols, drop = FALSE]
  write_csv(dependency, file.path(output_dir, "dependency_status.csv"))

  if (nrow(metrics)) {
    paired <- runner$build_paired_comparison(metrics)
    baseline <- runner$build_baseline_delta(metrics)
    rankings <- runner$build_landscape_rankings(metrics)
    pareto <- runner$build_pareto_front(metrics)
    summary_by <- c(
      "n_chr", "sample_depth", "grf_lambda", "package", "input_mode",
      "extrapolation_method", "minobs", "NN_prior_slot", "shell",
      "prediction_scale", "metric"
    )
    depth_lambda_summary <- runner$summary_stats(metrics, summary_by)
    overall_summary <- runner$summary_stats(
      metrics,
      c("n_chr", "package", "input_mode", "extrapolation_method", "minobs", "NN_prior_slot", "shell", "prediction_scale", "metric")
    )
    rank_summary <- runner$build_rank_summary(rankings)
    write_csv(paired, file.path(output_dir, "paired_landscape_comparison.csv"))
    write_csv(baseline, file.path(output_dir, "baseline_delta.csv"))
    write_csv(rankings, file.path(output_dir, "landscape_rankings.csv"))
    write_csv(pareto, file.path(output_dir, "pareto_front.csv"))
    write_csv(depth_lambda_summary, file.path(output_dir, "depth_lambda_summary.csv"))
    write_csv(depth_lambda_summary, file.path(output_dir, "lambda_summary.csv"))
    write_csv(overall_summary, file.path(output_dir, "overall_summary.csv"))
    write_csv(rank_summary, file.path(output_dir, "depth_lambda_method_rank_summary.csv"))
    write_csv(rank_summary, file.path(output_dir, "balanced_rank_summary.csv"))
    saveRDS(
      list(
        run_index = run_index_done,
        metrics_by_run = metrics,
        paired_landscape_comparison = paired,
        baseline_delta = baseline,
        landscape_rankings = rankings,
        depth_lambda_summary = depth_lambda_summary,
        overall_summary = overall_summary,
        method_rank_summary = rank_summary,
        pareto_front = pareto,
        failures = failures
      ),
      file.path(output_dir, "full_results.rds")
    )
  }

  report <- c(
    "# HPC Selected Depth/Lambda Comparison",
    "",
    sprintf("- expected_total_runs: %d", nrow(run_index_done)),
    sprintf("- expected_alfak2_runs: %d", sum(run_index_done$package == "alfak2")),
    sprintf("- expected_alfakR_runs: %d", sum(run_index_done$package == "alfakR")),
    sprintf("- actual_completed_runs: %d", sum(run_index_done$fit_status != "missing")),
    sprintf("- actual_successful_runs: %d", sum(run_index_done$fit_status == "ok")),
    sprintf("- actual_failed_or_missing_runs: %d", sum(run_index_done$fit_status != "ok")),
    sprintf("- n_chr_values: %s", paste(sort(unique(run_index_done$n_chr %||% NA_integer_)), collapse = ", ")),
    sprintf("- sample_depth_values: %s", paste(sort(unique(run_index_done$sample_depth)), collapse = ", ")),
    sprintf("- grf_lambda_values: %s", paste(sort(unique(run_index_done$grf_lambda)), collapse = ", ")),
    sprintf("- landscape_reps: %s", paste(sort(unique(run_index_done$landscape_rep)), collapse = ", ")),
    sprintf("- fit_repeats: %s", paste(sort(unique(run_index_done$fit_repeat)), collapse = ", ")),
    "",
    "## alfak2 combinations",
    "- All 9 extrapolation methods are evaluated under each of the three alfak2 input modes: `full`, `minobs_matched`, and `soft_minobs`.",
    paste(sprintf("- `%s:%s`", hpc_selected_alfak2_grid()$input_mode, hpc_selected_alfak2_grid()$extrapolation_method), collapse = "\n"),
    "",
    "## alfakR settings",
    "- All current 15 minobs / NN_prior slots are included.",
    "",
    "## Outputs",
    "- `metrics_by_run.csv`",
    "- `depth_lambda_summary.csv`",
    "- `paired_landscape_comparison.csv`",
    "- `baseline_delta.csv`",
    "- `landscape_rankings.csv`",
    "- `balanced_rank_summary.csv`",
    "- `pareto_front.csv`",
    "- `failures.csv`"
  )
  writeLines(report, file.path(output_dir, "report.md"))
  cat(
    "combined runs:", nrow(run_index_done), "\n",
    "successful:", sum(run_index_done$fit_status == "ok"), "\n",
    "failed_or_missing:", sum(run_index_done$fit_status != "ok"), "\n",
    "output_dir:", output_dir, "\n",
    sep = ""
  )
}

main <- function() {
  args <- parse_hpc_args(commandArgs(trailingOnly = TRUE))
  repo_dir <- repo_root()
  output_dir <- args$output_dir
  if (!grepl("^/", output_dir)) output_dir <- file.path(repo_dir, output_dir)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  alfakR_repo <- normalizePath(args$alfakR_repo %||% file.path(dirname(repo_dir), "alfakR"), mustWork = FALSE)
  runner <- load_runner_env(repo_dir)

  if (isTRUE(args$init)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package `pkgload` is required to initialize the benchmark index.", call. = FALSE)
    }
    pkgload::load_all(repo_dir, quiet = TRUE)
    run_index <- hpc_selected_build_run_index()
    write_csv(run_index, file.path(output_dir, "run_index.csv"))
    slurm_script <- write_slurm_script(
      repo_dir = repo_dir,
      output_dir = output_dir,
      alfakR_repo = alfakR_repo,
      n_tasks = nrow(run_index),
      qos = args$qos,
      mem = args$mem,
      cpus_per_task = args$cpus_per_task,
      time = args$time
    )
    summary_script <- write_summary_slurm_script(
      repo_dir = repo_dir,
      output_dir = output_dir,
      alfakR_repo = alfakR_repo,
      qos = args$qos,
      mem = args$summary_mem,
      cpus_per_task = args$summary_cpus_per_task,
      time = args$summary_time
    )
    cat(
      "initialized run_index:", file.path(output_dir, "run_index.csv"), "\n",
      "expected_total_runs:", nrow(run_index), "\n",
      "expected_alfak2_runs:", sum(run_index$package == "alfak2"), "\n",
      "expected_alfakR_runs:", sum(run_index$package == "alfakR"), "\n",
      "slurm_script:", slurm_script, "\n",
      "summary_slurm_script:", summary_script, "\n",
      sep = ""
    )
    if (isTRUE(args$submit)) {
      dependency_type <- match.arg(args$summary_dependency, c("afterany", "afterok"))
      array_submit <- submit_sbatch(slurm_script)
      summary_submit <- submit_sbatch(
        summary_script,
        extra_args = sprintf("--dependency=%s:%s", dependency_type, array_submit$job_id)
      )
      cat(
        "submitted_array_job_id:", array_submit$job_id, "\n",
        "submitted_summary_job_id:", summary_submit$job_id, "\n",
        "summary_dependency:", dependency_type, ":", array_submit$job_id, "\n",
        sep = ""
      )
    }
  }

  if (!is.null(args$task_id)) {
    run_task(args, runner, repo_dir, output_dir, alfakR_repo)
  }

  if (isTRUE(args$combine)) {
    combine_results(runner, output_dir, threads = hpc_worker_threads())
  }
}

if (identical(environment(), globalenv())) main()
