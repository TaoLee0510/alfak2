#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_grf_alfak2_parameter_calibration.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_grf_alfak2_vs_alfakR_benchmark.R"))

usage <- function() {
  cat(
    "Calibrate alfak2 parameters on sparse GRF two-timepoint data.\n\n",
    "This script is intentionally separate from the alfak2-vs-alfakR benchmark.\n",
    "It uses the same GRF/ABM input generator, fits only alfak2 over a parameter\n",
    "grid, and ranks parameter sets by sparse-data recovery against known GRF truth.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_grf_alfak2_parameter_calibration.R [options]\n\n",
    "Core options:\n",
    "  --mode=all|prepare|fit|fit-task|summarize\n",
    "  --task-index=1\n",
    "  --alfak2-repo=/Users/you/path/alfak2\n",
    "  --alfakR-repo=/Users/you/path/alfakR\n",
    "  --output-dir=benchmark/results/grf_alfak2_calibration\n",
    "  --source-input-dir=benchmark/results/grf_alfak2_vs_alfakR\n",
    "  --n-cores=1\n",
    "  --force-refit=false\n",
    "  --force-sim=false\n",
    "  --recompile-dll=false\n\n",
    "Sparse GRF input options:\n",
    "  --minobs=5,10,20\n",
    "  --n-sim=1\n",
    "  --lambdas=0.8\n",
    "  --time-gaps=2,4\n",
    "  --sample-depth=2000\n",
    "  --grf-centroid-mode=method_blind  # method_blind|nn_initial\n",
    "  --grf-centroid-min-cn=0\n",
    "  --grf-centroid-max-cn=4\n",
    "  --time-max=360\n",
    "  --passage-interval=45\n\n",
    "Calibration grid options:\n",
    "  --legacy-weights=pi0,directly_informed,uniform\n",
    "  --correct-efflux-values=true,false\n",
    "  --graph-edge-weights=normalized  # default normalized; add mutation for legacy baseline\n",
    "  --lambda-l-values=0.2,1,5\n",
    "  --lambda-e-values=0.05,0.25,1\n",
    "  --sigma-obs-values=0.02,0.05,0.1\n",
    "  --dm-concentrations=50\n",
    "  --effective-depth-modes=min\n",
    "  --local-shell-depths=0\n",
    "  --global-extra-shells=1\n\n",
    "Ranking options:\n",
    "  --objective-scope=nn|direct|all|holdout\n",
    "  --objective-metric=sparse_composite\n",
    "  --direct-weight=0.25\n",
    "  --holdout-weight=0.15\n",
    "  --holdout-fraction=0.25\n",
    "  --holdout-min-direct=6\n",
    "  --holdout-failure-penalty=1\n",
    "  --bias-weight=0.10\n",
    "  --spearman-weight=0.10\n",
    "  --false-high-weight=0.10\n",
    sep = ""
  )
}

default_alfakR_repo <- function(repo_dir) {
  sibling <- normalizePath(file.path(dirname(repo_dir), "alfakR"), winslash = "/", mustWork = FALSE)
  if (dir.exists(sibling)) sibling else "/share/lab_crd/lab_crd/taoli/Project/alfakR"
}

arg_logical_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  raw <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  raw <- raw[nzchar(raw)]
  if (!length(raw)) return(default)
  vapply(raw, function(x) {
    lx <- tolower(x)
    if (lx %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (lx %in% c("false", "f", "0", "no", "n")) return(FALSE)
    stop("Expected boolean values for --", gsub("_", "-", name), call. = FALSE)
  }, logical(1))
}

to_num <- function(x) suppressWarnings(as.numeric(x))

hash_r_object <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}

calibration_fit_config_hash <- function(task, cfg, repo_versions, actual_input_md5 = NA_character_) {
  task_for_hash <- task[setdiff(names(task), "outdir")]
  cfg_keys <- c(
    "pm", "n0", "nb", "drop_diploid", "input_policy",
    "alfak2_input_depth", "alfak2_effective_depth", "alfak2_observation_model",
    "alfak2_min_cn", "alfak2_max_cn", "max_nodes",
    "alfak2_anchor_count_reference", "alfak2_anchor_count_power",
    "eval_max", "iter_max", "retry_max",
    "holdout_fraction", "holdout_min_direct", "holdout_seed", "holdout_mode"
  )
  hash_r_object(list(
    cache_schema = 2L,
    task = task_for_hash,
    actual_input_md5 = as.character(actual_input_md5),
    config = cfg[intersect(cfg_keys, names(cfg))],
    repo_versions = repo_versions,
    implementation = list(
      holdout_mode = "zero_observation_weight",
      weighted_likelihood = "weighted_dirichlet_multinomial_v1"
    )
  ))
}

build_calibration_dirs <- function(output_dir) {
  dirs <- list(
    root = output_dir,
    cache = file.path(output_dir, "cache"),
    fits = file.path(output_dir, "fits"),
    tables = file.path(output_dir, "tables"),
    fit_parts = file.path(output_dir, "tables", "fit_results_parts")
  )
  for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dirs
}

build_calibration_config <- function(args, repo_dir) {
  legacy_weights <- arg_character_vec(args, "legacy_weights", c("pi0", "directly_informed", "uniform"))
  bad_weights <- setdiff(legacy_weights, c("pi0", "directly_informed", "uniform"))
  if (length(bad_weights)) stop("Unsupported legacy weights: ", paste(bad_weights, collapse = ", "), call. = FALSE)

  objective_metric <- as.character(arg_value(args, "objective_metric", "sparse_composite"))
  objective_metric <- match.arg(
    objective_metric,
    c("mae", "rmse", "centered_mae", "centered_rmse", "median_abs_error", "sparse_composite")
  )
  objective_scope <- as.character(arg_value(args, "objective_scope", "nn"))
  objective_scope <- match.arg(objective_scope, c("nn", "direct", "all", "holdout"))
  input_policy <- as.character(arg_value(args, "input_policy", "minobs_matched"))
  input_policy <- match.arg(input_policy, c("minobs_matched", "full", "soft_minobs"))
  input_depth <- as.character(arg_value(args, "alfak2_input_depth", "effective"))
  input_depth <- match.arg(input_depth, c("raw", "effective"))
  observation_model <- as.character(arg_value(args, "alfak2_observation_model", "dirichlet_multinomial"))
  if (!nzchar(observation_model)) observation_model <- NULL
  graph_edge_weights <- arg_character_vec(args, "graph_edge_weights", "normalized")
  bad_graph_weights <- setdiff(graph_edge_weights, c("mutation", "unit", "normalized"))
  if (length(bad_graph_weights)) {
    stop("Unsupported graph edge weights: ", paste(bad_graph_weights, collapse = ", "), call. = FALSE)
  }

  list(
    repo_dir = repo_dir,
    alfak2_repo = normalizePath(arg_value(args, "alfak2_repo", repo_dir), winslash = "/", mustWork = FALSE),
    alfakR_repo = normalizePath(arg_value(args, "alfakR_repo", default_alfakR_repo(repo_dir)), winslash = "/", mustWork = FALSE),
    output_dir = normalize_output_dir(repo_dir, arg_value(args, "output_dir", "benchmark/results/grf_alfak2_calibration")),
    source_input_dir = {
      x <- arg_value(args, "source_input_dir", NULL)
      if (is.null(x)) NULL else normalize_output_dir(repo_dir, x)
    },
    mode = tolower(as.character(arg_value(args, "mode", "all"))),
    task_index = arg_integer(args, "task_index", NA_integer_),
    n_cores = arg_integer(args, "n_cores", 1L),
    seed = arg_integer(args, "seed", 424242L),
    force_refit = arg_logical(args, "force_refit", FALSE),
    force_sim = arg_logical(args, "force_sim", FALSE),
    reuse_dirty_cache = arg_logical(args, "reuse_dirty_cache", FALSE),
    recompile_dll = arg_logical(args, "recompile_dll", FALSE),
    minobs = sort(unique(arg_integer_vec(args, "minobs", c(5L, 10L, 20L)))),
    n_sim = arg_integer(args, "n_sim", 1L),
    lambdas = arg_numeric_vec(args, "lambdas", 0.8),
    time_starts = arg_numeric_vec(args, "time_starts", 0),
    time_gaps = arg_numeric_vec(args, "time_gaps", c(2, 4)),
    pm = arg_numeric(args, "pm", 5e-05),
    n0 = arg_numeric(args, "n0", 100000),
    nb = arg_numeric(args, "nb", 10000000),
    drop_diploid = arg_logical(args, "drop_diploid", TRUE),
    k_dim = arg_integer(args, "k_dim", 22L),
    n_centroids = arg_integer(args, "n_centroids", 64L),
    grf_centroid_mode = normalize_grf_centroid_mode(arg_value(args, "grf_centroid_mode", "method_blind")),
    grf_centroid_min_cn = arg_integer(args, "grf_centroid_min_cn", 0L),
    grf_centroid_max_cn = arg_integer(args, "grf_centroid_max_cn", 4L),
    grf_centroid_jitter_sd = {
      x <- arg_value(args, "grf_centroid_jitter_sd", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    time_max = arg_numeric(args, "time_max", 360),
    passage_interval = arg_numeric(args, "passage_interval", 45),
    sample_depth = arg_integer(args, "sample_depth", 2000L),
    abm_pop_size = arg_numeric(args, "abm_pop_size", 50000),
    abm_delta_t = arg_numeric(args, "abm_delta_t", 1),
    abm_max_pop = arg_numeric(args, "abm_max_pop", 2000000),
    abm_culling_survival = arg_numeric(args, "abm_culling_survival", 0.01),
    input_policy = input_policy,
    alfak2_input_depth = input_depth,
    effective_depth_modes = arg_character_vec(args, "effective_depth_modes", "min"),
    alfak2_effective_depth = {
      x <- arg_value(args, "alfak2_effective_depth", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    alfak2_observation_model = observation_model,
    dm_concentrations = arg_numeric_vec(args, "dm_concentrations", 50),
    alfak2_min_cn = arg_integer(args, "alfak2_min_cn", 0L),
    alfak2_max_cn = {
      x <- arg_value(args, "alfak2_max_cn", NA_character_)
      y <- suppressWarnings(as.integer(x))
      if (is.finite(y)) y else NA_integer_
    },
    local_shell_depths = sort(unique(arg_integer_vec(args, "local_shell_depths", 0L))),
    global_extra_shells = sort(unique(arg_integer_vec(args, "global_extra_shells", 1L))),
    max_nodes = arg_integer(args, "alfak2_max_nodes", 150000L),
    lambda_l_values = arg_numeric_vec(args, "lambda_l_values", c(0.2, 1, 5)),
    lambda_e_values = arg_numeric_vec(args, "lambda_e_values", c(0.05, 0.25, 1)),
    sigma_obs_values = arg_numeric_vec(args, "sigma_obs_values", c(0.02, 0.05, 0.1)),
    graph_edge_weights = graph_edge_weights,
    alfak2_anchor_count_reference = arg_anchor_count_reference(args),
    alfak2_anchor_count_power = arg_numeric(args, "alfak2_anchor_count_power", 1),
    legacy_weights = legacy_weights,
    correct_efflux_values = arg_logical_vec(args, "correct_efflux_values", c(TRUE, FALSE)),
    eval_max = arg_integer(args, "alfak2_eval_max", 500L),
    iter_max = arg_integer(args, "alfak2_iter_max", 500L),
    retry_max = arg_integer(args, "alfak2_retry_max", 2000L),
    objective_scope = objective_scope,
    objective_metric = objective_metric,
    direct_weight = arg_numeric(args, "direct_weight", 0.25),
    holdout_weight = arg_numeric(args, "holdout_weight", 0.15),
    holdout_fraction = arg_numeric(args, "holdout_fraction", 0.25),
    holdout_min_direct = arg_integer(args, "holdout_min_direct", 6L),
    holdout_seed = arg_integer(args, "holdout_seed", 71011L),
    holdout_mode = "zero_observation_weight",
    holdout_failure_penalty = arg_numeric(args, "holdout_failure_penalty", 1),
    bias_weight = arg_numeric(args, "bias_weight", 0.10),
    spearman_weight = arg_numeric(args, "spearman_weight", 0.10),
    false_high_weight = arg_numeric(args, "false_high_weight", 0.10),
    write_node_table = arg_logical(args, "write_node_table", FALSE)
  )
}

make_param_grid <- function(cfg) {
  grid <- expand.grid(
    legacy_weight = cfg$legacy_weights,
    correct_efflux = cfg$correct_efflux_values,
    lambda_l = cfg$lambda_l_values,
    lambda_e = cfg$lambda_e_values,
    sigma_obs = cfg$sigma_obs_values,
    graph_edge_weight = cfg$graph_edge_weights,
    dm_concentration = cfg$dm_concentrations,
    effective_depth_mode = cfg$effective_depth_modes,
    local_shell_depth = cfg$local_shell_depths,
    global_extra_shell = cfg$global_extra_shells,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$param_id <- seq_len(nrow(grid))
  grid <- grid[c("param_id", setdiff(names(grid), "param_id"))]
  grid
}

numeric_in <- function(x, values, tol = 1e-8) {
  x <- as.numeric(x)
  values <- as.numeric(values)
  vapply(x, function(xx) any(abs(xx - values) <= tol), logical(1))
}

resolve_source_cache_path <- function(path, source_dir) {
  path <- as.character(path)
  if (!length(path) || !nzchar(path[[1L]])) return(NA_character_)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  fallback <- file.path(source_dir, "cache", basename(path))
  if (file.exists(fallback)) return(normalizePath(fallback, winslash = "/", mustWork = TRUE))
  path
}

filter_source_input_table <- function(source_tbl, cfg) {
  needed <- c(
    "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "patient_id", "grf_key", "grf_rds", "input_rds", "input_md5", "minobs"
  )
  missing <- setdiff(needed, names(source_tbl))
  if (length(missing)) {
    stop("Source input table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  source_tbl$simulation_id <- as.integer(source_tbl$simulation_id)
  source_tbl$lambda <- as.numeric(source_tbl$lambda)
  source_tbl$time_start <- as.numeric(source_tbl$time_start)
  source_tbl$time_gap <- as.numeric(source_tbl$time_gap)
  source_tbl$time_delta <- as.numeric(source_tbl$time_delta)
  source_tbl$minobs <- as.integer(source_tbl$minobs)
  keep <- source_tbl$simulation_id %in% seq_len(cfg$n_sim) &
    numeric_in(source_tbl$lambda, cfg$lambdas) &
    numeric_in(source_tbl$time_start, cfg$time_starts) &
    numeric_in(source_tbl$time_gap, cfg$time_gaps) &
    source_tbl$minobs %in% cfg$minobs
  if ("grf_centroid_mode" %in% names(source_tbl)) {
    keep <- keep & as.character(source_tbl$grf_centroid_mode) == as.character(cfg$grf_centroid_mode)
  }
  out <- source_tbl[keep, , drop = FALSE]
  out <- out[order(out$simulation_id, out$lambda, out$time_start, out$time_gap, out$minobs), , drop = FALSE]
  if (!nrow(out)) {
    stop(
      "No rows in source input table match requested calibration grid. ",
      "Check --n-sim, --lambdas, --time-gaps, --time-starts, and --minobs.",
      call. = FALSE
    )
  }
  out
}

prepare_calibration_inputs_from_source <- function(cfg, dirs) {
  source_dir <- normalizePath(cfg$source_input_dir, winslash = "/", mustWork = TRUE)
  source_input_path <- file.path(source_dir, "tables", "input_table.tsv")
  if (!file.exists(source_input_path)) {
    stop("Missing source input table: ", source_input_path, call. = FALSE)
  }
  param_grid <- make_param_grid(cfg)
  input_tbl <- filter_source_input_table(read_tsv(source_input_path), cfg)
  input_tbl$source_input_dir <- source_dir
  input_tbl$input_source <- "reused"
  input_tbl$grf_rds <- vapply(input_tbl$grf_rds, resolve_source_cache_path, character(1), source_dir = source_dir)
  input_tbl$input_rds <- vapply(input_tbl$input_rds, resolve_source_cache_path, character(1), source_dir = source_dir)
  if ("input_csv" %in% names(input_tbl)) {
    input_tbl$input_csv <- vapply(input_tbl$input_csv, resolve_source_cache_path, character(1), source_dir = source_dir)
  }
  missing_grf <- unique(input_tbl$grf_rds[!file.exists(input_tbl$grf_rds)])
  missing_input <- unique(input_tbl$input_rds[!file.exists(input_tbl$input_rds)])
  if (length(missing_grf) || length(missing_input)) {
    msg <- c(
      if (length(missing_grf)) paste0("missing grf_rds: ", paste(utils::head(missing_grf, 3), collapse = ", ")),
      if (length(missing_input)) paste0("missing input_rds: ", paste(utils::head(missing_input, 3), collapse = ", "))
    )
    stop("Source input files are not available: ", paste(msg, collapse = " | "), call. = FALSE)
  }
  if (!"input_md5" %in% names(input_tbl)) input_tbl$input_md5 <- unname(tools::md5sum(input_tbl$input_rds))

  tasks <- list()
  task_idx <- 0L
  for (i in seq_len(nrow(input_tbl))) {
    row <- input_tbl[i, , drop = FALSE]
    for (pg_idx in seq_len(nrow(param_grid))) {
      pg <- param_grid[pg_idx, , drop = FALSE]
      task_idx <- task_idx + 1L
      tasks[[task_idx]] <- data.frame(
        task_order = task_idx,
        engine = "alfak2",
        method = "alfak2_calibration",
        input_policy = cfg$input_policy,
        simulation_id = row$simulation_id,
        lambda = row$lambda,
        lambda_label = row$lambda_label,
        time_start = row$time_start,
        time_gap = row$time_gap,
        time_delta = row$time_delta,
        patient_id = row$patient_id,
        grf_key = row$grf_key,
        grf_rds = row$grf_rds,
        grf_centroid_mode = row_field(row, "grf_centroid_mode", NA_character_),
        grf_centroid_min_cn = row_field(row, "grf_centroid_min_cn", NA_integer_),
        grf_centroid_max_cn = row_field(row, "grf_centroid_max_cn", NA_integer_),
        input_rds = row$input_rds,
        input_md5 = row$input_md5,
        minobs = as.integer(row$minobs),
        pm = cfg$pm,
        benchmark_seed = as.integer(cfg$seed + row$simulation_id * 10000L + pg$param_id),
        alfak2_input_depth = cfg$alfak2_input_depth,
        alfak2_effective_depth = if (is.null(cfg$alfak2_effective_depth)) NA_real_ else cfg$alfak2_effective_depth,
        alfak2_observation_model = if (is.null(cfg$alfak2_observation_model)) "" else cfg$alfak2_observation_model,
        param_grid[pg_idx, , drop = FALSE],
        outdir = file.path(
          dirs$fits,
          paste0("lambda_", row$lambda_label),
          paste0("gap_", path_token(row$time_gap)),
          paste0("MINOBS_", row$minobs),
          sprintf("param_%04d", pg$param_id),
          row$patient_id
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  task_tbl <- do.call(rbind, tasks)
  write_tsv(param_grid, file.path(dirs$tables, "parameter_grid.tsv"))
  write_tsv(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  message(
    "Prepared calibration tasks from source inputs: ",
    nrow(input_tbl), " source input rows, ", nrow(task_tbl), " fit tasks."
  )
  invisible(list(input_table = input_tbl, task_table = task_tbl, parameter_grid = param_grid))
}

prepare_calibration_inputs <- function(cfg, dirs, repo_versions) {
  if (!is.null(cfg$source_input_dir)) {
    return(prepare_calibration_inputs_from_source(cfg, dirs))
  }

  input_rows <- list()
  tasks <- list()
  input_idx <- 0L
  task_idx <- 0L
  param_grid <- make_param_grid(cfg)
  time_axis_label <- paste0("tmax_", path_token(cfg$time_max), "_pint_", path_token(cfg$passage_interval))
  grf_landscape_label <- grf_landscape_token(cfg)

  for (sim_idx in seq_len(cfg$n_sim)) {
    for (lambda_idx in seq_along(cfg$lambdas)) {
      lambda <- cfg$lambdas[[lambda_idx]]
      lambda_label <- format_grf_label(lambda)
      abm_seed <- cfg$seed + sim_idx * 10000L + lambda_idx * 100L
      grf_key <- paste(sim_idx, lambda_label, grf_landscape_label, time_axis_label, sep = "__")
      grf_path <- file.path(dirs$cache, paste0("grf_sim_", grf_key, ".rds"))
      use_grf_cache <- !isTRUE(cfg$force_sim) &&
        file.exists(grf_path) &&
        (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfakR"]))
      grf_sim <- if (use_grf_cache) {
        readRDS(grf_path)
      } else {
        message("Simulating calibration GRF ABM: sim=", sim_idx, " lambda=", lambda)
        out <- simulate_nn_prior_grf_abm(
          seed = abm_seed,
          lambda = lambda,
          p = cfg$pm,
          k_dim = cfg$k_dim,
          n_centroids = cfg$n_centroids,
          time_max = cfg$time_max,
          passage_interval = cfg$passage_interval,
          abm_pop_size = cfg$abm_pop_size,
          abm_delta_t = cfg$abm_delta_t,
          abm_max_pop = cfg$abm_max_pop,
          abm_culling_survival = cfg$abm_culling_survival,
          centroid_mode = cfg$grf_centroid_mode,
          centroid_min_cn = cfg$grf_centroid_min_cn,
          centroid_max_cn = cfg$grf_centroid_max_cn,
          centroid_jitter_sd = cfg$grf_centroid_jitter_sd
        )
        saveRDS(out, grf_path)
        out
      }

      for (time_start in cfg$time_starts) {
        for (time_gap in cfg$time_gaps) {
          patient_id <- paste0(
            "grf_", sim_idx, "_lambda_", lambda_label, "_", grf_landscape_label, "_", time_axis_label,
            "_start_", path_token(time_start), "_gap_", path_token(time_gap)
          )
          input_rds <- file.path(dirs$cache, paste0("input_", patient_id, ".rds"))
          input_csv <- file.path(dirs$cache, paste0("input_", patient_id, ".csv"))
          input_seed <- abm_seed + as.integer(round(time_start * 10)) + as.integer(round(time_gap * 1000))
          yi <- if (!file.exists(input_rds) || isTRUE(cfg$force_sim)) {
            out <- build_two_timepoint_yi_from_abm(
              sim_wide = grf_sim$sim_wide,
              time_start = time_start,
              time_gap = time_gap,
              passage_interval = cfg$passage_interval,
              sample_depth = cfg$sample_depth,
              seed = input_seed
            )
            saveRDS(out, input_rds)
            utils::write.csv(data.frame(karyotype = rownames(out$x), out$x, check.names = FALSE), input_csv, row.names = FALSE)
            out
          } else {
            readRDS(input_rds)
          }
          input_md5 <- unname(tools::md5sum(input_rds))
          input_summary <- summarize_input_rows(yi, cfg$minobs, drop_diploid = cfg$drop_diploid)
          for (rr in seq_len(nrow(input_summary))) {
            input_idx <- input_idx + 1L
            input_rows[[input_idx]] <- data.frame(
              simulation_id = sim_idx,
              lambda = lambda,
              lambda_label = lambda_label,
              time_start = time_start,
              time_gap = time_gap,
              time_delta = as.numeric(yi$metadata$time_delta),
              patient_id = patient_id,
              grf_key = grf_key,
              grf_rds = grf_path,
              grf_centroid_mode = cfg$grf_centroid_mode,
              grf_centroid_min_cn = cfg$grf_centroid_min_cn,
              grf_centroid_max_cn = cfg$grf_centroid_max_cn,
              input_rds = input_rds,
              input_csv = input_csv,
              input_md5 = input_md5,
              input_summary[rr, , drop = FALSE],
              stringsAsFactors = FALSE
            )
          }

          for (minobs in cfg$minobs) {
            for (pg_idx in seq_len(nrow(param_grid))) {
              pg <- param_grid[pg_idx, , drop = FALSE]
              task_idx <- task_idx + 1L
              tasks[[task_idx]] <- data.frame(
                task_order = task_idx,
                engine = "alfak2",
                method = "alfak2_calibration",
                input_policy = cfg$input_policy,
                simulation_id = sim_idx,
                lambda = lambda,
                lambda_label = lambda_label,
                time_start = time_start,
                time_gap = time_gap,
                time_delta = as.numeric(yi$metadata$time_delta),
                patient_id = patient_id,
                grf_key = grf_key,
                grf_rds = grf_path,
                grf_centroid_mode = cfg$grf_centroid_mode,
                grf_centroid_min_cn = cfg$grf_centroid_min_cn,
                grf_centroid_max_cn = cfg$grf_centroid_max_cn,
                input_rds = input_rds,
                input_md5 = input_md5,
                minobs = as.integer(minobs),
                pm = cfg$pm,
                benchmark_seed = as.integer(abm_seed + time_gap * 1000L + minobs * 100L + pg$param_id),
                alfak2_input_depth = cfg$alfak2_input_depth,
                alfak2_effective_depth = if (is.null(cfg$alfak2_effective_depth)) NA_real_ else cfg$alfak2_effective_depth,
                alfak2_observation_model = if (is.null(cfg$alfak2_observation_model)) "" else cfg$alfak2_observation_model,
                param_grid[pg_idx, , drop = FALSE],
                outdir = file.path(
                  dirs$fits,
                  paste0("lambda_", lambda_label),
                  paste0("gap_", path_token(time_gap)),
                  paste0("MINOBS_", minobs),
                  sprintf("param_%04d", pg$param_id),
                  patient_id
                ),
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }

  input_tbl <- do.call(rbind, input_rows)
  task_tbl <- do.call(rbind, tasks)
  write_tsv(param_grid, file.path(dirs$tables, "parameter_grid.tsv"))
  write_tsv(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  message("Prepared calibration inputs: ", nrow(input_tbl), " input rows, ", nrow(task_tbl), " fit tasks.")
  invisible(list(input_table = input_tbl, task_table = task_tbl, parameter_grid = param_grid))
}

calibration_estimate_columns <- function(summary) {
  est_col <- if ("fitness_mean_alfakR_scale" %in% names(summary)) "fitness_mean_alfakR_scale" else "fitness_mean"
  sd_col <- if ("fitness_sd_alfakR_scale" %in% names(summary)) "fitness_sd_alfakR_scale" else "fitness_sd"
  list(mean = est_col, sd = sd_col)
}

select_calibration_holdout_labels <- function(fit, task, cfg) {
  direct <- fit$local$summary[
    as.character(fit$local$summary$support_tier) == "directly_informed" &
      is.finite(fit$local$summary$fitness_mean),
    ,
    drop = FALSE
  ]
  labels <- intersect(as.character(direct$karyotype), as.character(fit$global$summary$karyotype))
  labels <- labels[nzchar(labels)]
  min_direct <- max(3L, as.integer(cfg$holdout_min_direct))
  if (length(labels) < min_direct) return(character(0))
  set.seed(as.integer(cfg$holdout_seed) + as.integer(task$task_order[[1L]]) * 1009L + as.integer(task$param_id[[1L]]))
  n_holdout <- max(1L, as.integer(round(length(labels) * as.numeric(cfg$holdout_fraction))))
  n_holdout <- min(n_holdout, length(labels) - 3L)
  if (n_holdout < 1L) return(character(0))
  sample(labels, n_holdout)
}

zero_weight_holdout_counts <- function(counts, labels) {
  labels <- intersect(as.character(labels), rownames(counts))
  out <- counts
  observation_weights <- attr(counts, "observation_weights", exact = TRUE)
  soft_minobs <- attr(counts, "soft_minobs", exact = TRUE)
  if (is.null(observation_weights)) {
    observation_weights <- matrix(
      1,
      nrow = nrow(counts),
      ncol = 2L,
      dimnames = list(rownames(counts), colnames(counts))
    )
  }
  observation_weights[labels, ] <- 0
  attr(out, "observation_weights") <- observation_weights
  if (!is.null(soft_minobs)) attr(out, "soft_minobs") <- soft_minobs
  attr(out, "holdout_mode") <- list(mode = "zero_observation_weight", labels = labels)
  out
}

fit_calibration_holdout <- function(task, cfg, counts, dt, max_cn, fit, outdir) {
  holdout <- select_calibration_holdout_labels(fit, task, cfg)
  holdout_fit_path <- file.path(outdir, "alfak2_holdout_fit.rds")
  holdout_summary_path <- file.path(outdir, "holdout_landscape.rds")
  if (!length(holdout)) {
    return(list(
      holdout_status = "skipped_insufficient_direct",
      holdout_error_message = NA_character_,
      holdout_labels = NA_character_,
      holdout_n_labels = 0L,
      holdout_mode = cfg$holdout_mode,
      holdout_fit_path = holdout_fit_path,
      holdout_landscape_path = holdout_summary_path
    ))
  }
  holdout_counts <- zero_weight_holdout_counts(counts, holdout)
  if (!nrow(holdout_counts)) {
    return(list(
      holdout_status = "skipped_empty_input",
      holdout_error_message = NA_character_,
      holdout_labels = paste(holdout, collapse = ","),
      holdout_n_labels = length(holdout),
      holdout_mode = cfg$holdout_mode,
      holdout_fit_path = holdout_fit_path,
      holdout_landscape_path = holdout_summary_path
    ))
  }
  anchor_count_reference <- resolve_anchor_count_reference(cfg, task)
  tryCatch({
    set.seed(as.integer(cfg$holdout_seed) + as.integer(task$task_order[[1L]]) * 7919L + as.integer(task$param_id[[1L]]))
    hfit <- alfak2::fit_alfak2(
      holdout_counts,
      dt = dt,
      beta = cfg$pm,
      min_cn = cfg$alfak2_min_cn,
      max_cn = as.integer(max_cn),
      local_shell_depth = as.integer(task$local_shell_depth),
      global_extra_shell = as.integer(task$global_extra_shell),
      max_nodes = cfg$max_nodes,
      lambda_l_grid = as.numeric(task$lambda_l),
      lambda_e_grid = as.numeric(task$lambda_e),
      sigma_obs_grid = as.numeric(task$sigma_obs),
      graph_edge_weight = as.character(row_field(task, "graph_edge_weight", "mutation")),
      anchor_support_tiers = "directly_informed",
      anchor_exclude = holdout,
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = cfg$alfak2_anchor_count_power,
      input_depth = cfg$alfak2_input_depth,
      effective_depth = cfg$alfak2_effective_depth,
      effective_depth_mode = as.character(task$effective_depth_mode),
      observation_model = cfg$alfak2_observation_model,
      dm_concentration = as.numeric(task$dm_concentration),
      alfakR_scale = TRUE,
      n0 = cfg$n0,
      nb = cfg$nb,
      correct_efflux = as.logical(task$correct_efflux),
      legacy_weight = as.character(task$legacy_weight),
      control = list(eval.max = cfg$eval_max, iter.max = cfg$iter_max),
      retry_control = list(eval.max = cfg$retry_max, iter.max = cfg$retry_max)
    )
    saveRDS(hfit, holdout_fit_path)
    saveRDS(alfak2::summarize_alfak2(hfit, layer = "global"), holdout_summary_path)
    list(
      holdout_status = "ok",
      holdout_error_message = NA_character_,
      holdout_labels = paste(holdout, collapse = ","),
      holdout_n_labels = length(holdout),
      holdout_mode = cfg$holdout_mode,
      holdout_fit_path = holdout_fit_path,
      holdout_landscape_path = holdout_summary_path
    )
  }, error = function(e) {
    list(
      holdout_status = "error",
      holdout_error_message = conditionMessage(e),
      holdout_labels = paste(holdout, collapse = ","),
      holdout_n_labels = length(holdout),
      holdout_mode = cfg$holdout_mode,
      holdout_fit_path = holdout_fit_path,
      holdout_landscape_path = holdout_summary_path
    )
  })
}

run_calibration_fit_task <- function(task, cfg, repo_versions) {
  task <- as.list(task)
  outdir <- task$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  actual_md5 <- unname(tools::md5sum(as.character(task$input_rds)))
  declared_md5 <- as.character(task$input_md5)
  if (length(declared_md5) && !is.na(declared_md5[[1L]]) && nzchar(declared_md5[[1L]]) &&
      !identical(actual_md5, declared_md5[[1L]])) {
    stop("Prepared input_rds md5 changed for task ", task$task_order, ".", call. = FALSE)
  }
  expected_config_hash <- calibration_fit_config_hash(
    task,
    cfg,
    repo_versions,
    actual_input_md5 = actual_md5
  )
  allow_cache <- !isTRUE(cfg$force_refit) &&
    (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfak2"]))
  cached <- if (allow_cache) read_rds_if_exists(result_path) else NULL
  cached_holdout_ok <- is.list(cached) && "holdout_status" %in% names(cached) &&
    (!identical(cached$holdout_status, "ok") || file.exists(cached$holdout_fit_path))
  cached_hash <- if (is.list(cached) && "fit_config_hash" %in% names(cached)) {
    as.character(cached$fit_config_hash[[1L]])
  } else {
    NA_character_
  }
  cached_config_ok <- identical(cached_hash, expected_config_hash)
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$fit_path) &&
      isTRUE(cached_holdout_ok) && isTRUE(cached_config_ok)) {
    cached$cached <- TRUE
    saveRDS(cached, result_path)
    return(cached)
  }

  yi <- readRDS(task$input_rds)
  counts <- prepare_alfak2_counts(
    yi,
    minobs = task$minobs,
    input_policy = task$input_policy,
    drop_diploid = cfg$drop_diploid
  )
  selected_times <- suppressWarnings(as.numeric(colnames(counts)))
  dt <- diff(selected_times)
  if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
  max_cn <- cfg$alfak2_max_cn
  if (!is.finite(max_cn)) {
    k_mat <- parse_karyotype_ids_base(rownames(counts))
    max_cn <- max(k_mat, na.rm = TRUE) + as.integer(task$local_shell_depth) + as.integer(task$global_extra_shell)
  }
  anchor_count_reference <- resolve_anchor_count_reference(cfg, task)

  started <- Sys.time()
  res <- tryCatch({
    set.seed(as.integer(task$benchmark_seed))
    fit <- alfak2::fit_alfak2(
      counts,
      dt = dt,
      beta = cfg$pm,
      min_cn = cfg$alfak2_min_cn,
      max_cn = as.integer(max_cn),
      local_shell_depth = as.integer(task$local_shell_depth),
      global_extra_shell = as.integer(task$global_extra_shell),
      max_nodes = cfg$max_nodes,
      lambda_l_grid = as.numeric(task$lambda_l),
      lambda_e_grid = as.numeric(task$lambda_e),
      sigma_obs_grid = as.numeric(task$sigma_obs),
      graph_edge_weight = as.character(row_field(task, "graph_edge_weight", "mutation")),
      anchor_support_tiers = "directly_informed",
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = cfg$alfak2_anchor_count_power,
      input_depth = cfg$alfak2_input_depth,
      effective_depth = cfg$alfak2_effective_depth,
      effective_depth_mode = as.character(task$effective_depth_mode),
      observation_model = cfg$alfak2_observation_model,
      dm_concentration = as.numeric(task$dm_concentration),
      alfakR_scale = TRUE,
      n0 = cfg$n0,
      nb = cfg$nb,
      correct_efflux = as.logical(task$correct_efflux),
      legacy_weight = as.character(task$legacy_weight),
      control = list(eval.max = cfg$eval_max, iter.max = cfg$iter_max),
      retry_control = list(eval.max = cfg$retry_max, iter.max = cfg$retry_max)
    )
    fit_path <- file.path(outdir, "alfak2_fit.rds")
    summary_path <- file.path(outdir, "landscape.rds")
    saveRDS(fit, fit_path)
    saveRDS(alfak2::summarize_alfak2(fit, layer = "global"), summary_path)
    holdout <- fit_calibration_holdout(task, cfg, counts, dt, max_cn, fit, outdir)
    c(
      task,
      list(
        status = "ok",
        cached = FALSE,
        error_message = NA_character_,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = fit_path,
        landscape_path = summary_path,
        fit_config_hash = expected_config_hash,
        holdout_status = holdout$holdout_status,
        holdout_error_message = holdout$holdout_error_message,
        holdout_labels = holdout$holdout_labels,
        holdout_n_labels = holdout$holdout_n_labels,
        holdout_mode = holdout$holdout_mode,
        holdout_fit_path = holdout$holdout_fit_path,
        holdout_landscape_path = holdout$holdout_landscape_path,
        graph_edge_weight = as.character(row_field(task, "graph_edge_weight", "mutation")),
        anchor_count_reference = if (is.null(anchor_count_reference)) NA_real_ else as.numeric(anchor_count_reference),
        anchor_count_power = cfg$alfak2_anchor_count_power,
        local_convergence = fit$local$diagnostics$convergence,
        local_gradient_norm = fit$local$diagnostics$gradient_norm,
        local_covariance_status = fit$local$diagnostics$covariance_status,
        local_retry_attempted = fit$local$diagnostics$retry_attempted,
        global_factorization_status = fit$global$diagnostics$factorization_status,
        local_nodes = nrow(fit$local$summary),
        global_nodes = nrow(fit$global$summary)
      )
    )
  }, error = function(e) {
    c(
      task,
      list(
        status = "error",
        cached = FALSE,
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = file.path(outdir, "alfak2_fit.rds"),
        landscape_path = file.path(outdir, "landscape.rds"),
        fit_config_hash = expected_config_hash,
        holdout_status = NA_character_,
        holdout_error_message = NA_character_,
        holdout_labels = NA_character_,
        holdout_n_labels = NA_integer_,
        holdout_mode = cfg$holdout_mode,
        holdout_fit_path = file.path(outdir, "alfak2_holdout_fit.rds"),
        holdout_landscape_path = file.path(outdir, "holdout_landscape.rds"),
        graph_edge_weight = as.character(row_field(task, "graph_edge_weight", "mutation")),
        anchor_count_reference = if (is.null(anchor_count_reference)) NA_real_ else as.numeric(anchor_count_reference),
        anchor_count_power = cfg$alfak2_anchor_count_power,
        local_convergence = NA_integer_,
        local_gradient_norm = NA_real_,
        local_covariance_status = NA_character_,
        local_retry_attempted = NA,
        global_factorization_status = NA_character_,
        local_nodes = NA_integer_,
        global_nodes = NA_integer_
      )
    )
  })
  saveRDS(res, result_path)
  res
}

run_calibration_task_table <- function(task_tbl, cfg, repo_versions) {
  if (!nrow(task_tbl)) return(data.frame())
  n_cores <- max(1L, as.integer(cfg$n_cores))
  task_list <- split(task_tbl, seq_len(nrow(task_tbl)))
  message("Running calibration fits: ", length(task_list), " tasks with n_cores=", n_cores, ".")
  rows <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(task_list, run_calibration_fit_task, cfg = cfg, repo_versions = repo_versions, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(task_list, run_calibration_fit_task, cfg = cfg, repo_versions = repo_versions)
  }
  list_to_data_frame(rows)
}

run_calibration_single_task <- function(task_index, cfg, dirs, repo_versions) {
  task_tbl <- read_tsv(file.path(dirs$tables, "task_table.tsv"))
  task_index <- as.integer(task_index)
  if (!is.finite(task_index) || task_index < 1L || task_index > nrow(task_tbl)) {
    stop("`task_index` must be between 1 and ", nrow(task_tbl), ".", call. = FALSE)
  }
  result <- run_calibration_fit_task(task_tbl[task_index, , drop = FALSE], cfg, repo_versions)
  fit_tbl <- list_to_data_frame(list(result))
  part_path <- file.path(dirs$fit_parts, sprintf("task_%06d.tsv", task_index))
  write_tsv(fit_tbl, part_path)
  writeLines("ok", file.path(dirs$fit_parts, sprintf("task_%06d.done", task_index)))
  message("Wrote calibration fit part: ", part_path)
  fit_tbl
}

read_calibration_fit_parts <- function(dirs) {
  part_files <- list.files(dirs$fit_parts, pattern = "^task_[0-9]+\\.tsv$", full.names = TRUE)
  if (!length(part_files)) {
    fallback <- file.path(dirs$tables, "fit_results.tsv")
    if (file.exists(fallback)) return(read_tsv(fallback))
    return(data.frame())
  }
  parts <- lapply(sort(part_files), read_tsv)
  nms <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(nms, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, nms, drop = FALSE]
  })
  out <- do.call(rbind, parts)
  if ("task_order" %in% names(out)) out <- out[order(as.integer(out$task_order)), , drop = FALSE]
  out
}

accuracy_metric_row <- function(task_order, param_id, support_scope, est, tru, n_nodes) {
  ok <- is.finite(est) & is.finite(tru)
  if (!sum(ok)) {
    return(data.frame(
      task_order = as.integer(task_order),
      param_id = as.integer(param_id),
      support_scope = support_scope,
      n_nodes = as.integer(n_nodes),
      n_scored = 0L,
      mae = NA_real_,
      rmse = NA_real_,
      median_abs_error = NA_real_,
      signed_bias = NA_real_,
      centered_mae = NA_real_,
      centered_rmse = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_,
      sign_accuracy = NA_real_,
      false_high_rate = NA_real_,
      estimate_sd = NA_real_,
      truth_sd = NA_real_,
      estimate_sd_ratio = NA_real_,
      estimate_range_ratio = NA_real_,
      estimate_iqr_ratio = NA_real_,
      amplitude_collapse = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  est <- est[ok]
  tru <- tru[ok]
  err <- est - tru
  est_c <- est - mean(est)
  tru_c <- tru - mean(tru)
  centered_err <- est_c - tru_c
  estimate_sd <- if (length(est) >= 2L) stats::sd(est) else NA_real_
  truth_sd <- if (length(tru) >= 2L) stats::sd(tru) else NA_real_
  estimate_sd_ratio <- if (is.finite(truth_sd) && truth_sd > 0) estimate_sd / truth_sd else NA_real_
  estimate_range <- if (length(est) >= 2L) {
    stats::quantile(est, 0.95, na.rm = TRUE, names = FALSE) -
      stats::quantile(est, 0.05, na.rm = TRUE, names = FALSE)
  } else NA_real_
  truth_range <- if (length(tru) >= 2L) {
    stats::quantile(tru, 0.95, na.rm = TRUE, names = FALSE) -
      stats::quantile(tru, 0.05, na.rm = TRUE, names = FALSE)
  } else NA_real_
  estimate_iqr <- if (length(est) >= 2L) stats::IQR(est, na.rm = TRUE) else NA_real_
  truth_iqr <- if (length(tru) >= 2L) stats::IQR(tru, na.rm = TRUE) else NA_real_
  data.frame(
    task_order = as.integer(task_order),
    param_id = as.integer(param_id),
    support_scope = support_scope,
    n_nodes = as.integer(n_nodes),
    n_scored = as.integer(sum(ok)),
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2)),
    median_abs_error = stats::median(abs(err)),
    signed_bias = mean(err),
    centered_mae = mean(abs(centered_err)),
    centered_rmse = sqrt(mean(centered_err^2)),
    pearson = safe_cor(est, tru, method = "pearson"),
    spearman = safe_cor(est, tru, method = "spearman"),
    sign_accuracy = mean(sign(est_c) == sign(tru_c), na.rm = TRUE),
    false_high_rate = mean(est_c > 0 & tru_c <= 0, na.rm = TRUE),
    estimate_sd = estimate_sd,
    truth_sd = truth_sd,
    estimate_sd_ratio = estimate_sd_ratio,
    estimate_range_ratio = if (is.finite(truth_range) && truth_range > 0) estimate_range / truth_range else NA_real_,
    estimate_iqr_ratio = if (is.finite(truth_iqr) && truth_iqr > 0) estimate_iqr / truth_iqr else NA_real_,
    amplitude_collapse = !is.finite(estimate_sd_ratio) || estimate_sd_ratio < 0.02,
    stringsAsFactors = FALSE
  )
}

parse_holdout_labels <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(character(0))
  out <- trimws(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

score_graph_holdout <- function(fit, grf_sim, fr, cfg) {
  holdout <- parse_holdout_labels(fr$holdout_labels)
  if (!length(holdout)) {
    return(accuracy_metric_row(fr$task_order[[1L]], fr$param_id[[1L]], "holdout", numeric(0), numeric(0), 0L))
  }
  holdout_status <- as.character(fr$holdout_status[[1L]])
  if (!identical(holdout_status, "ok")) {
    return(accuracy_metric_row(fr$task_order[[1L]], fr$param_id[[1L]], "holdout", numeric(0), numeric(0), length(holdout)))
  }
  hfit <- read_rds_if_exists(as.character(fr$holdout_fit_path[[1L]]))
  if (is.null(hfit)) {
    return(accuracy_metric_row(fr$task_order[[1L]], fr$param_id[[1L]], "holdout", numeric(0), numeric(0), length(holdout)))
  }
  hs <- alfak2::summarize_alfak2(hfit, layer = "global")
  cols <- calibration_estimate_columns(hs)
  truth <- compute_grf_fitness_truth(holdout, grf_sim$centroids, as.numeric(fr$lambda[[1L]]))
  pred <- hs[[cols$mean]][match(holdout, as.character(hs$karyotype))]
  tru <- as.numeric(truth[holdout])
  accuracy_metric_row(fr$task_order[[1L]], fr$param_id[[1L]], "holdout", as.numeric(pred), tru, length(holdout))
}

score_one_fit <- function(fr, cfg) {
  fit <- read_rds_if_exists(as.character(fr$fit_path[[1L]]))
  grf_sim <- read_rds_if_exists(as.character(fr$grf_rds[[1L]]))
  if (is.null(fit) || is.null(grf_sim)) return(list(metrics = data.frame(), nodes = data.frame()))
  s <- alfak2::summarize_alfak2(fit, layer = "global")
  cols <- calibration_estimate_columns(s)
  truth <- compute_grf_fitness_truth(s$karyotype, grf_sim$centroids, as.numeric(fr$lambda[[1L]]))
  support_scope <- support_scope_for_node("alfak2", as.character(s$support_tier), FALSE, FALSE)
  node_tbl <- data.frame(
    task_order = as.integer(fr$task_order[[1L]]),
    param_id = as.integer(fr$param_id[[1L]]),
    simulation_id = as.integer(fr$simulation_id[[1L]]),
    lambda = as.numeric(fr$lambda[[1L]]),
    lambda_label = as.character(fr$lambda_label[[1L]]),
    time_start = as.numeric(fr$time_start[[1L]]),
    time_gap = as.numeric(fr$time_gap[[1L]]),
    time_delta = as.numeric(fr$time_delta[[1L]]),
    minobs = as.integer(fr$minobs[[1L]]),
    k = as.character(s$karyotype),
    support_tier = as.character(s$support_tier),
    support_scope = support_scope,
    estimated_fitness = as.numeric(s[[cols$mean]]),
    estimated_sd = as.numeric(s[[cols$sd]]),
    true_fitness = as.numeric(truth[match(as.character(s$karyotype), names(truth))]),
    stringsAsFactors = FALSE
  )
  node_tbl$estimation_error <- node_tbl$estimated_fitness - node_tbl$true_fitness

  scopes <- c("all", "direct", "nn")
  metrics <- lapply(scopes, function(scope) {
    x <- if (identical(scope, "all")) node_tbl else node_tbl[node_tbl$support_scope == scope, , drop = FALSE]
    accuracy_metric_row(
      fr$task_order[[1L]],
      fr$param_id[[1L]],
      scope,
      as.numeric(x$estimated_fitness),
      as.numeric(x$true_fitness),
      nrow(x)
    )
  })
  metrics <- do.call(rbind, metrics)
  metrics <- rbind(metrics, score_graph_holdout(fit, grf_sim, fr, cfg))
  meta_cols <- c(
    "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs", "param_id", "legacy_weight", "correct_efflux", "lambda_l", "lambda_e",
    "sigma_obs", "graph_edge_weight", "anchor_count_reference", "anchor_count_power",
    "dm_concentration", "effective_depth_mode", "local_shell_depth", "global_extra_shell",
    "local_covariance_status", "global_factorization_status",
    "holdout_status", "holdout_n_labels", "holdout_mode"
  )
  meta <- fr[rep(1L, nrow(metrics)), intersect(meta_cols, names(fr)), drop = FALSE]
  metrics <- cbind(meta, metrics[setdiff(names(metrics), names(meta))])
  list(metrics = metrics, nodes = node_tbl)
}

summarize_calibration_results <- function(cfg, dirs) {
  task_tbl <- read_tsv(file.path(dirs$tables, "task_table.tsv"))
  fit_tbl <- read_calibration_fit_parts(dirs)
  if (!nrow(fit_tbl)) stop("No calibration fit results found.", call. = FALSE)
  if ("task_order" %in% names(fit_tbl)) fit_tbl <- fit_tbl[order(as.integer(fit_tbl$task_order)), , drop = FALSE]
  write_tsv(fit_tbl, file.path(dirs$tables, "fit_results.tsv"))

  completed <- unique(as.integer(fit_tbl$task_order))
  missing <- task_tbl[!(as.integer(task_tbl$task_order) %in% completed), , drop = FALSE]
  write_tsv(missing, file.path(dirs$tables, "missing_fit_tasks.tsv"))

  ok_fit <- fit_tbl[as.character(fit_tbl$status) == "ok", , drop = FALSE]
  scored <- lapply(seq_len(nrow(ok_fit)), function(i) score_one_fit(ok_fit[i, , drop = FALSE], cfg))
  metric_tbl <- if (length(scored)) do.call(rbind, lapply(scored, `[[`, "metrics")) else data.frame()
  write_tsv(metric_tbl, file.path(dirs$tables, "calibration_metrics_by_fit.tsv"))
  if (isTRUE(cfg$write_node_table)) {
    node_tbl <- if (length(scored)) do.call(rbind, lapply(scored, `[[`, "nodes")) else data.frame()
    write_tsv(node_tbl, file.path(dirs$tables, "calibration_node_accuracy.tsv"))
  }

  ranked <- rank_calibration_parameters(metric_tbl, fit_tbl, cfg)
  write_tsv(ranked, file.path(dirs$tables, "calibration_ranked_params.tsv"))
  best <- ranked[ranked$rank == 1L, , drop = FALSE]
  write_tsv(best, file.path(dirs$tables, "best_params.tsv"))
  if (nrow(best)) {
    writeLines(best_param_cli(best), file.path(dirs$tables, "best_params_cli_args.txt"))
  }
  saveRDS(
    list(config = cfg, task_table = task_tbl, fit_results = fit_tbl, metrics = metric_tbl, ranked = ranked, best = best),
    file.path(dirs$root, "grf_alfak2_parameter_calibration.rds")
  )
  message("Wrote calibration summaries under: ", dirs$tables)
  invisible(list(fit_results = fit_tbl, metrics = metric_tbl, ranked = ranked, best = best, missing = missing))
}

rank_calibration_parameters <- function(metric_tbl, fit_tbl, cfg) {
  if (!nrow(metric_tbl)) return(data.frame())
  metric_names <- c(
    "mae", "rmse", "centered_mae", "centered_rmse", "median_abs_error",
    "signed_bias", "false_high_rate", "sign_accuracy", "spearman",
    "estimate_sd_ratio", "estimate_range_ratio", "estimate_iqr_ratio"
  )
  for (nm in intersect(metric_names, names(metric_tbl))) metric_tbl[[nm]] <- to_num(metric_tbl[[nm]])
  for (nm in c("estimate_sd_ratio", "estimate_range_ratio", "estimate_iqr_ratio")) {
    if (!nm %in% names(metric_tbl)) metric_tbl[[nm]] <- NA_real_
  }
  if (!"amplitude_collapse" %in% names(metric_tbl)) metric_tbl$amplitude_collapse <- NA
  metric_tbl$param_id <- as.integer(metric_tbl$param_id)
  fit_tbl$param_id <- as.integer(fit_tbl$param_id)

  param_cols <- c(
    "param_id", "input_policy", "alfak2_input_depth", "alfak2_effective_depth",
    "alfak2_observation_model", "legacy_weight", "correct_efflux", "lambda_l", "lambda_e", "sigma_obs",
    "graph_edge_weight", "anchor_count_reference", "anchor_count_power", "dm_concentration",
    "effective_depth_mode", "local_shell_depth", "global_extra_shell"
  )
  param_tbl <- unique(fit_tbl[, intersect(param_cols, names(fit_tbl)), drop = FALSE])
  status_tbl <- aggregate(
    list(
      n_fit_tasks = fit_tbl$param_id,
      n_ok = as.character(fit_tbl$status) == "ok",
      n_error = as.character(fit_tbl$status) != "ok"
    ),
    by = list(param_id = fit_tbl$param_id),
    FUN = function(x) if (is.logical(x)) sum(x, na.rm = TRUE) else length(x)
  )
  if ("holdout_status" %in% names(fit_tbl)) {
    holdout_status <- as.character(fit_tbl$holdout_status)
    holdout_status[is.na(holdout_status) | !nzchar(holdout_status)] <- "missing"
    holdout_status_tbl <- aggregate(
      list(
        n_holdout_ok = holdout_status == "ok",
        n_holdout_error = holdout_status == "error",
        n_holdout_skipped = grepl("^skipped", holdout_status),
        n_holdout_missing = holdout_status == "missing"
      ),
      by = list(param_id = fit_tbl$param_id),
      FUN = function(x) sum(x, na.rm = TRUE)
    )
  } else {
    holdout_status_tbl <- data.frame(param_id = integer())
  }

  summarize_scope <- function(scope, prefix) {
    x <- metric_tbl[metric_tbl$support_scope == scope, , drop = FALSE]
    if (!nrow(x)) return(data.frame(param_id = integer()))
    split_x <- split(x, x$param_id)
    rows <- lapply(split_x, function(z) {
      data.frame(
        param_id = as.integer(z$param_id[[1L]]),
        setNames(list(nrow(z)), paste0(prefix, "_n_conditions")),
        setNames(list(sum(z$n_scored > 0, na.rm = TRUE)), paste0(prefix, "_n_scored_conditions")),
        setNames(list(sum(z$n_scored, na.rm = TRUE)), paste0(prefix, "_total_scored")),
        setNames(list(if (cfg$objective_metric == "sparse_composite") NA_real_ else median(z[[cfg$objective_metric]], na.rm = TRUE)), paste0(prefix, "_median_", cfg$objective_metric)),
        setNames(list(if (cfg$objective_metric == "sparse_composite") NA_real_ else stats::quantile(z[[cfg$objective_metric]], 0.75, na.rm = TRUE, names = FALSE)), paste0(prefix, "_q75_", cfg$objective_metric)),
        setNames(list(median(z$centered_rmse, na.rm = TRUE)), paste0(prefix, "_median_centered_rmse")),
        setNames(list(median(z$centered_mae, na.rm = TRUE)), paste0(prefix, "_median_centered_mae")),
        setNames(list(median(abs(z$signed_bias), na.rm = TRUE)), paste0(prefix, "_median_abs_bias")),
        setNames(list(median(z$spearman, na.rm = TRUE)), paste0(prefix, "_median_spearman")),
        setNames(list(median(z$false_high_rate, na.rm = TRUE)), paste0(prefix, "_median_false_high_rate")),
        setNames(list(median(z$sign_accuracy, na.rm = TRUE)), paste0(prefix, "_median_sign_accuracy")),
        setNames(list(median(z$estimate_sd_ratio, na.rm = TRUE)), paste0(prefix, "_median_estimate_sd_ratio")),
        setNames(list(median(z$estimate_range_ratio, na.rm = TRUE)), paste0(prefix, "_median_estimate_range_ratio")),
        setNames(list(median(z$estimate_iqr_ratio, na.rm = TRUE)), paste0(prefix, "_median_estimate_iqr_ratio")),
        setNames(list(mean(z$amplitude_collapse, na.rm = TRUE)), paste0(prefix, "_amplitude_collapse_fraction")),
        check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }

  objective_tbl <- summarize_scope(cfg$objective_scope, "objective")
  direct_tbl <- summarize_scope("direct", "direct")
  holdout_tbl <- summarize_scope("holdout", "holdout")
  all_tbl <- summarize_scope("all", "all")
  out <- merge(param_tbl, status_tbl, by = "param_id", all.x = TRUE)
  out <- merge(out, holdout_status_tbl, by = "param_id", all.x = TRUE)
  out <- merge(out, objective_tbl, by = "param_id", all.x = TRUE)
  out <- merge(out, direct_tbl, by = "param_id", all.x = TRUE)
  out <- merge(out, holdout_tbl, by = "param_id", all.x = TRUE)
  out <- merge(out, all_tbl, by = "param_id", all.x = TRUE)

  objective_col <- paste0("objective_median_", cfg$objective_metric)
  direct_col <- paste0("direct_median_", cfg$objective_metric)
  bias_col <- "objective_median_abs_bias"
  if (!objective_col %in% names(out)) out[[objective_col]] <- NA_real_
  if (!direct_col %in% names(out)) out[[direct_col]] <- NA_real_
  if (!bias_col %in% names(out)) out[[bias_col]] <- NA_real_
  for (nm in c(
    "objective_median_centered_rmse", "direct_median_centered_rmse",
    "holdout_median_centered_rmse",
    "objective_median_spearman", "objective_median_false_high_rate",
    "objective_median_estimate_sd_ratio"
  )) {
    if (!nm %in% names(out)) out[[nm]] <- NA_real_
  }
  for (nm in c(
    "n_ok", "n_holdout_ok", "n_holdout_error", "n_holdout_skipped", "n_holdout_missing",
    "holdout_n_conditions", "holdout_n_scored_conditions", "holdout_total_scored"
  )) {
    if (!nm %in% names(out)) out[[nm]] <- 0
    out[[nm]][!is.finite(out[[nm]])] <- 0
  }
  holdout_expected <- pmax(as.numeric(out$n_ok), 1)
  holdout_unscored_fraction <- pmax(
    0,
    1 - as.numeric(out$holdout_n_scored_conditions) / holdout_expected
  )
  holdout_bad_status_fraction <- pmax(
    0,
    (as.numeric(out$n_holdout_error) + as.numeric(out$n_holdout_skipped) + as.numeric(out$n_holdout_missing)) /
      holdout_expected
  )
  holdout_no_scores <- !is.finite(out$holdout_total_scored) | as.numeric(out$holdout_total_scored) <= 0
  out$holdout_failure_fraction <- pmax(holdout_unscored_fraction, holdout_bad_status_fraction, as.numeric(holdout_no_scores))
  out$holdout_failure_penalty_value <- as.numeric(cfg$holdout_failure_penalty) * out$holdout_failure_fraction
  if (cfg$objective_metric == "sparse_composite") {
    objective_amplitude_penalty <- abs(log(pmax(
      ifelse(is.finite(out$objective_median_estimate_sd_ratio), out$objective_median_estimate_sd_ratio, 1e-4),
      1e-4
    )))
    out$calibration_score <- out$objective_median_centered_rmse +
      cfg$direct_weight * ifelse(is.finite(out$direct_median_centered_rmse), out$direct_median_centered_rmse, 0) +
      cfg$holdout_weight * ifelse(is.finite(out$holdout_median_centered_rmse), out$holdout_median_centered_rmse, 0) +
      cfg$bias_weight * ifelse(is.finite(out[[bias_col]]), out[[bias_col]], 0) -
      cfg$spearman_weight * ifelse(is.finite(out$objective_median_spearman), out$objective_median_spearman, 0) +
      cfg$false_high_weight * ifelse(is.finite(out$objective_median_false_high_rate), out$objective_median_false_high_rate, 0) +
      0.25 * objective_amplitude_penalty +
      out$holdout_failure_penalty_value
    out$objective_amplitude_penalty <- objective_amplitude_penalty
    objective_order_col <- "objective_median_centered_rmse"
  } else {
    out$calibration_score <- out[[objective_col]] +
      cfg$direct_weight * ifelse(is.finite(out[[direct_col]]), out[[direct_col]], 0) +
      cfg$holdout_weight * ifelse(is.finite(out$holdout_median_centered_rmse), out$holdout_median_centered_rmse, 0) +
      cfg$bias_weight * ifelse(is.finite(out[[bias_col]]), out[[bias_col]], 0) +
      out$holdout_failure_penalty_value
    objective_order_col <- objective_col
  }
  out$objective_scope <- cfg$objective_scope
  out$objective_metric <- cfg$objective_metric
  out$direct_weight <- cfg$direct_weight
  out$holdout_weight <- cfg$holdout_weight
  out$holdout_failure_penalty <- cfg$holdout_failure_penalty
  out$bias_weight <- cfg$bias_weight
  out$spearman_weight <- cfg$spearman_weight
  out$false_high_weight <- cfg$false_high_weight
  out <- out[order(out$calibration_score, out[[objective_order_col]], out$param_id), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  out
}

best_param_cli <- function(best) {
  if (!nrow(best)) return(character(0))
  b <- best[1L, , drop = FALSE]
  args <- c(
    sprintf("--alfak2-input-policies=%s", b$input_policy),
    sprintf("--alfak2-input-depth=%s", b$alfak2_input_depth),
    sprintf("--alfak2-legacy-weight=%s", b$legacy_weight),
    sprintf("--correct-efflux=%s", tolower(as.character(b$correct_efflux))),
    sprintf("--alfak2-lambda-l-grid=%s", b$lambda_l),
    sprintf("--alfak2-lambda-e-grid=%s", b$lambda_e),
    sprintf("--alfak2-sigma-obs-grid=%s", b$sigma_obs),
    sprintf("--alfak2-graph-edge-weight=%s", b$graph_edge_weight),
    sprintf("--alfak2-dm-concentration=%s", b$dm_concentration),
    sprintf("--alfak2-effective-depth-mode=%s", b$effective_depth_mode),
    sprintf("--alfak2-local-shell-depth=%s", b$local_shell_depth),
    sprintf("--alfak2-global-extra-shell=%s", b$global_extra_shell)
  )
  if ("alfak2_observation_model" %in% names(b) && nzchar(as.character(b$alfak2_observation_model))) {
    args <- c(args, sprintf("--alfak2-observation-model=%s", b$alfak2_observation_model))
  }
  if ("alfak2_effective_depth" %in% names(b) && is.finite(suppressWarnings(as.numeric(b$alfak2_effective_depth)))) {
    args <- c(args, sprintf("--alfak2-effective-depth=%s", b$alfak2_effective_depth))
  }
  if ("anchor_count_reference" %in% names(b) && is.finite(suppressWarnings(as.numeric(b$anchor_count_reference)))) {
    args <- c(args, sprintf("--alfak2-anchor-count-reference=%s", b$anchor_count_reference))
    args <- c(args, sprintf("--alfak2-anchor-count-power=%s", b$anchor_count_power))
  }
  paste(args, collapse = "\n")
}

run_mode_fit_all <- function(cfg, dirs, repo_versions) {
  task_tbl <- read_tsv(file.path(dirs$tables, "task_table.tsv"))
  fit_tbl <- run_calibration_task_table(task_tbl, cfg, repo_versions)
  write_tsv(fit_tbl, file.path(dirs$tables, "fit_results.tsv"))
  fit_tbl
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  repo_dir <- normalizePath(arg_value(args, "repo_dir", resolve_repo_dir()), winslash = "/", mustWork = FALSE)
  cfg <- build_calibration_config(args, repo_dir)
  dirs <- build_calibration_dirs(cfg$output_dir)
  valid_modes <- c("all", "prepare", "fit", "fit-task", "fit_task", "summarize")
  if (!cfg$mode %in% valid_modes) {
    stop("Unsupported --mode=", cfg$mode, ". Use all, prepare, fit, fit-task, or summarize.", call. = FALSE)
  }

  message("Loading current source trees with pkgload::load_all().")
  message("  alfakR: ", cfg$alfakR_repo)
  message("  alfak2: ", cfg$alfak2_repo)
  load_current_repos(cfg$alfakR_repo, cfg$alfak2_repo, recompile_dll = cfg$recompile_dll)
  repo_versions <- rbind(repo_state(cfg$alfakR_repo, "alfakR"), repo_state(cfg$alfak2_repo, "alfak2"))
  write_tsv(repo_versions, file.path(dirs$tables, "repo_versions.tsv"))

  if (cfg$mode %in% c("all", "prepare")) {
    prepare_calibration_inputs(cfg, dirs, repo_versions)
  }
  if (cfg$mode %in% c("all", "fit")) {
    run_mode_fit_all(cfg, dirs, repo_versions)
  }
  if (cfg$mode %in% c("fit-task", "fit_task")) {
    run_calibration_single_task(cfg$task_index, cfg, dirs, repo_versions)
  }
  if (cfg$mode %in% c("all", "summarize")) {
    summarize_calibration_results(cfg, dirs)
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  main()
}
