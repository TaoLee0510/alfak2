#!/usr/bin/env Rscript

find_repo_root <- function(start = getwd()) {
  cand <- normalizePath(file.path(start, c(".", "..", "../..", "../../..")), winslash = "/", mustWork = FALSE)
  for (x in unique(cand)) {
    if (file.exists(file.path(x, "DESCRIPTION")) && dir.exists(file.path(x, "benchmark", "scr"))) return(x)
  }
  stop("Could not locate repository root.", call. = FALSE)
}

repo_dir0 <- find_repo_root()
helper_path <- file.path(repo_dir0, "benchmark", "scr", "run_hybrid_alfak2_direct_alfakR_nn_benchmark.R")
if (!file.exists(helper_path)) stop("Missing helper benchmark script: ", helper_path, call. = FALSE)
suppressMessages(source(helper_path, local = .GlobalEnv))

usage <- function() {
  cat(
    "alfak2-native direct -> NN/Kriging GRF benchmark\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_alfak2_native_nn_kriging_benchmark.R --mode=all [options]\n\n",
    "Modes: prepare, fit, fit-task, summarize, all\n",
    "Core options mirror the hybrid runner, with defaults:\n",
    "  --output-dir=benchmark/results/alfak2_native_nn_kriging\n",
    "  --lambdas=0.2,0.6,0.8\n",
    "  --nn-priors=none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step\n",
    "  --n-repeats=3 --nboot=20 --n-cores=9\n",
    "Native options:\n",
    "  --nn-shell-depth=1\n",
    "  --kriging-max-anchors=300\n",
    "  --kriging-nugget=1e-4\n",
    sep = ""
  )
}

build_native_config <- function(args, repo_dir) {
  cfg <- build_hybrid_config(args, repo_dir)
  if (is.null(arg_value(args, "output_dir", NULL))) {
    cfg$output_dir <- normalize_output_dir(repo_dir, "benchmark/results/alfak2_native_nn_kriging")
  }
  cfg$nn_shell_depth <- arg_integer(args, "nn_shell_depth", 1L)
  cfg$kriging_max_anchors <- arg_integer(args, "kriging_max_anchors", 300L)
  cfg$kriging_nugget <- arg_numeric(args, "kriging_nugget", 1e-4)
  cfg$native_global_extra_shell <- arg_integer(args, "native_global_extra_shell", 2L)
  cfg$native_backend <- "alfak2_native_nn_kriging"
  cfg
}

build_native_dirs <- function(output_dir) {
  dirs <- list(root = output_dir, cache = file.path(output_dir, "cache"), fits = file.path(output_dir, "fits"),
               tables = file.path(output_dir, "tables"), results = file.path(output_dir, "results"),
               fit_parts = file.path(output_dir, "tables", "fit_results_parts"))
  for (p in dirs) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  dirs
}

build_native_task_table <- function(input_tbl, cfg, dirs) {
  rows <- list(); idx <- 0L
  base_rows <- input_tbl[!duplicated(input_tbl[c("simulation_id", "lambda", "minobs", "input_rds")]), , drop = FALSE]
  for (i in seq_len(nrow(base_rows))) {
    r <- base_rows[i, , drop = FALSE]
    for (repeat_id in cfg$repeat_ids) {
      repeat_dir <- paste0("repeat_", repeat_id)
      lambda_dir <- paste0("lambda_", as.character(r$lambda_label))
      base_seed <- as.integer(cfg$seed + as.integer(r$simulation_id) * 100000L +
                                match(as.numeric(r$lambda), cfg$lambdas) * 10000L +
                                as.integer(r$minobs) * 100L + as.integer(repeat_id) * 1000000L)
      for (prior in cfg$nn_priors) {
        idx <- idx + 1L
        method <- paste0("alfakR_baseline_", prior)
        rows[[idx]] <- data.frame(engine = "alfakR_baseline", method = method, input_policy = "alfakR_minobs_internal",
          nn_prior = prior, simulation_id = as.integer(r$simulation_id), lambda = as.numeric(r$lambda), lambda_label = as.character(r$lambda_label),
          repeat_id = as.integer(repeat_id), time_start = as.numeric(r$time_start), time_gap = as.numeric(r$time_gap),
          time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)), pm = cfg$pm,
          fit_beta_label = pm_to_label(cfg$pm), patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key),
          grf_rds = as.character(r$grf_rds), input_rds = as.character(r$input_rds),
          input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))),
          minobs = as.integer(r$minobs), benchmark_seed = as.integer(base_seed + match(prior, cfg$nn_priors)),
          outdir = file.path(dirs$fits, method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir),
          stringsAsFactors = FALSE)
      }
      for (policy in cfg$input_policies) {
        idx <- idx + 1L
        graph_method <- paste0("alfak2_graphgp_", method_policy_token(policy))
        rows[[idx]] <- data.frame(engine = "alfak2_graphgp", method = graph_method, input_policy = policy, nn_prior = NA_character_,
          simulation_id = as.integer(r$simulation_id), lambda = as.numeric(r$lambda), lambda_label = as.character(r$lambda_label),
          repeat_id = as.integer(repeat_id), time_start = as.numeric(r$time_start), time_gap = as.numeric(r$time_gap),
          time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)), pm = cfg$pm,
          fit_beta_label = pm_to_label(cfg$pm), patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key),
          grf_rds = as.character(r$grf_rds), input_rds = as.character(r$input_rds),
          input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))),
          minobs = as.integer(r$minobs), benchmark_seed = as.integer(base_seed + 500L + match(policy, cfg$input_policies)),
          outdir = file.path(dirs$fits, graph_method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir),
          stringsAsFactors = FALSE)
        for (prior in cfg$nn_priors) {
          idx <- idx + 1L
          method <- paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_", prior)
          rows[[idx]] <- data.frame(engine = "alfak2_native_nn_kriging", method = method, input_policy = policy, nn_prior = prior,
            simulation_id = as.integer(r$simulation_id), lambda = as.numeric(r$lambda), lambda_label = as.character(r$lambda_label),
            repeat_id = as.integer(repeat_id), time_start = as.numeric(r$time_start), time_gap = as.numeric(r$time_gap),
            time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)), pm = cfg$pm,
            fit_beta_label = pm_to_label(cfg$pm), patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key),
            grf_rds = as.character(r$grf_rds), input_rds = as.character(r$input_rds),
            input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))),
            minobs = as.integer(r$minobs),
            benchmark_seed = as.integer(base_seed + 700L + match(policy, cfg$input_policies) * 10L + match(prior, cfg$nn_priors)),
            outdir = file.path(dirs$fits, method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$task_order <- seq_len(nrow(out))
  out
}

native_load_or_prepare_inputs <- function(cfg, dirs) {
  input_tbl <- load_source_inputs(cfg)
  if (!nrow(input_tbl)) input_tbl <- prepare_generated_inputs(cfg, dirs)
  input_tbl <- input_tbl[input_tbl$simulation_id %in% cfg$simulation_ids & input_tbl$minobs %in% cfg$minobs &
                           numeric_in(input_tbl$lambda, cfg$lambdas), , drop = FALSE]
  input_tbl <- input_tbl[order(input_tbl$lambda, input_tbl$simulation_id, input_tbl$minobs), , drop = FALSE]
  task_tbl <- build_native_task_table(input_tbl, cfg, dirs)
  write_tsv_safe(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv_safe(task_tbl, file.path(dirs$tables, "native_task_table.tsv"))
  write_tsv_safe(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  saveRDS(cfg, file.path(dirs$root, "benchmark_config.rds"))
  list(input_table = input_tbl, task_table = task_tbl)
}

native_counts_for_task <- function(task, cfg) {
  yi <- readRDS(as.character(task$input_rds))
  prepare_counts_for_policy(yi, as.integer(task$minobs), as.character(task$input_policy), cfg$drop_diploid)
}

native_max_cn <- function(counts, cfg) {
  max_cn <- cfg$alfak2_max_cn
  if (!is.finite(max_cn)) {
    k_mat <- parse_karyotype_ids_base(rownames(counts))
    max_cn <- max(k_mat, na.rm = TRUE) + cfg$native_global_extra_shell
  }
  as.integer(max_cn)
}

native_dt <- function(counts, task) {
  selected_times <- suppressWarnings(as.numeric(colnames(counts)))
  dt <- if (length(selected_times) >= 2L && all(is.finite(selected_times))) diff(selected_times[1:2]) else as.numeric(task$time_delta)
  if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
  dt
}

graphgp_landscape <- function(fit) {
  s <- fit$global$summary
  data.frame(k = as.character(s$karyotype), mean = as.numeric(s$fitness_mean), sd = as.numeric(s$fitness_sd),
             fq = as.character(s$support_tier) == "directly_informed" | as.integer(s$support_distance) == 0L,
             nn = FALSE, support_scope_native = as.character(s$support_tier), stringsAsFactors = FALSE)
}

native_landscape <- function(fit) {
  s <- fit$summary
  data.frame(k = as.character(s$karyotype), mean = as.numeric(s$fitness_mean), sd = as.numeric(s$fitness_sd),
             fq = as.logical(s$fq), nn = as.logical(s$nn), support_scope_native = as.character(s$support_scope),
             stringsAsFactors = FALSE)
}

run_native_graphgp_fit <- function(task, cfg) {
  outdir <- task$outdir; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  cached <- if (!isTRUE(cfg$force_refit)) safe_read_rds(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) { cached$cached <- TRUE; return(cached) }
  counts <- native_counts_for_task(task, cfg)
  counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
  dt <- native_dt(counts, task)
  started <- Sys.time(); warnings <- character()
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    cap <- capture_warnings(alfak2::fit_alfak2(
      counts, dt = dt, beta = as.numeric(task$pm), min_cn = cfg$alfak2_min_cn, max_cn = native_max_cn(counts, cfg),
      local_shell_depth = 0, global_extra_shell = cfg$native_global_extra_shell, max_nodes = cfg$alfak2_max_nodes,
      lambda_l_grid = cfg$alfak2_lambda_l_grid, lambda_e_grid = cfg$alfak2_lambda_e_grid,
      sigma_obs_grid = cfg$alfak2_sigma_obs_grid, graph_edge_weight = cfg$graph_edge_weight,
      anchor_support_tiers = "directly_informed", input_depth = cfg$alfak2_input_depth,
      effective_depth_mode = cfg$alfak2_effective_depth_mode, observation_model = cfg$alfak2_observation_model,
      dm_concentration = cfg$alfak2_dm_concentration,
      control = list(eval.max = cfg$alfak2_eval_max, iter.max = cfg$alfak2_iter_max),
      retry_control = list(eval.max = cfg$alfak2_retry_max, iter.max = cfg$alfak2_retry_max)))
    fit <- cap$value; warnings <- cap$warnings
    saveRDS(fit, file.path(outdir, "alfak2_graphgp_fit.rds"))
    saveRDS(graphgp_landscape(fit), file.path(outdir, "landscape.Rds"))
    saveRDS(matrix(numeric(0), nrow = 0, ncol = 0), file.path(outdir, "landscape_posterior_samples.Rds"))
    list(status = "ok", cached = FALSE, error_message = NA_character_, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = NA_character_,
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_graphgp_fit.rds"),
      nboot_success = 0L, nboot_failed = 0L, local_convergence = fit$local$diagnostics$convergence,
      local_gradient_norm = fit$local$diagnostics$gradient_norm, local_covariance_status = fit$local$diagnostics$covariance_status,
      scale_status = "alfak2_native", bridge_status = "not_applicable", nn_status = "not_run", kriging_status = "graphgp")
  }, error = function(e) {
    list(status = "error", cached = FALSE, error_message = conditionMessage(e), elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = NA_character_,
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_graphgp_fit.rds"),
      nboot_success = 0L, nboot_failed = 0L, local_convergence = NA_integer_, local_gradient_norm = NA_real_,
      local_covariance_status = NA_character_, scale_status = "error", bridge_status = "not_applicable", nn_status = "not_run", kriging_status = "error")
  })
  out <- c(as.list(task), res); saveRDS(out, result_path); out
}

run_native_nn_kriging_fit <- function(task, cfg) {
  outdir <- task$outdir; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  cached <- if (!isTRUE(cfg$force_refit)) safe_read_rds(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) { cached$cached <- TRUE; return(cached) }
  counts <- native_counts_for_task(task, cfg)
  counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
  dt <- native_dt(counts, task)
  started <- Sys.time(); warnings <- character()
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    cap <- capture_warnings(alfak2:::fit_alfak2_nn_kriging_experimental(
      counts, dt = dt, beta = as.numeric(task$pm), min_cn = cfg$alfak2_min_cn, max_cn = native_max_cn(counts, cfg),
      local_shell_depth = 0, global_extra_shell = cfg$native_global_extra_shell, max_nodes = cfg$alfak2_max_nodes,
      lambda_l_grid = cfg$alfak2_lambda_l_grid, lambda_e_grid = cfg$alfak2_lambda_e_grid,
      sigma_obs_grid = cfg$alfak2_sigma_obs_grid, graph_edge_weight = cfg$graph_edge_weight,
      input_depth = cfg$alfak2_input_depth, effective_depth_mode = cfg$alfak2_effective_depth_mode,
      observation_model = cfg$alfak2_observation_model, dm_concentration = cfg$alfak2_dm_concentration,
      nn_prior = as.character(task$nn_prior), nn_shell_depth = cfg$nn_shell_depth,
      kriging_nugget = cfg$kriging_nugget, kriging_max_anchors = cfg$kriging_max_anchors,
      nboot = cfg$nboot, seed = task$benchmark_seed,
      control = list(eval.max = cfg$alfak2_eval_max, iter.max = cfg$alfak2_iter_max),
      retry_control = list(eval.max = cfg$alfak2_retry_max, iter.max = cfg$alfak2_retry_max)))
    fit <- cap$value; warnings <- cap$warnings
    saveRDS(fit, file.path(outdir, "alfak2_native_nn_kriging_fit.rds"))
    saveRDS(native_landscape(fit), file.path(outdir, "landscape.Rds"))
    saveRDS(fit$posterior_samples, file.path(outdir, "landscape_posterior_samples.Rds"))
    saveRDS(list(direct = fit$direct_state$fitness_boot, nn = fit$nn$fitness_boot, posterior = fit$posterior_samples),
            file.path(outdir, "bootstrap_res.Rds"))
    saveRDS(list(replicate = fit$direct_state$diagnostics$bootstrap, node = fit$nn$diagnostics),
            file.path(outdir, "nn_prior_diagnostics.Rds"))
    boot_diag <- fit$direct_state$diagnostics$bootstrap
    nboot_success <- if (is.data.frame(boot_diag) && nrow(boot_diag)) sum(boot_diag$status == "ok") else 0L
    list(status = "ok", cached = FALSE, error_message = NA_character_, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_native_nn_kriging_fit.rds"),
      nboot_success = nboot_success, nboot_failed = cfg$nboot - nboot_success,
      local_convergence = fit$local$diagnostics$convergence, local_gradient_norm = fit$local$diagnostics$gradient_norm,
      local_covariance_status = fit$local$diagnostics$covariance_status, scale_status = "alfak2_native",
      bridge_status = "not_applicable", nn_status = "ok", kriging_status = fit$kriging$diagnostics$solve_status)
  }, error = function(e) {
    list(status = "error", cached = FALSE, error_message = conditionMessage(e), elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_native_nn_kriging_fit.rds"),
      nboot_success = 0L, nboot_failed = cfg$nboot, local_convergence = NA_integer_, local_gradient_norm = NA_real_,
      local_covariance_status = NA_character_, scale_status = "error", bridge_status = "not_applicable", nn_status = "error", kriging_status = "error")
  })
  out <- c(as.list(task), res); saveRDS(out, result_path); out
}

run_native_one_task <- function(task, cfg) {
  task <- as.list(task)
  if (identical(task$engine, "alfakR_baseline")) return(run_alfakR_baseline_fit(task, cfg))
  if (identical(task$engine, "alfak2_graphgp")) return(run_native_graphgp_fit(task, cfg))
  if (identical(task$engine, "alfak2_native_nn_kriging")) return(run_native_nn_kriging_fit(task, cfg))
  stop("Unknown task engine: ", task$engine, call. = FALSE)
}

run_native_fit_tasks <- function(task_tbl, cfg) {
  n_cores <- max(1L, min(as.integer(cfg$n_cores), nrow(task_tbl)))
  lst <- lapply(seq_len(nrow(task_tbl)), function(i) task_tbl[i, , drop = FALSE])
  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(lst, run_native_one_task, cfg = cfg, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(lst, run_native_one_task, cfg = cfg)
  }
}

native_read_fit_results <- function(dirs) {
  parts <- list.files(dirs$fit_parts, pattern = "^task_[0-9]+\\.tsv$", full.names = TRUE)
  out <- list()
  p <- file.path(dirs$tables, "native_fit_results.tsv")
  if (file.exists(p)) out[[length(out) + 1L]] <- read_tsv(p)
  if (length(parts)) out[[length(out) + 1L]] <- read_fit_result_parts(dirs)
  out <- Filter(function(x) is.data.frame(x) && nrow(x), out)
  if (!length(out)) return(data.frame())
  x <- do.call(rbind, out)
  if ("task_order" %in% names(x)) {
    x$.row_priority <- ifelse(as.character(x$status) == "ok", 2L, 1L)
    x <- x[order(as.integer(x$task_order), x$.row_priority), , drop = FALSE]
    x <- x[!duplicated(x$task_order, fromLast = TRUE), , drop = FALSE]
    x$.row_priority <- NULL
  }
  x
}

native_landscape_to_nodes <- function(fr, cfg) {
  grf <- safe_read_rds(as.character(fr$grf_rds[[1]])); if (is.null(grf)) return(data.frame())
  land <- safe_read_rds(as.character(fr$landscape_path[[1]])); if (is.null(land)) return(data.frame())
  x <- as.data.frame(land, stringsAsFactors = FALSE)
  if (!"k" %in% names(x) && "karyotype" %in% names(x)) x$k <- x$karyotype
  if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x$mean <- x$fitness_mean
  if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x$sd <- x$fitness_sd
  if (!"sd" %in% names(x)) x$sd <- NA_real_
  if (!"fq" %in% names(x)) x$fq <- FALSE
  if (!"nn" %in% names(x)) x$nn <- FALSE
  truth <- compute_grf_fitness_truth(as.character(x$k), grf$centroids, as.numeric(fr$lambda[[1]]))
  scope <- ifelse(x$fq %in% TRUE, "direct_or_fq", ifelse(x$nn %in% TRUE, "nn", "farfield_kriging"))
  fine <- ifelse(fr$engine[[1]] == "alfakR_baseline" & x$fq %in% TRUE, "alfakR_fq",
          ifelse(fr$engine[[1]] %in% c("alfak2_graphgp", "alfak2_native_nn_kriging") & x$fq %in% TRUE, "alfak2_direct",
          ifelse(x$nn %in% TRUE & as.character(x$k) %in% rownames(readRDS(as.character(fr$input_rds[[1]]))$x), "nn_observed",
          ifelse(x$nn %in% TRUE, "nn_latent", "kriging_only"))))
  data.frame(simulation_id = as.integer(fr$simulation_id[[1]]), minobs = as.integer(fr$minobs[[1]]),
    method = as.character(fr$method[[1]]), engine = as.character(fr$engine[[1]]),
    input_policy = as.character(fr$input_policy[[1]]), nn_prior = as.character(fr$nn_prior[[1]]),
    lambda = as.numeric(fr$lambda[[1]]), lambda_label = as.character(fr$lambda_label[[1]]),
    repeat_id = as.integer(row_field(fr, "repeat_id", 1L)), scale = "evaluation",
    time_start = as.numeric(fr$time_start[[1]]), time_gap = as.numeric(fr$time_gap[[1]]), time_delta = as.numeric(fr$time_delta[[1]]),
    k = as.character(x$k), estimated_fitness = as.numeric(x$mean), estimated_sd = as.numeric(x$sd),
    true_fitness = as.numeric(truth[match(as.character(x$k), names(truth))]),
    support_scope = scope, support_scope_detail = fine, fq = as.logical(x$fq), nn = as.logical(x$nn),
    status = as.character(fr$status[[1]]), stringsAsFactors = FALSE)
}

native_make_comparisons <- function(nodes, cfg) {
  rbind_nonempty <- function(x) {
    x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
    if (length(x)) do.call(rbind, x) else data.frame()
  }
  baseline_prior <- cfg$nn_prior_baseline %||% "empirical"
  graphgp_vs_alfakR <- rbind_nonempty(lapply(cfg$input_policies, function(policy) {
    pair_delta(nodes, paste0("alfak2_graphgp_", method_policy_token(policy)), paste0("alfakR_baseline_", baseline_prior), "farfield_kriging")
  }))
  native_vs_graphgp <- rbind_nonempty(lapply(cfg$input_policies, function(policy) rbind_nonempty(lapply(cfg$nn_priors, function(pr) {
    pair_delta(nodes, paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_", pr),
               paste0("alfak2_graphgp_", method_policy_token(policy)), "farfield_kriging")
  }))))
  native_vs_alfakR_nn <- rbind_nonempty(lapply(cfg$input_policies, function(policy) rbind_nonempty(lapply(cfg$nn_priors, function(pr) {
    out <- pair_delta(nodes, paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_", pr),
                      paste0("alfakR_baseline_", pr), "nn")
    if (nrow(out)) names(out)[names(out) == "delta_rmse"] <- "delta_nn_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_centered_rmse"] <- "delta_nn_centered_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_pearson"] <- "delta_nn_pearson"
    if (nrow(out)) names(out)[names(out) == "delta_spearman"] <- "delta_nn_spearman"
    out
  }))))
  native_vs_alfakR_far <- rbind_nonempty(lapply(cfg$input_policies, function(policy) rbind_nonempty(lapply(cfg$nn_priors, function(pr) {
    out <- pair_delta(nodes, paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_", pr),
                      paste0("alfakR_baseline_", pr), "farfield_kriging")
    if (nrow(out)) names(out)[names(out) == "delta_rmse"] <- "delta_farfield_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_centered_rmse"] <- "delta_farfield_centered_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_pearson"] <- "delta_farfield_pearson"
    if (nrow(out)) names(out)[names(out) == "delta_spearman"] <- "delta_farfield_spearman"
    out
  }))))
  prior_delta <- if (all(c("empirical_two_step", "empirical_censored") %in% cfg$nn_priors)) {
    rbind_nonempty(lapply(cfg$input_policies, function(policy) {
      pair_delta(nodes, paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_empirical_two_step"),
                 paste0("alfak2_native_nn_kriging_", method_policy_token(policy), "_empirical_censored"),
                 "all_landscape")
    }))
  } else data.frame()
  list(graphgp_vs_alfakR = graphgp_vs_alfakR, native_vs_graphgp = native_vs_graphgp,
       native_vs_alfakR_nn = native_vs_alfakR_nn, native_vs_alfakR_farfield = native_vs_alfakR_far,
       prior_delta = prior_delta)
}

native_win_rates <- function(comp) {
  nn_win <- if (nrow(comp$native_vs_alfakR_nn)) with(comp$native_vs_alfakR_nn, delta_nn_centered_rmse <= 0 | delta_nn_spearman > 0) else logical()
  far_alfakR_win <- if (nrow(comp$native_vs_alfakR_farfield)) with(comp$native_vs_alfakR_farfield, delta_farfield_centered_rmse <= 0 | delta_farfield_spearman > 0) else logical()
  far_graphgp_win <- if (nrow(comp$native_vs_graphgp)) with(comp$native_vs_graphgp, delta_centered_rmse <= 0 | delta_spearman > 0) else logical()
  data.frame(native_nn_win_rate_vs_alfakR = mean(nn_win, na.rm = TRUE),
             native_farfield_win_rate_vs_alfakR = mean(far_alfakR_win, na.rm = TRUE),
             native_farfield_win_rate_vs_graphgp = mean(far_graphgp_win, na.rm = TRUE),
             n_nn_comparisons = length(nn_win),
             n_farfield_vs_alfakR = length(far_alfakR_win),
             n_farfield_vs_graphgp = length(far_graphgp_win),
             stringsAsFactors = FALSE)
}

native_convergence_diagnostics <- function(fit_tbl) {
  cols <- intersect(c("method", "engine", "simulation_id", "lambda", "lambda_label", "minobs", "repeat_id",
                      "input_policy", "nn_prior", "status", "nboot_success", "nboot_failed", "warning_count",
                      "elapsed_sec", "local_convergence", "local_gradient_norm", "local_covariance_status",
                      "nn_status", "kriging_status", "error_message"), names(fit_tbl))
  fit_tbl[, cols, drop = FALSE]
}

write_native_report <- function(cfg, dirs, fit_tbl, metrics, comp, win, rec) {
  path <- file.path(dirs$root, "alfak2_native_nn_kriging_report.md")
  fmt_mean <- function(x, cols) {
    if (!nrow(x)) return("- no rows")
    paste(vapply(cols, function(cl) paste0(cl, "=", signif(mean(x[[cl]], na.rm = TRUE), 4)), character(1)), collapse = ", ")
  }
  lines <- c(
    "# alfak2 native NN/Kriging benchmark report",
    "",
    "## Setup",
    paste0("- alfak2 branch: ", repo_branch(cfg$alfak2_repo), " commit ", repo_head(cfg$alfak2_repo)),
    paste0("- alfakR branch: ", repo_branch(cfg$alfakR_repo), " commit ", repo_head(cfg$alfakR_repo)),
    paste0("- simulations: ", paste(cfg$simulation_ids, collapse = ",")),
    paste0("- lambdas: ", paste(cfg$lambdas, collapse = ",")),
    paste0("- minobs: ", paste(cfg$minobs, collapse = ",")),
    paste0("- repeats: ", paste(cfg$repeat_ids, collapse = ",")),
    paste0("- nn priors: ", paste(cfg$nn_priors, collapse = ",")),
    paste0("- nboot: ", cfg$nboot),
    paste0("- native_global_extra_shell: ", cfg$native_global_extra_shell),
    paste0("- nn_shell_depth: ", cfg$nn_shell_depth),
    paste0("- kriging_max_anchors: ", cfg$kriging_max_anchors),
    paste0("- quick: ", cfg$quick),
    "",
    "## Implementation",
    "- This runner uses alfak2 TMB direct-informed nodes as the native direct layer.",
    "- The current graph Gaussian backend is retained as `alfak2_graphgp_*`.",
    "- The experimental branch runs alfak2-native NN and C++/RcppEigen graph-distance kernel Kriging without converting into alfakR fq/parent state.",
    "- `fit_alfak2()` default behavior is unchanged.",
    "",
    "## Fit status",
    paste0("- rows: ", nrow(fit_tbl), ", ok: ", sum(fit_tbl$status == "ok", na.rm = TRUE), ", failed: ", sum(fit_tbl$status != "ok", na.rm = TRUE)),
    "",
    "## Win rates",
    paste(capture.output(print(win)), collapse = "\n"),
    "",
    "## Native NN vs alfakR NN",
    paste0("- ", fmt_mean(comp$native_vs_alfakR_nn, c("delta_nn_centered_rmse", "delta_nn_spearman", "delta_false_high_rate"))),
    "",
    "## Native Kriging vs alfakR Kriging",
    paste0("- ", fmt_mean(comp$native_vs_alfakR_farfield, c("delta_farfield_centered_rmse", "delta_farfield_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))),
    if (!nrow(comp$native_vs_alfakR_farfield)) "- No farfield delta rows were available for this run; use native_global_extra_shell > nn_shell_depth for Kriging-only nodes." else character(),
    "",
    "## Native Kriging vs alfak2 graphgp",
    paste0("- ", fmt_mean(comp$native_vs_graphgp, c("delta_centered_rmse", "delta_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))),
    if (!nrow(comp$native_vs_graphgp)) "- No graphgp farfield delta rows were available for this run; the fast smoke uses a shallow graph." else character(),
    "",
    "## Recommendation",
    paste0("- ", rec$recommendation[1]),
    paste0("- next_priority: ", rec$next_priority[1])
  )
  writeLines(lines, path)
  invisible(path)
}

native_recommendation <- function(win, comp) {
  far_sd_drop <- if (nrow(comp$native_vs_graphgp)) mean(comp$native_vs_graphgp$delta_estimate_sd_ratio, na.rm = TRUE) else NA_real_
  ok <- is.finite(win$native_farfield_win_rate_vs_graphgp[1]) &&
    win$native_farfield_win_rate_vs_graphgp[1] > 0.5 &&
    (!is.finite(far_sd_drop) || far_sd_drop > -0.25)
  data.frame(
    recommendation = if (ok) "native NN/Kriging is promising enough for broader benchmarking" else
      "keep native NN/Kriging experimental until shape/amplitude diagnostics improve",
    next_priority = "Run shell-2 GRF comparisons on the full lambda/repeat grid, then tune Kriging amplitude and range diagnostics against graphgp.",
    stringsAsFactors = FALSE
  )
}

run_native_fit_mode <- function(cfg, dirs) {
  task_tbl <- read_tsv(file.path(dirs$tables, "native_task_table.tsv"))
  rows <- run_native_fit_tasks(task_tbl, cfg)
  fit_tbl <- list_to_data_frame(rows)
  write_tsv_safe(fit_tbl, file.path(dirs$tables, "native_fit_results.tsv"))
  fit_tbl
}

run_native_fit_task_mode <- function(cfg, dirs, args) {
  task_tbl <- read_tsv(file.path(dirs$tables, "native_task_table.tsv"))
  slurm_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
  task_index <- arg_integer(args, "task_index", if (is.finite(slurm_idx)) slurm_idx else 1L)
  task <- select_task_row(task_tbl, task_index)
  res <- run_native_one_task(task, cfg)
  fit_tbl <- list_to_data_frame(list(res))
  write_tsv_safe(fit_tbl, file.path(dirs$fit_parts, sprintf("task_%06d.tsv", task_index)))
  writeLines("ok", file.path(dirs$fit_parts, sprintf("task_%06d.done", task_index)))
  fit_tbl
}

summarize_native_all <- function(cfg, dirs) {
  task_tbl <- read_tsv(file.path(dirs$tables, "native_task_table.tsv"))
  fit_tbl <- native_read_fit_results(dirs)
  if (!nrow(fit_tbl)) stop("No native fit results found.", call. = FALSE)
  ok <- fit_tbl[fit_tbl$status == "ok", , drop = FALSE]
  node_rows <- lapply(seq_len(nrow(ok)), function(i) native_landscape_to_nodes(ok[i, , drop = FALSE], cfg))
  nodes <- Filter(function(x) is.data.frame(x) && nrow(x), node_rows)
  nodes <- if (length(nodes)) do.call(rbind, nodes) else data.frame()
  if (nrow(nodes)) nodes$estimation_error <- nodes$estimated_fitness - nodes$true_fitness
  metrics <- compute_metrics(nodes)
  comp <- native_make_comparisons(nodes, cfg)
  win <- native_win_rates(comp)
  rec <- native_recommendation(win, comp)
  conv <- native_convergence_diagnostics(fit_tbl)
  amp <- metrics[metrics$support_scope %in% c("farfield_kriging", "all_landscape"), , drop = FALSE]
  write_tsv_safe(fit_tbl, file.path(dirs$tables, "native_fit_results.tsv"))
  write_tsv_safe(nodes, file.path(dirs$tables, "native_node_accuracy.tsv"))
  saveRDS(nodes, file.path(dirs$tables, "native_node_accuracy.rds"))
  write_tsv_safe(metrics, file.path(dirs$tables, "native_all_metrics_long.tsv"))
  write_tsv_safe(metrics[metrics$support_scope %in% c("direct_or_fq", "direct_only", "alfak2_direct", "alfakR_fq"), , drop = FALSE],
                 file.path(dirs$tables, "native_direct_metrics.tsv"))
  write_tsv_safe(metrics[metrics$support_scope %in% c("nn", "nn_observed", "nn_latent"), , drop = FALSE],
                 file.path(dirs$tables, "native_nn_metrics.tsv"))
  write_tsv_safe(metrics[metrics$support_scope %in% c("farfield_kriging", "kriging_only", "all_landscape"), , drop = FALSE],
                 file.path(dirs$tables, "native_kriging_metrics.tsv"))
  write_tsv_safe(comp$graphgp_vs_alfakR, file.path(dirs$tables, "native_graphgp_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$native_vs_graphgp, file.path(dirs$tables, "native_nn_kriging_vs_graphgp_delta.tsv"))
  write_tsv_safe(comp$native_vs_alfakR_nn, file.path(dirs$tables, "native_nn_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$native_vs_alfakR_farfield, file.path(dirs$tables, "native_kriging_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$prior_delta, file.path(dirs$tables, "native_empirical_two_step_vs_censored.tsv"))
  write_tsv_safe(amp, file.path(dirs$tables, "native_amplitude_diagnostics.tsv"))
  write_tsv_safe(conv, file.path(dirs$tables, "native_convergence_diagnostics.tsv"))
  write_tsv_safe(win, file.path(dirs$tables, "native_win_rate_summary.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "native_recommendation.tsv"))
  saveRDS(list(config = cfg, task_table = task_tbl, fit_results = fit_tbl, node_accuracy = nodes,
               metrics = metrics, comparisons = comp, convergence = conv, amplitude = amp,
               win_rate = win, recommendation = rec),
          file.path(dirs$results, "alfak2_native_nn_kriging_all_results.rds"))
  write_native_report(cfg, dirs, fit_tbl, metrics, comp, win, rec)
  invisible(list(fit_results = fit_tbl, nodes = nodes, metrics = metrics, comparisons = comp, win = win, recommendation = rec))
}

validate_native_outputs <- function(dirs) {
  required <- c("native_task_table.tsv", "native_fit_results.tsv", "native_direct_metrics.tsv",
                "native_nn_metrics.tsv", "native_kriging_metrics.tsv", "native_all_metrics_long.tsv",
                "native_nn_kriging_vs_graphgp_delta.tsv", "native_nn_vs_alfakR_delta.tsv",
                "native_kriging_vs_alfakR_delta.tsv", "native_amplitude_diagnostics.tsv",
                "native_convergence_diagnostics.tsv", "native_win_rate_summary.tsv", "native_recommendation.tsv")
  paths <- file.path(dirs$tables, required)
  missing <- paths[!file.exists(paths) | file.info(paths)$size <= 0]
  rds <- file.path(dirs$results, "alfak2_native_nn_kriging_all_results.rds")
  report <- file.path(dirs$root, "alfak2_native_nn_kriging_report.md")
  c(missing, if (!file.exists(rds)) rds else character(), if (!file.exists(report) || file.info(report)$size <= 0) report else character())
}

main_native <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) { usage(); return(invisible(NULL)) }
  repo_dir <- normalizePath(arg_value(args, "repo_dir", find_repo_root()), winslash = "/", mustWork = FALSE)
  cfg <- build_native_config(args, repo_dir)
  dirs <- build_native_dirs(cfg$output_dir)
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  if (identical(mode, "fit_task")) mode <- "fit-task"
  if (!mode %in% c("prepare", "fit", "fit-task", "summarize", "all")) stop("Unsupported --mode=", mode, call. = FALSE)
  if (mode %in% c("fit", "fit-task", "summarize")) {
    cfg0 <- safe_read_rds(file.path(dirs$root, "benchmark_config.rds"))
    if (is.list(cfg0)) {
      cfg_from_args <- cfg
      override_names <- unique(c(intersect(names(args), names(cfg_from_args)),
                                 "repo_dir", "alfak2_repo", "alfakR_repo", "output_dir",
                                 "n_cores", "force_refit", "force_sim", "reuse_dirty_cache", "recompile_dll"))
      cfg <- cfg0
      for (nm in intersect(override_names, names(cfg_from_args))) cfg[[nm]] <- cfg_from_args[[nm]]
    }
  }
  repo_versions <- rbind(repo_state(cfg$alfak2_repo, "alfak2"), repo_state(cfg$alfakR_repo, "alfakR"))
  write_tsv_safe(repo_versions, file.path(dirs$tables, "repo_versions.tsv"))
  message("Loading current source trees.")
  message("  alfak2: ", cfg$alfak2_repo)
  message("  alfakR: ", cfg$alfakR_repo)
  if (!isTRUE(cfg$recompile_dll)) ensure_alfak2_runtime_dll(cfg)
  load_current_repos(cfg$alfakR_repo, cfg$alfak2_repo, recompile_dll = isTRUE(cfg$recompile_dll))
  if (mode == "prepare") return(invisible(native_load_or_prepare_inputs(cfg, dirs)))
  if (mode == "fit") return(invisible(run_native_fit_mode(cfg, dirs)))
  if (mode == "fit-task") return(invisible(run_native_fit_task_mode(cfg, dirs, args)))
  if (mode == "summarize") return(invisible(summarize_native_all(cfg, dirs)))
  native_load_or_prepare_inputs(cfg, dirs)
  run_native_fit_mode(cfg, dirs)
  out <- summarize_native_all(cfg, dirs)
  missing <- validate_native_outputs(dirs)
  if (length(missing)) warning("Missing or empty required outputs: ", paste(missing, collapse = ", "))
  invisible(out)
}

if (sys.nframe() == 0L) main_native()
