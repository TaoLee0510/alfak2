#!/usr/bin/env Rscript

find_repo_root <- function(start = getwd()) {
  cand <- normalizePath(file.path(start, c(".", "..", "../..", "../../..")), winslash = "/", mustWork = FALSE)
  for (x in unique(cand)) {
    if (file.exists(file.path(x, "DESCRIPTION")) && dir.exists(file.path(x, "benchmark", "scr"))) return(x)
  }
  stop("Could not locate repository root.", call. = FALSE)
}

repo_dir0 <- find_repo_root()
helper_path <- file.path(repo_dir0, "benchmark", "scr", "run_grf_alfak2_vs_alfakR_benchmark.R")
if (!file.exists(helper_path)) stop("Missing helper benchmark script: ", helper_path, call. = FALSE)
suppressMessages(source(helper_path, local = .GlobalEnv))

usage <- function() {
  cat(
    "Hybrid alfak2-direct -> alfakR NN/Kriging GRF benchmark\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_hybrid_alfak2_direct_alfakR_nn_benchmark.R --mode=all [options]\n\n",
    "Modes: prepare, fit, fit-task, summarize, all\n",
    "Core options:\n",
    "  --alfak2-repo=.\n",
    "  --alfakR-repo=../alfakR\n",
    "  --output-dir=benchmark/results/hybrid_alfak2_direct_alfakR_nn\n",
    "  --source-input-dir=benchmark/results/grf_alfak2_vs_alfakR_shared_inputs\n",
    "  --simulation-ids=1,2,3,4,5,6,7,8,9,10\n",
    "  --lambdas=0.2,0.6,0.8\n",
    "  --minobs=5,10,20\n",
    "  --input-policies=full,minobs_matched\n",
    "  --nn-priors=none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step\n",
    "  --nn-prior-baseline=empirical\n",
    "  --n-repeats=3\n",
    "  --nboot=20\n",
    "  --sample-depth=2000\n",
    "  --graph-edge-weight=normalized\n",
    "  --n-cores=9\n",
    "  --quick=false\n",
    sep = ""
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

normalize_repo_path <- function(repo_dir, x) {
  if (is.null(x) || !nzchar(as.character(x))) return(NULL)
  x <- as.character(x)
  if (grepl("^/", x)) normalizePath(x, winslash = "/", mustWork = FALSE) else
    normalizePath(file.path(repo_dir, x), winslash = "/", mustWork = FALSE)
}

resolve_alfakR_repo <- function(repo_dir, requested) {
  req <- normalize_repo_path(repo_dir, requested)
  candidates <- unique(c(
    if (!is.null(req) && dir.exists(req)) req else character(),
    normalizePath(file.path(repo_dir, "..", "alfakR"), winslash = "/", mustWork = FALSE),
    "/share/lab_crd/lab_crd/taoli/Project/alfakR",
    req
  ))
  for (cand in candidates) {
    if (!is.na(cand) && nzchar(cand) && dir.exists(cand) && file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
  }
  req %||% candidates[[1L]]
}

repo_branch <- function(path) system_text("git", c("branch", "--show-current"), cwd = path)
repo_head <- function(path) system_text("git", c("rev-parse", "HEAD"), cwd = path)
repo_dirty <- function(path) nzchar(system_text("git", c("status", "--short"), cwd = path))

ensure_alfak2_runtime_dll <- function(cfg) {
  own <- file.path(cfg$alfak2_repo, "src", .Platform$dynlib.ext %||% "alfak2.so")
  own <- file.path(cfg$alfak2_repo, "src", paste0("alfak2", .Platform$dynlib.ext))
  candidates <- unique(c(own, file.path(dirname(cfg$alfak2_repo), "alfak2", "src", paste0("alfak2", .Platform$dynlib.ext))))
  if ("alfak2" %in% names(getLoadedDLLs())) return(invisible(TRUE))
  for (cand in candidates) {
    if (file.exists(cand)) {
      dyn.load(cand)
      message("Loaded alfak2 runtime DLL: ", cand)
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

build_hybrid_config <- function(args, repo_dir) {
  quick <- arg_logical(args, "quick", FALSE)
  sim_ids <- arg_integer_vec(args, "simulation_ids", 1:10)
  minobs <- arg_integer_vec(args, "minobs", c(5L, 10L, 20L))
  nboot <- arg_integer(args, "nboot", 20L)
  lambdas <- if (!is.null(arg_value(args, "lambdas", NULL))) {
    arg_numeric_vec(args, "lambdas", c(0.2, 0.6, 0.8))
  } else {
    arg_numeric_vec(args, "lambda", c(0.2, 0.6, 0.8))
  }
  lambdas <- sort(unique(lambdas[is.finite(lambdas) & lambdas > 0]))
  if (!length(lambdas)) stop("No positive lambdas requested.", call. = FALSE)
  n_repeats <- arg_integer(args, "n_repeats", arg_integer(args, "fit_repeats", 3L))
  repeat_ids <- if (!is.null(arg_value(args, "repeats", NULL))) {
    arg_integer_vec(args, "repeats", seq_len(max(1L, n_repeats)))
  } else {
    seq_len(max(1L, n_repeats))
  }
  repeat_ids <- sort(unique(as.integer(repeat_ids[is.finite(repeat_ids) & repeat_ids > 0L])))
  if (!length(repeat_ids)) repeat_ids <- 1L
  if (isTRUE(quick)) {
    sim_ids <- intersect(sim_ids, 1:3)
    if (!length(sim_ids)) sim_ids <- 1:3
    minobs <- min(minobs)
    nboot <- min(nboot, 5L)
  }
  allowed_priors <- c("none", "empirical", "empirical_censored", "empirical_censored_weighted", "empirical_two_step")
  priors <- arg_character_vec(args, "nn_priors", allowed_priors)
  priors <- sub("^nn_prior_", "", priors)
  priors <- priors[priors %in% allowed_priors]
  baseline_prior <- sub("^nn_prior_", "", as.character(arg_value(args, "nn_prior_baseline", "empirical")))
  if (!baseline_prior %in% allowed_priors) stop("Unsupported --nn-prior-baseline=", baseline_prior, call. = FALSE)
  if (!baseline_prior %in% priors) priors <- c(baseline_prior, priors)
  policies <- arg_character_vec(args, "input_policies", c("full", "minobs_matched"))
  graph_edge_weight <- as.character(arg_value(args, "graph_edge_weight", arg_value(args, "alfak2_graph_edge_weight", "normalized")))
  graph_edge_weight <- match.arg(graph_edge_weight, c("mutation", "unit", "normalized"))
  alfakR_repo <- resolve_alfakR_repo(repo_dir, arg_value(args, "alfakR_repo", "../alfakR"))
  list(
    repo_dir = repo_dir,
    alfak2_repo = normalize_repo_path(repo_dir, arg_value(args, "alfak2_repo", repo_dir)),
    alfakR_repo = alfakR_repo,
    output_dir = normalize_output_dir(repo_dir, arg_value(args, "output_dir", "benchmark/results/hybrid_alfak2_direct_alfakR_nn")),
    source_input_dir = normalize_repo_path(repo_dir, arg_value(args, "source_input_dir", NULL)),
    simulation_ids = sort(unique(as.integer(sim_ids))),
    minobs = sort(unique(as.integer(minobs))),
    input_policies = intersect(policies, c("full", "minobs_matched")),
    nn_priors = unique(priors),
    nn_prior_baseline = baseline_prior,
    lambdas = lambdas,
    repeat_ids = repeat_ids,
    n_repeats = length(repeat_ids),
    nboot = as.integer(nboot),
    quick = quick,
    sample_depth = arg_integer(args, "sample_depth", 2000L),
    graph_edge_weight = graph_edge_weight,
    seed = arg_integer(args, "seed", 424242L),
    pm = arg_numeric(args, "pm", 5e-05),
    k_dim = arg_integer(args, "k_dim", 22L),
    n_centroids = arg_integer(args, "n_centroids", 64L),
    grf_centroid_mode = normalize_grf_centroid_mode(arg_value(args, "grf_centroid_mode", "method_blind")),
    grf_centroid_min_cn = arg_integer(args, "grf_centroid_min_cn", 0L),
    grf_centroid_max_cn = arg_integer(args, "grf_centroid_max_cn", 4L),
    grf_centroid_jitter_sd = {
      y <- suppressWarnings(as.numeric(arg_value(args, "grf_centroid_jitter_sd", NA_character_)))
      if (is.finite(y)) y else NULL
    },
    time_start = arg_numeric(args, "time_start", 0),
    time_gap = arg_numeric(args, "time_gap", 2),
    time_max = arg_numeric(args, "time_max", 360),
    passage_interval = arg_numeric(args, "passage_interval", 45),
    abm_pop_size = arg_numeric(args, "abm_pop_size", 50000),
    abm_delta_t = arg_numeric(args, "abm_delta_t", 1),
    abm_max_pop = arg_numeric(args, "abm_max_pop", 2000000),
    abm_culling_survival = arg_numeric(args, "abm_culling_survival", 0.01),
    n0 = arg_numeric(args, "n0", 100000),
    nb = arg_numeric(args, "nb", 10000000),
    correct_efflux = arg_logical(args, "correct_efflux", TRUE),
    drop_diploid = arg_logical(args, "drop_diploid", TRUE),
    n_cores = arg_integer(args, "n_cores", 9L),
    force_refit = arg_logical(args, "force_refit", FALSE),
    force_sim = arg_logical(args, "force_sim", FALSE),
    reuse_dirty_cache = arg_logical(args, "reuse_dirty_cache", TRUE),
    recompile_dll = arg_logical(args, "recompile_dll", FALSE),
    grid_n = arg_integer(args, "grid_n", 81L),
    nn_prior_fit_subset = as.character(arg_value(args, "nn_prior_fit_subset", "hybrid")),
    nn_prior_zero_exposure_quantile = arg_numeric(args, "nn_prior_zero_exposure_quantile", 0.10),
    nn_prior_zero_weight_scale = arg_numeric(args, "nn_prior_zero_weight_scale", 0.50),
    nn_prior_zero_weight_cap_ratio = { y <- suppressWarnings(as.numeric(arg_value(args, "nn_prior_zero_weight_cap_ratio", NA_character_))); if (is.finite(y)) y else NULL },
    nn_prior_zero_birth_fallback_weight = { y <- suppressWarnings(as.numeric(arg_value(args, "nn_prior_zero_birth_fallback_weight", NA_character_))); if (is.finite(y)) y else NULL },
    nn_prior_zero_birth_child_floor = arg_numeric(args, "nn_prior_zero_birth_child_floor", 0.25),
    nn_prior_zero_birth_child_shape = arg_numeric(args, "nn_prior_zero_birth_child_shape", 1),
    nn_prior_zero_birth_replicate_floor = arg_numeric(args, "nn_prior_zero_birth_replicate_floor", 0.50),
    nn_prior_zero_birth_replicate_shape = arg_numeric(args, "nn_prior_zero_birth_replicate_shape", 1),
    nn_prior_two_step_support = as.character(arg_value(args, "nn_prior_two_step_support", "rescue")),
    nn_prior_two_step_support_min = arg_numeric(args, "nn_prior_two_step_support_min", 0.15),
    nn_prior_two_step_cap_floor = arg_numeric(args, "nn_prior_two_step_cap_floor", 0.30),
    alfak2_input_depth = as.character(arg_value(args, "alfak2_input_depth", "effective")),
    alfak2_effective_depth_mode = as.character(arg_value(args, "alfak2_effective_depth_mode", "min")),
    alfak2_observation_model = { x <- as.character(arg_value(args, "alfak2_observation_model", "")); if (nzchar(x)) x else NULL },
    alfak2_dm_concentration = { y <- suppressWarnings(as.numeric(arg_value(args, "alfak2_dm_concentration", NA_character_))); if (is.finite(y)) y else NULL },
    alfak2_min_cn = arg_integer(args, "alfak2_min_cn", 0L),
    alfak2_max_cn = { y <- suppressWarnings(as.integer(arg_value(args, "alfak2_max_cn", NA_character_))); if (is.finite(y)) y else NA_integer_ },
    alfak2_max_nodes = arg_integer(args, "alfak2_max_nodes", 150000L),
    alfak2_lambda_l_grid = arg_numeric_vec(args, "alfak2_lambda_l_grid", 1),
    alfak2_lambda_e_grid = arg_numeric_vec(args, "alfak2_lambda_e_grid", 0.25),
    alfak2_sigma_obs_grid = arg_numeric_vec(args, "alfak2_sigma_obs_grid", 0.05),
    alfak2_legacy_weight = as.character(arg_value(args, "alfak2_legacy_weight", "pi0")),
    alfak2_eval_max = arg_integer(args, "alfak2_eval_max", 500L),
    alfak2_iter_max = arg_integer(args, "alfak2_iter_max", 500L),
    alfak2_retry_max = arg_integer(args, "alfak2_retry_max", 2000L)
  )
}

build_hybrid_dirs <- function(output_dir) {
  dirs <- list(root = output_dir, cache = file.path(output_dir, "cache"), fits = file.path(output_dir, "fits"),
               tables = file.path(output_dir, "tables"), results = file.path(output_dir, "results"),
               fit_parts = file.path(output_dir, "tables", "fit_results_parts"))
  for (p in dirs) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  dirs
}

candidate_source_dirs <- function(cfg) {
  unique(na.omit(c(
    cfg$source_input_dir,
    file.path(cfg$repo_dir, "benchmark/results/grf_alfak2_vs_alfakR_shared_inputs"),
    file.path(cfg$repo_dir, "benchmark/results/grf_alfak2_vs_alfakR"),
    file.path(cfg$repo_dir, "benchmark/results/farfield_shape_probe_default"),
    file.path(cfg$repo_dir, "benchmark/results/farfield_core_fix_probe"),
    file.path(cfg$repo_dir, "benchmark/results/farfield_shell1_empirical_uncertainty")
  )))
}

load_source_inputs <- function(cfg) {
  for (src in candidate_source_dirs(cfg)) {
    tab <- file.path(src, "tables", "input_table.tsv")
    if (!file.exists(tab)) next
    x <- tryCatch(read_tsv(tab), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) next
    needed <- c("simulation_id", "minobs", "input_rds", "grf_rds", "lambda", "lambda_label", "time_start", "time_gap", "time_delta", "patient_id", "grf_key")
    if (!all(needed %in% names(x))) next
    x$simulation_id <- as.integer(x$simulation_id)
    x$minobs <- as.integer(x$minobs)
    x$lambda <- as.numeric(x$lambda)
    x$time_start <- as.numeric(x$time_start)
    x$time_gap <- as.numeric(x$time_gap)
    if ("sim_pm" %in% names(x)) x$sim_pm <- as.numeric(x$sim_pm)
    keep <- x$simulation_id %in% cfg$simulation_ids &
      x$minobs %in% cfg$minobs &
      numeric_in(x$lambda, cfg$lambdas) &
      numeric_in(x$time_start, cfg$time_start) &
      numeric_in(x$time_gap, cfg$time_gap)
    if ("sim_pm" %in% names(x)) keep <- keep & numeric_in(x$sim_pm, cfg$pm)
    if ("grf_centroid_mode" %in% names(x)) {
      mode_vec <- vapply(as.character(x$grf_centroid_mode), normalize_grf_centroid_mode, character(1))
      keep <- keep & mode_vec == cfg$grf_centroid_mode
    }
    if ("grf_centroid_min_cn" %in% names(x)) keep <- keep & numeric_in(x$grf_centroid_min_cn, cfg$grf_centroid_min_cn)
    if ("grf_centroid_max_cn" %in% names(x)) keep <- keep & numeric_in(x$grf_centroid_max_cn, cfg$grf_centroid_max_cn)
    x <- x[keep, , drop = FALSE]
    if (!nrow(x)) next
    x$grf_rds <- vapply(x$grf_rds, resolve_source_cache_path, character(1), source_dir = src)
    x$input_rds <- vapply(x$input_rds, resolve_source_cache_path, character(1), source_dir = src)
    ok <- file.exists(x$grf_rds) & file.exists(x$input_rds)
    x <- x[ok, , drop = FALSE]
    if (nrow(x)) {
      x$input_source <- "reused"
      x$source_input_dir <- normalizePath(src, winslash = "/", mustWork = FALSE)
      return(x)
    }
  }
  data.frame()
}

prepare_generated_inputs <- function(cfg, dirs) {
  rows <- list(); idx <- 0L
  sim_pm_label <- pm_to_label(cfg$pm)
  time_axis_label <- paste0("tmax_", path_token(cfg$time_max), "_pint_", path_token(cfg$passage_interval))
  grf_landscape_label <- grf_landscape_token(list(grf_centroid_mode = cfg$grf_centroid_mode,
                                                  grf_centroid_min_cn = cfg$grf_centroid_min_cn,
                                                  grf_centroid_max_cn = cfg$grf_centroid_max_cn))
  for (sim_id in cfg$simulation_ids) {
    for (lambda_idx in seq_along(cfg$lambdas)) {
      lambda <- cfg$lambdas[[lambda_idx]]
      lambda_label <- format_grf_label(lambda)
      abm_seed <- cfg$seed + as.integer(sim_id) * 10000L + as.integer(lambda_idx) * 100L
      grf_key <- paste(sim_id, lambda_label, paste0("simpm_", sim_pm_label), grf_landscape_label, time_axis_label, sep = "__")
      grf_path <- file.path(dirs$cache, paste0("grf_sim_", grf_key, ".rds"))
      grf_sim <- if (!isTRUE(cfg$force_sim) && file.exists(grf_path)) {
        readRDS(grf_path)
      } else {
        message("Simulating GRF ABM for hybrid input: simulation_id=", sim_id, " lambda=", lambda)
        out <- simulate_nn_prior_grf_abm(seed = abm_seed, lambda = lambda, p = cfg$pm,
          k_dim = cfg$k_dim, n_centroids = cfg$n_centroids, time_max = cfg$time_max,
          passage_interval = cfg$passage_interval, abm_pop_size = cfg$abm_pop_size,
          abm_delta_t = cfg$abm_delta_t, abm_max_pop = cfg$abm_max_pop,
          abm_culling_survival = cfg$abm_culling_survival, centroid_mode = cfg$grf_centroid_mode,
          centroid_min_cn = cfg$grf_centroid_min_cn, centroid_max_cn = cfg$grf_centroid_max_cn,
          centroid_jitter_sd = cfg$grf_centroid_jitter_sd)
        saveRDS(out, grf_path)
        out
      }
      patient_id <- paste0("grf_", sim_id, "_lambda_", lambda_label, "_simpm_", sim_pm_label,
                           "_", grf_landscape_label, "_", time_axis_label, "_start_", path_token(cfg$time_start),
                           "_gap_", path_token(cfg$time_gap))
      input_rds <- file.path(dirs$cache, paste0("input_", patient_id, ".rds"))
      input_csv <- file.path(dirs$cache, paste0("input_", patient_id, ".csv"))
      yi <- if (!file.exists(input_rds) || isTRUE(cfg$force_sim)) {
        out <- build_two_timepoint_yi_from_abm(grf_sim$sim_wide, cfg$time_start, cfg$time_gap,
                                              cfg$passage_interval, cfg$sample_depth, abm_seed + 2000L)
        saveRDS(out, input_rds)
        utils::write.csv(data.frame(karyotype = rownames(out$x), out$x, check.names = FALSE), input_csv, row.names = FALSE)
        out
      } else readRDS(input_rds)
      input_summary <- summarize_input_rows(yi, cfg$minobs, drop_diploid = cfg$drop_diploid)
      for (i in seq_len(nrow(input_summary))) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(simulation_id = as.integer(sim_id), lambda = lambda, lambda_label = lambda_label,
          time_start = cfg$time_start, time_gap = cfg$time_gap, time_delta = as.numeric(yi$metadata$time_delta),
          sim_pm = cfg$pm, patient_id = patient_id, grf_key = grf_key, grf_rds = grf_path,
          grf_centroid_mode = cfg$grf_centroid_mode, grf_centroid_min_cn = cfg$grf_centroid_min_cn,
          grf_centroid_max_cn = cfg$grf_centroid_max_cn,
          grf_centroid_jitter_sd = if (is.null(cfg$grf_centroid_jitter_sd)) NA_real_ else as.numeric(cfg$grf_centroid_jitter_sd),
          input_rds = input_rds, input_csv = input_csv,
          input_md5 = unname(tools::md5sum(input_rds)), input_source = "generated", source_input_dir = dirs$cache,
          input_summary[i, , drop = FALSE], stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

method_policy_token <- function(policy) if (identical(policy, "minobs_matched")) "minobs" else "full"

build_hybrid_task_table <- function(input_tbl, cfg, dirs) {
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
          repeat_id = as.integer(repeat_id),
          time_start = as.numeric(r$time_start), time_gap = as.numeric(r$time_gap), time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)),
          pm = cfg$pm, fit_beta_label = pm_to_label(cfg$pm), patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key), grf_rds = as.character(r$grf_rds),
          input_rds = as.character(r$input_rds), input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))), minobs = as.integer(r$minobs),
          benchmark_seed = as.integer(base_seed + match(prior, cfg$nn_priors)),
          outdir = file.path(dirs$fits, method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir), stringsAsFactors = FALSE)
      }
      for (policy in cfg$input_policies) {
        idx <- idx + 1L
        method <- paste0("alfak2_direct_only_", if (policy == "full") "full" else "minobs_matched")
        rows[[idx]] <- data.frame(engine = "alfak2_direct_only", method = method, input_policy = policy, nn_prior = NA_character_,
          simulation_id = as.integer(r$simulation_id), lambda = as.numeric(r$lambda), lambda_label = as.character(r$lambda_label), time_start = as.numeric(r$time_start),
          repeat_id = as.integer(repeat_id),
          time_gap = as.numeric(r$time_gap), time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)), pm = cfg$pm, fit_beta_label = pm_to_label(cfg$pm),
          patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key), grf_rds = as.character(r$grf_rds), input_rds = as.character(r$input_rds),
          input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))), minobs = as.integer(r$minobs),
          benchmark_seed = as.integer(base_seed + 500L + match(policy, cfg$input_policies)),
          outdir = file.path(dirs$fits, method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir), stringsAsFactors = FALSE)
        for (prior in cfg$nn_priors) {
          idx <- idx + 1L
          method <- paste0("hybrid_alfak2_direct_", method_policy_token(policy), "_", prior)
          rows[[idx]] <- data.frame(engine = "hybrid", method = method, input_policy = policy, nn_prior = prior,
            simulation_id = as.integer(r$simulation_id), lambda = as.numeric(r$lambda), lambda_label = as.character(r$lambda_label), time_start = as.numeric(r$time_start),
            repeat_id = as.integer(repeat_id),
            time_gap = as.numeric(r$time_gap), time_delta = as.numeric(r$time_delta), sim_pm = as.numeric(row_field(r, "sim_pm", cfg$pm)), pm = cfg$pm, fit_beta_label = pm_to_label(cfg$pm),
            patient_id = as.character(r$patient_id), grf_key = as.character(r$grf_key), grf_rds = as.character(r$grf_rds), input_rds = as.character(r$input_rds),
            input_md5 = as.character(row_field(r, "input_md5", unname(tools::md5sum(as.character(r$input_rds))))), minobs = as.integer(r$minobs),
            benchmark_seed = as.integer(base_seed + 700L + match(policy, cfg$input_policies) * 10L + match(prior, cfg$nn_priors)),
            outdir = file.path(dirs$fits, method, lambda_dir, paste0("sim_", r$simulation_id), paste0("MINOBS_", r$minobs), repeat_dir), stringsAsFactors = FALSE)
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$task_order <- seq_len(nrow(out))
  out
}

load_or_prepare_grf_inputs <- function(cfg, dirs) {
  input_tbl <- load_source_inputs(cfg)
  if (!nrow(input_tbl)) input_tbl <- prepare_generated_inputs(cfg, dirs)
  input_tbl <- input_tbl[input_tbl$simulation_id %in% cfg$simulation_ids & input_tbl$minobs %in% cfg$minobs &
                           numeric_in(input_tbl$lambda, cfg$lambdas), , drop = FALSE]
  input_tbl <- input_tbl[order(input_tbl$lambda, input_tbl$simulation_id, input_tbl$minobs), , drop = FALSE]
  task_tbl <- build_hybrid_task_table(input_tbl, cfg, dirs)
  write_tsv_safe(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv_safe(task_tbl, file.path(dirs$tables, "hybrid_task_table.tsv"))
  write_tsv_safe(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  saveRDS(cfg, file.path(dirs$root, "benchmark_config.rds"))
  list(input_table = input_tbl, task_table = task_tbl)
}

prepare_counts_for_policy <- function(yi, minobs, policy, drop_diploid = TRUE) {
  counts <- as.matrix(yi$x); storage.mode(counts) <- "integer"
  if (isTRUE(drop_diploid)) counts <- drop_diploid_counts(counts)
  counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
  if (identical(policy, "minobs_matched")) counts <- counts[rowSums(counts, na.rm = TRUE) >= as.integer(minobs), , drop = FALSE]
  if (!nrow(counts)) stop("No rows remain for input policy ", policy, call. = FALSE)
  counts
}

fit_alfak2_direct_once <- function(counts, task, cfg) {
  counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
  selected_times <- suppressWarnings(as.numeric(colnames(counts)))
  dt <- if (length(selected_times) >= 2L && all(is.finite(selected_times))) diff(selected_times[1:2]) else as.numeric(task$time_delta)
  if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
  max_cn <- cfg$alfak2_max_cn
  if (!is.finite(max_cn)) {
    k_mat <- parse_karyotype_ids_base(rownames(counts))
    max_cn <- max(k_mat, na.rm = TRUE)
  }
  direct_graph_edge_weight <- if (identical(cfg$graph_edge_weight, "normalized")) "unit" else cfg$graph_edge_weight
  alfak2::fit_alfak2(counts, dt = dt, beta = as.numeric(task$pm), min_cn = cfg$alfak2_min_cn, max_cn = as.integer(max_cn),
    local_shell_depth = 0, global_extra_shell = 0, max_nodes = cfg$alfak2_max_nodes,
    lambda_l_grid = cfg$alfak2_lambda_l_grid, lambda_e_grid = cfg$alfak2_lambda_e_grid, sigma_obs_grid = cfg$alfak2_sigma_obs_grid,
    graph_edge_weight = direct_graph_edge_weight, anchor_support_tiers = "directly_informed", input_depth = cfg$alfak2_input_depth,
    effective_depth_mode = cfg$alfak2_effective_depth_mode, observation_model = cfg$alfak2_observation_model, dm_concentration = cfg$alfak2_dm_concentration,
    alfakR_scale = TRUE, n0 = cfg$n0, nb = cfg$nb, correct_efflux = cfg$correct_efflux, legacy_weight = cfg$alfak2_legacy_weight,
    control = list(eval.max = cfg$alfak2_eval_max, iter.max = cfg$alfak2_iter_max), retry_control = list(eval.max = cfg$alfak2_retry_max, iter.max = cfg$alfak2_retry_max))
}

extract_direct_state <- function(fit, parent_ids) {
  s <- alfak2::summarize_alfak2(fit, layer = "local")
  s <- s[as.character(s$support_tier) == "directly_informed", , drop = FALSE]
  ids <- as.character(parent_ids)
  f <- stats::setNames(rep(NA_real_, length(ids)), ids)
  fn <- stats::setNames(rep(NA_real_, length(ids)), ids)
  x0 <- stats::setNames(rep(NA_real_, length(ids)), ids)
  x0n <- stats::setNames(rep(NA_real_, length(ids)), ids)
  sd <- stats::setNames(rep(NA_real_, length(ids)), ids)
  m <- match(ids, as.character(s$karyotype))
  ok <- !is.na(m)
  f[ok] <- as.numeric(s$fitness_mean_alfakR_scale[m[ok]])
  fn[ok] <- as.numeric(s$fitness_mean[m[ok]])
  sd[ok] <- as.numeric(s$fitness_sd_alfakR_scale[m[ok]])
  if ("pi0" %in% names(s)) x0[ok] <- as.numeric(s$pi0[m[ok]])
  if ("pi1" %in% names(s)) x0n[ok] <- as.numeric(s$pi1[m[ok]])
  x0[!is.finite(x0) | x0 < 0] <- 0
  if (sum(x0, na.rm = TRUE) > 0) x0 <- x0 / sum(x0, na.rm = TRUE)
  list(f = f, f_native = fn, x0 = x0, sd = sd, summary = s)
}

build_alfakR_parent_state_from_alfak2 <- function(count_data, current_fq, current_timepoints, current_n0, current_nb,
                                                  current_nn_info, correct_efflux, task, cfg, diag_env, context) {
  direct_counts <- count_data[current_fq, , drop = FALSE]
  direct_counts <- direct_counts[rowSums(direct_counts, na.rm = TRUE) > 0, , drop = FALSE]
  if (!nrow(direct_counts)) stop("Hybrid alfak2 direct bootstrap has no nonzero parent rows.", call. = FALSE)
  fit <- fit_alfak2_direct_once(direct_counts, task, cfg)
  state <- extract_direct_state(fit, current_fq)
  f_final <- state$f
  x0_final <- state$x0
  finite_parent <- is.finite(f_final) & is.finite(x0_final)
  fpar <- f_final[finite_parent]
  x0par <- x0_final[finite_parent]
  if (length(fpar) < 2L) stop("Hybrid alfak2 direct produced fewer than two finite parents.", call. = FALSE)
  if (sum(x0par, na.rm = TRUE) > 0) x0par <- x0par / sum(x0par, na.rm = TRUE)
  opt_res <- list(x0 = unname(x0par), f = unname(fpar))
  names(opt_res$x0) <- names(x0par); names(opt_res$f) <- names(fpar)
  x_norm <- alfakR:::normalize_columns(count_data[names(fpar), , drop = FALSE])
  peak_times <- current_timepoints[apply(x_norm, 1, which.max)]
  bt <- alfakR:::sanitize_birth_times(
    alfakR:::find_birth_times(opt_res, time_range = c(-1000, max(current_timepoints)), minF = 1 / current_n0),
    peak_times = peak_times, timepoints = current_timepoints)
  birth_times <- bt$birth_times; names(birth_times) <- names(fpar)
  birth_fallback <- bt$fallback_mask; names(birth_fallback) <- names(fpar)
  xfit <- alfakR:::project_forward_log(x0par, fpar, current_timepoints)
  rownames(xfit) <- names(fpar)
  ntot <- round(colSums(count_data))
  nn_child_contexts <- lapply(current_nn_info, function(nni_item) {
    alfakR:::prepare_nn_child_context(nni_item, count_data, fpar, birth_times, birth_fallback, xfit, current_timepoints, ntot)
  })
  names(nn_child_contexts) <- names(current_nn_info)
  row <- data.frame(context = context, n_input_rows = nrow(count_data), n_parent_requested = length(current_fq),
    n_parent_nonzero = nrow(direct_counts), n_parent_finite = length(fpar), direct_native_mean = mean(state$f_native, na.rm = TRUE),
    direct_alfakR_scale_mean = mean(state$f, na.rm = TRUE), local_convergence = fit$local$diagnostics$convergence,
    local_gradient_norm = fit$local$diagnostics$gradient_norm, local_covariance_status = fit$local$diagnostics$covariance_status,
    stringsAsFactors = FALSE)
  diag_env$rows[[length(diag_env$rows) + 1L]] <- row
  list(f_initial = f_final, f_final = f_final, x0_initial = x0_final, x0_final = x0_final,
       fpar = fpar, x0par = x0par, ntot_rounded = ntot, nn_child_contexts = nn_child_contexts,
       search_interval = alfakR:::expand_nn_fitness_search_interval(fpar))
}

with_hybrid_parent_state <- function(task, cfg, diag_env, expr) {
  ns <- asNamespace("alfakR")
  name <- "prepare_bootstrap_nn_dataset_state"
  old <- get(name, envir = ns)
  replacement <- function(count_data, current_fq, current_timepoints, current_epsilon, current_n0, current_nb,
                          current_viability, current_nn_info, correct_efflux = FALSE, context = "bootstrap replicate") {
    build_alfakR_parent_state_from_alfak2(count_data, current_fq, current_timepoints, current_n0, current_nb,
      current_nn_info, correct_efflux, task, cfg, diag_env, context)
  }
  was_locked <- bindingIsLocked(name, ns)
  if (was_locked) unlockBinding(name, ns)
  assign(name, replacement, envir = ns)
  if (was_locked) lockBinding(name, ns)
  on.exit({
    if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
    assign(name, old, envir = ns)
    if (was_locked) lockBinding(name, ns)
  }, add = TRUE)
  force(expr)
}

run_alfakR_baseline_fit <- function(task, cfg) {
  outdir <- task$outdir; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  cached <- if (!isTRUE(cfg$force_refit)) safe_read_rds(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) { cached$cached <- TRUE; return(cached) }
  yi <- prepare_alfakR_yi(readRDS(task$input_rds), drop_diploid = cfg$drop_diploid)
  started <- Sys.time(); warnings <- character()
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    cap <- capture_warnings(alfakR::alfak(yi = yi, outdir = outdir, passage_times = NULL, minobs = task$minobs,
      nboot = cfg$nboot, n0 = cfg$n0, nb = cfg$nb, pm = task$pm, correct_efflux = cfg$correct_efflux,
      nn_prior = task$nn_prior, nn_prior_grid_n = cfg$grid_n, nn_prior_fit_subset = cfg$nn_prior_fit_subset,
      nn_prior_zero_exposure_quantile = cfg$nn_prior_zero_exposure_quantile, nn_prior_zero_weight_scale = cfg$nn_prior_zero_weight_scale,
      nn_prior_zero_weight_cap_ratio = cfg$nn_prior_zero_weight_cap_ratio, nn_prior_zero_birth_fallback_weight = cfg$nn_prior_zero_birth_fallback_weight,
      nn_prior_zero_birth_child_floor = cfg$nn_prior_zero_birth_child_floor, nn_prior_zero_birth_child_shape = cfg$nn_prior_zero_birth_child_shape,
      nn_prior_zero_birth_replicate_floor = cfg$nn_prior_zero_birth_replicate_floor, nn_prior_zero_birth_replicate_shape = cfg$nn_prior_zero_birth_replicate_shape,
      nn_prior_two_step_support = cfg$nn_prior_two_step_support, nn_prior_two_step_support_min = cfg$nn_prior_two_step_support_min,
      nn_prior_two_step_cap_floor = cfg$nn_prior_two_step_cap_floor))
    warnings <- cap$warnings
    if (length(warnings)) writeLines(warnings, file.path(outdir, "fit_warnings.log"))
    list(status = "ok", cached = FALSE, error_message = NA_character_, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = result_path,
      nboot_success = cfg$nboot, nboot_failed = 0L, bridge_status = "not_applicable", nn_status = "ok", kriging_status = "ok")
  }, error = function(e) {
    list(status = "error", cached = FALSE, error_message = conditionMessage(e), elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"), posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"),
      fit_path = result_path, nboot_success = 0L, nboot_failed = cfg$nboot, bridge_status = "not_applicable", nn_status = "error", kriging_status = "error")
  })
  out <- c(as.list(task), res)
  saveRDS(out, result_path)
  out
}

run_alfak2_direct_only_fit <- function(task, cfg) {
  outdir <- task$outdir; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  cached <- if (!isTRUE(cfg$force_refit)) safe_read_rds(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) { cached$cached <- TRUE; return(cached) }
  yi <- readRDS(task$input_rds)
  counts0 <- prepare_counts_for_policy(yi, task$minobs, task$input_policy, cfg$drop_diploid)
  started <- Sys.time(); boot_rows <- list(); boot_mat <- NULL
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    fit <- fit_alfak2_direct_once(counts0, task, cfg)
    state <- extract_direct_state(fit, rownames(counts0))
    boot <- lapply(seq_len(cfg$nboot), function(b) {
      bc <- alfakR:::bootstrap_counts(counts0)
      ff <- fit_alfak2_direct_once(bc, task, cfg)
      extract_direct_state(ff, rownames(counts0))$f
    })
    boot_mat <- do.call(rbind, boot); colnames(boot_mat) <- rownames(counts0)
    s <- state$summary
    sd_boot <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
    landscape <- data.frame(k = names(state$f), mean = as.numeric(state$f), median = as.numeric(apply(boot_mat, 2, stats::median, na.rm = TRUE)),
      sd = as.numeric(sd_boot[names(state$f)]), fq = TRUE, nn = FALSE, native_mean = as.numeric(state$f_native), stringsAsFactors = FALSE)
    saveRDS(fit, file.path(outdir, "alfak2_direct_fit.rds"))
    saveRDS(list(final_fitness = boot_mat, initial_fitness = boot_mat, initial_frequencies = matrix(rep(state$x0, each = cfg$nboot), nrow = cfg$nboot),
                 final_frequencies = matrix(rep(state$x0, each = cfg$nboot), nrow = cfg$nboot), nn_fitness = matrix(numeric(0), nrow = cfg$nboot, ncol = 0)),
            file.path(outdir, "bootstrap_res.Rds"))
    saveRDS(landscape, file.path(outdir, "landscape.Rds"))
    saveRDS(matrix(numeric(0), nrow = 0, ncol = 0), file.path(outdir, "landscape_posterior_samples.Rds"))
    list(status = "ok", cached = FALSE, error_message = NA_character_, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = 0L, warning_messages = NA_character_, landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_direct_fit.rds"),
      nboot_success = nrow(boot_mat), nboot_failed = cfg$nboot - nrow(boot_mat), local_convergence = fit$local$diagnostics$convergence,
      local_gradient_norm = fit$local$diagnostics$gradient_norm, local_covariance_status = fit$local$diagnostics$covariance_status,
      scale_status = "alfakR_scale", bridge_status = "direct_only", nn_status = "not_run", kriging_status = "not_run")
  }, error = function(e) {
    list(status = "error", cached = FALSE, error_message = conditionMessage(e), elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = 0L, warning_messages = NA_character_, landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
      posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = file.path(outdir, "alfak2_direct_fit.rds"), nboot_success = 0L, nboot_failed = cfg$nboot,
      local_convergence = NA_integer_, local_gradient_norm = NA_real_, local_covariance_status = NA_character_, scale_status = "error", bridge_status = "error", nn_status = "not_run", kriging_status = "not_run")
  })
  out <- c(as.list(task), res); saveRDS(out, result_path); out
}

run_hybrid_fit <- function(task, cfg) {
  outdir <- task$outdir; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  cached <- if (!isTRUE(cfg$force_refit)) safe_read_rds(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) { cached$cached <- TRUE; return(cached) }
  yi <- prepare_alfakR_yi(readRDS(task$input_rds), drop_diploid = cfg$drop_diploid)
  internal_minobs <- if (identical(task$input_policy, "full")) 1L else as.integer(task$minobs)
  diag_env <- new.env(parent = emptyenv()); diag_env$rows <- list()
  started <- Sys.time(); warnings <- character()
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    cap <- capture_warnings(with_hybrid_parent_state(task, cfg, diag_env, {
      alfakR:::solve_fitness_bootstrap(yi, minobs = internal_minobs, nboot = cfg$nboot, n0 = cfg$n0, nb = cfg$nb, pm = task$pm,
        correct_efflux = cfg$correct_efflux, nn_prior = task$nn_prior, nn_prior_grid_n = cfg$grid_n, nn_prior_fit_subset = cfg$nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = cfg$nn_prior_zero_exposure_quantile, nn_prior_zero_weight_scale = cfg$nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = cfg$nn_prior_zero_weight_cap_ratio, nn_prior_zero_birth_fallback_weight = cfg$nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = cfg$nn_prior_zero_birth_child_floor, nn_prior_zero_birth_child_shape = cfg$nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = cfg$nn_prior_zero_birth_replicate_floor, nn_prior_zero_birth_replicate_shape = cfg$nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = cfg$nn_prior_two_step_support, nn_prior_two_step_support_min = cfg$nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = cfg$nn_prior_two_step_cap_floor)
    }))
    fq_boot <- cap$value; warnings <- cap$warnings
    if (length(warnings)) writeLines(warnings, file.path(outdir, "fit_warnings.log"))
    saveRDS(fq_boot, file.path(outdir, "bootstrap_res.Rds"))
    if (!is.null(fq_boot$nn_prior_diagnostics)) saveRDS(list(replicate = fq_boot$nn_prior_diagnostics, node = fq_boot$nn_two_step_node_diagnostics %||% data.frame()), file.path(outdir, "nn_prior_diagnostics.Rds"))
    landscape_data <- alfakR:::fitKrig(fq_boot, cfg$nboot)
    saveRDS(landscape_data$summary_stats, file.path(outdir, "landscape.Rds"))
    saveRDS(landscape_data$posterior_samples, file.path(outdir, "landscape_posterior_samples.Rds"))
    bridge_diag <- if (length(diag_env$rows)) do.call(rbind, diag_env$rows) else data.frame()
    saveRDS(bridge_diag, file.path(outdir, "bridge_diagnostics.Rds"))
    list(status = "ok", cached = FALSE, error_message = NA_character_, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_,
      landscape_path = file.path(outdir, "landscape.Rds"), bootstrap_path = file.path(outdir, "bootstrap_res.Rds"), posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"),
      fit_path = result_path, nboot_success = nrow(fq_boot$final_fitness), nboot_failed = cfg$nboot - nrow(fq_boot$final_fitness),
      local_convergence = paste(unique(bridge_diag$local_convergence), collapse = ","), local_gradient_norm = suppressWarnings(max(bridge_diag$local_gradient_norm, na.rm = TRUE)),
      local_covariance_status = paste(unique(bridge_diag$local_covariance_status), collapse = ","), scale_status = "alfakR_scale", bridge_status = "ok", nn_status = "ok", kriging_status = "ok")
  }, error = function(e) {
    bridge_diag <- if (length(diag_env$rows)) do.call(rbind, diag_env$rows) else data.frame()
    saveRDS(bridge_diag, file.path(outdir, "bridge_diagnostics.Rds"))
    list(status = "error", cached = FALSE, error_message = conditionMessage(e), elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
      warning_count = length(warnings), warning_messages = if (length(warnings)) paste(warnings, collapse = " || ") else NA_character_, landscape_path = file.path(outdir, "landscape.Rds"),
      bootstrap_path = file.path(outdir, "bootstrap_res.Rds"), posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"), fit_path = result_path,
      nboot_success = 0L, nboot_failed = cfg$nboot, local_convergence = NA_character_, local_gradient_norm = NA_real_, local_covariance_status = NA_character_,
      scale_status = "error", bridge_status = "error", nn_status = "error", kriging_status = "error")
  })
  out <- c(as.list(task), res); saveRDS(out, result_path); out
}

run_one_task <- function(task, cfg) {
  task <- as.list(task)
  if (identical(task$engine, "alfakR_baseline")) return(run_alfakR_baseline_fit(task, cfg))
  if (identical(task$engine, "alfak2_direct_only")) return(run_alfak2_direct_only_fit(task, cfg))
  if (identical(task$engine, "hybrid")) return(run_hybrid_fit(task, cfg))
  stop("Unknown task engine: ", task$engine, call. = FALSE)
}

run_fit_tasks <- function(task_tbl, cfg) {
  n_cores <- max(1L, min(as.integer(cfg$n_cores), nrow(task_tbl)))
  lst <- lapply(seq_len(nrow(task_tbl)), function(i) task_tbl[i, , drop = FALSE])
  if (.Platform$OS.type == "unix" && n_cores > 1L) parallel::mclapply(lst, run_one_task, cfg = cfg, mc.cores = n_cores, mc.preschedule = FALSE) else lapply(lst, run_one_task, cfg = cfg)
}

read_hybrid_fit_results <- function(dirs) {
  parts <- list.files(dirs$fit_parts, pattern = "^task_[0-9]+\\.tsv$", full.names = TRUE)
  p <- file.path(dirs$tables, "hybrid_fit_results.tsv")
  out <- list()
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

landscape_to_nodes <- function(fr, cfg) {
  grf <- safe_read_rds(as.character(fr$grf_rds[[1]])); if (is.null(grf)) return(data.frame())
  land <- safe_read_rds(as.character(fr$landscape_path[[1]])); if (is.null(land)) return(data.frame())
  x <- as.data.frame(land, stringsAsFactors = FALSE)
  if (!"k" %in% names(x) && "karyotype" %in% names(x)) x$k <- x$karyotype
  if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x$mean <- x$fitness_mean
  if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x$sd <- x$fitness_sd
  if (!"sd" %in% names(x)) x$sd <- NA_real_
  if (!"fq" %in% names(x)) x$fq <- grepl("direct_only", as.character(fr$method[[1]]))
  if (!"nn" %in% names(x)) x$nn <- FALSE
  truth <- compute_grf_fitness_truth(as.character(x$k), grf$centroids, as.numeric(fr$lambda[[1]]))
  scope <- ifelse(x$fq %in% TRUE, "direct_or_fq", ifelse(x$nn %in% TRUE, "nn", "farfield_kriging"))
  fine <- ifelse(fr$engine[[1]] == "alfakR_baseline" & x$fq %in% TRUE, "alfakR_fq",
          ifelse(fr$engine[[1]] %in% c("hybrid", "alfak2_direct_only") & x$fq %in% TRUE, "alfak2_direct",
          ifelse(x$nn %in% TRUE & as.character(x$k) %in% rownames(readRDS(as.character(fr$input_rds[[1]]))$x), "nn_observed",
          ifelse(x$nn %in% TRUE, "nn_latent", "kriging_only"))))
  data.frame(simulation_id = as.integer(fr$simulation_id[[1]]), minobs = as.integer(fr$minobs[[1]]), method = as.character(fr$method[[1]]), engine = as.character(fr$engine[[1]]),
    input_policy = as.character(fr$input_policy[[1]]), nn_prior = as.character(fr$nn_prior[[1]]), lambda = as.numeric(fr$lambda[[1]]), lambda_label = as.character(fr$lambda_label[[1]]),
    repeat_id = as.integer(row_field(fr, "repeat_id", 1L)), scale = "alfakR_scale",
    time_start = as.numeric(fr$time_start[[1]]), time_gap = as.numeric(fr$time_gap[[1]]), time_delta = as.numeric(fr$time_delta[[1]]), k = as.character(x$k),
    estimated_fitness = as.numeric(x$mean), estimated_sd = as.numeric(x$sd), true_fitness = as.numeric(truth[match(as.character(x$k), names(truth))]),
    support_scope = scope, support_scope_detail = fine, fq = as.logical(x$fq), nn = as.logical(x$nn), status = as.character(fr$status[[1]]), stringsAsFactors = FALSE)
}

safe_metric_cor <- function(x, y, method) safe_cor(as.numeric(x), as.numeric(y), method)

topk_overlap_value <- function(est, truth) {
  ok <- is.finite(est) & is.finite(truth); est <- est[ok]; truth <- truth[ok]
  if (length(est) < 2L) return(NA_real_)
  k <- max(1L, min(10L, ceiling(0.1 * length(est))))
  mean(order(est, decreasing = TRUE)[seq_len(k)] %in% order(truth, decreasing = TRUE)[seq_len(k)])
}

metric_row <- function(df) {
  ok <- is.finite(df$estimated_fitness) & is.finite(df$true_fitness)
  est <- df$estimated_fitness[ok]; truth <- df$true_fitness[ok]
  err <- est - truth
  ec <- est - mean(est); tc <- truth - mean(truth); cerr <- ec - tc
  esd <- if (length(est) >= 2L) stats::sd(est) else NA_real_; tsd <- if (length(truth) >= 2L) stats::sd(truth) else NA_real_
  erange <- if (length(est) >= 2L) diff(stats::quantile(est, c(.05, .95), na.rm = TRUE, names = FALSE)) else NA_real_
  trange <- if (length(truth) >= 2L) diff(stats::quantile(truth, c(.05, .95), na.rm = TRUE, names = FALSE)) else NA_real_
  eiqr <- if (length(est) >= 2L) stats::IQR(est, na.rm = TRUE) else NA_real_; tiqr <- if (length(truth) >= 2L) stats::IQR(truth, na.rm = TRUE) else NA_real_
  ratio <- if (is.finite(tsd) && tsd > 0) esd / tsd else NA_real_
  data.frame(n_nodes = nrow(df), n_scored = sum(ok), rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_, mae = if (length(err)) mean(abs(err)) else NA_real_,
    centered_rmse = if (length(cerr)) sqrt(mean(cerr^2)) else NA_real_, centered_mae = if (length(cerr)) mean(abs(cerr)) else NA_real_,
    pearson = safe_metric_cor(est, truth, "pearson"), spearman = safe_metric_cor(est, truth, "spearman"), signed_bias = if (length(err)) mean(err) else NA_real_,
    sign_accuracy = if (length(ec)) mean(sign(ec) == sign(tc), na.rm = TRUE) else NA_real_, false_high_rate = if (length(ec)) mean(ec > 0 & tc <= 0, na.rm = TRUE) else NA_real_,
    estimate_sd = esd, truth_sd = tsd, estimate_sd_ratio = ratio, estimate_range_ratio = if (is.finite(trange) && trange > 0) erange / trange else NA_real_,
    estimate_iqr_ratio = if (is.finite(tiqr) && tiqr > 0) eiqr / tiqr else NA_real_, topk_overlap = topk_overlap_value(est, truth),
    coverage95 = if (any(ok) && "estimated_sd" %in% names(df)) mean(abs(err) <= 1.96 * df$estimated_sd[ok], na.rm = TRUE) else NA_real_,
    mean_estimated_sd = mean(df$estimated_sd, na.rm = TRUE), amplitude_collapse = is.finite(ratio) && ratio < 0.02,
    shape_classification = if (!is.finite(ratio) || ratio < 0.02) "amplitude_collapse" else if (is.finite(safe_metric_cor(est, truth, "spearman")) && safe_metric_cor(est, truth, "spearman") < 0) "wrong_direction" else "noncollapsed", stringsAsFactors = FALSE)
}

edge_delta_metrics <- function(nodes, method, simulation_id, minobs, input_policy, nn_prior,
                               lambda, repeat_id, scale = "alfakR_scale") {
  df <- nodes[nodes$method == method & nodes$simulation_id == simulation_id & nodes$minobs == minobs &
                nodes$input_policy == input_policy & nodes$nn_prior == nn_prior &
                numeric_in(nodes$lambda, lambda) & nodes$repeat_id == repeat_id & nodes$scale == scale, , drop = FALSE]
  parents <- df[df$support_scope == "direct_or_fq", , drop = FALSE]
  kids <- df[df$support_scope == "nn", , drop = FALSE]
  if (!nrow(parents) || !nrow(kids)) return(data.frame(edge_delta_sign_agreement = NA_real_, edge_delta_spearman = NA_real_, edge_delta_pearson = NA_real_, edge_delta_sd_ratio = NA_real_))
  pmat <- parse_karyotype_ids_base(parents$k); kmat <- parse_karyotype_ids_base(kids$k)
  rows <- list(); idx <- 0L
  for (i in seq_len(nrow(kmat))) {
    d <- rowSums(abs(sweep(pmat, 2, kmat[i, ], "-")))
    hit <- which(d == 1)
    for (j in hit) { idx <- idx + 1L; rows[[idx]] <- c(est = kids$estimated_fitness[i] - parents$estimated_fitness[j], truth = kids$true_fitness[i] - parents$true_fitness[j]) }
  }
  if (!length(rows)) return(data.frame(edge_delta_sign_agreement = NA_real_, edge_delta_spearman = NA_real_, edge_delta_pearson = NA_real_, edge_delta_sd_ratio = NA_real_))
  m <- do.call(rbind, rows); est <- as.numeric(m[, "est"]); truth <- as.numeric(m[, "truth"]); ok <- is.finite(est) & is.finite(truth)
  data.frame(edge_delta_sign_agreement = if (any(ok)) mean(sign(est[ok]) == sign(truth[ok]), na.rm = TRUE) else NA_real_,
    edge_delta_spearman = safe_metric_cor(est, truth, "spearman"), edge_delta_pearson = safe_metric_cor(est, truth, "pearson"),
    edge_delta_sd_ratio = if (sum(ok) >= 2L && stats::sd(truth[ok]) > 0) stats::sd(est[ok]) / stats::sd(truth[ok]) else NA_real_)
}

compute_metrics <- function(nodes) {
  if (!nrow(nodes)) return(data.frame())
  keys <- c("method", "engine", "simulation_id", "minobs", "input_policy", "nn_prior", "lambda", "lambda_label", "repeat_id", "scale", "time_gap")
  scopes <- c("direct_or_fq", "nn", "farfield_kriging", "all_landscape", "direct_only", "alfakR_fq", "alfak2_direct", "nn_observed", "nn_latent", "kriging_only", "common_nodes", "alfakR_scope_common", "alfak2_scope_common")
  groups <- split(nodes, interaction(nodes[keys], drop = TRUE, lex.order = TRUE))
  out <- list(); idx <- 0L
  for (g in groups) {
    for (scope in scopes) {
      if (scope == "all_landscape") df <- g else if (scope %in% c("alfakR_fq", "alfak2_direct", "nn_observed", "nn_latent", "kriging_only")) df <- g[g$support_scope_detail == scope, , drop = FALSE] else if (scope == "direct_only") df <- g[g$support_scope == "direct_or_fq", , drop = FALSE] else if (scope %in% c("common_nodes", "alfakR_scope_common", "alfak2_scope_common")) df <- g[FALSE, , drop = FALSE] else df <- g[g$support_scope == scope, , drop = FALSE]
      base <- g[1L, keys, drop = FALSE]
      idx <- idx + 1L
      m <- metric_row(df)
      ed <- if (scope == "nn") {
        edge_delta_metrics(nodes, base$method, base$simulation_id, base$minobs, base$input_policy, base$nn_prior,
                           base$lambda, base$repeat_id, base$scale)
      } else data.frame(edge_delta_sign_agreement = NA_real_, edge_delta_spearman = NA_real_, edge_delta_pearson = NA_real_, edge_delta_sd_ratio = NA_real_)
      out[[idx]] <- cbind(base, support_scope = scope, m, ed, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

pair_delta <- function(nodes, lhs_method, rhs_method, scope, extra_filter = function(x) rep(TRUE, nrow(x))) {
  lhs <- nodes[nodes$method == lhs_method & extra_filter(nodes), , drop = FALSE]
  rhs <- nodes[nodes$method == rhs_method & extra_filter(nodes), , drop = FALSE]
  keys <- c("simulation_id", "minobs", "lambda", "lambda_label", "repeat_id", "scale", "time_gap")
  if (!nrow(lhs) || !nrow(rhs)) return(data.frame())
  rows <- list(); idx <- 0L
  for (key in unique(interaction(lhs[keys], drop = TRUE))) {
    l0 <- lhs[interaction(lhs[keys], drop = TRUE) == key, , drop = FALSE]
    r0 <- rhs[interaction(rhs[keys], drop = TRUE) == key, , drop = FALSE]
    if (!nrow(l0) || !nrow(r0)) next
    l0 <- if (scope == "all_landscape") l0 else l0[l0$support_scope == scope, , drop = FALSE]
    r0 <- if (scope == "all_landscape") r0 else r0[r0$support_scope == scope, , drop = FALSE]
    common <- intersect(l0$k, r0$k)
    if (!length(common)) next
    l <- l0[match(common, l0$k), , drop = FALSE]; r <- r0[match(common, r0$k), , drop = FALSE]
    lm <- metric_row(l); rm <- metric_row(r)
    led <- if (identical(scope, "nn")) {
      edge_delta_metrics(nodes, lhs_method, l$simulation_id[1], l$minobs[1], l$input_policy[1], l$nn_prior[1],
                         l$lambda[1], l$repeat_id[1], l$scale[1])
    } else data.frame(edge_delta_sign_agreement = NA_real_, edge_delta_spearman = NA_real_, edge_delta_pearson = NA_real_, edge_delta_sd_ratio = NA_real_)
    red <- if (identical(scope, "nn")) {
      edge_delta_metrics(nodes, rhs_method, r$simulation_id[1], r$minobs[1], r$input_policy[1], r$nn_prior[1],
                         r$lambda[1], r$repeat_id[1], r$scale[1])
    } else data.frame(edge_delta_sign_agreement = NA_real_, edge_delta_spearman = NA_real_, edge_delta_pearson = NA_real_, edge_delta_sd_ratio = NA_real_)
    idx <- idx + 1L
    rows[[idx]] <- data.frame(simulation_id = l$simulation_id[1], minobs = l$minobs[1],
      lambda = l$lambda[1], lambda_label = l$lambda_label[1], repeat_id = l$repeat_id[1],
      scale = l$scale[1], time_gap = l$time_gap[1], lhs_input_policy = l$input_policy[1],
      rhs_input_policy = r$input_policy[1], lhs_nn_prior = l$nn_prior[1], rhs_nn_prior = r$nn_prior[1],
      support_scope = scope, lhs_method = lhs_method, rhs_method = rhs_method,
      n_common = length(common), delta_rmse = lm$rmse - rm$rmse, delta_centered_rmse = lm$centered_rmse - rm$centered_rmse,
      delta_pearson = lm$pearson - rm$pearson, delta_spearman = lm$spearman - rm$spearman, delta_signed_bias = lm$signed_bias - rm$signed_bias,
      delta_estimate_sd_ratio = lm$estimate_sd_ratio - rm$estimate_sd_ratio, delta_topk_overlap = lm$topk_overlap - rm$topk_overlap,
      delta_false_high_rate = lm$false_high_rate - rm$false_high_rate,
      delta_edge_delta_sign_agreement = led$edge_delta_sign_agreement - red$edge_delta_sign_agreement,
      delta_edge_delta_spearman = led$edge_delta_spearman - red$edge_delta_spearman,
      delta_edge_delta_pearson = led$edge_delta_pearson - red$edge_delta_pearson,
      delta_edge_delta_sd_ratio = led$edge_delta_sd_ratio - red$edge_delta_sd_ratio,
      stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

make_comparison_tables <- function(nodes, cfg) {
  rbind_nonempty <- function(x) {
    x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
    if (length(x)) do.call(rbind, x) else data.frame()
  }
  baseline_prior <- cfg$nn_prior_baseline %||% "empirical"
  direct <- rbind_nonempty(lapply(cfg$input_policies, function(policy) {
    dm <- paste0("alfak2_direct_only_", if (policy == "full") "full" else "minobs_matched")
    pair_delta(nodes, dm, paste0("alfakR_baseline_", baseline_prior), "direct_or_fq")
  }))
  nn <- rbind_nonempty(lapply(cfg$input_policies, function(policy) rbind_nonempty(lapply(cfg$nn_priors, function(pr) {
    out <- pair_delta(nodes, paste0("hybrid_alfak2_direct_", method_policy_token(policy), "_", pr), paste0("alfakR_baseline_", pr), "nn")
    if (nrow(out)) names(out)[names(out) == "delta_rmse"] <- "delta_nn_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_centered_rmse"] <- "delta_nn_centered_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_pearson"] <- "delta_nn_pearson"
    if (nrow(out)) names(out)[names(out) == "delta_spearman"] <- "delta_nn_spearman"
    out
  }))))
  krig <- rbind_nonempty(lapply(cfg$input_policies, function(policy) rbind_nonempty(lapply(cfg$nn_priors, function(pr) {
    out <- pair_delta(nodes, paste0("hybrid_alfak2_direct_", method_policy_token(policy), "_", pr), paste0("alfakR_baseline_", pr), "farfield_kriging")
    if (nrow(out)) names(out)[names(out) == "delta_rmse"] <- "delta_farfield_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_centered_rmse"] <- "delta_farfield_centered_rmse"
    if (nrow(out)) names(out)[names(out) == "delta_pearson"] <- "delta_farfield_pearson"
    if (nrow(out)) names(out)[names(out) == "delta_spearman"] <- "delta_farfield_spearman"
    out
  }))))
  full_vs <- rbind_nonempty(lapply(cfg$nn_priors, function(pr) pair_delta(nodes, paste0("hybrid_alfak2_direct_full_", pr), paste0("hybrid_alfak2_direct_minobs_", pr), "all_landscape")))
  compare_priors <- setdiff(cfg$nn_priors, baseline_prior)
  prior_cmp <- rbind_nonempty(lapply(c("alfakR_baseline", paste0("hybrid_alfak2_direct_", sapply(cfg$input_policies, method_policy_token))), function(prefix) {
    rbind_nonempty(lapply(compare_priors, function(pr) {
      pair_delta(nodes, paste0(prefix, "_", pr), paste0(prefix, "_", baseline_prior), "all_landscape")
    }))
  }))
  two_step_vs_censored <- if (all(c("empirical_two_step", "empirical_censored") %in% cfg$nn_priors)) {
    rbind_nonempty(lapply(c("alfakR_baseline", paste0("hybrid_alfak2_direct_", sapply(cfg$input_policies, method_policy_token))), function(prefix) {
      pair_delta(nodes, paste0(prefix, "_empirical_two_step"), paste0(prefix, "_empirical_censored"), "all_landscape")
    }))
  } else data.frame()
  list(direct = direct, nn = nn, kriging = krig, full_vs_minobs = full_vs,
       prior_cmp = prior_cmp, two_step_vs_censored = two_step_vs_censored)
}

scale_diagnostics <- function(nodes, cfg) {
  rows <- list(); idx <- 0L
  for (policy in cfg$input_policies) for (pr in cfg$nn_priors) {
    dm <- paste0("alfak2_direct_only_", if (policy == "full") "full" else "minobs_matched")
    ar <- paste0("alfakR_baseline_", pr)
    for (lambda in unique(nodes$lambda)) for (sid in unique(nodes$simulation_id)) for (mo in unique(nodes$minobs)) for (rid in unique(nodes$repeat_id)) {
      a <- nodes[nodes$method == dm & nodes$simulation_id == sid & nodes$minobs == mo &
                   numeric_in(nodes$lambda, lambda) & nodes$repeat_id == rid &
                   nodes$support_scope == "direct_or_fq", , drop = FALSE]
      b <- nodes[nodes$method == ar & nodes$simulation_id == sid & nodes$minobs == mo &
                   numeric_in(nodes$lambda, lambda) & nodes$repeat_id == rid &
                   nodes$support_scope == "direct_or_fq", , drop = FALSE]
      common <- intersect(a$k, b$k); if (!length(common)) next
      aa <- a[match(common, a$k), , drop = FALSE]; bb <- b[match(common, b$k), , drop = FALSE]
      idx <- idx + 1L
      rows[[idx]] <- data.frame(simulation_id = sid, minobs = mo, lambda = lambda,
        lambda_label = aa$lambda_label[1], repeat_id = rid, input_policy = policy, nn_prior = pr,
        direct_native_mean = NA_real_, direct_alfakR_scale_mean = mean(aa$estimated_fitness, na.rm = TRUE), alfakR_baseline_fq_mean = mean(bb$estimated_fitness, na.rm = TRUE),
        mean_shift = mean(aa$estimated_fitness - bb$estimated_fitness, na.rm = TRUE), sd_ratio = stats::sd(aa$estimated_fitness, na.rm = TRUE) / stats::sd(bb$estimated_fitness, na.rm = TRUE),
        correlation = safe_metric_cor(aa$estimated_fitness, bb$estimated_fitness, "pearson"), bias_vs_truth = mean(aa$estimated_fitness - aa$true_fitness, na.rm = TRUE),
        bias_vs_alfakR = mean(aa$estimated_fitness - bb$estimated_fitness, na.rm = TRUE), scale_status = "alfakR_scale", stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

bootstrap_diagnostics <- function(fit_tbl) {
  if (!nrow(fit_tbl)) return(data.frame())
  out <- fit_tbl[, intersect(c("method", "engine", "simulation_id", "lambda", "lambda_label", "minobs", "repeat_id", "input_policy", "nn_prior", "status", "nboot_success", "nboot_failed", "warning_count", "elapsed_sec", "bridge_status", "nn_status", "kriging_status", "local_convergence", "local_gradient_norm", "local_covariance_status", "error_message"), names(fit_tbl)), drop = FALSE]
  out$fallback_rate <- NA_real_
  for (i in seq_len(nrow(out))) {
    diag_path <- file.path(as.character(fit_tbl$outdir[i]), "nn_prior_diagnostics.Rds")
    d <- safe_read_rds(diag_path)
    rep <- if (is.list(d) && !is.null(d$replicate)) as.data.frame(d$replicate) else data.frame()
    if (nrow(rep) && "nn_prior_mode_used" %in% names(rep)) out$fallback_rate[i] <- mean(rep$nn_prior_mode_used %in% c("none", "empirical_censored_weighted"), na.rm = TRUE)
  }
  out
}

win_rate_summary <- function(comp) {
  direct_win <- if (nrow(comp$direct)) with(comp$direct, delta_centered_rmse <= 0 & delta_spearman >= 0) else logical()
  nn_win <- if (nrow(comp$nn)) with(comp$nn, delta_nn_centered_rmse <= 0 | delta_nn_spearman > 0) else logical()
  far_win <- if (nrow(comp$kriging)) with(comp$kriging, delta_farfield_centered_rmse <= 0 | delta_farfield_spearman > 0) else logical()
  all_three <- logical()
  if (nrow(comp$direct) && nrow(comp$nn) && nrow(comp$kriging)) {
    d <- comp$direct
    d$direct_win <- direct_win
    n <- comp$nn
    n$nn_win <- nn_win
    f <- comp$kriging
    f$farfield_win <- far_win
    dn_keys <- intersect(c("simulation_id", "minobs", "lambda", "lambda_label", "repeat_id", "lhs_input_policy"), names(d))
    nf_keys <- intersect(c("simulation_id", "minobs", "lambda", "lambda_label", "repeat_id", "lhs_input_policy", "lhs_nn_prior"), names(n))
    nf <- merge(n[, c(nf_keys, "nn_win"), drop = FALSE], f[, c(nf_keys, "farfield_win"), drop = FALSE], by = nf_keys)
    dnf <- merge(nf, d[, c(dn_keys, "direct_win"), drop = FALSE], by = dn_keys)
    all_three <- with(dnf, direct_win & nn_win & farfield_win)
  }
  data.frame(direct_win_rate = mean(direct_win, na.rm = TRUE), nn_win_rate = mean(nn_win, na.rm = TRUE), farfield_win_rate = mean(far_win, na.rm = TRUE),
    all_three_layers_win_rate = mean(all_three, na.rm = TRUE), n_direct_comparisons = length(direct_win), n_nn_comparisons = length(nn_win), n_farfield_comparisons = length(far_win), stringsAsFactors = FALSE)
}

recommendation_table <- function(win, comp) {
  dw <- win$direct_win_rate[1]; nw <- win$nn_win_rate[1]; fw <- win$farfield_win_rate[1]
  conclusion <- if (is.finite(dw) && is.finite(nw) && is.finite(fw) && dw >= .5 && nw >= .5 && fw >= .5) {
    "hybrid is worth upstreaming / modularizing into alfak2"
  } else if (is.finite(dw) && dw >= .5) {
    "alfak2 direct improves direct estimation but does not improve alfakR extrapolation"
  } else {
    "benchmark-only bridge should be retained until direct/NN/farfield stability improves"
  }
  data.frame(recommendation = conclusion, upstream_alfakR_nn_kriging = if (grepl("worth upstreaming", conclusion)) "yes" else "no", next_priority = "Run the full 10-simulation nboot=20 grid after addressing full-input amplitude collapse and checking direct TMB convergence diagnostics.", stringsAsFactors = FALSE)
}

write_report <- function(cfg, dirs, repo_versions, fit_tbl, metrics, comp, scale_tbl, boot_tbl, win, rec) {
  path <- file.path(dirs$root, "hybrid_alfak2_direct_alfakR_nn_report.md")
  ok <- fit_tbl[fit_tbl$status == "ok", , drop = FALSE]
  fmt <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!length(x) || !is.finite(x[1])) return("NA")
    format(signif(x[1], 4), scientific = FALSE, trim = TRUE)
  }
  method_line <- function(x, prefix, cols) {
    if (!nrow(x)) return(paste0("- ", prefix, ": no comparable rows."))
    apply(x, 1, function(r) {
      vals <- vapply(cols, function(cc) paste0(cc, "=", fmt(r[[cc]])), character(1))
      paste0("- ", prefix, ": ", r[["lhs_method"]], " vs ", r[["rhs_method"]], " (n_common=", r[["n_common"]], "): ", paste(vals, collapse = ", "))
    })
  }
  direct_lines <- method_line(comp$direct, "Direct/fq", c("delta_centered_rmse", "delta_spearman", "delta_signed_bias", "delta_estimate_sd_ratio"))
  nn_lines <- method_line(comp$nn, "NN", c("delta_nn_centered_rmse", "delta_nn_spearman", "delta_edge_delta_sign_agreement", "delta_edge_delta_spearman", "delta_false_high_rate"))
  kriging_lines <- method_line(comp$kriging, "Kriging/farfield", c("delta_farfield_centered_rmse", "delta_farfield_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))
  policy_lines <- method_line(comp$full_vs_minobs, "Full vs minobs_matched", c("delta_centered_rmse", "delta_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))
  prior_lines <- method_line(comp$prior_cmp, paste0("NN prior vs ", cfg$nn_prior_baseline %||% "empirical"), c("delta_centered_rmse", "delta_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))
  two_step_lines <- method_line(comp$two_step_vs_censored, "empirical_two_step vs empirical_censored", c("delta_centered_rmse", "delta_spearman", "delta_estimate_sd_ratio", "delta_false_high_rate"))
  scale_lines <- if (nrow(scale_tbl)) {
    apply(scale_tbl, 1, function(r) paste0("- ", r[["input_policy"]], " / ", r[["nn_prior"]], ": mean_shift=", fmt(r[["mean_shift"]]),
      ", sd_ratio=", fmt(r[["sd_ratio"]]), ", correlation=", fmt(r[["correlation"]]), ", bias_vs_truth=", fmt(r[["bias_vs_truth"]]),
      ", status=", r[["scale_status"]]))
  } else "- No scale diagnostics rows were available."
  boot_lines <- if (nrow(boot_tbl)) {
    apply(boot_tbl, 1, function(r) paste0("- ", r[["method"]], ": status=", r[["status"]], ", nboot_success=", r[["nboot_success"]],
      ", nboot_failed=", r[["nboot_failed"]], ", bridge=", r[["bridge_status"]], ", nn=", r[["nn_status"]],
      ", kriging=", r[["kriging_status"]], ", warnings=", r[["warning_count"]]))
  } else "- No bootstrap diagnostics rows were available."
  downsample_lines <- if (isTRUE(cfg$quick)) c(
    "## Downsampling note",
    "- This is a quick/downsampled GRF benchmark run, not the full requested 10-simulation, minobs 5/10/20, nboot=20 grid.",
    paste0("- Executed: simulation_ids=", paste(cfg$simulation_ids, collapse = ","), "; lambdas=", paste(cfg$lambdas, collapse = ","),
           "; minobs=", paste(cfg$minobs, collapse = ","), "; repeats=", paste(cfg$repeat_ids, collapse = ","), "; nboot=", cfg$nboot, "."),
    "- Omitted: remaining simulation/minobs/nboot combinations; reason: local runtime was dominated by alfakR two-step and full-input hybrid Kriging, while this run was sufficient to validate the bridge and expose stability issues."
  ) else character()
  lines <- c(
    "# Hybrid alfak2 direct -> alfakR NN/Kriging GRF benchmark", "",
    "## Data source and settings",
    paste0("- Output directory: `", dirs$root, "`"),
    paste0("- Simulation ids: ", paste(cfg$simulation_ids, collapse = ",")),
    paste0("- Lambdas: ", paste(cfg$lambdas, collapse = ",")),
    paste0("- Minobs: ", paste(cfg$minobs, collapse = ",")),
    paste0("- Repeat ids: ", paste(cfg$repeat_ids, collapse = ",")),
    paste0("- Input policies: ", paste(cfg$input_policies, collapse = ",")),
    paste0("- NN priors: ", paste(cfg$nn_priors, collapse = ",")),
    paste0("- NN prior comparison baseline: ", cfg$nn_prior_baseline),
    paste0("- nboot: ", cfg$nboot, if (isTRUE(cfg$quick)) " (quick mode)" else ""),
    paste0("- sample_depth: ", cfg$sample_depth),
    paste0("- graph_edge_weight: ", cfg$graph_edge_weight), "",
    downsample_lines, if (length(downsample_lines)) "" else character(),
    "## Repository state",
    paste0("- alfak2 branch = ", repo_versions$branch[repo_versions$package == "alfak2"], "; commit = ", repo_versions$head[repo_versions$package == "alfak2"]),
    paste0("- alfakR branch = ", repo_versions$branch[repo_versions$package == "alfakR"], "; commit = ", repo_versions$head[repo_versions$package == "alfakR"]),
    paste0("- alfak2 dirty = ", repo_versions$dirty[repo_versions$package == "alfak2"], "; alfakR dirty = ", repo_versions$dirty[repo_versions$package == "alfakR"]), "",
    "## Bridge implementation",
    "- Existing GRF input generation and alfak2 scale helpers were reused from `run_grf_alfak2_vs_alfakR_benchmark.R`.",
    "- No package public API was changed. The bridge is benchmark-only in this runner.",
    "- During hybrid fits, alfakR's internal `prepare_bootstrap_nn_dataset_state()` is temporarily replaced so each bootstrap replicate fits alfak2 direct TMB (`local_shell_depth = 0`, `global_extra_shell = 0`) and returns alfakR-compatible `fpar`, `x0par`, birth times, parent trajectories, and NN child contexts.",
    "- Scale uses alfak2 `alfakR_scale = TRUE` with `n0`, `nb`, `correct_efflux`, and `legacy_weight`; diagnostics are written to `hybrid_scale_diagnostics.tsv`.",
    "- Bootstrap uses independent count bootstraps; the hybrid does not duplicate one MAP row across bootstrap replicates.", "",
    "## Methods",
    paste0("- alfakR baseline priors: ", paste(cfg$nn_priors, collapse = ","), "."),
    "- `alfak2_direct_only_full` and `alfak2_direct_only_minobs_matched`.",
    paste0("- hybrid priors by input policy: ", paste(cfg$nn_priors, collapse = ","), "."), "",
    "## Direct/fq layer results",
    direct_lines, "",
    "## NN layer results",
    nn_lines, "",
    "## Kriging/farfield layer results",
    kriging_lines, "",
    "## Full vs minobs_matched input policy results",
    policy_lines, "",
    "## NN prior vs empirical baseline results",
    prior_lines, "",
    "## empirical_two_step vs empirical_censored results",
    two_step_lines, "",
    "## Scale diagnostics",
    scale_lines, "",
    "## Bootstrap diagnostics",
    boot_lines, "",
    "## Runtime diagnostics",
    paste0("- Successful tasks: ", nrow(ok), " / ", nrow(fit_tbl), "."),
    paste0("- Failed tasks: ", sum(fit_tbl$status != "ok", na.rm = TRUE), ". See `hybrid_fit_results.tsv`."), "",
    "## Win-rate summary",
    paste0("- direct_win_rate: ", signif(win$direct_win_rate[1], 4)),
    paste0("- nn_win_rate: ", signif(win$nn_win_rate[1], 4)),
    paste0("- farfield_win_rate: ", signif(win$farfield_win_rate[1], 4)),
    paste0("- all_three_layers_win_rate: ", signif(win$all_three_layers_win_rate[1], 4)), "",
    "## Final conclusion",
    paste0("- ", rec$recommendation[1]),
    paste0("- Module migration recommendation: ", rec$upstream_alfakR_nn_kriging[1], "."),
    paste0("- Highest priority next step: ", rec$next_priority[1])
  )
  writeLines(lines, path)
  invisible(path)
}

summarize_all <- function(cfg, dirs) {
  task_tbl <- read_tsv(file.path(dirs$tables, "hybrid_task_table.tsv"))
  fit_tbl <- read_hybrid_fit_results(dirs)
  if (!nrow(fit_tbl)) {
    res_files <- list.files(dirs$fits, pattern = "fit_result\\.rds$", recursive = TRUE, full.names = TRUE)
    fits <- lapply(res_files, safe_read_rds); fits <- fits[vapply(fits, is.list, logical(1))]
    fit_tbl <- list_to_data_frame(fits)
  }
  if ("task_order" %in% names(fit_tbl)) fit_tbl <- fit_tbl[order(as.integer(fit_tbl$task_order)), , drop = FALSE]
  write_tsv_safe(fit_tbl, file.path(dirs$tables, "hybrid_fit_results.tsv"))
  node_rows <- list(); idx <- 0L
  for (i in seq_len(nrow(fit_tbl))) if (identical(as.character(fit_tbl$status[i]), "ok")) {
    nd <- landscape_to_nodes(fit_tbl[i, , drop = FALSE], cfg)
    if (nrow(nd)) { idx <- idx + 1L; node_rows[[idx]] <- nd }
  }
  nodes <- if (length(node_rows)) do.call(rbind, node_rows) else data.frame()
  if (nrow(nodes)) nodes$estimation_error <- nodes$estimated_fitness - nodes$true_fitness
  write_tsv_safe(nodes, file.path(dirs$tables, "node_accuracy.tsv")); saveRDS(nodes, file.path(dirs$tables, "node_accuracy.rds"))
  metrics <- compute_metrics(nodes)
  direct_m <- metrics[metrics$support_scope %in% c("direct_or_fq", "direct_only", "alfakR_fq", "alfak2_direct"), , drop = FALSE]
  nn_m <- metrics[metrics$support_scope %in% c("nn", "nn_observed", "nn_latent"), , drop = FALSE]
  far_m <- metrics[metrics$support_scope %in% c("farfield_kriging", "kriging_only", "all_landscape"), , drop = FALSE]
  comp <- make_comparison_tables(nodes, cfg)
  scale_tbl <- scale_diagnostics(nodes, cfg)
  boot_tbl <- bootstrap_diagnostics(fit_tbl)
  win <- win_rate_summary(comp)
  rec <- recommendation_table(win, comp)
  write_tsv_safe(direct_m, file.path(dirs$tables, "hybrid_direct_metrics.tsv"))
  write_tsv_safe(nn_m, file.path(dirs$tables, "hybrid_nn_metrics.tsv"))
  write_tsv_safe(far_m, file.path(dirs$tables, "hybrid_farfield_kriging_metrics.tsv"))
  write_tsv_safe(metrics, file.path(dirs$tables, "hybrid_all_metrics_long.tsv"))
  write_tsv_safe(comp$direct, file.path(dirs$tables, "hybrid_direct_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$nn, file.path(dirs$tables, "hybrid_nn_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$kriging, file.path(dirs$tables, "hybrid_kriging_vs_alfakR_delta.tsv"))
  write_tsv_safe(comp$full_vs_minobs, file.path(dirs$tables, "hybrid_full_vs_minobs_policy.tsv"))
  write_tsv_safe(comp$two_step_vs_censored, file.path(dirs$tables, "hybrid_empirical_censored_vs_two_step.tsv"))
  write_tsv_safe(comp$prior_cmp, file.path(dirs$tables, "hybrid_nn_prior_vs_empirical.tsv"))
  write_tsv_safe(scale_tbl, file.path(dirs$tables, "hybrid_scale_diagnostics.tsv"))
  write_tsv_safe(boot_tbl, file.path(dirs$tables, "hybrid_bootstrap_diagnostics.tsv"))
  write_tsv_safe(win, file.path(dirs$tables, "hybrid_win_rate_summary.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "hybrid_recommendation.tsv"))
  repo_versions <- if (file.exists(file.path(dirs$tables, "repo_versions.tsv"))) read_tsv(file.path(dirs$tables, "repo_versions.tsv")) else data.frame()
  saveRDS(list(config = cfg, task_table = task_tbl, fit_results = fit_tbl, node_accuracy = nodes, metrics = metrics, comparisons = comp, scale_diagnostics = scale_tbl, bootstrap_diagnostics = boot_tbl, win_rate = win, recommendation = rec),
          file.path(dirs$results, "hybrid_alfak2_direct_alfakR_nn_all_results.rds"))
  write_report(cfg, dirs, repo_versions, fit_tbl, metrics, comp, scale_tbl, boot_tbl, win, rec)
  invisible(list(fit_results = fit_tbl, nodes = nodes, metrics = metrics, comparisons = comp, win = win, recommendation = rec))
}

run_fit_mode <- function(cfg, dirs) {
  task_tbl <- read_tsv(file.path(dirs$tables, "hybrid_task_table.tsv"))
  rows <- run_fit_tasks(task_tbl, cfg)
  fit_tbl <- list_to_data_frame(rows)
  write_tsv_safe(fit_tbl, file.path(dirs$tables, "hybrid_fit_results.tsv"))
  fit_tbl
}

run_fit_task_mode <- function(cfg, dirs, args) {
  task_tbl <- read_tsv(file.path(dirs$tables, "hybrid_task_table.tsv"))
  slurm_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
  task_index <- arg_integer(args, "task_index", if (is.finite(slurm_idx)) slurm_idx else 1L)
  task <- select_task_row(task_tbl, task_index)
  res <- run_one_task(task, cfg)
  fit_tbl <- list_to_data_frame(list(res))
  write_tsv_safe(fit_tbl, file.path(dirs$fit_parts, sprintf("task_%06d.tsv", task_index)))
  writeLines("ok", file.path(dirs$fit_parts, sprintf("task_%06d.done", task_index)))
  fit_tbl
}

validate_required_outputs <- function(dirs) {
  required <- c("hybrid_task_table.tsv", "hybrid_fit_results.tsv", "hybrid_direct_metrics.tsv", "hybrid_nn_metrics.tsv", "hybrid_farfield_kriging_metrics.tsv", "hybrid_all_metrics_long.tsv", "hybrid_direct_vs_alfakR_delta.tsv", "hybrid_nn_vs_alfakR_delta.tsv", "hybrid_kriging_vs_alfakR_delta.tsv", "hybrid_full_vs_minobs_policy.tsv", "hybrid_empirical_censored_vs_two_step.tsv", "hybrid_scale_diagnostics.tsv", "hybrid_bootstrap_diagnostics.tsv", "hybrid_win_rate_summary.tsv", "hybrid_recommendation.tsv")
  paths <- file.path(dirs$tables, required)
  missing <- paths[!file.exists(paths) | file.info(paths)$size <= 0]
  rds <- file.path(dirs$results, "hybrid_alfak2_direct_alfakR_nn_all_results.rds")
  report <- file.path(dirs$root, "hybrid_alfak2_direct_alfakR_nn_report.md")
  c(missing, if (!file.exists(rds)) rds else character(), if (!file.exists(report) || file.info(report)$size <= 0) report else character())
}

main_hybrid <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) { usage(); return(invisible(NULL)) }
  repo_dir <- normalizePath(arg_value(args, "repo_dir", find_repo_root()), winslash = "/", mustWork = FALSE)
  cfg <- build_hybrid_config(args, repo_dir)
  dirs <- build_hybrid_dirs(cfg$output_dir)
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  if (identical(mode, "fit_task")) mode <- "fit-task"
  if (!mode %in% c("prepare", "fit", "fit-task", "summarize", "all")) stop("Unsupported --mode=", mode, call. = FALSE)
  if (mode %in% c("fit", "fit-task", "summarize")) {
    cfg0 <- safe_read_rds(file.path(dirs$root, "benchmark_config.rds"))
    if (is.list(cfg0)) {
      cfg_from_args <- cfg
      override_names <- unique(c(
        intersect(names(args), names(cfg_from_args)),
        "repo_dir", "alfak2_repo", "alfakR_repo", "output_dir",
        "n_cores", "force_refit", "force_sim", "reuse_dirty_cache", "recompile_dll"
      ))
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
  if (mode == "prepare") return(invisible(load_or_prepare_grf_inputs(cfg, dirs)))
  if (mode == "fit") return(invisible(run_fit_mode(cfg, dirs)))
  if (mode == "fit-task") return(invisible(run_fit_task_mode(cfg, dirs, args)))
  if (mode == "summarize") return(invisible(summarize_all(cfg, dirs)))
  load_or_prepare_grf_inputs(cfg, dirs)
  run_fit_mode(cfg, dirs)
  out <- summarize_all(cfg, dirs)
  missing <- validate_required_outputs(dirs)
  if (length(missing)) warning("Missing or empty required outputs: ", paste(missing, collapse = ", "))
  invisible(out)
}

if (sys.nframe() == 0L) main_hybrid()
