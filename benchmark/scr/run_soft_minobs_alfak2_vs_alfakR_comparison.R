#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Run a compact soft_minobs comparison for alfak2 input modes vs alfakR.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_soft_minobs_alfak2_vs_alfakR_comparison.R \\\n",
    "    --source-input-dir=benchmark/results/local_grf_method_blind_lambda0p2_full_minobs_local1_global1_smoke_pm_5e_05 \\\n",
    "    --output-dir=benchmark/results/soft_minobs_alfak2_vs_alfakR_pm_5e_05\n",
    sep = ""
  )
}

source(file.path("benchmark", "scr", "run_grf_alfak2_vs_alfakR_benchmark.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0)
}

repo_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
source_input_dir <- normalize_output_dir(repo_dir, arg_value(
  args,
  "source_input_dir",
  "benchmark/results/local_grf_method_blind_lambda0p2_full_minobs_local1_global1_smoke_pm_5e_05"
))
output_dir <- normalize_output_dir(repo_dir, arg_value(
  args,
  "output_dir",
  "benchmark/results/soft_minobs_alfak2_vs_alfakR_pm_5e_05"
))
alfakR_repo <- normalizePath(arg_value(args, "alfakR_repo", "/Users/4482173/Documents/GitHub/alfakR"),
                             winslash = "/", mustWork = TRUE)
minobs <- arg_integer(args, "minobs", 5L)
nboot <- arg_integer(args, "nboot", 5L)
force <- arg_logical(args, "force", FALSE)
local_shell_depth <- arg_integer(args, "local_shell_depth", 0L)
global_extra_shell <- arg_integer(args, "global_extra_shell", 1L)
compute_sd <- arg_logical(args, "compute_sd", FALSE)

dirs <- list(
  root = output_dir,
  cache = file.path(output_dir, "cache"),
  fits = file.path(output_dir, "fits"),
  tables = file.path(output_dir, "tables"),
  results = file.path(output_dir, "results")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

load_current_repos(alfakR_repo, repo_dir, recompile_dll = FALSE)
repo_versions <- rbind(repo_state(alfakR_repo, "alfakR"), repo_state(repo_dir, "alfak2"))
write_tsv(repo_versions, file.path(dirs$tables, "repo_versions.tsv"))

input_tbl <- read_tsv(file.path(source_input_dir, "tables", "input_table.tsv"))
input_tbl <- input_tbl[input_tbl$simulation_id == 1L &
                         abs(input_tbl$lambda - 0.2) < 1e-12 &
                         input_tbl$time_start == 0 &
                         input_tbl$time_gap == 2 &
                         input_tbl$minobs == minobs, , drop = FALSE]
if (!nrow(input_tbl)) stop("No matching source input row.", call. = FALSE)
input_row <- input_tbl[1, , drop = FALSE]
write_tsv(input_row, file.path(dirs$tables, "input_table.tsv"))

yi <- readRDS(as.character(input_row$input_rds[[1L]]))
grf <- readRDS(as.character(input_row$grf_rds[[1L]]))
fit_beta <- as.numeric(row_field(input_row, "pm", row_field(input_row, "sim_pm", 5e-05)))
if (!is.finite(fit_beta)) fit_beta <- 5e-05
lambda <- as.numeric(input_row$lambda[[1L]])
time_delta <- as.numeric(input_row$time_delta[[1L]])
patient_id <- as.character(input_row$patient_id[[1L]])
lambda_label <- as.character(input_row$lambda_label[[1L]])

fit_rows <- list()
node_rows <- list()
fit_idx <- 0L
node_idx <- 0L

bind_df_fill <- function(dfs) {
  dfs <- Filter(function(x) !is.null(x) && nrow(x), dfs)
  if (!length(dfs)) return(data.frame())
  nms <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    miss <- setdiff(nms, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, nms, drop = FALSE]
  })
  do.call(rbind, dfs)
}

make_fit_row <- function(engine, method, input_policy, nn_prior, outdir, extra = list()) {
  base <- data.frame(
    engine = engine,
    method = method,
    input_policy = input_policy,
    nn_prior = nn_prior,
    simulation_id = as.integer(input_row$simulation_id[[1L]]),
    lambda = lambda,
    lambda_label = lambda_label,
    time_start = as.numeric(input_row$time_start[[1L]]),
    time_gap = as.numeric(input_row$time_gap[[1L]]),
    time_delta = time_delta,
    sim_pm = as.numeric(row_field(input_row, "sim_pm", fit_beta)),
    pm = fit_beta,
    fit_beta_label = pm_to_label(fit_beta),
    patient_id = patient_id,
    grf_key = as.character(input_row$grf_key[[1L]]),
    grf_rds = as.character(input_row$grf_rds[[1L]]),
    input_rds = as.character(input_row$input_rds[[1L]]),
    input_md5 = as.character(input_row$input_md5[[1L]]),
    minobs = minobs,
    outdir = outdir,
    stringsAsFactors = FALSE
  )
  if (length(extra)) {
    for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  }
  base
}

run_one_alfak2 <- function(mode) {
  method <- paste0("alfak2_effective_", mode$name)
  outdir <- file.path(dirs$fits, "alfak2", method)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fit_path <- file.path(outdir, "alfak2_fit.rds")
  summary_path <- file.path(outdir, "landscape.rds")
  result_path <- file.path(outdir, "fit_result.rds")
  if (!force && file.exists(result_path) && file.exists(fit_path)) {
    row <- readRDS(result_path)
    fit <- readRDS(fit_path)
    if (!inherits(fit, "alfak2_fit") && is.list(fit) && !is.null(fit$global$summary)) {
      class(fit) <- "alfak2_fit"
      saveRDS(fit, fit_path)
    }
    if (!identical(as.character(row$status[[1L]]), "ok") &&
        inherits(fit, "alfak2_fit") && !is.null(fit$global$summary)) {
      row$status <- "ok"
      row$error_message <- NA_character_
      row$local_shell_depth <- local_shell_depth
      row$global_extra_shell <- global_extra_shell
      row$compute_sd <- compute_sd
      saveRDS(row, result_path)
      saveRDS(fit$global$summary, summary_path)
    }
    nodes <- coerce_alfak2_nodes(
      fit = fit,
      method = method,
      input_policy = mode$name,
      fit_row = row,
      centroids = grf$centroids,
      lambda = lambda,
      use_legacy_scale = TRUE
    )
    if (nrow(nodes)) {
      node_idx <<- node_idx + 1L
      node_rows[[node_idx]] <<- nodes
    }
    return(row)
  }

  counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = mode$counts_policy, drop_diploid = TRUE)
  selected_times <- suppressWarnings(as.numeric(colnames(counts)))
  dt <- diff(selected_times)
  if (!is.finite(dt) || dt <= 0) dt <- time_delta
  k_mat <- parse_karyotype_ids_base(rownames(counts))
  max_cn <- max(k_mat, na.rm = TRUE) + local_shell_depth + global_extra_shell
  started <- Sys.time()
  res <- tryCatch({
    data <- alfak2:::prepare_counts_for_input_depth(
      counts,
      dt = dt,
      beta = fit_beta,
      input_depth = "effective",
      effective_depth_mode = "min"
    )
    local_graph <- alfak2::build_karyotype_graph(
      data,
      shell_depth = local_shell_depth,
      min_cn = 0,
      max_cn = as.integer(max_cn),
      max_nodes = 150000
    )
    local <- alfak2::fit_local_posterior(
      data,
      local_graph,
      control = list(eval.max = 500, iter.max = 500),
      retry_control = list(eval.max = 2000, iter.max = 2000)
    )
    global_graph <- alfak2::build_karyotype_graph(
      data,
      shell_depth = local_shell_depth + global_extra_shell,
      min_cn = 0,
      max_cn = as.integer(max_cn),
      max_nodes = 150000
    )
    global <- alfak2::fit_graph_posterior(
      local,
      global_graph,
      lambda_l_grid = 0.2,
      lambda_e_grid = 0.01,
      sigma_obs_grid = 0.05,
      graph_edge_weight = "normalized",
      anchor_support_tiers = "directly_informed",
      anchor_count_reference = mode$anchor_count_reference,
      anchor_count_power = 1,
      compute_sd = compute_sd
    )
    fit <- list(data = data, local = local, global = global,
                diagnostics = list(local = local$diagnostics, global = global$diagnostics,
                                   graph_edge_weight = "normalized",
                                   local_shell_depth = local_shell_depth,
                                   global_extra_shell = global_extra_shell,
                                   compute_sd = compute_sd))
    class(fit) <- "alfak2_fit"
    fit <- alfak2:::add_alfakR_scale_to_fit(
      fit,
      n0 = 100000,
      nb = 10000000,
      correct_efflux = TRUE,
      legacy_weight = "pi0"
    )
    saveRDS(fit, fit_path)
    saveRDS(alfak2::summarize_alfak2(fit, layer = "global"), summary_path)
    row <- make_fit_row(
      "alfak2", method, mode$name, NA_character_, outdir,
      list(
        status = "ok",
        error_message = NA_character_,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = fit_path,
        landscape_path = summary_path,
        graph_edge_weight = "normalized",
        anchor_count_reference = if (is.null(mode$anchor_count_reference)) NA_real_ else mode$anchor_count_reference,
        anchor_count_power = 1,
        local_shell_depth = local_shell_depth,
        global_extra_shell = global_extra_shell,
        compute_sd = compute_sd,
        local_convergence = fit$local$diagnostics$convergence,
        local_gradient_norm = fit$local$diagnostics$gradient_norm,
        local_covariance_status = fit$local$diagnostics$covariance_status,
        local_retry_attempted = fit$local$diagnostics$retry_attempted,
        global_factorization_status = fit$global$diagnostics$factorization_status,
        local_nodes = nrow(fit$local$summary),
        global_nodes = nrow(fit$global$summary)
      )
    )
    nodes <- coerce_alfak2_nodes(
      fit = fit,
      method = method,
      input_policy = mode$name,
      fit_row = row,
      centroids = grf$centroids,
      lambda = lambda,
      use_legacy_scale = TRUE
    )
    if (nrow(nodes)) {
      node_idx <<- node_idx + 1L
      node_rows[[node_idx]] <<- nodes
    }
    row
  }, error = function(e) {
    make_fit_row(
      "alfak2", method, mode$name, NA_character_, outdir,
      list(
        status = "error",
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = fit_path,
        landscape_path = summary_path,
        graph_edge_weight = "normalized",
        anchor_count_reference = if (is.null(mode$anchor_count_reference)) NA_real_ else mode$anchor_count_reference,
        anchor_count_power = 1
        ,
        local_shell_depth = local_shell_depth,
        global_extra_shell = global_extra_shell,
        compute_sd = compute_sd
      )
    )
  })
  saveRDS(res, result_path)
  res
}

run_one_alfakR <- function(method_name) {
  method <- paste0("alfakR_", method_name)
  outdir <- file.path(dirs$fits, "alfakR", paste0("nn_prior_", method_name))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  landscape_path <- file.path(outdir, "landscape.Rds")
  if (!force && file.exists(result_path) && file.exists(landscape_path)) {
    row <- readRDS(result_path)
    landscape <- readRDS(landscape_path)
    nodes <- coerce_alfakR_landscape_nodes(
      landscape = landscape,
      method = method,
      nn_prior = method_name,
      fit_row = row,
      centroids = grf$centroids,
      lambda = lambda
    )
    if (nrow(nodes)) {
      node_idx <<- node_idx + 1L
      node_rows[[node_idx]] <<- nodes
    }
    return(row)
  }

  yi_r <- prepare_alfakR_yi(yi, drop_diploid = TRUE)
  started <- Sys.time()
  res <- tryCatch({
    set.seed(1000 + match(method_name, alfakR_methods))
    cap <- capture_warnings(
      alfakR::alfak(
        yi = yi_r,
        outdir = outdir,
        passage_times = NULL,
        minobs = minobs,
        nboot = nboot,
        n0 = 100000,
        nb = 10000000,
        pm = fit_beta,
        correct_efflux = TRUE,
        nn_prior = method_name,
        nn_prior_grid_n = 81,
        nn_prior_fit_subset = "hybrid",
        nn_prior_zero_exposure_quantile = 0.10,
        nn_prior_zero_weight_scale = 0.50,
        nn_prior_zero_weight_cap_ratio = NULL,
        nn_prior_zero_birth_fallback_weight = NULL,
        nn_prior_zero_birth_child_floor = 0.25,
        nn_prior_zero_birth_child_shape = 1,
        nn_prior_zero_birth_replicate_floor = 0.50,
        nn_prior_zero_birth_replicate_shape = 1,
        nn_prior_two_step_support = "rescue",
        nn_prior_two_step_support_min = 0.15,
        nn_prior_two_step_cap_floor = 0.30
      )
    )
    row <- make_fit_row(
      "alfakR", method, "alfakR_minobs_internal", method_name, outdir,
      list(
        status = "ok",
        error_message = NA_character_,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        warning_count = length(cap$warnings),
        warning_messages = if (length(cap$warnings)) paste(cap$warnings, collapse = " || ") else NA_character_,
        landscape_path = landscape_path,
        bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
        posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"),
        xval_path = file.path(outdir, "xval.Rds")
      )
    )
    landscape <- readRDS(landscape_path)
    nodes <- coerce_alfakR_landscape_nodes(
      landscape = landscape,
      method = method,
      nn_prior = method_name,
      fit_row = row,
      centroids = grf$centroids,
      lambda = lambda
    )
    if (nrow(nodes)) {
      node_idx <<- node_idx + 1L
      node_rows[[node_idx]] <<- nodes
    }
    row
  }, error = function(e) {
    make_fit_row(
      "alfakR", method, "alfakR_minobs_internal", method_name, outdir,
      list(
        status = "error",
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        warning_count = 0L,
        warning_messages = NA_character_,
        landscape_path = landscape_path
      )
    )
  })
  saveRDS(res, result_path)
  res
}

alfak2_modes <- list(
  list(name = "full", counts_policy = "full", anchor_count_reference = NULL),
  list(name = "minobs_matched", counts_policy = "minobs_matched", anchor_count_reference = NULL),
  list(name = "soft_minobs", counts_policy = "soft_minobs", anchor_count_reference = as.numeric(minobs)),
  list(name = "full_count_weighted_anchors", counts_policy = "full", anchor_count_reference = as.numeric(minobs))
)
alfakR_methods <- c("none", "empirical", "empirical_censored", "empirical_censored_weighted", "empirical_two_step")

for (mode in alfak2_modes) {
  message("Fitting alfak2 mode: ", mode$name)
  fit_idx <- fit_idx + 1L
  fit_rows[[fit_idx]] <- run_one_alfak2(mode)
}
for (method in alfakR_methods) {
  message("Fitting alfakR nn_prior: ", method)
  fit_idx <- fit_idx + 1L
  fit_rows[[fit_idx]] <- run_one_alfakR(method)
}

fit_tbl <- list_to_data_frame(fit_rows)
write_tsv(fit_tbl, file.path(dirs$tables, "fit_results.tsv"))
node_tbl <- bind_df_fill(node_rows)
if (nrow(node_tbl)) node_tbl$estimation_error <- node_tbl$estimated_fitness - node_tbl$true_fitness
write_tsv(node_tbl, file.path(dirs$tables, "node_accuracy.tsv"))
summary_tbl <- summarize_accuracy(node_tbl)
write_tsv(summary_tbl, file.path(dirs$tables, "summary_by_lambda_time_minobs_method.tsv"))
delta_tbl <- make_delta_vs_alfakR(summary_tbl)
write_tsv(delta_tbl, file.path(dirs$tables, "alfak2_delta_vs_alfakR.tsv"))

whole <- summary_tbl[summary_tbl$support_scope %in% c("whole", "all"), , drop = FALSE]
whole <- whole[order(whole$mae, whole$rmse, -whole$pearson, -whole$spearman), , drop = FALSE]
write_tsv(whole, file.path(dirs$tables, "whole_scope_accuracy_ranking.tsv"))

metric_row <- function(df, scope_name) {
  ok <- is.finite(df$estimated_fitness) & is.finite(df$true_fitness)
  df <- df[ok, , drop = FALSE]
  if (!nrow(df)) {
    return(data.frame(
      n_nodes = 0L, mae = NA_real_, rmse = NA_real_, centered_rmse = NA_real_,
      pearson = NA_real_, spearman = NA_real_, stringsAsFactors = FALSE
    ))
  }
  err <- df$estimated_fitness - df$true_fitness
  ce <- df$estimated_fitness - mean(df$estimated_fitness)
  ct <- df$true_fitness - mean(df$true_fitness)
  data.frame(
    n_nodes = nrow(df),
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2)),
    centered_rmse = sqrt(mean((ce - ct)^2)),
    pearson = if (nrow(df) >= 3L && stats::sd(df$estimated_fitness) > 0 && stats::sd(df$true_fitness) > 0) {
      stats::cor(df$estimated_fitness, df$true_fitness, method = "pearson")
    } else NA_real_,
    spearman = if (nrow(df) >= 3L && stats::sd(df$estimated_fitness) > 0 && stats::sd(df$true_fitness) > 0) {
      suppressWarnings(stats::cor(df$estimated_fitness, df$true_fitness, method = "spearman"))
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

method_keys <- unique(node_tbl$method)
common_k <- Reduce(intersect, split(as.character(node_tbl$k), node_tbl$method)[method_keys])
strict_node_tbl <- node_tbl[node_tbl$k %in% common_k, , drop = FALSE]
strict_rows <- lapply(split(strict_node_tbl, strict_node_tbl$method), function(df) {
  meta <- df[1, c("engine", "method", "input_policy", "nn_prior"), drop = FALSE]
  cbind(meta, support_scope = "strict_common_all_methods", metric_row(df, "strict_common_all_methods"))
})
strict_common <- bind_df_fill(strict_rows)
strict_common <- strict_common[order(strict_common$mae, strict_common$rmse, -strict_common$pearson, -strict_common$spearman), , drop = FALSE]
write_tsv(strict_common, file.path(dirs$tables, "strict_common_accuracy_ranking.tsv"))

support_tier_rows <- lapply(split(node_tbl, node_tbl$method), function(df0) {
  meta <- df0[1, c("engine", "method", "input_policy", "nn_prior"), drop = FALSE]
  tiers <- unique(as.character(df0$support_tier))
  rows <- lapply(tiers, function(tier) {
    cbind(meta, support_tier = tier, metric_row(df0[df0$support_tier == tier, , drop = FALSE], tier))
  })
  rows[[length(rows) + 1L]] <- cbind(meta, support_tier = "all", metric_row(df0, "all"))
  bind_df_fill(rows)
})
support_tier_accuracy <- bind_df_fill(support_tier_rows)
support_tier_accuracy <- support_tier_accuracy[
  order(support_tier_accuracy$engine, support_tier_accuracy$method, support_tier_accuracy$support_tier),
  , drop = FALSE
]
write_tsv(support_tier_accuracy, file.path(dirs$tables, "support_tier_accuracy.tsv"))

compact <- summary_tbl[summary_tbl$support_scope %in% c("whole", "all", "direct", "nn", "other"), , drop = FALSE]
compact <- compact[order(compact$support_scope, compact$mae), , drop = FALSE]
write_tsv(compact, file.path(dirs$tables, "own_scope_accuracy_ranking.tsv"))

strict_cols <- intersect(c("engine", "method", "input_policy", "support_scope", "n_nodes", "mae", "rmse", "pearson", "spearman"), names(summary_tbl))
report <- c(
  "# soft_minobs alfak2 vs alfakR comparison",
  "",
  paste0("- source-input-dir: `", source_input_dir, "`"),
  paste0("- output-dir: `", output_dir, "`"),
  paste0("- condition: simulation_id=1, lambda=0.2, time_gap=2, minobs=", minobs, ", pm=", fit_beta),
  paste0("- alfak2 shell setting: local_shell_depth=", local_shell_depth,
         ", global_extra_shell=", global_extra_shell,
         ", compute_sd=", compute_sd),
  "- alfak2 modes: full, minobs_matched, soft_minobs, full_count_weighted_anchors",
  "- alfak2 global config: graph_edge_weight=normalized, lambda_l=0.2, lambda_e=0.01, sigma_obs=0.05",
  "- alfakR modes: none, empirical, empirical_censored, empirical_censored_weighted, empirical_two_step",
  "",
  "## Whole-scope ranking by MAE",
  paste(utils::capture.output(print(whole[, strict_cols, drop = FALSE], row.names = FALSE)), collapse = "\n"),
  "",
  "## Strict-common ranking by MAE",
  paste(utils::capture.output(print(strict_common[, intersect(c(strict_cols, "centered_rmse"), names(strict_common)), drop = FALSE], row.names = FALSE)), collapse = "\n"),
  "",
  "## Support-tier ranking",
  paste(utils::capture.output(print(support_tier_accuracy[, intersect(c("engine", "method", "input_policy", "support_tier", "n_nodes", "mae", "rmse", "centered_rmse", "pearson", "spearman"), names(support_tier_accuracy)), drop = FALSE], row.names = FALSE)), collapse = "\n"),
  "",
  "## Notes",
  "- `soft_minobs` keeps all observed karyotypes and downweights low-count rows in the local likelihood.",
  "- `full_count_weighted_anchors` is included as a diagnostic: full input plus anchor_count_reference=minobs at the global anchor stage.",
  "- This is a compact single-condition run intended to answer whether soft_minobs changes the mode comparison; broader multi-sim validation should be run before changing defaults."
)
writeLines(report, file.path(output_dir, "soft_minobs_alfak2_vs_alfakR_report.md"))

saveRDS(
  list(input_table = input_row, fit_results = fit_tbl, node_accuracy = node_tbl,
       summary = summary_tbl, delta = delta_tbl, whole_ranking = whole,
       strict_common = strict_common, support_tier_accuracy = support_tier_accuracy),
  file.path(dirs$results, "soft_minobs_alfak2_vs_alfakR_results.rds")
)

message("Wrote soft_minobs comparison under: ", output_dir)
