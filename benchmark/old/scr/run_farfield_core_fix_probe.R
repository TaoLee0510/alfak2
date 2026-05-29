#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_core_fix_probe.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_local_calibration_patch.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

usage <- function() {
  cat(
    "Run farfield core fix probe I1-I6.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_core_fix_probe.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --delta-debug-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --local-calibration-dir=benchmark/results/farfield_local_calibration_patch \\\n",
    "    --output-dir=benchmark/results/farfield_core_fix_probe \\\n",
    "    --simulation-ids=1,2,3,4,5 --minobs=5 --input-policy=full\n",
    sep = ""
  )
}

make_i_dirs <- function(output_dir) make_probe_dirs(output_dir)

read_existing_tsv <- function(path) {
  if (file.exists(path)) read_tsv_safe(path) else data.frame()
}

shape_status_i <- function(x, delta_based = FALSE, edge_gate = TRUE, deploy_gate = TRUE) {
  amp <- is.finite(x$estimate_sd_ratio) & x$estimate_sd_ratio >= 0.02
  rank <- is.finite(x$pearson) & is.finite(x$spearman) & x$pearson > 0 & x$spearman > 0
  if (isTRUE(delta_based) && (!isTRUE(edge_gate) || !isTRUE(deploy_gate))) return("delta_untrusted")
  if (!isTRUE(amp)) return("amplitude_collapse")
  if (!isTRUE(rank)) return("wrong_direction")
  "valid_shape_config"
}

local_shrink_controls_i <- function(shrink) {
  switch(
    as.character(shrink),
    strong = list(eta_borrowed_prior_mean = -8, eta_borrowed_prior_sd = 1.0, eta_distance_penalty = 1.5),
    very_strong = list(eta_borrowed_prior_mean = -10, eta_borrowed_prior_sd = 0.5, eta_distance_penalty = 2.5),
    list(eta_borrowed_prior_mean = -6, eta_borrowed_prior_sd = 1.5, eta_distance_penalty = 0.75)
  )
}

fit_local_variant_i <- function(data, graph, variant = "baseline", shrink = "current",
                                local_parameterization = "f", eval_max = 500,
                                support_tier_f_sd_multiplier = NULL,
                                borrowed_residual_sd = NULL,
                                weakly_supported_residual_sd = NULL) {
  ctrl <- local_shrink_controls_i(shrink)
  do.call(
    alfak2::fit_local_posterior,
    c(
      list(
        data = data,
        graph = graph,
        observation_model = "dirichlet_multinomial",
        dm_concentration = 50,
        control = list(eval.max = eval_max, iter.max = eval_max),
        retry_on_untrusted_covariance = FALSE,
        local_parameterization = local_parameterization,
        return_optimizer_diagnostics = TRUE,
        support_tier_f_sd_multiplier = support_tier_f_sd_multiplier,
        borrowed_residual_sd = borrowed_residual_sd,
        weakly_supported_residual_sd = weakly_supported_residual_sd
      ),
      ctrl
    )
  )
}

local_diag_rows_i <- function(fit, variant, shrink = NA_character_, local_parameterization = NA_character_,
                              eval_max = NA_integer_) {
  opt <- fit$diagnostics$optimizer
  block <- opt$gradient_block_summary
  cbind(
    data.frame(
      variant = variant,
      shrink = shrink,
      local_parameterization = local_parameterization,
      eval_max = eval_max,
      convergence = fit$diagnostics$convergence,
      message = fit$diagnostics$message,
      objective = fit$diagnostics$objective,
      gradient_norm = fit$diagnostics$gradient_norm,
      covariance_status = fit$diagnostics$covariance_status,
      covariance_fallback = fit$diagnostics$covariance_fallback,
      fitness_sd_source = fit$diagnostics$fitness_sd_source,
      stringsAsFactors = FALSE
    ),
    block
  )
}

edge_alignment_i <- function(fit, grf, lambda, variant = NA_character_) {
  summary <- fit$summary
  truth <- compute_grf_truth(summary$karyotype, grf$centroids, lambda)
  truth <- as.numeric(truth[as.character(summary$karyotype)])
  m <- setNames(as.numeric(summary$fitness_mean), as.character(summary$karyotype))
  t <- setNames(truth, as.character(summary$karyotype))
  graph <- fit$graph
  parent <- as.character(graph$labels)[as.integer(unlist(graph$parent_from0)) + 1L]
  child <- as.character(graph$labels)[as.integer(unlist(graph$parent_to0)) + 1L]
  keep <- is.finite(m[parent]) & is.finite(m[child]) & is.finite(t[parent]) & is.finite(t[child])
  if (!any(keep)) {
    return(data.frame(variant = variant, n_edges = 0L, delta_spearman = NA_real_,
                      delta_sign_agreement = NA_real_, estimated_delta_sd_ratio = NA_real_,
                      stringsAsFactors = FALSE))
  }
  ed <- as.numeric(m[child[keep]] - m[parent[keep]])
  td <- as.numeric(t[child[keep]] - t[parent[keep]])
  data.frame(
    variant = variant,
    n_edges = sum(keep),
    delta_spearman = safe_cor2(ed, td, "spearman"),
    delta_sign_agreement = mean(sign(ed) == sign(td), na.rm = TRUE),
    estimated_delta_sd_ratio = stats::sd(ed, na.rm = TRUE) / stats::sd(td, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

run_i1 <- function(local_calibration_dir, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i1_calibration_failure_state_patch.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  h1_root <- normalizePath(local_calibration_dir, winslash = "/", mustWork = TRUE)
  demo <- read_existing_tsv(file.path(h1_root, "tables", "h1_calibration_metrics_gate_demo.tsv"))
  ranked <- read_existing_tsv(file.path(h1_root, "tables", "h1_calibration_ranked_params_demo.tsv"))
  status <- read_existing_tsv(file.path(h1_root, "tables", "h1_calibration_recommended_status.tsv"))
  best_numeric <- read_existing_tsv(file.path(h1_root, "tables", "h1_best_numeric_only_params.tsv"))
  if (!nrow(ranked)) ranked <- demo
  gate_summary <- data.frame(
    n_total_configs = nrow(ranked),
    n_nonoracle_configs = if ("is_oracle" %in% names(ranked)) sum(!ranked$is_oracle, na.rm = TRUE) else NA_integer_,
    n_valid_shape_configs = sum(ranked$recommended_status == "valid_shape_config", na.rm = TRUE),
    n_nonoracle_valid_shape_configs = if ("is_oracle" %in% names(ranked)) sum(ranked$recommended_status == "valid_shape_config" & !ranked$is_oracle, na.rm = TRUE) else 0L,
    n_amplitude_collapse = sum(ranked$recommended_status == "amplitude_collapse", na.rm = TRUE),
    n_wrong_direction = sum(ranked$recommended_status == "wrong_direction", na.rm = TRUE),
    n_delta_untrusted = sum(ranked$recommended_status == "delta_untrusted", na.rm = TRUE),
    recommended_status = if (nrow(status)) status$recommended_status[[1L]] else "no_valid_shape_configuration",
    stringsAsFactors = FALSE
  )
  formal_patch <- data.frame(
    file = "benchmark/scr/run_grf_alfak2_parameter_calibration.R",
    modified_functions = "rank_calibration_parameters,summarize_calibration_results,classify_calibration_shape,calibration_gate_summary",
    compatibility = "Keeps fit_results.tsv, calibration_metrics_by_fit.tsv, calibration_ranked_params.tsv, best_params.tsv, best_params_cli_args.txt; adds gate/status/numeric-only outputs.",
    best_params_cli_args_when_no_valid = "no_valid_shape_configuration",
    stringsAsFactors = FALSE
  )
  write_tsv_safe(demo, file.path(dirs$tables, "i1_calibration_metrics_gate_demo.tsv"))
  write_tsv_safe(ranked, file.path(dirs$tables, "i1_calibration_ranked_params_demo.tsv"))
  write_tsv_safe(gate_summary, file.path(dirs$tables, "i1_calibration_gate_summary.tsv"))
  write_tsv_safe(status, file.path(dirs$tables, "i1_calibration_recommended_status.tsv"))
  write_tsv_safe(best_numeric, file.path(dirs$tables, "i1_best_numeric_only_params.tsv"))
  write_tsv_safe(formal_patch, file.path(dirs$tables, "i1_formal_calibration_patch_summary.tsv"))
  out <- list(demo = demo, ranked = ranked, gate_summary = gate_summary, status = status,
              best_numeric = best_numeric, formal_patch = formal_patch)
  saveRDS(out, rds)
  out
}

run_i2 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i2_global_mean_only_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  local <- bundle$local
  graph <- local$graph
  started <- Sys.time()
  fit_sd <- alfak2::fit_graph_posterior(local, graph, lambda_l_grid = 0.2, lambda_e_grid = 0.01,
                                        sigma_obs_grid = 0.05, graph_edge_weight = "normalized",
                                        compute_sd = TRUE)
  t_sd <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  started <- Sys.time()
  fit_mean <- alfak2::fit_graph_posterior(local, graph, lambda_l_grid = 0.2, lambda_e_grid = 0.01,
                                          sigma_obs_grid = 0.05, graph_edge_weight = "normalized",
                                          compute_sd = FALSE)
  t_mean <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  max_diff <- max(abs(fit_sd$summary$fitness_mean - fit_mean$summary$fitness_mean), na.rm = TRUE)
  equivalence <- data.frame(
    n_nodes = nrow(fit_sd$summary),
    max_abs_mean_diff = max_diff,
    sd_true_all_finite = all(is.finite(fit_sd$summary$fitness_sd)),
    sd_false_all_na = all(is.na(fit_mean$summary$fitness_sd)),
    mean_equivalent = is.finite(max_diff) && max_diff < 1e-10,
    stringsAsFactors = FALSE
  )
  timing <- data.frame(
    mode = c("compute_sd_TRUE", "compute_sd_FALSE"),
    n_nodes = nrow(fit_sd$summary),
    runtime_sec = c(t_sd, t_mean),
    speedup = c(1, if (is.finite(t_mean) && t_mean > 0) t_sd / t_mean else NA_real_),
    stringsAsFactors = FALSE
  )
  cfg <- data.frame(experiment = "I2", candidate_id = "normalized_mean_only", graph_edge_weight = "normalized",
                    lambda_l = 0.2, lambda_e = 0.01, sigma_obs = 0.05, anchor_var_mode = "current",
                    prior_mean_mode = "zero", prior_mean_scale = 0, anchor_count_reference_mode = "none",
                    stringsAsFactors = FALSE)
  metrics <- score_summary_abcd(fit_mean$summary, graph, grf, task_info$lambda, task_info, cfg, "mean_only_no_sd")
  write_tsv_safe(equivalence, file.path(dirs$tables, "i2_global_mean_only_equivalence.tsv"))
  write_tsv_safe(timing, file.path(dirs$tables, "i2_global_mean_only_timing.tsv"))
  write_tsv_safe(metrics, file.path(dirs$tables, "i2_global_mean_only_shape_metrics.tsv"))
  out <- list(equivalence = equivalence, timing = timing, metrics = metrics,
              fit_sd = fit_sd, fit_mean = fit_mean)
  saveRDS(out, rds)
  out
}

run_i3 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i3_local_tmb_diagnostics.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  shell0 <- alfak2::build_karyotype_graph(data, shell_depth = 0, max_nodes = 30000)
  shell1 <- bundle$local$graph
  variants <- data.frame(
    variant = c("shell0_control", "shell1_baseline", "shell1_strong_shrink", "shell1_very_strong_shrink"),
    shell_depth = c(0, 1, 1, 1),
    shrink = c("current", "current", "strong", "very_strong"),
    stringsAsFactors = FALSE
  )
  rows <- list(); tiers <- list(); tops <- list(); aligns <- list()
  for (i in seq_len(nrow(variants))) {
    graph <- if (variants$shell_depth[[i]] == 0) shell0 else shell1
    fit <- fit_local_variant_i(data, graph, variants$variant[[i]], variants$shrink[[i]], eval_max = 500)
    rows[[i]] <- local_diag_rows_i(fit, variants$variant[[i]], variants$shrink[[i]], "f", 500)
    tiers[[i]] <- cbind(data.frame(variant = variants$variant[[i]], stringsAsFactors = FALSE),
                        fit$diagnostics$optimizer$grad_f_by_support_tier)
    tops[[i]] <- cbind(data.frame(variant = variants$variant[[i]], stringsAsFactors = FALSE),
                       fit$diagnostics$optimizer$top_gradient_nodes)
    aligns[[i]] <- edge_alignment_i(fit, grf, task_info$lambda, variants$variant[[i]])
  }
  block <- bind_rows_fill(rows)
  tier_tbl <- bind_rows_fill(tiers)
  top_tbl <- bind_rows_fill(tops)
  align_tbl <- bind_rows_fill(aligns)
  write_tsv_safe(block, file.path(dirs$tables, "i3_local_gradient_block_summary.tsv"))
  write_tsv_safe(tier_tbl, file.path(dirs$tables, "i3_local_f_gradient_by_support_tier.tsv"))
  write_tsv_safe(top_tbl, file.path(dirs$tables, "i3_local_top_gradient_nodes.tsv"))
  write_tsv_safe(block, file.path(dirs$tables, "i3_local_variant_diagnostics.tsv"))
  out <- list(block = block, by_tier = tier_tbl, top_nodes = top_tbl, edge_alignment = align_tbl)
  saveRDS(out, rds)
  out
}

run_i4 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i4_g_parameterization_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  grid <- expand.grid(
    local_parameterization = c("f", "g_equivalent"),
    shrink = c("current", "strong", "very_strong"),
    eval_max = 500,
    stringsAsFactors = FALSE
  )
  rows <- list(); aligns <- list()
  for (i in seq_len(nrow(grid))) {
    variant <- paste("I4", grid$local_parameterization[[i]], grid$shrink[[i]], sep = "_")
    fit <- fit_local_variant_i(data, graph, variant, grid$shrink[[i]], grid$local_parameterization[[i]], grid$eval_max[[i]])
    rows[[i]] <- local_diag_rows_i(fit, variant, grid$shrink[[i]], grid$local_parameterization[[i]], grid$eval_max[[i]])
    aligns[[i]] <- cbind(data.frame(variant = variant, local_parameterization = grid$local_parameterization[[i]],
                                    shrink = grid$shrink[[i]], stringsAsFactors = FALSE),
                         edge_alignment_i(fit, grf, task_info$lambda, variant))
  }
  summary <- bind_rows_fill(rows)
  align <- bind_rows_fill(aligns)
  sdstatus <- summary[, intersect(names(summary), c("variant", "local_parameterization", "shrink", "covariance_status", "covariance_fallback", "fitness_sd_source")), drop = FALSE]
  write_tsv_safe(summary, file.path(dirs$tables, "i4_local_f_vs_g_parameterization.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "i4_local_g_gradient_by_block.tsv"))
  write_tsv_safe(align, file.path(dirs$tables, "i4_local_g_edge_alignment.tsv"))
  write_tsv_safe(sdstatus, file.path(dirs$tables, "i4_local_g_sdreport_status.tsv"))
  out <- list(summary = summary, edge_alignment = align, sdstatus = sdstatus)
  saveRDS(out, rds)
  out
}

run_i5 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i5_borrowed_f_variants.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  variants <- list(
    list(variant = "B0_baseline", shrink = "current"),
    list(variant = "B1_support_tier_prior", shrink = "current",
         support_tier_f_sd_multiplier = c(local_borrowed = 0.5, weakly_supported = 0.25)),
    list(variant = "B2_borrowed_residual_0p10", shrink = "current", borrowed_residual_sd = 0.10),
    list(variant = "B3_deterministic_borrowed_0p01", shrink = "current", borrowed_residual_sd = 0.01),
    list(variant = "B4_weakly_supported_fixed", shrink = "current", weakly_supported_residual_sd = 0.01),
    list(variant = "B5_g_strong_residual_0p10", shrink = "strong", local_parameterization = "g_equivalent", borrowed_residual_sd = 0.10)
  )
  rows <- list(); tiers <- list(); aligns <- list()
  for (i in seq_along(variants)) {
    v <- variants[[i]]
    fit <- fit_local_variant_i(
      data, graph, v$variant, v$shrink %||% "current", v$local_parameterization %||% "f", 500,
      support_tier_f_sd_multiplier = v$support_tier_f_sd_multiplier %||% NULL,
      borrowed_residual_sd = v$borrowed_residual_sd %||% NULL,
      weakly_supported_residual_sd = v$weakly_supported_residual_sd %||% NULL
    )
    rows[[i]] <- local_diag_rows_i(fit, v$variant, v$shrink %||% "current", v$local_parameterization %||% "f", 500)
    tiers[[i]] <- cbind(data.frame(variant = v$variant, stringsAsFactors = FALSE),
                        fit$diagnostics$optimizer$grad_f_by_support_tier)
    aligns[[i]] <- cbind(data.frame(variant = v$variant, stringsAsFactors = FALSE),
                         edge_alignment_i(fit, grf, task_info$lambda, v$variant))
  }
  result <- bind_rows_fill(rows)
  tier_tbl <- bind_rows_fill(tiers)
  align <- bind_rows_fill(aligns)
  rec <- data.frame(
    recommendation = if (any(result$covariance_status == "TMB_sdreport")) "borrowed_variant_candidate_found" else "no_trusted_shell_depth1_variant",
    best_gradient_variant = result$variant[which.min(result$gradient_norm)],
    best_gradient_norm = min(result$gradient_norm, na.rm = TRUE),
    best_edge_sign_agreement = max(align$delta_sign_agreement, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(result, file.path(dirs$tables, "i5_borrowed_f_variant_results.tsv"))
  write_tsv_safe(tier_tbl, file.path(dirs$tables, "i5_borrowed_f_gradient_by_tier.tsv"))
  write_tsv_safe(align, file.path(dirs$tables, "i5_borrowed_f_edge_alignment.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "i5_borrowed_f_variant_recommendation.tsv"))
  out <- list(result = result, by_tier = tier_tbl, edge_alignment = align, recommendation = rec)
  saveRDS(out, rds)
  out
}

fit_multisim_local_i <- function(counts, dt, local_config = "current") {
  data <- alfak2::prepare_alfak2_data(counts, dt = dt)
  graph <- alfak2::build_karyotype_graph(data, shell_depth = 1, max_nodes = 30000)
  if (identical(local_config, "best_fix")) {
    fit_local_variant_i(data, graph, "best_fix", "strong", "g_equivalent", 500, borrowed_residual_sd = 0.10)
  } else {
    fit_local_variant_i(data, graph, "current", "current", "f", 500)
  }
}

run_i6 <- function(source_input_dir, simulation_ids, minobs, input_policy, i5, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "i6_multisim_validation.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  global_rows <- list(); local_rows <- list(); norm_rows <- list(); runtime_rows <- list()
  gidx <- 0L; lidx <- 0L; ridx <- 0L
  cfgs <- data.frame(
    experiment = "I6",
    candidate_id = c("mutation_baseline", "normalized_default", "unit_stress"),
    graph_edge_weight = c("mutation", "normalized", "unit"),
    lambda_l = c(0.2, 0.2, 0.2),
    lambda_e = c(1, 0.01, 0.01),
    sigma_obs = c(0.05, 0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    stringsAsFactors = FALSE
  )
  for (sim in simulation_ids) {
    row <- resolve_shared_input_row(source_input_dir, sim, minobs)
    if (!nrow(row)) next
    yi <- readRDS(row$input_rds[[1L]])
    counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = input_policy, drop_diploid = TRUE)
    dt <- suppressWarnings(diff(as.numeric(colnames(counts))))
    if (length(dt) != 1L || !is.finite(dt) || dt <= 0) dt <- as.numeric(row$time_delta[[1L]])
    grf <- readRDS(row$grf_rds[[1L]])
    locals <- list()
    for (lc in c("current", "best_fix")) {
      started <- Sys.time()
      lf <- tryCatch(fit_multisim_local_i(counts, dt, lc), error = function(e) e)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      lidx <- lidx + 1L
      if (inherits(lf, "error")) {
        local_rows[[lidx]] <- data.frame(simulation_id = sim, local_config = lc, status = "error",
                                         error_message = conditionMessage(lf), elapsed_sec = elapsed,
                                         stringsAsFactors = FALSE)
      } else {
        locals[[lc]] <- lf
        local_rows[[lidx]] <- data.frame(
          simulation_id = sim, local_config = lc, status = "ok",
          convergence = lf$diagnostics$convergence,
          gradient_norm = lf$diagnostics$gradient_norm,
          covariance_status = lf$diagnostics$covariance_status,
          covariance_fallback = lf$diagnostics$covariance_fallback,
          edge_sign_agreement = edge_alignment_i(lf, grf, as.numeric(row$lambda[[1L]]), lc)$delta_sign_agreement,
          elapsed_sec = elapsed,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(locals$current)) next
    graph <- tryCatch(alfak2::build_karyotype_graph(locals$current$data, shell_depth = 2, max_nodes = 30000),
                      error = function(e) e)
    if (inherits(graph, "error")) next
    for (j in seq_len(nrow(cfgs))) {
      started <- Sys.time()
      fit <- tryCatch(
        alfak2::fit_graph_posterior(
          locals$current, graph,
          lambda_l_grid = cfgs$lambda_l[[j]],
          lambda_e_grid = cfgs$lambda_e[[j]],
          sigma_obs_grid = cfgs$sigma_obs[[j]],
          graph_edge_weight = cfgs$graph_edge_weight[[j]],
          compute_sd = FALSE
        ),
        error = function(e) e
      )
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      ridx <- ridx + 1L
      runtime_rows[[ridx]] <- data.frame(simulation_id = sim, candidate_id = cfgs$candidate_id[[j]],
                                         runtime_sec = elapsed, compute_sd = FALSE,
                                         status = if (inherits(fit, "error")) "error" else "ok",
                                         stringsAsFactors = FALSE)
      gidx <- gidx + 1L
      if (inherits(fit, "error")) {
        global_rows[[gidx]] <- cbind(data.frame(simulation_id = sim, support_scope = "farfield",
                                                metric_scale = "native", recommended_status = "global_error",
                                                failure_reason = conditionMessage(fit), stringsAsFactors = FALSE),
                                     cfgs[j, , drop = FALSE])
        next
      }
      task_info <- list(simulation_id = sim, minobs = minobs, input_policy = input_policy,
                        lambda = as.numeric(row$lambda[[1L]]), dt = dt,
                        beta = if ("sim_beta" %in% names(row)) as.numeric(row$sim_beta[[1L]]) else 0.00005)
      m <- score_summary_abcd(fit$summary, graph, grf, as.numeric(row$lambda[[1L]]), task_info, cfgs[j, , drop = FALSE], "mean_only_no_sd")
      far <- m[m$support_scope == "farfield" & m$metric_scale == "native", , drop = FALSE]
      far$shape_classification <- metric_shape_class(far)
      far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_i(far[k, , drop = FALSE]), character(1))
      far$failure_reason <- ifelse(far$recommended_status == "valid_shape_config", "passed_shape_gates", far$recommended_status)
      far$compute_sd <- FALSE
      far$runtime_sec <- elapsed
      global_rows[[gidx]] <- far
      if (cfgs$graph_edge_weight[[j]] == "normalized") norm_rows[[length(norm_rows) + 1L]] <- far
    }
  }
  global_tbl <- bind_rows_fill(global_rows)
  local_tbl <- bind_rows_fill(local_rows)
  norm_tbl <- bind_rows_fill(norm_rows)
  runtime_tbl <- bind_rows_fill(runtime_rows)
  write_tsv_safe(global_tbl, file.path(dirs$tables, "i6_multi_sim_failure_state_validation.tsv"))
  write_tsv_safe(local_tbl, file.path(dirs$tables, "i6_multi_sim_local_fix_validation.tsv"))
  write_tsv_safe(norm_tbl, file.path(dirs$tables, "i6_multi_sim_normalized_default_validation.tsv"))
  write_tsv_safe(runtime_tbl, file.path(dirs$tables, "i6_multi_sim_runtime.tsv"))
  out <- list(global = global_tbl, local = local_tbl, normalized = norm_tbl, runtime = runtime_tbl)
  saveRDS(out, rds)
  out
}

make_i_recommendations <- function(i1, i2, i3, i4, i5, i6) {
  data.frame(
    table = c("calibration_failure_state_recommendation", "global_mean_only_recommendation",
              "local_tmb_diagnostics_recommendation", "g_parameterization_recommendation",
              "borrowed_f_variant_recommendation", "multisim_validation_recommendation",
              "recommended_next_steps"),
    recommendation = c(
      "Keep formal no_valid_shape_configuration gate.",
      "Keep compute_sd=FALSE / posterior_sd_mode=none for benchmark-scale mean scoring.",
      "Keep return_optimizer_diagnostics for local TMB debugging.",
      "Keep g_equivalent as experimental; do not make default.",
      "Borrowed f variants did not establish trusted shell_depth=1 local fit.",
      "Keep normalized as default candidate with failure-state gate; do not promote collapsed/wrong-direction fits.",
      "Prioritize local f-block identifiability and calibration failure-state before edge-gradient."
    ),
    evidence = c(
      i1$gate_summary$recommended_status[[1L]],
      paste0("mean_diff=", fmt_metric(i2$equivalence$max_abs_mean_diff[[1L]]), "; speedup=", fmt_metric(i2$timing$speedup[i2$timing$mode == "compute_sd_FALSE"][[1L]])),
      paste0("max_block=", i3$block$max_gradient_block_name[which.max(i3$block$global_gradient_norm)]),
      paste0("best_g_gradient=", fmt_metric(min(i4$summary$gradient_norm[i4$summary$local_parameterization == "g_equivalent"], na.rm = TRUE))),
      paste0("best_variant=", i5$recommendation$best_gradient_variant[[1L]], "; status=", i5$recommendation$recommendation[[1L]]),
      paste0("normalized_collapse_fraction=", fmt_metric(mean(i6$normalized$amplitude_collapse, na.rm = TRUE))),
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_i_report <- function(dirs, args_info, ctx, i1, i2, i3, i4, i5, i6, recs) {
  all_long <- bind_rows_fill(list(
    transform(i1$ranked, experiment = "I1_calibration_failure_state"),
    transform(i2$metrics, experiment = "I2_global_mean_only"),
    transform(i3$block, experiment = "I3_local_tmb_diagnostics"),
    transform(i4$summary, experiment = "I4_g_parameterization"),
    transform(i5$result, experiment = "I5_borrowed_f_variants"),
    transform(i6$global, experiment = "I6_multisim_global"),
    transform(i6$local, experiment = "I6_multisim_local")
  ))
  summary <- data.frame(
    experiment = paste0("I", 1:6),
    key_result = c(
      i1$gate_summary$recommended_status[[1L]],
      paste0("mean_equivalent=", i2$equivalence$mean_equivalent[[1L]], "; speedup=", fmt_metric(i2$timing$speedup[i2$timing$mode == "compute_sd_FALSE"][[1L]])),
      paste0("max_gradient_block=", i3$block$max_gradient_block_name[which.max(i3$block$global_gradient_norm)]),
      paste0("best_g_gradient=", fmt_metric(min(i4$summary$gradient_norm[i4$summary$local_parameterization == "g_equivalent"], na.rm = TRUE))),
      i5$recommendation$recommendation[[1L]],
      paste0("n_full_global_sims=", length(unique(i6$global$simulation_id)))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_i_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "i_experiment_summary.tsv"))
  for (nm in recs$table) write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  i3_tier <- i3$by_tier[order(-i3$by_tier$grad_f_max_abs), , drop = FALSE][1L, , drop = FALSE]
  i4_best <- i4$summary[order(i4$summary$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  i5_best <- i5$result[order(i5$result$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  lines <- c(
    "# Farfield Core Fix Probe Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", args_info$source_input_dir, "`"),
    paste0("- abcd-dir: `", args_info$abcd_dir, "`"),
    paste0("- diagnostics-dir: `", args_info$diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", args_info$delta_probe_dir, "`"),
    paste0("- delta-debug-dir: `", args_info$delta_debug_dir, "`"),
    paste0("- local-calibration-dir: `", args_info$local_calibration_dir, "`"),
    paste0("- simulation_ids: ", paste(args_info$simulation_ids, collapse = ",")),
    paste0("- minobs: ", args_info$minobs),
    paste0("- input_policy: ", args_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## Prior results summary",
    "- ABCD/E/F/G/H show oracle edge information can recover farfield shape, but non-oracle delta deployment and local shell_depth=1 fits are not yet trustworthy.",
    "",
    "## I1 formal calibration failure-state patch",
    paste0("- recommended status: ", i1$gate_summary$recommended_status[[1L]], "."),
    "- Formal calibration keeps legacy outputs and adds gate summary, recommended status, and best_numeric_only_params.",
    "- `best_params_cli_args.txt` uses `no_valid_shape_configuration` when no non-oracle config passes shape gates.",
    "",
    "## I2 global mean-only",
    paste0("- compute_sd=FALSE mean max absolute difference: ", fmt_metric(i2$equivalence$max_abs_mean_diff[[1L]]), "."),
    paste0("- mean-only speedup on the smoke graph: ", fmt_metric(i2$timing$speedup[i2$timing$mode == "compute_sd_FALSE"][[1L]]), "."),
    "- When posterior sd is disabled, sd/confidence intervals are NA and prior_dominated tiers are skipped; shape metrics still score fitness_mean.",
    "",
    "## I3 local TMB diagnostics",
    paste0("- block diagnostics are available from package `fit_local_posterior()`: top block=", i3$block$max_gradient_block_name[which.max(i3$block$global_gradient_norm)], "."),
    paste0("- largest f-gradient tier: ", i3_tier$support_tier[[1L]], " / ", i3_tier$support_scope[[1L]], "."),
    "- Gradients remain concentrated in direct/local-borrowed f nodes; shrinkage alone does not make shell_depth=1 trustworthy.",
    "",
    "## I4 g = dt * f parameterization",
    paste0("- best variant: ", i4_best$variant[[1L]], " with gradient=", fmt_metric(i4_best$gradient_norm[[1L]]), "."),
    "- g_equivalent is now a real TMB parameterization and reports f-scale fitness, but covariance remains untrusted in this probe.",
    "",
    "## I5 borrowed f variants",
    paste0("- best borrowed variant: ", i5_best$variant[[1L]], " with gradient=", fmt_metric(i5_best$gradient_norm[[1L]]), "."),
    paste0("- recommendation: ", i5$recommendation$recommendation[[1L]], "."),
    "- Residual/deterministic borrowed penalties did not establish a trusted shell_depth=1 local fit.",
    "",
    "## I6 mean-only multi-sim validation",
    paste0("- full global scoring simulations: ", paste(sort(unique(i6$global$simulation_id)), collapse = ","), "."),
    paste0("- normalized collapse fraction: ", fmt_metric(mean(i6$normalized$amplitude_collapse, na.rm = TRUE)), "."),
    paste0("- any non-oracle valid shape: ", any(i6$global$recommended_status == "valid_shape_config", na.rm = TRUE), "."),
    "- Failure-state behavior remains necessary after mean-only full scoring.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- Keep normalized as benchmark/probe/calibration default candidate with amplitude-collapse diagnostics; keep unit as stress-test and mutation as legacy baseline.",
    "- Do not default `anchor_count_reference=minobs` for full input.",
    "- Keep `compute_sd=FALSE` and `return_optimizer_diagnostics=TRUE` as formal diagnostic/performance hooks.",
    "- Highest priorities: formal calibration gate tests, local f-block identifiability, g/residual local variants on more simulations, and then non-oracle delta estimator redesign."
  )
  writeLines(lines, file.path(dirs$root, "farfield_core_fix_probe_report.md"))
  saveRDS(list(i1 = i1, i2 = i2, i3 = i3, i4 = i4, i5 = i5, i6 = i6, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_core_fix_probe_all_results.rds"))
}

main_i <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage(); return(invisible(NULL))
  }
  mode <- match.arg(tolower(as.character(arg_value(args, "mode", "all"))),
                    c("prepare", "i1-calibration-failure-state", "i2-global-mean-only",
                      "i3-local-tmb-diagnostics", "i4-g-parameterization",
                      "i5-borrowed-f-variants", "i6-multisim-validation",
                      "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  delta_debug_dir <- as.character(arg_value(args, "delta_debug_dir", "benchmark/results/farfield_delta_debug"))
  local_calibration_dir <- as.character(arg_value(args, "local_calibration_dir", "benchmark/results/farfield_local_calibration_patch"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_core_fix_probe"))
  simulation_ids <- arg_integer_csv(args, "simulation_ids", 1:5)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_i_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, 1, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, 1, minobs, input_policy, force = FALSE)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(simulation_id = 1, minobs = minobs, input_policy = input_policy,
                    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
                    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
                    beta = if ("sim_beta" %in% names(ctx$input_table)) as.numeric(ctx$input_table$sim_beta[[1L]]) else 0.00005)
  saveRDS(list(context = ctx, simulation_ids = simulation_ids), file.path(dirs$results, "prepare_context.rds"))
  if (mode == "prepare") return(invisible(dirs$root))
  i1 <- if (mode %in% c("all", "i1-calibration-failure-state")) run_i1(local_calibration_dir, dirs, force = force) else readRDS(file.path(dirs$results, "i1_calibration_failure_state_patch.rds"))
  if (mode == "i1-calibration-failure-state") return(invisible(i1))
  i2 <- if (mode %in% c("all", "i2-global-mean-only")) run_i2(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "i2_global_mean_only_probe.rds"))
  if (mode == "i2-global-mean-only") return(invisible(i2))
  i3 <- if (mode %in% c("all", "i3-local-tmb-diagnostics")) run_i3(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "i3_local_tmb_diagnostics.rds"))
  if (mode == "i3-local-tmb-diagnostics") return(invisible(i3))
  i4 <- if (mode %in% c("all", "i4-g-parameterization")) run_i4(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "i4_g_parameterization_probe.rds"))
  if (mode == "i4-g-parameterization") return(invisible(i4))
  i5 <- if (mode %in% c("all", "i5-borrowed-f-variants")) run_i5(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "i5_borrowed_f_variants.rds"))
  if (mode == "i5-borrowed-f-variants") return(invisible(i5))
  i6 <- if (mode %in% c("all", "i6-multisim-validation")) run_i6(source_input_dir, simulation_ids, minobs, input_policy, i5, dirs, force = force) else readRDS(file.path(dirs$results, "i6_multisim_validation.rds"))
  if (mode == "i6-multisim-validation") return(invisible(i6))
  if (mode %in% c("all", "summarize")) {
    recs <- make_i_recommendations(i1, i2, i3, i4, i5, i6)
    args_info <- list(source_input_dir = source_input_dir, abcd_dir = abcd_dir,
                      diagnostics_dir = diagnostics_dir, delta_probe_dir = delta_probe_dir,
                      delta_debug_dir = delta_debug_dir, local_calibration_dir = local_calibration_dir,
                      simulation_ids = simulation_ids, minobs = minobs, input_policy = input_policy)
    write_i_report(dirs, args_info, ctx, i1, i2, i3, i4, i5, i6, recs)
  }
  message("Wrote farfield core fix probe under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) main_i()
