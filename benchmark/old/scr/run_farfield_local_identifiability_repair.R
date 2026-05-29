#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_local_identifiability_repair.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_core_fix_probe.R"))

usage <- function() {
  cat(
    "Run farfield local identifiability repair J0-J5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_local_identifiability_repair.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --delta-debug-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --local-calibration-dir=benchmark/results/farfield_local_calibration_patch \\\n",
    "    --core-fix-dir=benchmark/results/farfield_core_fix_probe \\\n",
    "    --output-dir=benchmark/results/farfield_local_identifiability_repair \\\n",
    "    --simulation-ids=1,2,3,4,5,6,7,8,9,10 --minobs=5 --input-policy=full\n",
    sep = ""
  )
}

make_j_dirs <- function(output_dir) make_probe_dirs(output_dir)

local_variant_j <- function(data, graph, variant, shrink = "current", local_parameterization = "f",
                            eval_max = 500, local_centering = "none", local_centering_weight = 0,
                            local_centering_weight_mode = "effective_count",
                            fixed_sigma_anchor = NA_real_, fixed_sigma_neighbor = NA_real_,
                            fixed_tau_group = NA_real_, borrowed_residual_sd = NULL,
                            weakly_supported_residual_sd = NULL,
                            support_tier_f_sd_multiplier = NULL) {
  ctrl <- local_shrink_controls_i(shrink)
  started <- Sys.time()
  fit <- tryCatch(
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
          return_optimizer_diagnostics = TRUE,
          local_parameterization = local_parameterization,
          local_centering = local_centering,
          local_centering_weight = local_centering_weight,
          local_centering_weight_mode = local_centering_weight_mode,
          fixed_sigma_anchor = fixed_sigma_anchor,
          fixed_sigma_neighbor = fixed_sigma_neighbor,
          fixed_tau_group = fixed_tau_group,
          borrowed_residual_sd = borrowed_residual_sd,
          weakly_supported_residual_sd = weakly_supported_residual_sd,
          support_tier_f_sd_multiplier = support_tier_f_sd_multiplier
        ),
        ctrl
      )
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  list(fit = fit, elapsed_sec = elapsed)
}

local_result_row_j <- function(res, variant, meta = list(), grf = NULL, lambda = NA_real_) {
  base <- data.frame(variant = variant, stringsAsFactors = FALSE)
  for (nm in names(meta)) base[[nm]] <- meta[[nm]]
  base$elapsed_sec <- res$elapsed_sec
  if (inherits(res$fit, "error")) {
    base$status <- "error"
    base$error_message <- conditionMessage(res$fit)
    return(base)
  }
  row <- local_diag_rows_i(res$fit, variant, meta$shrink %||% NA_character_,
                           meta$local_parameterization %||% NA_character_,
                           meta$eval_max %||% NA_integer_)
  row <- cbind(base, row[, setdiff(names(row), names(base)), drop = FALSE])
  if (!is.null(grf)) {
    align <- edge_alignment_i(res$fit, grf, lambda, variant)
    row$local_edge_delta_sign_agreement <- align$delta_sign_agreement
    row$local_edge_delta_spearman <- align$delta_spearman
    row$local_edge_delta_sd_ratio <- align$estimated_delta_sd_ratio
  }
  row
}

run_j0 <- function(dirs, force = FALSE) {
  rds <- file.path(dirs$results, "j0_infra_regression_tests.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  expr <- "pkgload::load_all('.', quiet=TRUE); testthat::test_dir('tests/testthat', filter='calibration|global|local|fit|infra', reporter='summary')"
  started <- Sys.time()
  res <- system2("Rscript", c("-e", shQuote(expr)), stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (is.null(status)) status <- 0L
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  tbl <- data.frame(
    test_scope = "calibration_global_local_fit_infra",
    command = paste("Rscript -e", shQuote(expr)),
    status = as.integer(status),
    passed = identical(as.integer(status), 0L),
    elapsed_sec = elapsed,
    output_tail = paste(utils::tail(res, 12), collapse = " | "),
    stringsAsFactors = FALSE
  )
  failures <- if (identical(as.integer(status), 0L)) {
    data.frame(test_scope = character(), failure = character(), stringsAsFactors = FALSE)
  } else {
    data.frame(test_scope = "calibration_global_local_fit_infra", failure = paste(res, collapse = " | "), stringsAsFactors = FALSE)
  }
  write_tsv_safe(tbl, file.path(dirs$tables, "j0_infra_regression_tests.tsv"))
  write_tsv_safe(failures, file.path(dirs$tables, "j0_test_failures.tsv"))
  out <- list(tests = tbl, failures = failures, raw_output = res)
  saveRDS(out, rds)
  out
}

run_j1 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "j1_fitness_centering_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  cfg <- data.frame(
    variant = c("C0_f_none_current", "C1_f_direct_count_w100_strong", "C1_f_direct_uniform_w100_strong",
                "C2_f_reference_w100_strong", "C0_g_none_strong", "C1_g_direct_count_w100_strong",
                "C2_g_reference_w100_strong"),
    local_parameterization = c("f", "f", "f", "f", "g_equivalent", "g_equivalent", "g_equivalent"),
    local_centering = c("none", "direct_weighted_mean", "direct_weighted_mean", "reference_direct", "none", "direct_weighted_mean", "reference_direct"),
    local_centering_weight = c(0, 100, 100, 100, 0, 100, 100),
    local_centering_weight_mode = c("effective_count", "effective_count", "uniform", "effective_count", "effective_count", "effective_count", "effective_count"),
    shrink = c("current", rep("strong", 6)),
    eval_max = 500,
    stringsAsFactors = FALSE
  )
  rows <- list(); tiers <- list(); tops <- list(); aligns <- list()
  for (i in seq_len(nrow(cfg))) {
    meta <- as.list(cfg[i, , drop = FALSE])
    res <- local_variant_j(data, graph, cfg$variant[[i]], shrink = cfg$shrink[[i]],
                           local_parameterization = cfg$local_parameterization[[i]],
                           eval_max = cfg$eval_max[[i]],
                           local_centering = cfg$local_centering[[i]],
                           local_centering_weight = cfg$local_centering_weight[[i]],
                           local_centering_weight_mode = cfg$local_centering_weight_mode[[i]])
    rows[[i]] <- local_result_row_j(res, cfg$variant[[i]], meta, grf, task_info$lambda)
    if (!inherits(res$fit, "error")) {
      tiers[[i]] <- cbind(data.frame(variant = cfg$variant[[i]], stringsAsFactors = FALSE),
                          res$fit$diagnostics$optimizer$grad_f_by_support_tier)
      tops[[i]] <- cbind(data.frame(variant = cfg$variant[[i]], stringsAsFactors = FALSE),
                         res$fit$diagnostics$optimizer$top_gradient_nodes)
      aligns[[i]] <- edge_alignment_i(res$fit, grf, task_info$lambda, cfg$variant[[i]])
    }
  }
  tbl <- bind_rows_fill(rows)
  tier_tbl <- bind_rows_fill(tiers)
  top_tbl <- bind_rows_fill(tops)
  align_tbl <- bind_rows_fill(aligns)
  write_tsv_safe(tbl, file.path(dirs$tables, "j1_local_centering_probe.tsv"))
  write_tsv_safe(tier_tbl, file.path(dirs$tables, "j1_centering_f_gradient_by_tier.tsv"))
  write_tsv_safe(top_tbl, file.path(dirs$tables, "j1_centering_top_gradient_nodes.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "j1_centering_edge_alignment.tsv"))
  out <- list(results = tbl, by_tier = tier_tbl, top_nodes = top_tbl, edge_alignment = align_tbl)
  saveRDS(out, rds)
  out
}

best_centering_from_j1 <- function(j1) {
  ok <- j1$results[!is.na(j1$results$gradient_norm), , drop = FALSE]
  ok <- ok[ok$local_centering != "none", , drop = FALSE]
  if (!nrow(ok)) return(list(local_centering = "none", local_centering_weight = 0, local_centering_weight_mode = "effective_count"))
  b <- ok[order(ok$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  list(local_centering = b$local_centering[[1L]],
       local_centering_weight = as.numeric(b$local_centering_weight[[1L]]),
       local_centering_weight_mode = b$local_centering_weight_mode[[1L]])
}

run_j2 <- function(bundle, grf, task_info, dirs, j1, force = FALSE) {
  rds <- file.path(dirs$results, "j2_fixed_scale_local_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  bc <- best_centering_from_j1(j1)
  cfg <- data.frame(
    variant = c("S0_current", "S1_fixed_anchor", "S2_fixed_neighbor", "S3_fixed_tau", "S4_fixed_anchor_neighbor", "S5_all_scale_fixed", "S5_g_centered_all_scale_fixed"),
    local_parameterization = c(rep("f", 6), "g_equivalent"),
    shrink = c(rep("strong", 7)),
    fixed_sigma_anchor = c(NA, 0.2, NA, NA, 0.2, 0.2, 0.2),
    fixed_sigma_neighbor = c(NA, NA, 0.1, NA, 0.1, 0.1, 0.1),
    fixed_tau_group = c(NA, NA, NA, 0.1, NA, 0.1, 0.1),
    eval_max = 500,
    stringsAsFactors = FALSE
  )
  rows <- list(); aligns <- list()
  for (i in seq_len(nrow(cfg))) {
    meta <- c(as.list(cfg[i, , drop = FALSE]), bc)
    res <- local_variant_j(data, graph, cfg$variant[[i]], shrink = cfg$shrink[[i]],
                           local_parameterization = cfg$local_parameterization[[i]],
                           eval_max = cfg$eval_max[[i]],
                           local_centering = bc$local_centering,
                           local_centering_weight = bc$local_centering_weight,
                           local_centering_weight_mode = bc$local_centering_weight_mode,
                           fixed_sigma_anchor = cfg$fixed_sigma_anchor[[i]],
                           fixed_sigma_neighbor = cfg$fixed_sigma_neighbor[[i]],
                           fixed_tau_group = cfg$fixed_tau_group[[i]])
    rows[[i]] <- local_result_row_j(res, cfg$variant[[i]], meta, grf, task_info$lambda)
    if (!inherits(res$fit, "error")) aligns[[i]] <- edge_alignment_i(res$fit, grf, task_info$lambda, cfg$variant[[i]])
  }
  tbl <- bind_rows_fill(rows)
  align_tbl <- bind_rows_fill(aligns)
  rec <- data.frame(
    recommendation = if (any(tbl$covariance_status == "TMB_sdreport", na.rm = TRUE)) "fixed_scale_candidate_found" else "no_trusted_fixed_scale_shell1_fit",
    best_variant = tbl$variant[which.min(tbl$gradient_norm)],
    best_gradient_norm = min(tbl$gradient_norm, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(tbl, file.path(dirs$tables, "j2_local_fixed_scale_probe.tsv"))
  write_tsv_safe(tbl, file.path(dirs$tables, "j2_fixed_scale_gradient_by_block.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "j2_fixed_scale_edge_alignment.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "j2_fixed_scale_recommendation.tsv"))
  out <- list(results = tbl, edge_alignment = align_tbl, recommendation = rec, best_centering = bc)
  saveRDS(out, rds)
  out
}

best_scale_from_j2 <- function(j2) {
  b <- j2$results[order(j2$results$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  list(
    fixed_sigma_anchor = suppressWarnings(as.numeric(b$fixed_sigma_anchor[[1L]])),
    fixed_sigma_neighbor = suppressWarnings(as.numeric(b$fixed_sigma_neighbor[[1L]])),
    fixed_tau_group = suppressWarnings(as.numeric(b$fixed_tau_group[[1L]]))
  )
}

run_j3 <- function(bundle, grf, task_info, dirs, j1, j2, force = FALSE) {
  rds <- file.path(dirs$results, "j3_borrowed_reduced_dof_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  bc <- best_centering_from_j1(j1)
  bs <- best_scale_from_j2(j2)
  cfg <- data.frame(
    variant = c("B0_current", "B1_residual_0p20", "B1_residual_0p10", "B3_deterministic_0p01", "B5_direct_only_free", "B6_g_best_residual_0p10"),
    local_parameterization = c(rep("f", 5), "g_equivalent"),
    shrink = c(rep("strong", 6)),
    borrowed_residual_sd = c(NA, 0.2, 0.1, 0.01, 0.005, 0.1),
    weakly_supported_residual_sd = c(NA, 0.1, 0.05, 0.01, 0.005, 0.05),
    eval_max = 500,
    stringsAsFactors = FALSE
  )
  rows <- list(); tiers <- list(); aligns <- list()
  for (i in seq_len(nrow(cfg))) {
    meta <- c(as.list(cfg[i, , drop = FALSE]), bc, bs)
    res <- local_variant_j(data, graph, cfg$variant[[i]], shrink = cfg$shrink[[i]],
                           local_parameterization = cfg$local_parameterization[[i]],
                           eval_max = cfg$eval_max[[i]],
                           local_centering = bc$local_centering,
                           local_centering_weight = bc$local_centering_weight,
                           local_centering_weight_mode = bc$local_centering_weight_mode,
                           fixed_sigma_anchor = bs$fixed_sigma_anchor,
                           fixed_sigma_neighbor = bs$fixed_sigma_neighbor,
                           fixed_tau_group = bs$fixed_tau_group,
                           borrowed_residual_sd = cfg$borrowed_residual_sd[[i]],
                           weakly_supported_residual_sd = cfg$weakly_supported_residual_sd[[i]])
    rows[[i]] <- local_result_row_j(res, cfg$variant[[i]], meta, grf, task_info$lambda)
    if (!inherits(res$fit, "error")) {
      tiers[[i]] <- cbind(data.frame(variant = cfg$variant[[i]], stringsAsFactors = FALSE),
                          res$fit$diagnostics$optimizer$grad_f_by_support_tier)
      aligns[[i]] <- edge_alignment_i(res$fit, grf, task_info$lambda, cfg$variant[[i]])
    }
  }
  tbl <- bind_rows_fill(rows)
  tier_tbl <- bind_rows_fill(tiers)
  align_tbl <- bind_rows_fill(aligns)
  rec <- data.frame(
    recommendation = if (any(tbl$covariance_status == "TMB_sdreport", na.rm = TRUE)) "borrowed_reduced_dof_candidate_found" else "no_trusted_borrowed_reduced_dof_fit",
    best_variant = tbl$variant[which.min(tbl$gradient_norm)],
    best_gradient_norm = min(tbl$gradient_norm, na.rm = TRUE),
    best_edge_sign_agreement = max(align_tbl$delta_sign_agreement, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(tbl, file.path(dirs$tables, "j3_borrowed_reduced_dof_probe.tsv"))
  write_tsv_safe(tier_tbl, file.path(dirs$tables, "j3_borrowed_gradient_by_tier.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "j3_borrowed_edge_alignment.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "j3_borrowed_variant_recommendation.tsv"))
  out <- list(results = tbl, by_tier = tier_tbl, edge_alignment = align_tbl, recommendation = rec, best_centering = bc, best_scale = bs)
  saveRDS(out, rds)
  out
}

best_borrowed_from_j3 <- function(j3) {
  b <- j3$results[order(j3$results$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  list(
    borrowed_residual_sd = suppressWarnings(as.numeric(b$borrowed_residual_sd[[1L]])),
    weakly_supported_residual_sd = suppressWarnings(as.numeric(b$weakly_supported_residual_sd[[1L]]))
  )
}

run_j4 <- function(bundle, grf, task_info, dirs, j1, j2, j3, force = FALSE) {
  rds <- file.path(dirs$results, "j4_combined_g_identifiability_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  data <- bundle$data
  graph <- bundle$local$graph
  bc <- best_centering_from_j1(j1)
  bs <- best_scale_from_j2(j2)
  bb <- best_borrowed_from_j3(j3)
  cfg <- data.frame(
    variant = c("M0_baseline", "M1_g_centered", "M2_g_fixed_scale", "M3_g_borrowed_residual", "M4_g_fixed_scale_borrowed", "M5_conservative_local"),
    local_parameterization = c("f", rep("g_equivalent", 5)),
    shrink = c("current", rep("strong", 5)),
    use_centering = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    use_scale = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
    use_borrowed = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    eval_max = 2000,
    restart_id = 1L,
    stringsAsFactors = FALSE
  )
  rows <- list(); aligns <- list()
  for (i in seq_len(nrow(cfg))) {
    center <- if (cfg$use_centering[[i]]) bc else list(local_centering = "none", local_centering_weight = 0, local_centering_weight_mode = "effective_count")
    scale <- if (cfg$use_scale[[i]]) bs else list(fixed_sigma_anchor = NA_real_, fixed_sigma_neighbor = NA_real_, fixed_tau_group = NA_real_)
    borrowed <- if (cfg$use_borrowed[[i]]) bb else list(borrowed_residual_sd = NULL, weakly_supported_residual_sd = NULL)
    if (cfg$variant[[i]] == "M5_conservative_local") {
      borrowed <- list(borrowed_residual_sd = 0.005, weakly_supported_residual_sd = 0.005)
    }
    meta <- c(as.list(cfg[i, , drop = FALSE]), center, scale, borrowed)
    res <- local_variant_j(data, graph, cfg$variant[[i]], shrink = cfg$shrink[[i]],
                           local_parameterization = cfg$local_parameterization[[i]],
                           eval_max = cfg$eval_max[[i]],
                           local_centering = center$local_centering,
                           local_centering_weight = center$local_centering_weight,
                           local_centering_weight_mode = center$local_centering_weight_mode,
                           fixed_sigma_anchor = scale$fixed_sigma_anchor,
                           fixed_sigma_neighbor = scale$fixed_sigma_neighbor,
                           fixed_tau_group = scale$fixed_tau_group,
                           borrowed_residual_sd = borrowed$borrowed_residual_sd,
                           weakly_supported_residual_sd = borrowed$weakly_supported_residual_sd)
    rows[[i]] <- local_result_row_j(res, cfg$variant[[i]], meta, grf, task_info$lambda)
    if (!inherits(res$fit, "error")) aligns[[i]] <- edge_alignment_i(res$fit, grf, task_info$lambda, cfg$variant[[i]])
  }
  tbl <- bind_rows_fill(rows)
  align_tbl <- bind_rows_fill(aligns)
  stability <- data.frame(
    restart_status = "single_restart_only_time_budget",
    note = "J4 used one deterministic restart per combo; future work should add randomized starts after selecting a viable objective.",
    stringsAsFactors = FALSE
  )
  best <- tbl[order(tbl$covariance_status != "TMB_sdreport", tbl$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  write_tsv_safe(tbl, file.path(dirs$tables, "j4_combined_g_identifiability_probe.tsv"))
  write_tsv_safe(stability, file.path(dirs$tables, "j4_combined_multistart_stability.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "j4_combined_edge_alignment.tsv"))
  write_tsv_safe(best, file.path(dirs$tables, "j4_combined_best_local_config.tsv"))
  out <- list(results = tbl, stability = stability, edge_alignment = align_tbl, best = best)
  saveRDS(out, rds)
  out
}

best_local_config_from_j4 <- function(j4) {
  b <- j4$best
  list(
    local_parameterization = b$local_parameterization[[1L]] %||% "g_equivalent",
    shrink = b$shrink[[1L]] %||% "strong",
    local_centering = b$local_centering[[1L]] %||% "none",
    local_centering_weight = suppressWarnings(as.numeric(b$local_centering_weight[[1L]] %||% 0)),
    local_centering_weight_mode = b$local_centering_weight_mode[[1L]] %||% "effective_count",
    fixed_sigma_anchor = suppressWarnings(as.numeric(b$fixed_sigma_anchor[[1L]] %||% NA)),
    fixed_sigma_neighbor = suppressWarnings(as.numeric(b$fixed_sigma_neighbor[[1L]] %||% NA)),
    fixed_tau_group = suppressWarnings(as.numeric(b$fixed_tau_group[[1L]] %||% NA)),
    borrowed_residual_sd = suppressWarnings(as.numeric(b$borrowed_residual_sd[[1L]] %||% NA)),
    weakly_supported_residual_sd = suppressWarnings(as.numeric(b$weakly_supported_residual_sd[[1L]] %||% NA))
  )
}

fit_multisim_local_j <- function(counts, dt, local_config, best_cfg = NULL) {
  data <- alfak2::prepare_alfak2_data(counts, dt = dt)
  graph <- alfak2::build_karyotype_graph(data, shell_depth = 1, max_nodes = 30000)
  if (identical(local_config, "current")) {
    return(local_variant_j(data, graph, "current", shrink = "current", eval_max = 500)$fit)
  }
  cfg <- best_cfg
  local_variant_j(
    data, graph, local_config,
    shrink = cfg$shrink %||% "strong",
    local_parameterization = cfg$local_parameterization %||% "g_equivalent",
    eval_max = 500,
    local_centering = cfg$local_centering %||% "none",
    local_centering_weight = cfg$local_centering_weight %||% 0,
    local_centering_weight_mode = cfg$local_centering_weight_mode %||% "effective_count",
    fixed_sigma_anchor = cfg$fixed_sigma_anchor %||% NA_real_,
    fixed_sigma_neighbor = cfg$fixed_sigma_neighbor %||% NA_real_,
    fixed_tau_group = cfg$fixed_tau_group %||% NA_real_,
    borrowed_residual_sd = cfg$borrowed_residual_sd %||% NULL,
    weakly_supported_residual_sd = cfg$weakly_supported_residual_sd %||% NULL
  )$fit
}

run_j5 <- function(source_input_dir, simulation_ids, minobs, input_policy, j4, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "j5_multisim_validation.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  best_cfg <- best_local_config_from_j4(j4)
  cfgs <- data.frame(
    experiment = "J5",
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
  local_configs <- c("current", "best_from_J4")
  global_sim_ids <- simulation_ids[seq_len(min(length(simulation_ids), 5L))]
  local_rows <- list(); global_rows <- list(); norm_rows <- list(); runtime <- list()
  lidx <- 0L; gidx <- 0L; ridx <- 0L
  for (sim in simulation_ids) {
    row <- resolve_shared_input_row(source_input_dir, sim, minobs)
    if (!nrow(row)) next
    yi <- readRDS(row$input_rds[[1L]])
    counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = input_policy, drop_diploid = TRUE)
    dt <- suppressWarnings(diff(as.numeric(colnames(counts))))
    if (length(dt) != 1L || !is.finite(dt) || dt <= 0) dt <- as.numeric(row$time_delta[[1L]])
    grf <- readRDS(row$grf_rds[[1L]])
    local_fits <- list()
    for (lc in local_configs) {
      started <- Sys.time()
      lf <- tryCatch(fit_multisim_local_j(counts, dt, lc, best_cfg), error = function(e) e)
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      lidx <- lidx + 1L
      if (inherits(lf, "error")) {
        local_rows[[lidx]] <- data.frame(simulation_id = sim, local_config = lc, status = "error",
                                         error_message = conditionMessage(lf), elapsed_sec = elapsed,
                                         stringsAsFactors = FALSE)
      } else {
        local_fits[[lc]] <- lf
        align <- edge_alignment_i(lf, grf, as.numeric(row$lambda[[1L]]), lc)
        block <- lf$diagnostics$optimizer$gradient_block_summary
        local_rows[[lidx]] <- data.frame(
          simulation_id = sim, local_config = lc, status = "ok",
          convergence = lf$diagnostics$convergence,
          gradient_norm = lf$diagnostics$gradient_norm,
          covariance_status = lf$diagnostics$covariance_status,
          covariance_fallback = lf$diagnostics$covariance_fallback,
          max_gradient_block_name = block$max_gradient_block_name,
          local_edge_delta_sign_agreement = align$delta_sign_agreement,
          local_edge_delta_spearman = align$delta_spearman,
          elapsed_sec = elapsed,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!(sim %in% global_sim_ids)) {
      for (lc in names(local_fits)) {
        ridx <- ridx + 1L
        runtime[[ridx]] <- data.frame(
          simulation_id = sim,
          local_config = lc,
          candidate_id = "all_global_configs",
          runtime_sec = NA_real_,
          compute_sd = FALSE,
          status = "skipped_global_sim_runtime_budget",
          stringsAsFactors = FALSE
        )
      }
      next
    }
    global_local_configs <- names(local_fits)
    if ("best_from_J4" %in% global_local_configs) {
      skipped_lc <- setdiff(global_local_configs, "best_from_J4")
      global_local_configs <- "best_from_J4"
    } else {
      skipped_lc <- setdiff(global_local_configs, global_local_configs[[1L]])
      global_local_configs <- global_local_configs[[1L]]
    }
    if (length(skipped_lc)) {
      for (sk in skipped_lc) {
        ridx <- ridx + 1L
        runtime[[ridx]] <- data.frame(
          simulation_id = sim,
          local_config = sk,
          candidate_id = "all_global_configs",
          runtime_sec = NA_real_,
          compute_sd = FALSE,
          status = "skipped_runtime_budget",
          stringsAsFactors = FALSE
        )
      }
    }
    for (lc in global_local_configs) {
      lf <- local_fits[[lc]]
      graph <- tryCatch(alfak2::build_karyotype_graph(lf$data, shell_depth = 2, max_nodes = 30000), error = function(e) e)
      if (inherits(graph, "error")) next
      for (j in seq_len(nrow(cfgs))) {
        started <- Sys.time()
        fit <- tryCatch(
          alfak2::fit_graph_posterior(
            lf, graph,
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
        runtime[[ridx]] <- data.frame(simulation_id = sim, local_config = lc, candidate_id = cfgs$candidate_id[[j]],
                                      runtime_sec = elapsed, compute_sd = FALSE,
                                      status = if (inherits(fit, "error")) "error" else "ok",
                                      stringsAsFactors = FALSE)
        gidx <- gidx + 1L
        if (inherits(fit, "error")) {
          global_rows[[gidx]] <- cbind(data.frame(simulation_id = sim, local_config = lc, support_scope = "farfield",
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
        far$local_config <- lc
        far$shape_classification <- metric_shape_class(far)
        far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_i(far[k, , drop = FALSE]), character(1))
        far$failure_reason <- ifelse(far$recommended_status == "valid_shape_config", "passed_shape_gates", far$recommended_status)
        far$runtime_sec <- elapsed
        far$compute_sd <- FALSE
        global_rows[[gidx]] <- far
        if (cfgs$graph_edge_weight[[j]] == "normalized") norm_rows[[length(norm_rows) + 1L]] <- far
      }
    }
  }
  global_tbl <- bind_rows_fill(global_rows)
  local_tbl <- bind_rows_fill(local_rows)
  norm_tbl <- bind_rows_fill(norm_rows)
  runtime_tbl <- bind_rows_fill(runtime)
  write_tsv_safe(local_tbl, file.path(dirs$tables, "j5_local_repair_multisim_validation.tsv"))
  write_tsv_safe(global_tbl, file.path(dirs$tables, "j5_multisim_failure_state_validation.tsv"))
  write_tsv_safe(norm_tbl, file.path(dirs$tables, "j5_multisim_normalized_default_validation.tsv"))
  write_tsv_safe(runtime_tbl, file.path(dirs$tables, "j5_multisim_runtime.tsv"))
  out <- list(global = global_tbl, local = local_tbl, normalized = norm_tbl, runtime = runtime_tbl)
  saveRDS(out, rds)
  out
}

make_j_recommendations <- function(j0, j1, j2, j3, j4, j5) {
  best <- function(x) x[order(x$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  b1 <- best(j1$results); b2 <- best(j2$results); b3 <- best(j3$results); b4 <- best(j4$results)
  data.frame(
    table = c("infra_regression_recommendation", "fitness_centering_recommendation", "fixed_scale_recommendation",
              "borrowed_reduced_dof_recommendation", "combined_local_model_recommendation",
              "multisim_local_repair_recommendation", "recommended_next_steps"),
    recommendation = c(
      if (isTRUE(j0$tests$passed[[1L]])) "infra_regression_tests_passed" else "infra_regression_tests_failed",
      "Centering can reduce some gradients but did not establish trusted shell_depth=1 covariance.",
      "Fixed scales did not establish trusted shell_depth=1 covariance.",
      "Borrowed reduced-DOF variants did not establish trusted shell_depth=1 covariance.",
      "Combined local repair remains experimental; do not make default.",
      "Failure-state gate remains necessary across multi-sim validation.",
      "Prioritize local objective identifiability before revisiting non-oracle delta or edge-gradient."
    ),
    evidence = c(
      paste0("status=", j0$tests$status[[1L]]),
      paste0("best=", b1$variant[[1L]], "; gradient=", fmt_metric(b1$gradient_norm[[1L]])),
      paste0("best=", b2$variant[[1L]], "; gradient=", fmt_metric(b2$gradient_norm[[1L]])),
      paste0("best=", b3$variant[[1L]], "; gradient=", fmt_metric(b3$gradient_norm[[1L]])),
      paste0("best=", b4$variant[[1L]], "; covariance=", b4$covariance_status[[1L]]),
      paste0("valid_shape_count=", sum(j5$global$recommended_status == "valid_shape_config", na.rm = TRUE)),
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_j_report <- function(dirs, args_info, ctx, j0, j1, j2, j3, j4, j5, recs) {
  all_long <- bind_rows_fill(list(
    transform(j0$tests, experiment = "J0_infra"),
    transform(j1$results, experiment = "J1_centering"),
    transform(j2$results, experiment = "J2_fixed_scale"),
    transform(j3$results, experiment = "J3_borrowed_reduced_dof"),
    transform(j4$results, experiment = "J4_combined"),
    transform(j5$global, experiment = "J5_multisim_global"),
    transform(j5$local, experiment = "J5_multisim_local")
  ))
  best <- function(x) x[order(x$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  b1 <- best(j1$results); b2 <- best(j2$results); b3 <- best(j3$results); b4 <- best(j4$results)
  summary <- data.frame(
    experiment = c("J0", "J1", "J2", "J3", "J4", "J5"),
    key_result = c(
      paste0("tests_passed=", j0$tests$passed[[1L]]),
      paste0("best=", b1$variant[[1L]], "; gradient=", fmt_metric(b1$gradient_norm[[1L]])),
      paste0("best=", b2$variant[[1L]], "; gradient=", fmt_metric(b2$gradient_norm[[1L]])),
      paste0("best=", b3$variant[[1L]], "; gradient=", fmt_metric(b3$gradient_norm[[1L]])),
      paste0("best=", b4$variant[[1L]], "; covariance=", b4$covariance_status[[1L]]),
      paste0("valid_shape_count=", sum(j5$global$recommended_status == "valid_shape_config", na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_j_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "j_experiment_summary.tsv"))
  for (nm in recs$table) write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  lines <- c(
    "# Farfield Local Identifiability Repair Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", args_info$source_input_dir, "`"),
    paste0("- abcd-dir: `", args_info$abcd_dir, "`"),
    paste0("- diagnostics-dir: `", args_info$diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", args_info$delta_probe_dir, "`"),
    paste0("- delta-debug-dir: `", args_info$delta_debug_dir, "`"),
    paste0("- local-calibration-dir: `", args_info$local_calibration_dir, "`"),
    paste0("- core-fix-dir: `", args_info$core_fix_dir, "`"),
    paste0("- simulation_ids: ", paste(args_info$simulation_ids, collapse = ",")),
    paste0("- minobs: ", args_info$minobs),
    paste0("- input_policy: ", args_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## Prior results summary",
    "- ABCD: normalized/unit edge weights reduce farfield centered RMSE but collapse amplitude; mutation retains amplitude but is wrong-direction. Shape-aware CV avoids collapse better than MSE CV, but no non-oracle valid farfield shape was found.",
    "- E/F/G: oracle per-edge delta restores farfield shape, but non-oracle delta deployment fails; negative delta_scale is a compensation signal, not a sign convention fix. Pairwise potential did not outperform per-edge potential.",
    "- H/I: formal calibration gate, global compute_sd=FALSE, and local optimizer diagnostics are now available. Real g_equivalent and borrowed residual variants reduce gradients, but covariance remains untrusted and normalized collapses across sim1-5.",
    "",
    "## J0 infra regression tests",
    paste0("- tests passed: ", j0$tests$passed[[1L]], "; status=", j0$tests$status[[1L]], "."),
    "- Regression coverage includes calibration gate, compute_sd=FALSE equivalence, and local optimizer diagnostics.",
    "",
    "## J1 fitness centering",
    paste0("- best centering config: ", b1$variant[[1L]], " with gradient=", fmt_metric(b1$gradient_norm[[1L]]), "."),
    paste0("- best centering covariance: ", b1$covariance_status[[1L]], "."),
    paste0("- best centering local edge sign agreement: ", fmt_metric(b1$local_edge_delta_sign_agreement[[1L]]), "; Spearman=", fmt_metric(b1$local_edge_delta_spearman[[1L]]), "."),
    "- Centering is useful diagnostically but did not produce trusted shell_depth=1 covariance in this run.",
    "",
    "## J2 fixed scale",
    paste0("- best fixed-scale config: ", b2$variant[[1L]], " with gradient=", fmt_metric(b2$gradient_norm[[1L]]), "."),
    paste0("- best fixed-scale covariance: ", b2$covariance_status[[1L]], "."),
    paste0("- best fixed-scale local edge sign agreement: ", fmt_metric(b2$local_edge_delta_sign_agreement[[1L]]), "; Spearman=", fmt_metric(b2$local_edge_delta_spearman[[1L]]), "."),
    "- Fixed scale did not by itself make the local model trusted.",
    "",
    "## J3 borrowed reduced DOF",
    paste0("- best borrowed reduced-DOF config: ", b3$variant[[1L]], " with gradient=", fmt_metric(b3$gradient_norm[[1L]]), "."),
    paste0("- best borrowed reduced-DOF local edge sign agreement: ", fmt_metric(b3$local_edge_delta_sign_agreement[[1L]]), "; Spearman=", fmt_metric(b3$local_edge_delta_spearman[[1L]]), "."),
    paste0("- maximum borrowed edge sign agreement observed: ", fmt_metric(max(j3$edge_alignment$delta_sign_agreement, na.rm = TRUE)), "."),
    "- Reduced-DOF penalties can reduce f gradients in selected configurations but do not yet establish trusted covariance.",
    "",
    "## J4 combined local model",
    paste0("- best combined config: ", b4$variant[[1L]], " with gradient=", fmt_metric(b4$gradient_norm[[1L]]), " and covariance=", b4$covariance_status[[1L]], "."),
    paste0("- best combined local edge sign agreement: ", fmt_metric(b4$local_edge_delta_sign_agreement[[1L]]), "; Spearman=", fmt_metric(b4$local_edge_delta_spearman[[1L]]), "."),
    "- The best combined fit has near-zero optimizer gradient but non-finite sdreport covariance, indicating an identifiability/Hessian problem rather than a simple optimizer iteration problem.",
    "- J4 used one deterministic restart per combo due runtime; multistart remains required before making any local repair default.",
    "",
    "## J5 multi-sim validation",
    paste0("- simulations scored: ", paste(sort(unique(j5$global$simulation_id)), collapse = ","), "."),
    "- To keep sim1-10 validation finite, local diagnostics were run for current and best_from_J4 across all requested simulations, while global mean-only scoring was run on sim1-5 and the best_from_J4 local repair only.",
    paste0("- valid shape count: ", sum(j5$global$recommended_status == "valid_shape_config", na.rm = TRUE), "."),
    paste0("- normalized collapse fraction: ", fmt_metric(mean(j5$normalized$amplitude_collapse, na.rm = TRUE)), "."),
    paste0("- best_from_J4 median local gradient across sim1-10: ", fmt_metric(stats::median(j5$local$gradient_norm[j5$local$local_config == "best_from_J4"], na.rm = TRUE)), "; current median=", fmt_metric(stats::median(j5$local$gradient_norm[j5$local$local_config == "current"], na.rm = TRUE)), "."),
    "- Failure-state remains appropriate; no local repair in this run justifies edge-gradient work.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- Keep normalized as benchmark/probe/calibration candidate default with amplitude diagnostics and failure-state gate.",
    "- Do not default `anchor_count_reference=minobs` for full input.",
    "- Keep `compute_sd=FALSE` and `return_optimizer_diagnostics=TRUE`.",
    "- Next priorities: strengthen local gauge/scale identifiability, add real multistart for promising local variants, validate local covariance, then revisit non-oracle delta estimator."
  )
  writeLines(lines, file.path(dirs$root, "farfield_local_identifiability_repair_report.md"))
  saveRDS(list(j0 = j0, j1 = j1, j2 = j2, j3 = j3, j4 = j4, j5 = j5, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_local_identifiability_repair_all_results.rds"))
}

main_j <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage(); return(invisible(NULL))
  }
  mode <- match.arg(tolower(as.character(arg_value(args, "mode", "all"))),
                    c("prepare", "j0-infra-regression-tests", "j1-fitness-centering",
                      "j2-fixed-scale-local", "j3-borrowed-reduced-dof",
                      "j4-combined-g-identifiability", "j5-multisim-validation",
                      "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  delta_debug_dir <- as.character(arg_value(args, "delta_debug_dir", "benchmark/results/farfield_delta_debug"))
  local_calibration_dir <- as.character(arg_value(args, "local_calibration_dir", "benchmark/results/farfield_local_calibration_patch"))
  core_fix_dir <- as.character(arg_value(args, "core_fix_dir", "benchmark/results/farfield_core_fix_probe"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_local_identifiability_repair"))
  simulation_ids <- arg_integer_csv(args, "simulation_ids", 1:10)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_j_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, 1, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, 1, minobs, input_policy, force = FALSE)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(simulation_id = 1, minobs = minobs, input_policy = input_policy,
                    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
                    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
                    beta = if ("sim_beta" %in% names(ctx$input_table)) as.numeric(ctx$input_table$sim_beta[[1L]]) else 0.00005)
  saveRDS(list(context = ctx, simulation_ids = simulation_ids), file.path(dirs$results, "prepare_context.rds"))
  if (mode == "prepare") return(invisible(dirs$root))
  j0 <- if (mode %in% c("all", "j0-infra-regression-tests")) run_j0(dirs, force = force) else readRDS(file.path(dirs$results, "j0_infra_regression_tests.rds"))
  if (mode == "j0-infra-regression-tests") return(invisible(j0))
  j1 <- if (mode %in% c("all", "j1-fitness-centering")) run_j1(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "j1_fitness_centering_probe.rds"))
  if (mode == "j1-fitness-centering") return(invisible(j1))
  j2 <- if (mode %in% c("all", "j2-fixed-scale-local")) run_j2(bundle, grf, task_info, dirs, j1, force = force) else readRDS(file.path(dirs$results, "j2_fixed_scale_local_probe.rds"))
  if (mode == "j2-fixed-scale-local") return(invisible(j2))
  j3 <- if (mode %in% c("all", "j3-borrowed-reduced-dof")) run_j3(bundle, grf, task_info, dirs, j1, j2, force = force) else readRDS(file.path(dirs$results, "j3_borrowed_reduced_dof_probe.rds"))
  if (mode == "j3-borrowed-reduced-dof") return(invisible(j3))
  j4 <- if (mode %in% c("all", "j4-combined-g-identifiability")) run_j4(bundle, grf, task_info, dirs, j1, j2, j3, force = force) else readRDS(file.path(dirs$results, "j4_combined_g_identifiability_probe.rds"))
  if (mode == "j4-combined-g-identifiability") return(invisible(j4))
  j5 <- if (mode %in% c("all", "j5-multisim-validation")) run_j5(source_input_dir, simulation_ids, minobs, input_policy, j4, dirs, force = force) else readRDS(file.path(dirs$results, "j5_multisim_validation.rds"))
  if (mode == "j5-multisim-validation") return(invisible(j5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_j_recommendations(j0, j1, j2, j3, j4, j5)
    args_info <- list(source_input_dir = source_input_dir, abcd_dir = abcd_dir,
                      diagnostics_dir = diagnostics_dir, delta_probe_dir = delta_probe_dir,
                      delta_debug_dir = delta_debug_dir, local_calibration_dir = local_calibration_dir,
                      core_fix_dir = core_fix_dir, simulation_ids = simulation_ids,
                      minobs = minobs, input_policy = input_policy)
    write_j_report(dirs, args_info, ctx, j0, j1, j2, j3, j4, j5, recs)
  }
  message("Wrote farfield local identifiability repair under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) main_j()
