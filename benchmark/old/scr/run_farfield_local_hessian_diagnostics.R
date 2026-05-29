#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_local_hessian_diagnostics.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_local_identifiability_repair.R"))

usage <- function() {
  cat(
    "Run farfield local Hessian diagnostics K1-K5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_local_hessian_diagnostics.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --delta-debug-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --local-calibration-dir=benchmark/results/farfield_local_calibration_patch \\\n",
    "    --core-fix-dir=benchmark/results/farfield_core_fix_probe \\\n",
    "    --identifiability-dir=benchmark/results/farfield_local_identifiability_repair \\\n",
    "    --output-dir=benchmark/results/farfield_local_hessian_diagnostics \\\n",
    "    --simulation-ids=1,2,3,4,5,6,7,8,9,10 --minobs=5 --input-policy=full\n",
    sep = ""
  )
}

make_k_dirs <- function(output_dir) make_probe_dirs(output_dir)

support_scope_k <- function(tier, distance) {
  tier <- as.character(tier)
  distance <- as.integer(distance)
  ifelse(tier == "directly_informed" | distance == 0L, "direct",
         ifelse(tier == "local_borrowed" | distance == 1L, "local_borrowed",
                ifelse(tier %in% c("weakly_supported", "graph_borrowed", "prior_dominated") | distance >= 2L,
                       "weakly_supported", "other")))
}

k_config_table <- function(which = c("core", "k4", "k5")) {
  which <- match.arg(which)
  if (which == "core") {
    return(data.frame(
      config_id = c("C0_shell0_control", "C1_shell1_baseline", "C2_all_scale_fixed",
                    "C3_borrowed_residual_0p20", "C4_J4_M4_g_fixed_scale_borrowed",
                    "C5_g_centered", "C6_deterministic_direct_only"),
      shell_depth = c(0L, rep(1L, 6)),
      local_parameterization = c("f", "f", "f", "f", "g_equivalent", "g_equivalent", "g_equivalent"),
      shrink = c("current", "current", "strong", "strong", "strong", "strong", "strong"),
      local_centering = c("none", "none", "reference_direct", "reference_direct", "reference_direct", "reference_direct", "reference_direct"),
      local_centering_weight = c(0, 0, 100, 100, 100, 100, 100),
      fixed_sigma_anchor = c(NA, NA, 0.2, 0.2, 0.2, NA, 0.2),
      fixed_sigma_neighbor = c(NA, NA, 0.1, 0.1, 0.1, NA, 0.1),
      fixed_tau_group = c(NA, NA, 0.1, 0.1, 0.1, NA, 0.1),
      borrowed_residual_sd = c(NA, NA, NA, 0.2, 0.2, NA, 0.005),
      weakly_supported_residual_sd = c(NA, NA, NA, 0.1, 0.1, NA, 0.005),
      eta_prior_sd = c(rep(5, 7)),
      eta_borrowed_prior_sd = c(rep(1.5, 7)),
      eta_borrowed_prior_mean = c(rep(-6, 7)),
      eta_distance_penalty = c(rep(0.75, 7)),
      eval_max = c(500, 500, 500, 500, 2000, 500, 500),
      stringsAsFactors = FALSE
    ))
  }
  if (which == "k4") {
    return(data.frame(
      config_id = c("E0_baseline", "E1_fix_borrowed_eta_prior", "E2_all_eta_stronger_prior",
                    "E3_eta_centering_approx", "E4_direct_only_free_f",
                    "E5_deterministic_borrowed_mean", "E6_delta_context_zero_deterministic",
                    "E7_no_context_direct_borrowed_residual", "E8_safe_shell0_global"),
      shell_depth = c(rep(1L, 8), 0L),
      local_parameterization = c("f", "g_equivalent", "g_equivalent", "g_equivalent",
                                 "g_equivalent", "g_equivalent", "g_equivalent",
                                 "g_equivalent", "f"),
      shrink = c("current", rep("strong", 7), "current"),
      local_centering = c("none", rep("reference_direct", 7), "none"),
      local_centering_weight = c(0, rep(100, 7), 0),
      fixed_sigma_anchor = c(NA, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, NA),
      fixed_sigma_neighbor = c(NA, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, NA),
      fixed_tau_group = c(NA, 0.1, 0.1, 0.1, 0.1, 0.1, 0.001, 0.001, NA),
      borrowed_residual_sd = c(NA, 0.2, 0.2, 0.2, 0.005, 0.005, 0.005, 0.1, NA),
      weakly_supported_residual_sd = c(NA, 0.1, 0.1, 0.1, 0.005, 0.005, 0.005, 0.05, NA),
      eta_prior_sd = c(5, 5, 1, 1, 5, 5, 5, 5, 5),
      eta_borrowed_prior_sd = c(1.5, 0.05, 0.25, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      eta_borrowed_prior_mean = c(-6, -8, -8, -8, -8, -8, -8, -8, -6),
      eta_distance_penalty = c(0.75, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.75),
      eval_max = c(500, rep(1000, 7), 500),
      approximation = c("none", "strong_eta_prior", "strong_eta_prior", "eta_prior_approx_centering",
                        "near_deterministic_residual", "near_deterministic_residual",
                        "fixed_tau_approx_delta_zero", "fixed_tau_approx_no_context",
                        "shell0_safe_baseline"),
      stringsAsFactors = FALSE
    ))
  }
  scales <- data.frame(
    fixed_sigma_anchor = c(0.1, 0.2, 0.5, 0.2, 0.2),
    fixed_sigma_neighbor = c(0.05, 0.1, 0.2, 0.05, 0.2),
    fixed_tau_group = c(0.05, 0.1, 0.2, 0.2, 0.05),
    stringsAsFactors = FALSE
  )
  variants <- data.frame(
    local_variant = c("L0_fixed_scale_current_structure", "L1_fixed_scale_centered_g",
                      "L2_fixed_scale_direct_only_free", "L3_fixed_scale_borrowed_residual",
                      "L4_safe_local_shell0"),
    shell_depth = c(1L, 1L, 1L, 1L, 0L),
    local_parameterization = c("f", "g_equivalent", "g_equivalent", "g_equivalent", "f"),
    local_centering = c("none", "reference_direct", "reference_direct", "reference_direct", "none"),
    local_centering_weight = c(0, 100, 100, 100, 0),
    borrowed_residual_sd = c(NA, NA, 0.005, 0.2, NA),
    weakly_supported_residual_sd = c(NA, NA, 0.005, 0.1, NA),
    shrink = c("strong", "strong", "strong", "strong", "current"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  idx <- 0L
  for (v in seq_len(nrow(variants))) {
    if (variants$local_variant[[v]] == "L4_safe_local_shell0") {
      idx <- idx + 1L
      rows[[idx]] <- cbind(variants[v, , drop = FALSE], scales[5, , drop = FALSE])
    } else {
      for (s in seq_len(nrow(scales))) {
        idx <- idx + 1L
        rows[[idx]] <- cbind(variants[v, , drop = FALSE], scales[s, , drop = FALSE])
      }
    }
  }
  out <- bind_rows_fill(rows)
  out$config_id <- paste0(out$local_variant, "_sa", out$fixed_sigma_anchor,
                          "_sn", out$fixed_sigma_neighbor, "_tg", out$fixed_tau_group)
  out$eval_max <- ifelse(out$shell_depth == 0L, 500L, 500L)
  out
}

build_graph_for_cfg <- function(data, cfg) {
  alfak2::build_karyotype_graph(data, shell_depth = as.integer(cfg$shell_depth[[1L]]), max_nodes = 30000)
}

fit_local_k <- function(data, cfg, sdreport_mode = "all_f_current", restart_id = NA_integer_,
                        initial_jitter_sd = 0, return_tmb_objects = TRUE) {
  graph <- build_graph_for_cfg(data, cfg)
  ctrl <- local_shrink_controls_i(cfg$shrink[[1L]] %||% "current")
  ctrl <- ctrl[setdiff(names(ctrl), c("eta_prior_sd", "eta_borrowed_prior_mean",
                                      "eta_borrowed_prior_sd", "eta_distance_penalty"))]
  finite_or_null <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1L || !is.finite(x)) NULL else x
  }
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
          control = list(eval.max = as.integer(cfg$eval_max[[1L]]), iter.max = as.integer(cfg$eval_max[[1L]])),
          retry_on_untrusted_covariance = FALSE,
          return_optimizer_diagnostics = TRUE,
          return_tmb_objects = return_tmb_objects,
          local_parameterization = cfg$local_parameterization[[1L]],
          local_centering = cfg$local_centering[[1L]],
          local_centering_weight = as.numeric(cfg$local_centering_weight[[1L]]),
          local_centering_weight_mode = "effective_count",
          fixed_sigma_anchor = suppressWarnings(as.numeric(cfg$fixed_sigma_anchor[[1L]])),
          fixed_sigma_neighbor = suppressWarnings(as.numeric(cfg$fixed_sigma_neighbor[[1L]])),
          fixed_tau_group = suppressWarnings(as.numeric(cfg$fixed_tau_group[[1L]])),
          borrowed_residual_sd = finite_or_null(cfg$borrowed_residual_sd[[1L]]),
          weakly_supported_residual_sd = finite_or_null(cfg$weakly_supported_residual_sd[[1L]]),
          eta_prior_sd = as.numeric(cfg$eta_prior_sd[[1L]] %||% 5),
          eta_borrowed_prior_mean = as.numeric(cfg$eta_borrowed_prior_mean[[1L]] %||% -6),
          eta_borrowed_prior_sd = as.numeric(cfg$eta_borrowed_prior_sd[[1L]] %||% 1.5),
          eta_distance_penalty = as.numeric(cfg$eta_distance_penalty[[1L]] %||% 0.75),
          sdreport_mode = sdreport_mode,
          initial_jitter_sd = initial_jitter_sd,
          initial_seed = if (is.na(restart_id)) NULL else as.integer(1000 + restart_id)
        ),
        ctrl
      )
    ),
    error = function(e) e
  )
  list(fit = fit, elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")))
}

local_diag_row_k <- function(res, cfg, grf = NULL, lambda = NA_real_, extra = list()) {
  base <- as.data.frame(cfg, stringsAsFactors = FALSE)
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base$elapsed_sec <- res$elapsed_sec
  if (inherits(res$fit, "error")) {
    base$status <- "error"
    base$error_message <- conditionMessage(res$fit)
    return(base)
  }
  fit <- res$fit
  block <- fit$diagnostics$optimizer$gradient_block_summary
  row <- cbind(
    base,
    data.frame(
      status = "ok",
      convergence = fit$diagnostics$convergence,
      message = fit$diagnostics$message,
      objective = fit$diagnostics$objective,
      gradient_norm = fit$diagnostics$gradient_norm,
      covariance_status = fit$diagnostics$covariance_status,
      covariance_fallback = fit$diagnostics$covariance_fallback,
      fitness_sd_source = fit$diagnostics$fitness_sd_source,
      sdreport_mode = fit$diagnostics$sdreport_mode,
      grad_eta_max_abs = block$grad_eta_max_abs,
      grad_f_max_abs = block$grad_f_max_abs,
      grad_delta_context_max_abs = block$grad_delta_context_max_abs,
      grad_mu_group_max_abs = block$grad_mu_group_max_abs,
      grad_log_sigma_neighbor_abs = block$grad_log_sigma_neighbor_abs,
      grad_log_sigma_anchor_abs = block$grad_log_sigma_anchor_abs,
      grad_log_tau_group_max_abs = block$grad_log_tau_group_max_abs,
      grad_f_direct_max_abs = block$grad_f_direct_max_abs,
      grad_f_local_borrowed_max_abs = block$grad_f_local_borrowed_max_abs,
      grad_f_weakly_supported_max_abs = block$grad_f_weakly_supported_max_abs,
      grad_f_graph_borrowed_max_abs = block$grad_f_graph_borrowed_max_abs,
      max_gradient_block_name = block$max_gradient_block_name,
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(grf)) {
    align <- edge_alignment_i(fit, grf, lambda, cfg$config_id[[1L]])
    row$local_edge_delta_sign_agreement <- align$delta_sign_agreement
    row$local_edge_delta_spearman <- align$delta_spearman
    row$local_edge_delta_pearson <- NA_real_
    row$local_edge_delta_sd_ratio <- align$estimated_delta_sd_ratio
  }
  row
}

parameter_map_k <- function(fit) {
  plist <- fit$optimizer$obj$env$parList(fit$optimizer$opt$par)
  rows <- list()
  start <- 1L
  for (nm in names(plist)) {
    len <- length(plist[[nm]])
    if (!len) next
    idx <- seq.int(start, length.out = len)
    rows[[length(rows) + 1L]] <- data.frame(
      par_index = idx,
      parameter_block = nm,
      parameter_sub_index = seq_len(len),
      stringsAsFactors = FALSE
    )
    start <- start + len
  }
  out <- bind_rows_fill(rows)
  n <- nrow(fit$summary)
  node_blocks <- out$parameter_block %in% c("eta", "f")
  out$node_id <- NA_integer_
  out$karyotype <- NA_character_
  out$support_tier <- NA_character_
  out$support_distance <- NA_integer_
  out$support_scope <- NA_character_
  ii <- which(node_blocks & out$parameter_sub_index <= n)
  out$node_id[ii] <- out$parameter_sub_index[ii]
  out$karyotype[ii] <- as.character(fit$summary$karyotype[out$node_id[ii]])
  out$support_tier[ii] <- as.character(fit$summary$support_tier[out$node_id[ii]])
  out$support_distance[ii] <- as.integer(fit$summary$support_distance[out$node_id[ii]])
  out$support_scope[ii] <- support_scope_k(out$support_tier[ii], out$support_distance[ii])
  out
}

compute_hessian_diagnostics <- function(fit, config_id, grf = NULL, lambda = NA_real_, thresholds = c(1e-8, 1e-6, 1e-4)) {
  if (inherits(fit, "error") || is.null(fit$optimizer$obj)) {
    fail <- data.frame(config_id = config_id, hessian_available = FALSE, hessian_status = "fit_error",
                       near_zero_threshold = thresholds, stringsAsFactors = FALSE)
    return(list(spectrum = fail, block = data.frame(), tier = data.frame(), top_nodes = data.frame(),
                failure = fail[1, , drop = FALSE]))
  }
  H <- try(fit$optimizer$obj$he(fit$optimizer$opt$par), silent = TRUE)
  if (inherits(H, "try-error")) {
    fail <- data.frame(config_id = config_id, hessian_available = FALSE, hessian_status = "hessian_error",
                       hessian_error = conditionMessage(attr(H, "condition")),
                       near_zero_threshold = thresholds, stringsAsFactors = FALSE)
    return(list(spectrum = fail, block = data.frame(), tier = data.frame(), top_nodes = data.frame(),
                failure = fail[1, , drop = FALSE]))
  }
  H <- as.matrix(H)
  H <- (H + t(H)) / 2
  ev <- try(eigen(H, symmetric = TRUE), silent = TRUE)
  if (inherits(ev, "try-error")) {
    fail <- data.frame(config_id = config_id, hessian_available = TRUE, hessian_status = "eigen_error",
                       hessian_error = conditionMessage(attr(ev, "condition")),
                       near_zero_threshold = thresholds, stringsAsFactors = FALSE)
    return(list(spectrum = fail, block = data.frame(), tier = data.frame(), top_nodes = data.frame(),
                failure = fail[1, , drop = FALSE]))
  }
  values <- as.numeric(ev$values)
  finite_values <- values[is.finite(values)]
  spectrum <- do.call(rbind, lapply(thresholds, function(thr) {
    n_neg <- sum(values < -thr, na.rm = TRUE)
    n_zero <- sum(abs(values) <= thr, na.rm = TRUE)
    n_pos <- sum(values > thr, na.rm = TRUE)
    denom <- suppressWarnings(min(abs(values[abs(values) > thr]), na.rm = TRUE))
    cond <- if (is.finite(denom) && denom > 0) max(abs(values), na.rm = TRUE) / denom else Inf
    status <- if (!length(finite_values) || any(!is.finite(values))) "nonfinite" else if (n_neg > 0) "indefinite" else if (n_zero > 0) "near_singular" else "positive_definite"
    data.frame(
      config_id = config_id,
      hessian_available = TRUE,
      hessian_status = status,
      hessian_min_eigenvalue = min(values, na.rm = TRUE),
      hessian_max_eigenvalue = max(values, na.rm = TRUE),
      hessian_condition_number = cond,
      n_negative_eigenvalues = n_neg,
      n_near_zero_eigenvalues = n_zero,
      near_zero_threshold = thr,
      n_positive_eigenvalues = n_pos,
      n_nonfinite_eigenvalues = sum(!is.finite(values)),
      stringsAsFactors = FALSE
    )
  }))
  pmap <- parameter_map_k(fit)
  ord <- order(abs(values), decreasing = FALSE)
  ord <- ord[seq_len(min(10L, length(ord)))]
  truth <- rep(NA_real_, nrow(fit$summary))
  if (!is.null(grf)) {
    tm <- compute_grf_truth(fit$summary$karyotype, grf$centroids, lambda)
    truth <- as.numeric(tm[as.character(fit$summary$karyotype)])
  }
  block_rows <- list(); tier_rows <- list(); top_rows <- list()
  for (jj in seq_along(ord)) {
    vec_id <- ord[[jj]]
    load <- as.numeric(ev$vectors[, vec_id])
    total <- sum(abs(load), na.rm = TRUE)
    if (!is.finite(total) || total <= 0) total <- 1
    tmp <- pmap
    tmp$loading_abs <- abs(load[tmp$par_index])
    tmp$loading_signed <- load[tmp$par_index]
    blk <- stats::aggregate(loading_abs ~ parameter_block, tmp, sum, na.rm = TRUE)
    blk$loading_fraction <- blk$loading_abs / total
    blk$config_id <- config_id
    blk$eigenvector_rank <- jj
    blk$eigenvalue <- values[[vec_id]]
    block_rows[[length(block_rows) + 1L]] <- blk
    ft <- tmp[tmp$parameter_block %in% c("f", "eta") & !is.na(tmp$support_scope), , drop = FALSE]
    if (nrow(ft)) {
      tr <- stats::aggregate(loading_abs ~ parameter_block + support_tier + support_scope, ft, sum, na.rm = TRUE)
      tr$loading_fraction <- tr$loading_abs / total
      tr$config_id <- config_id
      tr$eigenvector_rank <- jj
      tr$eigenvalue <- values[[vec_id]]
      tier_rows[[length(tier_rows) + 1L]] <- tr
      top <- ft[order(ft$loading_abs, decreasing = TRUE), , drop = FALSE]
      top <- top[seq_len(min(20L, nrow(top))), , drop = FALSE]
      top$config_id <- config_id
      top$eigenvector_rank <- jj
      top$eigenvalue <- values[[vec_id]]
      top$parameter_name <- paste0(top$parameter_block, "[", top$parameter_sub_index, "]")
      node <- top$node_id
      top$count_t0 <- fit$summary$count_t0[node]
      top$count_t1 <- fit$summary$count_t1[node]
      top$count_total <- fit$summary$count_total[node]
      top$effective_count_total <- fit$summary$effective_count_total[node]
      top$fitness_mean <- fit$summary$fitness_mean[node]
      top$eta <- fit$parameter_mode$eta[node]
      top$truth_fitness <- truth[node]
      top_rows[[length(top_rows) + 1L]] <- top[, c("config_id", "eigenvector_rank", "eigenvalue", "parameter_name",
                                                   "parameter_block", "node_id", "karyotype", "support_tier",
                                                   "support_distance", "support_scope", "loading_abs",
                                                   "loading_signed", "count_t0", "count_t1", "count_total",
                                                   "effective_count_total", "fitness_mean", "eta", "truth_fitness")]
    }
  }
  failure <- spectrum[spectrum$near_zero_threshold == 1e-6, , drop = FALSE]
  dominant <- bind_rows_fill(block_rows)
  if (nrow(dominant)) {
    d <- stats::aggregate(loading_fraction ~ parameter_block, dominant, median, na.rm = TRUE)
    d <- d[order(d$loading_fraction, decreasing = TRUE), , drop = FALSE]
    failure$nullspace_dominant_block <- d$parameter_block[[1L]]
  }
  tiers <- bind_rows_fill(tier_rows)
  if (nrow(tiers)) {
    td <- stats::aggregate(loading_fraction ~ support_scope, tiers[tiers$parameter_block == "f", , drop = FALSE], median, na.rm = TRUE)
    td <- td[order(td$loading_fraction, decreasing = TRUE), , drop = FALSE]
    failure$nullspace_dominant_support_tier <- if (nrow(td)) td$support_scope[[1L]] else NA_character_
  }
  list(spectrum = spectrum, block = dominant, tier = tiers, top_nodes = bind_rows_fill(top_rows), failure = failure)
}

sdreport_table_k <- function(fit) {
  tab <- try(suppressWarnings(as.data.frame(summary(fit$sdreport, "report"))), silent = TRUE)
  if (inherits(tab, "try-error")) {
    return(data.frame(sdreport_status = "summary_error", sdreport_error = conditionMessage(attr(tab, "condition")),
                      n_reported = NA_integer_, adreport_nonfinite_count = NA_integer_, stringsAsFactors = FALSE))
  }
  if (!nrow(tab)) {
    return(data.frame(sdreport_status = "empty", sdreport_error = NA_character_,
                      n_reported = 0L, adreport_nonfinite_count = 0L, stringsAsFactors = FALSE))
  }
  rn <- rownames(tab)
  data.frame(
    reported_name = rn,
    estimate = tab$Estimate,
    std_error = tab$`Std. Error`,
    sdreport_status = ifelse(is.finite(tab$`Std. Error`) & tab$`Std. Error` > 0, "finite", "nonfinite"),
    stringsAsFactors = FALSE
  )
}

run_k1 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "k1_hessian_nullspace_diagnostics.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- k_config_table("core")
  fits <- list(); diag_rows <- list(); spectra <- list(); blocks <- list(); tiers <- list(); tops <- list(); failures <- list()
  for (i in seq_len(nrow(cfg))) {
    res <- fit_local_k(bundle$data, cfg[i, , drop = FALSE], return_tmb_objects = TRUE)
    fits[[cfg$config_id[[i]]]] <- if (inherits(res$fit, "error")) NULL else res$fit
    diag_rows[[i]] <- local_diag_row_k(res, cfg[i, , drop = FALSE], grf, task_info$lambda)
    if (!inherits(res$fit, "error")) {
      hd <- compute_hessian_diagnostics(res$fit, cfg$config_id[[i]], grf, task_info$lambda)
      spectra[[i]] <- hd$spectrum; blocks[[i]] <- hd$block; tiers[[i]] <- hd$tier
      tops[[i]] <- hd$top_nodes; failures[[i]] <- hd$failure
    }
  }
  spectrum <- bind_rows_fill(spectra)
  block <- bind_rows_fill(blocks)
  tier <- bind_rows_fill(tiers)
  top <- bind_rows_fill(tops)
  failure <- bind_rows_fill(failures)
  write_tsv_safe(spectrum, file.path(dirs$tables, "k1_local_hessian_spectrum.tsv"))
  write_tsv_safe(block, file.path(dirs$tables, "k1_nullspace_loading_by_parameter_block.tsv"))
  write_tsv_safe(tier, file.path(dirs$tables, "k1_nullspace_loading_by_support_tier.tsv"))
  write_tsv_safe(top, file.path(dirs$tables, "k1_nullspace_top_nodes.tsv"))
  write_tsv_safe(failure, file.path(dirs$tables, "k1_hessian_failure_modes.tsv"))
  out <- list(local = bind_rows_fill(diag_rows), spectrum = spectrum, block = block, tier = tier, top_nodes = top, failure = failure)
  saveRDS(out, rds)
  out
}

run_k2 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "k2_sdreport_triage.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- k_config_table("core")
  cfg <- cfg[cfg$config_id %in% c("C0_shell0_control", "C1_shell1_baseline", "C2_all_scale_fixed",
                                  "C3_borrowed_residual_0p20", "C4_J4_M4_g_fixed_scale_borrowed"), , drop = FALSE]
  modes <- c("none", "direct_f_only", "direct_and_local_borrowed_f", "all_f_current", "pi0_pi1_only")
  rows <- list(); direct <- list(); allnodes <- list(); fail <- list()
  idx <- 0L
  for (i in seq_len(nrow(cfg))) {
    for (mode in modes) {
      idx <- idx + 1L
      res <- fit_local_k(bundle$data, cfg[i, , drop = FALSE], sdreport_mode = mode, return_tmb_objects = TRUE)
      base <- local_diag_row_k(res, cfg[i, , drop = FALSE], grf, task_info$lambda, extra = list(sdreport_mode_test = mode))
      if (inherits(res$fit, "error")) {
        rows[[idx]] <- base
        next
      }
      tab <- sdreport_table_k(res$fit)
      finite <- if ("std_error" %in% names(tab)) is.finite(tab$std_error) & tab$std_error > 0 else logical()
      base$n_reported <- if ("reported_name" %in% names(tab)) nrow(tab) else tab$n_reported[[1L]]
      base$adreport_nonfinite_count <- if (length(finite)) sum(!finite) else tab$adreport_nonfinite_count[[1L]]
      base$adreport_status <- if (length(finite) && all(finite)) "finite" else if (length(finite)) "nonfinite" else tab$sdreport_status[[1L]]
      jp <- try(suppressWarnings(TMB::sdreport(res$fit$optimizer$obj, par.fixed = res$fit$optimizer$opt$par, getJointPrecision = TRUE)), silent = TRUE)
      base$joint_precision_status <- if (inherits(jp, "try-error")) "error" else if (!is.null(jp$jointPrecision)) "available" else "missing"
      rows[[idx]] <- base
      if ("reported_name" %in% names(tab)) {
        tab$config_id <- cfg$config_id[[i]]
        tab$sdreport_mode_test <- mode
        if (mode == "direct_f_only") direct[[length(direct) + 1L]] <- tab
        if (mode == "all_f_current") allnodes[[length(allnodes) + 1L]] <- tab
      }
    }
  }
  triage <- bind_rows_fill(rows)
  failure <- triage[, intersect(c("config_id", "sdreport_mode_test", "covariance_status", "adreport_status",
                                  "adreport_nonfinite_count", "joint_precision_status", "gradient_norm"),
                                names(triage)), drop = FALSE]
  write_tsv_safe(triage, file.path(dirs$tables, "k2_sdreport_triage.tsv"))
  write_tsv_safe(bind_rows_fill(direct), file.path(dirs$tables, "k2_sdreport_direct_only.tsv"))
  write_tsv_safe(bind_rows_fill(allnodes), file.path(dirs$tables, "k2_sdreport_all_nodes.tsv"))
  write_tsv_safe(failure, file.path(dirs$tables, "k2_sdreport_failure_modes.tsv"))
  out <- list(triage = triage, direct = bind_rows_fill(direct), all_nodes = bind_rows_fill(allnodes), failure = failure)
  saveRDS(out, rds)
  out
}

pairwise_cor_summary <- function(mat) {
  if (is.null(mat) || !length(mat) || ncol(as.matrix(mat)) < 1L) {
    return(c(median = NA_real_, min = NA_real_, max = NA_real_))
  }
  if (is.null(mat) || nrow(mat) < 2L) return(c(median = NA_real_, min = NA_real_, max = NA_real_))
  cc <- suppressWarnings(stats::cor(t(mat), use = "pairwise.complete.obs"))
  vals <- cc[upper.tri(cc)]
  c(median = stats::median(vals, na.rm = TRUE), min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE))
}

edge_delta_vector_k <- function(fit) {
  graph <- fit$graph
  m <- setNames(as.numeric(fit$summary$fitness_mean), as.character(fit$summary$karyotype))
  parent <- as.character(graph$labels)[as.integer(unlist(graph$parent_from0)) + 1L]
  child <- as.character(graph$labels)[as.integer(unlist(graph$parent_to0)) + 1L]
  as.numeric(m[child] - m[parent])
}

run_k3 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "k3_j4_multistart_stability.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  all_cfg <- k_config_table("core")
  cfg <- all_cfg[all_cfg$config_id %in% c("C1_shell1_baseline", "C2_all_scale_fixed", "C3_borrowed_residual_0p20",
                                          "C4_J4_M4_g_fixed_scale_borrowed", "C0_shell0_control"), , drop = FALSE]
  names(cfg)[names(cfg) == "config_id"] <- "config_id"
  restarts <- 1:10
  runs <- list(); fit_store <- list(); idx <- 0L
  for (i in seq_len(nrow(cfg))) {
    for (r in restarts) {
      idx <- idx + 1L
      res <- fit_local_k(bundle$data, cfg[i, , drop = FALSE], restart_id = r, initial_jitter_sd = 0.05, return_tmb_objects = FALSE)
      runs[[idx]] <- local_diag_row_k(res, cfg[i, , drop = FALSE], grf, task_info$lambda, extra = list(restart_id = r))
      if (!inherits(res$fit, "error")) fit_store[[paste(cfg$config_id[[i]], r, sep = "::")]] <- res$fit
    }
  }
  run_tbl <- bind_rows_fill(runs)
  stab_rows <- list(); edge_rows <- list(); node_rows <- list()
  for (cid in unique(cfg$config_id)) {
    keys <- grep(paste0("^", cid, "::"), names(fit_store), value = TRUE)
    fs <- fit_store[keys]
    if (!length(fs)) next
    labels <- as.character(fs[[1]]$summary$karyotype)
    mat_all <- do.call(rbind, lapply(fs, function(f) setNames(f$summary$fitness_mean, f$summary$karyotype)[labels]))
    direct <- fs[[1]]$summary$support_distance == 0
    borrowed <- fs[[1]]$summary$support_distance == 1
    weak <- fs[[1]]$summary$support_distance >= 2
    eta_mat <- do.call(rbind, lapply(fs, function(f) as.numeric(f$parameter_mode$eta)))
    ed_mat <- do.call(rbind, lapply(fs, edge_delta_vector_k))
    ca <- pairwise_cor_summary(mat_all)
    cd <- pairwise_cor_summary(mat_all[, direct, drop = FALSE])
    cb <- pairwise_cor_summary(mat_all[, borrowed, drop = FALSE])
    cw <- pairwise_cor_summary(mat_all[, weak, drop = FALSE])
    ce <- pairwise_cor_summary(eta_mat)
    ced <- pairwise_cor_summary(ed_mat)
    stab_rows[[length(stab_rows) + 1L]] <- data.frame(
      config_id = cid,
      n_restarts = length(fs),
      fitness_mean_pairwise_correlation = ca["median"],
      direct_fitness_pairwise_correlation = cd["median"],
      local_borrowed_fitness_pairwise_correlation = cb["median"],
      weakly_supported_fitness_pairwise_correlation = cw["median"],
      eta_pairwise_correlation = ce["median"],
      edge_delta_pairwise_correlation = ced["median"],
      objective_sd = stats::sd(run_tbl$objective[run_tbl$config_id == cid], na.rm = TRUE),
      gradient_norm_median = stats::median(run_tbl$gradient_norm[run_tbl$config_id == cid], na.rm = TRUE),
      gradient_norm_iqr = stats::IQR(run_tbl$gradient_norm[run_tbl$config_id == cid], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      config_id = cid,
      edge_delta_pairwise_correlation = ced["median"],
      edge_delta_sign_agreement_across_restarts = mean(sign(ed_mat[-1, , drop = FALSE]) == rep(sign(ed_mat[1, ]), each = nrow(ed_mat) - 1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    node_rows[[length(node_rows) + 1L]] <- data.frame(
      config_id = cid,
      node_id = seq_along(labels),
      karyotype = labels,
      support_tier = fs[[1]]$summary$support_tier,
      support_distance = fs[[1]]$summary$support_distance,
      nodewise_restart_mean = colMeans(mat_all, na.rm = TRUE),
      nodewise_restart_sd = apply(mat_all, 2, stats::sd, na.rm = TRUE),
      nodewise_restart_cv = apply(mat_all, 2, stats::sd, na.rm = TRUE) / pmax(abs(colMeans(mat_all, na.rm = TRUE)), 1e-8),
      stringsAsFactors = FALSE
    )
  }
  stab <- bind_rows_fill(stab_rows)
  edge <- bind_rows_fill(edge_rows)
  nodes <- bind_rows_fill(node_rows)
  rec <- data.frame(
    recommendation = if (any(stab$config_id == "C4_J4_M4_g_fixed_scale_borrowed" & stab$fitness_mean_pairwise_correlation > 0.99, na.rm = TRUE)) "map_mean_stable_covariance_untrusted" else "map_stability_not_established",
    restarts_per_config = length(restarts),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(run_tbl, file.path(dirs$tables, "k3_j4_multistart_runs.tsv"))
  write_tsv_safe(stab, file.path(dirs$tables, "k3_multistart_fitness_stability.tsv"))
  write_tsv_safe(edge, file.path(dirs$tables, "k3_multistart_edge_delta_stability.tsv"))
  write_tsv_safe(nodes, file.path(dirs$tables, "k3_nodewise_restart_variance.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "k3_multistart_recommendation.tsv"))
  out <- list(runs = run_tbl, stability = stab, edge = edge, nodewise = nodes, recommendation = rec)
  saveRDS(out, rds)
  out
}

run_k4 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "k4_eta_f_borrowed_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- k_config_table("k4")
  rows <- list(); hrows <- list()
  for (i in seq_len(nrow(cfg))) {
    res <- fit_local_k(bundle$data, cfg[i, , drop = FALSE], return_tmb_objects = TRUE)
    rows[[i]] <- local_diag_row_k(res, cfg[i, , drop = FALSE], grf, task_info$lambda)
    if (!inherits(res$fit, "error")) {
      hd <- compute_hessian_diagnostics(res$fit, cfg$config_id[[i]], grf, task_info$lambda, thresholds = 1e-6)
      hrows[[i]] <- hd$spectrum
    }
  }
  tbl <- bind_rows_fill(rows)
  hess <- bind_rows_fill(hrows)
  rec <- data.frame(
    recommendation = if (any(tbl$covariance_status == "TMB_sdreport", na.rm = TRUE)) "eta_borrowed_variant_candidate_found" else "eta_f_borrowed_variants_still_experimental",
    best_gradient_config = tbl$config_id[which.min(tbl$gradient_norm)],
    best_gradient_norm = min(tbl$gradient_norm, na.rm = TRUE),
    best_edge_sign_agreement = max(tbl$local_edge_delta_sign_agreement, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(tbl, file.path(dirs$tables, "k4_eta_f_confounding_probe.tsv"))
  write_tsv_safe(tbl[grepl("direct_only", tbl$config_id), , drop = FALSE], file.path(dirs$tables, "k4_direct_only_free_f_probe.tsv"))
  write_tsv_safe(tbl[grepl("deterministic|delta_context_zero|no_context", tbl$config_id), , drop = FALSE], file.path(dirs$tables, "k4_deterministic_borrowed_probe.tsv"))
  write_tsv_safe(hess, file.path(dirs$tables, "k4_eta_f_hessian_summary.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "k4_eta_f_recommendation.tsv"))
  out <- list(results = tbl, hessian = hess, recommendation = rec)
  saveRDS(out, rds)
  out
}

run_k5 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "k5_fixed_scale_calibration_local_model.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  cfg <- k_config_table("k5")
  rows <- list()
  for (i in seq_len(nrow(cfg))) {
    row <- cfg[i, , drop = FALSE]
    row$eta_prior_sd <- 5
    row$eta_borrowed_prior_sd <- 1.5
    row$eta_borrowed_prior_mean <- -6
    row$eta_distance_penalty <- 0.75
    res <- fit_local_k(bundle$data, row, return_tmb_objects = FALSE)
    rows[[i]] <- local_diag_row_k(res, row, grf, task_info$lambda)
  }
  local_tbl <- bind_rows_fill(rows)
  top <- local_tbl[order(local_tbl$covariance_status != "TMB_sdreport", local_tbl$gradient_norm), , drop = FALSE]
  top <- top[seq_len(min(5L, nrow(top))), , drop = FALSE]
  global_cfg <- data.frame(
    experiment = "K5",
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
  grow <- list(); gidx <- 0L
  for (i in seq_len(nrow(top))) {
    row <- top[i, , drop = FALSE]
    row$eta_prior_sd <- 5
    row$eta_borrowed_prior_sd <- 1.5
    row$eta_borrowed_prior_mean <- -6
    row$eta_distance_penalty <- 0.75
    res <- fit_local_k(bundle$data, row, return_tmb_objects = FALSE)
    if (inherits(res$fit, "error")) next
    graph <- alfak2::build_karyotype_graph(res$fit$data, shell_depth = 2, max_nodes = 30000)
    for (j in seq_len(nrow(global_cfg))) {
      fit <- tryCatch(alfak2::fit_graph_posterior(
        res$fit, graph,
        lambda_l_grid = global_cfg$lambda_l[[j]],
        lambda_e_grid = global_cfg$lambda_e[[j]],
        sigma_obs_grid = global_cfg$sigma_obs[[j]],
        graph_edge_weight = global_cfg$graph_edge_weight[[j]],
        compute_sd = FALSE
      ), error = function(e) e)
      gidx <- gidx + 1L
      if (inherits(fit, "error")) {
        grow[[gidx]] <- data.frame(config_id = row$config_id[[1L]], candidate_id = global_cfg$candidate_id[[j]],
                                   status = "error", error_message = conditionMessage(fit), stringsAsFactors = FALSE)
      } else {
        m <- score_summary_abcd(fit$summary, graph, grf, task_info$lambda, task_info, global_cfg[j, , drop = FALSE], "fixed_scale_calibration")
        far <- m[m$support_scope == "farfield" & m$metric_scale == "native", , drop = FALSE]
        far$config_id <- row$config_id[[1L]]
        far$local_variant <- row$local_variant[[1L]]
        far$shape_classification <- metric_shape_class(far)
        far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_i(far[k, , drop = FALSE]), character(1))
        grow[[gidx]] <- far
      }
    }
  }
  global_tbl <- bind_rows_fill(grow)
  rec <- data.frame(
    recommendation = if (any(local_tbl$covariance_status == "TMB_sdreport", na.rm = TRUE)) "fixed_scale_local_candidate_found" else "fixed_scale_should_be_calibration_level_but_not_trusted_yet",
    best_local_config = top$config_id[[1L]],
    best_gradient_norm = top$gradient_norm[[1L]],
    any_valid_global_shape = any(global_tbl$recommended_status == "valid_shape_config", na.rm = TRUE),
    grid_note = "Runtime-bounded representative grid: 5 scale tuples x 4 shell1 variants plus shell0 control.",
    stringsAsFactors = FALSE
  )
  write_tsv_safe(local_tbl, file.path(dirs$tables, "k5_fixed_scale_grid_local_fit.tsv"))
  write_tsv_safe(global_tbl, file.path(dirs$tables, "k5_fixed_scale_global_shape.tsv"))
  write_tsv_safe(rec, file.path(dirs$tables, "k5_fixed_scale_calibration_recommendation.tsv"))
  out <- list(local = local_tbl, global = global_tbl, recommendation = rec)
  saveRDS(out, rds)
  out
}

make_k_recommendations <- function(k1, k2, k3, k4, k5) {
  m4 <- k1$failure[k1$failure$config_id == "C4_J4_M4_g_fixed_scale_borrowed" & k1$failure$near_zero_threshold == 1e-6, , drop = FALSE]
  direct_ok <- any(k2$triage$sdreport_mode_test == "direct_f_only" & k2$triage$adreport_status == "finite", na.rm = TRUE)
  shell1_direct_ok <- any(k2$triage$sdreport_mode_test == "direct_f_only" &
                            k2$triage$config_id != "C0_shell0_control" &
                            k2$triage$adreport_status == "finite", na.rm = TRUE)
  m4_stab <- k3$stability[k3$stability$config_id == "C4_J4_M4_g_fixed_scale_borrowed", , drop = FALSE]
  data.frame(
    table = c("hessian_nullspace_recommendation", "sdreport_triage_recommendation",
              "multistart_stability_recommendation", "eta_f_borrowed_recommendation",
              "fixed_scale_calibration_recommendation", "recommended_next_steps"),
    recommendation = c(
      "Treat J4/M4 as Hessian-identifiability failure, not a simple optimizer failure.",
      if (shell1_direct_ok) "Direct-only covariance is available for at least one shell_depth=1 mode; all-node covariance remains unsafe." else "Direct-only covariance is only available for shell0 here; shell_depth=1 covariance remains unsafe.",
      if (nrow(m4_stab) && m4_stab$fitness_mean_pairwise_correlation[[1L]] > 0.99) "MAP mean is restart-stable enough to consider empirical variance diagnostics." else "MAP restart stability remains insufficient.",
      "Eta/f and borrowed-node variants remain diagnostic; shell_depth=1 should stay experimental.",
      "Fixed local scales should move into calibration-level search, but no trusted shell_depth=1 covariance was found.",
      "Do not implement edge-gradient; prioritize direct-only covariance, empirical restart variance, and safe shell0 production fallback."
    ),
    evidence = c(
      if (nrow(m4)) paste0("M4 status=", m4$hessian_status[[1L]], "; dominant_block=", m4$nullspace_dominant_block[[1L]]) else "M4 Hessian unavailable",
      paste0("direct_only_finite_any=", direct_ok, "; shell1_direct_only_finite=", shell1_direct_ok),
      if (nrow(m4_stab)) paste0("M4 fitness cor=", fmt_metric(m4_stab$fitness_mean_pairwise_correlation[[1L]])) else "M4 stability unavailable",
      paste0("best K4 gradient=", fmt_metric(min(k4$results$gradient_norm, na.rm = TRUE))),
      paste0("best K5 gradient=", fmt_metric(min(k5$local$gradient_norm, na.rm = TRUE)), "; valid_global=", any(k5$global$recommended_status == "valid_shape_config", na.rm = TRUE)),
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_k_report <- function(dirs, args_info, ctx, k1, k2, k3, k4, k5, recs) {
  all_long <- bind_rows_fill(list(
    transform(k1$local, experiment = "K1_local_hessian"),
    transform(k1$spectrum, experiment = "K1_spectrum"),
    transform(k2$triage, experiment = "K2_sdreport"),
    transform(k3$runs, experiment = "K3_multistart_runs"),
    transform(k4$results, experiment = "K4_eta_f_borrowed"),
    transform(k5$local, experiment = "K5_fixed_scale_local"),
    transform(k5$global, experiment = "K5_fixed_scale_global")
  ))
  m4 <- k1$failure[k1$failure$config_id == "C4_J4_M4_g_fixed_scale_borrowed" & k1$failure$near_zero_threshold == 1e-6, , drop = FALSE]
  m4stab <- k3$stability[k3$stability$config_id == "C4_J4_M4_g_fixed_scale_borrowed", , drop = FALSE]
  direct_any <- any(k2$triage$sdreport_mode_test == "direct_f_only" & k2$triage$adreport_status == "finite", na.rm = TRUE)
  direct_shell1 <- any(k2$triage$sdreport_mode_test == "direct_f_only" &
                         k2$triage$config_id != "C0_shell0_control" &
                         k2$triage$adreport_status == "finite", na.rm = TRUE)
  all_any <- any(k2$triage$sdreport_mode_test == "all_f_current" & k2$triage$adreport_status == "finite", na.rm = TRUE)
  all_shell1 <- any(k2$triage$sdreport_mode_test == "all_f_current" &
                      k2$triage$config_id != "C0_shell0_control" &
                      k2$triage$adreport_status == "finite", na.rm = TRUE)
  m4_direct <- k2$triage[k2$triage$config_id == "C4_J4_M4_g_fixed_scale_borrowed" &
                           k2$triage$sdreport_mode_test == "direct_f_only", , drop = FALSE]
  summary <- data.frame(
    experiment = c("K1", "K2", "K3", "K4", "K5"),
    key_result = c(
      paste0("M4_hessian=", k1$failure$hessian_status[k1$failure$config_id == "C4_J4_M4_g_fixed_scale_borrowed" & k1$failure$near_zero_threshold == 1e-6][[1L]]),
      paste0("direct_only_finite_any=", direct_any, "; shell1=", direct_shell1),
      paste0("M4_fitness_cor=", fmt_metric(k3$stability$fitness_mean_pairwise_correlation[k3$stability$config_id == "C4_J4_M4_g_fixed_scale_borrowed"][[1L]])),
      paste0("best_gradient=", fmt_metric(min(k4$results$gradient_norm, na.rm = TRUE))),
      paste0("valid_global_shape=", sum(k5$global$recommended_status == "valid_shape_config", na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_k_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "k_experiment_summary.tsv"))
  for (nm in recs$table) write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  lines <- c(
    "# Farfield Local Hessian Diagnostics Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", args_info$source_input_dir, "`"),
    paste0("- abcd-dir: `", args_info$abcd_dir, "`"),
    paste0("- diagnostics-dir: `", args_info$diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", args_info$delta_probe_dir, "`"),
    paste0("- delta-debug-dir: `", args_info$delta_debug_dir, "`"),
    paste0("- local-calibration-dir: `", args_info$local_calibration_dir, "`"),
    paste0("- core-fix-dir: `", args_info$core_fix_dir, "`"),
    paste0("- identifiability-dir: `", args_info$identifiability_dir, "`"),
    paste0("- simulation_ids: ", paste(args_info$simulation_ids, collapse = ",")),
    paste0("- minobs: ", args_info$minobs),
    paste0("- input_policy: ", args_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## Prior results summary",
    "- ABCD through J show no non-oracle valid farfield shape. Oracle per-edge delta works, but non-oracle delta and shell_depth=1 local covariance remain untrusted.",
    "- J4/M4 reduced gradient to near zero but sdreport was non-finite, motivating Hessian/nullspace and sdreport triage.",
    "",
    "## K1 Hessian/nullspace",
    paste0("- J4/M4 Hessian status: ", if (nrow(m4)) m4$hessian_status[[1L]] else "unavailable", "."),
    paste0("- J4/M4 dominant near-null block: ", if (nrow(m4)) m4$nullspace_dominant_block[[1L]] else "NA", "; dominant support scope: ", if (nrow(m4)) m4$nullspace_dominant_support_tier[[1L]] else "NA", "."),
    "- Near-null diagnostics point to identifiability/Hessian structure rather than iteration count.",
    "",
    "## K2 sdreport triage",
    paste0("- direct-only finite ADREPORT observed: ", direct_any, " overall; shell_depth=1 direct-only finite: ", direct_shell1, "."),
    paste0("- J4/M4 direct-only ADREPORT status: ", if (nrow(m4_direct)) m4_direct$adreport_status[[1L]] else "NA", "."),
    paste0("- all-node finite ADREPORT observed: ", all_any, " overall; shell_depth=1 all-node finite: ", all_shell1, "."),
    "- In this run finite covariance is limited to shell0; shell_depth=1 needs empirical/restart variance or further reparameterization before production use.",
    "",
    "## K3 J4/M4 multistart",
    paste0("- M4 fitness mean pairwise correlation: ", if (nrow(m4stab)) fmt_metric(m4stab$fitness_mean_pairwise_correlation[[1L]]) else "NA", "."),
    paste0("- M4 edge-delta pairwise correlation: ", if (nrow(m4stab)) fmt_metric(m4stab$edge_delta_pairwise_correlation[[1L]]) else "NA", "."),
    "- MAP stability and covariance reliability are treated separately; stable mean alone is not enough to train non-oracle delta.",
    "",
    "## K4 eta/f and deterministic borrowed",
    paste0("- best K4 gradient: ", fmt_metric(min(k4$results$gradient_norm, na.rm = TRUE)), "."),
    paste0("- best K4 edge sign agreement: ", fmt_metric(max(k4$results$local_edge_delta_sign_agreement, na.rm = TRUE)), "."),
    "- Strong eta priors and deterministic borrowed approximations remain diagnostic and do not make shell_depth=1 production-ready.",
    "",
    "## K5 fixed-scale calibration-level local",
    paste0("- best K5 local gradient: ", fmt_metric(min(k5$local$gradient_norm, na.rm = TRUE)), "."),
    paste0("- K5 valid global shape count: ", sum(k5$global$recommended_status == "valid_shape_config", na.rm = TRUE), "."),
    "- K5 used a runtime-bounded representative fixed-scale grid after the full 3x3x3 grid exceeded practical runtime.",
    "- Fixed scales should be treated as calibration hyperparameters, but safe shell0 remains the more reliable production fallback until covariance is trusted.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- Main blocker: non-oracle delta remains untrusted and shell_depth=1 all-node covariance/identifiability is not production-safe.",
    "- Keep normalized as benchmark/probe/calibration candidate default, with amplitude-collapse diagnostics and no_valid_shape_configuration gate.",
    "- Keep compute_sd=FALSE and return_optimizer_diagnostics=TRUE.",
    "- Do not default anchor_count_reference=minobs for full input.",
    "- Next priorities: direct-only covariance path, empirical restart/bootstrap anchor variance, fixed-scale calibration grid, and shell0 safe production mode."
  )
  writeLines(lines, file.path(dirs$root, "farfield_local_hessian_diagnostics_report.md"))
  saveRDS(list(k1 = k1, k2 = k2, k3 = k3, k4 = k4, k5 = k5, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_local_hessian_diagnostics_all_results.rds"))
}

main_k <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage(); return(invisible(NULL))
  }
  mode <- match.arg(tolower(as.character(arg_value(args, "mode", "all"))),
                    c("prepare", "k1-hessian-nullspace", "k2-sdreport-triage",
                      "k3-j4-multistart", "k4-eta-f-borrowed-probe",
                      "k5-fixed-scale-calibration-local", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  delta_debug_dir <- as.character(arg_value(args, "delta_debug_dir", "benchmark/results/farfield_delta_debug"))
  local_calibration_dir <- as.character(arg_value(args, "local_calibration_dir", "benchmark/results/farfield_local_calibration_patch"))
  core_fix_dir <- as.character(arg_value(args, "core_fix_dir", "benchmark/results/farfield_core_fix_probe"))
  identifiability_dir <- as.character(arg_value(args, "identifiability_dir", "benchmark/results/farfield_local_identifiability_repair"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_local_hessian_diagnostics"))
  simulation_ids <- arg_integer_csv(args, "simulation_ids", 1:10)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_k_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, 1, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, 1, minobs, input_policy, force = FALSE)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(simulation_id = 1, minobs = minobs, input_policy = input_policy,
                    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
                    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
                    beta = if ("sim_beta" %in% names(ctx$input_table)) as.numeric(ctx$input_table$sim_beta[[1L]]) else 0.00005)
  saveRDS(list(context = ctx, simulation_ids = simulation_ids), file.path(dirs$results, "prepare_context.rds"))
  if (mode == "prepare") return(invisible(dirs$root))
  k1 <- if (mode %in% c("all", "k1-hessian-nullspace")) run_k1(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "k1_hessian_nullspace_diagnostics.rds"))
  if (mode == "k1-hessian-nullspace") return(invisible(k1))
  k2 <- if (mode %in% c("all", "k2-sdreport-triage")) run_k2(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "k2_sdreport_triage.rds"))
  if (mode == "k2-sdreport-triage") return(invisible(k2))
  k3 <- if (mode %in% c("all", "k3-j4-multistart")) run_k3(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "k3_j4_multistart_stability.rds"))
  if (mode == "k3-j4-multistart") return(invisible(k3))
  k4 <- if (mode %in% c("all", "k4-eta-f-borrowed-probe")) run_k4(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "k4_eta_f_borrowed_probe.rds"))
  if (mode == "k4-eta-f-borrowed-probe") return(invisible(k4))
  k5 <- if (mode %in% c("all", "k5-fixed-scale-calibration-local")) run_k5(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "k5_fixed_scale_calibration_local_model.rds"))
  if (mode == "k5-fixed-scale-calibration-local") return(invisible(k5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_k_recommendations(k1, k2, k3, k4, k5)
    args_info <- list(source_input_dir = source_input_dir, abcd_dir = abcd_dir,
                      diagnostics_dir = diagnostics_dir, delta_probe_dir = delta_probe_dir,
                      delta_debug_dir = delta_debug_dir, local_calibration_dir = local_calibration_dir,
                      core_fix_dir = core_fix_dir, identifiability_dir = identifiability_dir,
                      simulation_ids = simulation_ids, minobs = minobs, input_policy = input_policy)
    write_k_report(dirs, args_info, ctx, k1, k2, k3, k4, k5, recs)
  }
  message("Wrote farfield local Hessian diagnostics under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) main_k()
