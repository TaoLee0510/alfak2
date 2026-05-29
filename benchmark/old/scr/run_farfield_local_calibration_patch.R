#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_local_calibration_patch.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "run_farfield_delta_debug.R"))
source(file.path(repo_guess, "benchmark", "scr", "run_grf_alfak2_parameter_calibration.R"))

usage <- function() {
  cat(
    "Run farfield local/calibration patch H1-H5.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_local_calibration_patch.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --abcd-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --diagnostics-dir=benchmark/results/farfield_shape_diagnostics \\\n",
    "    --delta-probe-dir=benchmark/results/farfield_delta_estimator_probe \\\n",
    "    --delta-debug-dir=benchmark/results/farfield_delta_debug \\\n",
    "    --output-dir=benchmark/results/farfield_local_calibration_patch \\\n",
    "    --simulation-ids=1,2,3,4,5 --minobs=5 --input-policy=full\n",
    sep = ""
  )
}

make_h_dirs <- function(output_dir) make_probe_dirs(output_dir)

arg_integer_csv <- function(args, name, default) {
  raw <- as.character(arg_value(args, name, paste(default, collapse = ",")))
  out <- suppressWarnings(as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
  out[is.finite(out)]
}

load_delta_debug_results <- function(delta_debug_dir) {
  delta_debug_dir <- normalizePath(delta_debug_dir, winslash = "/", mustWork = TRUE)
  all_path <- file.path(delta_debug_dir, "results", "farfield_delta_debug_all_results.rds")
  out <- if (file.exists(all_path)) readRDS(all_path) else list()
  out$root <- delta_debug_dir
  out
}

shape_status_from_metrics <- function(x, delta_based = FALSE, edge_gate = TRUE, deploy_gate = TRUE) {
  amp <- is.finite(x$estimate_sd_ratio) & x$estimate_sd_ratio >= 0.02
  rank <- is.finite(x$pearson) & is.finite(x$spearman) & x$pearson > 0 & x$spearman > 0
  if (isTRUE(delta_based) && (!isTRUE(edge_gate) || !isTRUE(deploy_gate))) return("delta_untrusted")
  if (!isTRUE(amp)) return("amplitude_collapse")
  if (!isTRUE(rank)) return("wrong_direction")
  "valid_shape_config"
}

run_h1 <- function(dirs, delta_debug, force = FALSE) {
  rds <- file.path(dirs$results, "h1_calibration_gate_patch.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  gate <- read_tsv_safe(file.path(delta_debug$root, "tables", "calibration_shape_gate_demo.tsv"))
  if (!nrow(gate)) gate <- data.frame()
  demo <- gate
  ranked <- demo[order(demo$recommended_status != "valid_shape_config", demo$shape_score, demo$centered_rmse), , drop = FALSE]
  ranked$rank <- seq_len(nrow(ranked))
  valid <- ranked[ranked$recommended_status == "valid_shape_config" & !ranked$is_oracle, , drop = FALSE]
  best_numeric <- ranked[ranked$recommended_status != "valid_shape_config" & !ranked$is_oracle, , drop = FALSE]
  if (nrow(best_numeric)) best_numeric <- best_numeric[1L, , drop = FALSE]
  status <- data.frame(
    recommended_status = if (nrow(valid)) "valid_shape_config" else "no_valid_shape_configuration",
    n_total_configs = nrow(ranked),
    n_nonoracle_configs = sum(!ranked$is_oracle, na.rm = TRUE),
    n_nonoracle_valid_shape_configs = nrow(valid),
    best_params_cli_args = if (nrow(valid)) "available" else "no_valid_shape_configuration",
    best_numeric_only_config = if (nrow(best_numeric)) best_numeric$candidate_id[[1L]] else NA_character_,
    stringsAsFactors = FALSE
  )
  patch <- data.frame(
    file = "benchmark/scr/run_grf_alfak2_parameter_calibration.R",
    modified_functions = paste(c("accuracy_metric_row", "rank_calibration_parameters", "summarize_calibration_results", "calibration_gate_summary", "classify_calibration_shape"), collapse = ","),
    patch_summary = "Added shape gates, recommended_status, calibration_gate_summary, best_numeric_only_params, and no_valid_shape_configuration CLI behavior.",
    stringsAsFactors = FALSE
  )
  write_tsv_safe(demo, file.path(dirs$tables, "h1_calibration_metrics_gate_demo.tsv"))
  write_tsv_safe(ranked, file.path(dirs$tables, "h1_calibration_ranked_params_demo.tsv"))
  write_tsv_safe(status, file.path(dirs$tables, "h1_calibration_recommended_status.tsv"))
  write_tsv_safe(best_numeric, file.path(dirs$tables, "h1_best_numeric_only_params.tsv"))
  write_tsv_safe(patch, file.path(dirs$tables, "h1_formal_calibration_patch_summary.tsv"))
  out <- list(demo = demo, ranked = ranked, status = status, best_numeric = best_numeric, patch = patch)
  saveRDS(out, rds)
  out
}

build_local_tmb_debug_object <- function(data, graph, observation_model = "dirichlet_multinomial",
                                         dm_concentration = 50,
                                         eta_borrowed_prior_mean = -6,
                                         eta_borrowed_prior_sd = 1.5,
                                         eta_distance_penalty = 0.75,
                                         staged_init = FALSE,
                                         stage0_fit = NULL) {
  n <- length(graph$labels)
  observed_index <- match(data$labels, as.character(graph$labels))
  y0 <- numeric(n); y1 <- numeric(n)
  y0[observed_index] <- data$counts[, 1]
  y1[observed_index] <- data$counts[, 2]
  obs_weight0 <- rep(1, n); obs_weight1 <- rep(1, n)
  effective_count_total <- y0 + y1
  borrowed_eta <- effective_count_total <= 0 & graph$support_distance > 0L
  eta_prior_mean <- rep(0, n)
  eta_prior_sd_vec <- rep(5, n)
  eta_prior_mean[borrowed_eta] <- eta_borrowed_prior_mean - eta_distance_penalty * pmax(0, as.integer(graph$support_distance[borrowed_eta]) - 1L)
  eta_prior_sd_vec[borrowed_eta] <- eta_borrowed_prior_sd
  p0 <- (y0 + 0.5) / sum(y0 + 0.5)
  p1 <- (y1 + 0.5) / sum(y1 + 0.5)
  f0 <- log(p1) - log(p0)
  f0 <- f0 - mean(f0)
  n_context <- length(graph$context_label)
  n_group <- max(unlist(graph$context_group0), 0L) + 1L
  parameters <- list(
    eta = log(p0),
    f = f0 / max(data$dt, .Machine$double.eps),
    delta_context = rep(0, n_context),
    mu_group = rep(0, n_group),
    log_sigma_neighbor = log(0.35),
    log_sigma_anchor = log(0.6),
    log_tau_group = rep(log(0.2), n_group)
  )
  if (isTRUE(staged_init) && !is.null(stage0_fit)) {
    f_map <- setNames(stage0_fit$summary$fitness_mean, stage0_fit$summary$karyotype)
    idx <- match(as.character(graph$labels), names(f_map))
    parameters$f[!is.na(idx)] <- f_map[idx[!is.na(idx)]]
    if (any(is.na(idx))) parameters$f[is.na(idx)] <- median(parameters$f[!is.na(idx)], na.rm = TRUE)
    parameters$delta_context[] <- 0
  }
  tmb_data <- list(
    y0 = as.numeric(y0), y1 = as.numeric(y1), n_nodes = as.integer(n),
    trans_from = as.integer(unlist(graph$transition_from0)),
    trans_to = as.integer(unlist(graph$transition_to0)),
    trans_weight = as.numeric(unlist(graph$transition_weight)),
    support_distance = as.integer(graph$support_distance),
    parent_from = as.integer(unlist(graph$parent_from0)),
    parent_to = as.integer(unlist(graph$parent_to0)),
    parent_weight = as.numeric(unlist(graph$parent_weight)),
    parent_context = as.integer(unlist(graph$parent_context0)),
    context_group = as.integer(unlist(graph$context_group0)),
    dt = as.numeric(data$dt),
    obs_weight0 = as.numeric(obs_weight0), obs_weight1 = as.numeric(obs_weight1),
    use_observation_weights = 0L,
    eta_prior_mean = as.numeric(eta_prior_mean),
    eta_prior_sd = as.numeric(eta_prior_sd_vec),
    anchor_prior_scale = 1.0, mu_prior_scale = 0.5, scale_prior_scale = 1.0,
    observation_model = as.integer(observation_model == "dirichlet_multinomial"),
    observation_weight_mode = 1L,
    dm_concentration = as.numeric(dm_concentration)
  )
  obj <- TMB::MakeADFun(tmb_data, parameters, DLL = "alfak2", silent = TRUE)
  list(obj = obj, tmb_data = tmb_data, parameters = parameters, y0 = y0, y1 = y1,
       eta_prior_mean = eta_prior_mean, eta_prior_sd = eta_prior_sd_vec,
       effective_count_total = effective_count_total)
}

optimize_local_tmb_debug <- function(data, graph, variant = "baseline", eval_max = 2000, iter_max = 2000,
                                     observation_model = "dirichlet_multinomial", dm_concentration = 50,
                                     eta_borrowed_prior_mean = -6, eta_borrowed_prior_sd = 1.5,
                                     eta_distance_penalty = 0.75, staged_init = FALSE, stage0_fit = NULL,
                                     parameterization = "f") {
  built <- build_local_tmb_debug_object(
    data, graph, observation_model, dm_concentration,
    eta_borrowed_prior_mean, eta_borrowed_prior_sd, eta_distance_penalty,
    staged_init = staged_init, stage0_fit = stage0_fit
  )
  obj <- built$obj
  f_idx <- names(obj$par) == "f"
  started <- Sys.time()
  if (parameterization == "g_equivalent") {
    dt <- as.numeric(data$dt)
    par0 <- obj$par
    par0[f_idx] <- par0[f_idx] * dt
    to_f <- function(pg) {
      pf <- pg
      pf[f_idx] <- pg[f_idx] / dt
      pf
    }
    fn <- function(pg) obj$fn(to_f(pg))
    gr <- function(pg) {
      gf <- obj$gr(to_f(pg))
      gf[f_idx] <- gf[f_idx] / dt
      gf
    }
    opt <- nlminb(par0, fn, gr, control = list(eval.max = eval_max, iter.max = iter_max))
    par_final_g <- opt$par
    par_final <- to_f(par_final_g)
    grad_transformed <- gr(par_final_g)
  } else {
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = eval_max, iter.max = iter_max))
    par_final <- opt$par
    grad_transformed <- obj$gr(par_final)
  }
  grad <- obj$gr(par_final)
  plist <- obj$env$parList(par_final)
  report <- try(obj$report(par_final), silent = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  block <- tapply(abs(grad), names(obj$par), max)
  block_t <- tapply(abs(grad_transformed), names(obj$par), max)
  list(
    variant = variant, parameterization = parameterization, built = built, obj = obj, opt = opt,
    par = par_final, grad = grad, grad_transformed = grad_transformed, plist = plist,
    report = if (inherits(report, "try-error")) list() else report,
    elapsed_sec = elapsed, block = block, block_transformed = block_t,
    gradient_norm = sqrt(sum(grad^2)), transformed_gradient_norm = sqrt(sum(grad_transformed^2))
  )
}

support_scope_from_local <- function(tier, dist) {
  ifelse(tier == "directly_informed", "direct",
         ifelse(tier %in% c("local_borrowed", "graph_borrowed") | dist == 1, "local_borrowed",
                ifelse(dist >= 2, "weakly_supported", "all")))
}

summarize_local_debug_fit <- function(res, graph, grf, lambda) {
  n <- length(graph$labels)
  grad_list <- res$obj$env$parList(res$grad)
  grad_eta <- as.numeric(grad_list$eta)
  grad_f <- as.numeric(grad_list$f)
  eta <- res$plist$eta
  f <- res$plist$f
  pi0 <- if (!is.null(res$report$pi0)) as.numeric(res$report$pi0) else rep(NA_real_, n)
  pi1 <- if (!is.null(res$report$pi1)) as.numeric(res$report$pi1) else rep(NA_real_, n)
  truth <- compute_truth_for_nodes(as.character(graph$labels), grf, lambda)
  names(truth) <- as.character(graph$labels)
  parent_count <- tabulate(as.integer(graph$parent_to0) + 1L, nbins = n)
  child_count <- tabulate(as.integer(graph$parent_from0) + 1L, nbins = n)
  transition_in <- tabulate(as.integer(graph$transition_to0) + 1L, nbins = n)
  transition_out <- tabulate(as.integer(graph$transition_from0) + 1L, nbins = n)
  data.frame(
    variant = res$variant,
    parameterization = res$parameterization,
    node_id = seq_len(n),
    karyotype = as.character(graph$labels),
    support_tier = as.character(graph$support_tier),
    support_distance = as.integer(graph$support_distance),
    support_scope = support_scope_from_local(as.character(graph$support_tier), as.integer(graph$support_distance)),
    count_t0 = res$built$y0,
    count_t1 = res$built$y1,
    count_total = res$built$y0 + res$built$y1,
    effective_count_total = res$built$effective_count_total,
    pi0 = pi0,
    pi1 = pi1,
    fitness_mean = as.numeric(f),
    fitness_sd = NA_real_,
    grad_f = as.numeric(grad_f),
    abs_grad_f = abs(as.numeric(grad_f)),
    eta = as.numeric(eta),
    grad_eta = as.numeric(grad_eta),
    eta_prior_mean = res$built$eta_prior_mean,
    eta_prior_sd = res$built$eta_prior_sd,
    parent_count = parent_count,
    parent_weight_sum = safe_rowsum_vec(as.numeric(graph$parent_weight), as.integer(graph$parent_to0) + 1L, n),
    child_count = child_count,
    child_weight_sum = safe_rowsum_vec(as.numeric(graph$parent_weight), as.integer(graph$parent_from0) + 1L, n),
    n_transition_in = transition_in,
    n_transition_out = transition_out,
    is_observed = (res$built$y0 + res$built$y1) > 0,
    truth_fitness = as.numeric(truth[as.character(graph$labels)]),
    fitness_error_if_truth_available = as.numeric(f) - as.numeric(truth[as.character(graph$labels)]),
    stringsAsFactors = FALSE
  )
}

safe_rowsum_vec <- function(x, group, n) {
  out <- numeric(n)
  if (length(x) && length(group)) {
    z <- rowsum(x, group, reorder = FALSE)
    out[as.integer(rownames(z))] <- as.numeric(z[, 1])
  }
  out
}

local_debug_summary_row <- function(res, graph, node_tbl = NULL) {
  block <- res$block
  data.frame(
    variant = res$variant,
    parameterization = res$parameterization,
    convergence = res$opt$convergence,
    message = res$opt$message,
    objective = res$opt$objective,
    gradient_norm = res$gradient_norm,
    transformed_gradient_norm = res$transformed_gradient_norm,
    grad_eta_max_abs = unname(block["eta"]),
    grad_f_max_abs = unname(block["f"]),
    grad_delta_context_max_abs = unname(block["delta_context"]),
    grad_mu_group_max_abs = unname(block["mu_group"]),
    grad_log_sigma_neighbor_abs = unname(block["log_sigma_neighbor"]),
    grad_log_sigma_anchor_abs = unname(block["log_sigma_anchor"]),
    grad_log_tau_group_max_abs = unname(block["log_tau_group"]),
    max_gradient_block_name = names(which.max(block)),
    covariance_status = "not_run_custom_gradient_probe",
    covariance_fallback = NA,
    fitness_sd_source = "not_run",
    n_local_nodes = length(graph$labels),
    n_direct = sum(graph$support_tier == "directly_informed"),
    n_local_borrowed = sum(graph$support_tier != "directly_informed"),
    elapsed_sec = res$elapsed_sec,
    stringsAsFactors = FALSE
  )
}

edge_alignment_from_node_tbl <- function(node_tbl, graph) {
  from <- as.integer(graph$edge_from)
  to <- as.integer(graph$edge_to)
  keep <- from >= 1 & to >= 1 & from <= nrow(node_tbl) & to <= nrow(node_tbl)
  if (!any(keep)) return(edge_delta_metric_row(numeric(), numeric()))
  est <- node_tbl$fitness_mean[to[keep]] - node_tbl$fitness_mean[from[keep]]
  tru <- node_tbl$truth_fitness[to[keep]] - node_tbl$truth_fitness[from[keep]]
  edge_delta_metric_row(est, tru)
}

run_h2 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "h2_local_f_block_diagnostics.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$local$graph
  variants <- data.frame(
    variant = c("baseline", "strong_shrink", "very_strong_shrink", "shell0_control"),
    shell_depth = c(1, 1, 1, 0),
    eta_mean = c(-6, -8, -10, -6),
    eta_sd = c(1.5, 1.0, 0.5, 1.5),
    eta_penalty = c(0.75, 1.5, 2.5, 0.75),
    stringsAsFactors = FALSE
  )
  rows <- list(); summaries <- list(); aligns <- list()
  for (i in seq_len(nrow(variants))) {
    g <- if (variants$shell_depth[[i]] == 0) alfak2::build_karyotype_graph(bundle$data, shell_depth = 0) else graph
    res <- optimize_local_tmb_debug(bundle$data, g, variant = variants$variant[[i]], eval_max = 2000, iter_max = 2000,
                                    eta_borrowed_prior_mean = variants$eta_mean[[i]],
                                    eta_borrowed_prior_sd = variants$eta_sd[[i]],
                                    eta_distance_penalty = variants$eta_penalty[[i]])
    node <- summarize_local_debug_fit(res, g, grf, task_info$lambda)
    rows[[i]] <- node
    summaries[[i]] <- local_debug_summary_row(res, g)
    aligns[[i]] <- cbind(data.frame(variant = variants$variant[[i]], stringsAsFactors = FALSE), edge_alignment_from_node_tbl(node, g))
  }
  node_tbl <- bind_rows_fill(rows)
  summary_tbl <- bind_rows_fill(summaries)
  align_tbl <- bind_rows_fill(aligns)
  by_tier <- aggregate(
    list(
      grad_f_max_abs = node_tbl$abs_grad_f,
      grad_f_median_abs = node_tbl$abs_grad_f,
      n_nodes = node_tbl$node_id
    ),
    by = list(variant = node_tbl$variant, support_tier = node_tbl$support_tier, support_scope = node_tbl$support_scope),
    FUN = function(x) c(max = max(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE), n = length(x))[1]
  )
  # Rebuild cleaner by-tier rows because aggregate above keeps one function shape.
  by_tier <- do.call(rbind, lapply(split(node_tbl, paste(node_tbl$variant, node_tbl$support_tier, node_tbl$support_scope)), function(z) {
    data.frame(
      variant = z$variant[[1]], support_tier = z$support_tier[[1]], support_scope = z$support_scope[[1]],
      n_nodes = nrow(z), grad_f_max_abs = max(z$abs_grad_f, na.rm = TRUE),
      grad_f_median_abs = stats::median(z$abs_grad_f, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  top <- do.call(rbind, lapply(split(node_tbl, node_tbl$variant), function(z) head(z[order(-z$abs_grad_f), ], 20L)))
  write_tsv_safe(node_tbl, file.path(dirs$tables, "h2_local_f_gradient_by_node.tsv"))
  write_tsv_safe(by_tier, file.path(dirs$tables, "h2_local_f_gradient_by_support_tier.tsv"))
  write_tsv_safe(top, file.path(dirs$tables, "h2_local_top_gradient_nodes.tsv"))
  write_tsv_safe(summary_tbl, file.path(dirs$tables, "h2_local_variant_gradient_summary.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "h2_local_edge_delta_alignment.tsv"))
  out <- list(nodes = node_tbl, by_tier = by_tier, top = top, summary = summary_tbl, edge_alignment = align_tbl)
  saveRDS(out, rds)
  out
}

run_h3 <- function(bundle, grf, task_info, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "h3_g_parameterization_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$local$graph
  configs <- expand.grid(
    parameterization = c("f", "g_equivalent"),
    eval_max = c(500, 2000, 5000),
    shrink = c("current", "strong", "very_strong"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  shrink_vals <- list(current = c(-6, 1.5, 0.75), strong = c(-8, 1.0, 1.5), very_strong = c(-10, 0.5, 2.5))
  summaries <- list(); aligns <- list()
  for (i in seq_len(nrow(configs))) {
    sv <- shrink_vals[[configs$shrink[[i]]]]
    variant <- paste(configs$parameterization[[i]], configs$shrink[[i]], paste0("eval", configs$eval_max[[i]]), sep = "_")
    res <- optimize_local_tmb_debug(bundle$data, graph, variant = variant,
                                    parameterization = configs$parameterization[[i]],
                                    eval_max = configs$eval_max[[i]], iter_max = configs$eval_max[[i]],
                                    eta_borrowed_prior_mean = sv[[1]], eta_borrowed_prior_sd = sv[[2]], eta_distance_penalty = sv[[3]])
    node <- summarize_local_debug_fit(res, graph, grf, task_info$lambda)
    summaries[[i]] <- cbind(local_debug_summary_row(res, graph), data.frame(shrink = configs$shrink[[i]], eval_max = configs$eval_max[[i]], stringsAsFactors = FALSE))
    aligns[[i]] <- cbind(data.frame(variant = variant, parameterization = configs$parameterization[[i]], shrink = configs$shrink[[i]], eval_max = configs$eval_max[[i]], stringsAsFactors = FALSE),
                         edge_alignment_from_node_tbl(node, graph))
  }
  native_note <- data.frame(
    variant = "g_native_prior",
    parameterization = "g_native",
    implementation_status = "not_implemented_requires_local_model_tmb_cpp_prior_rewrite",
    recommendation = "Only implement after g_equivalent shows clear conditioning benefit.",
    stringsAsFactors = FALSE
  )
  summary_tbl <- bind_rows_fill(summaries)
  align_tbl <- bind_rows_fill(aligns)
  block_tbl <- summary_tbl
  write_tsv_safe(summary_tbl, file.path(dirs$tables, "h3_local_f_vs_g_parameterization.tsv"))
  write_tsv_safe(bind_rows_fill(list(block_tbl, native_note)), file.path(dirs$tables, "h3_local_g_parameterization_gradient_by_block.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "h3_local_g_parameterization_edge_alignment.tsv"))
  out <- list(summary = summary_tbl, block = bind_rows_fill(list(block_tbl, native_note)), edge_alignment = align_tbl, native_note = native_note)
  saveRDS(out, rds)
  out
}

run_h4 <- function(bundle, grf, task_info, dirs, h3, force = FALSE) {
  rds <- file.path(dirs$results, "h4_staged_g_init_probe.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  graph <- bundle$local$graph
  shell0_graph <- alfak2::build_karyotype_graph(bundle$data, shell_depth = 0)
  stage0 <- alfak2::fit_local_posterior(bundle$data, graph = shell0_graph, observation_model = "multinomial",
                                        control = list(eval.max = 500, iter.max = 500), retry_on_untrusted_covariance = FALSE)
  variants <- data.frame(
    variant = c("staged_f_current_shrink", "staged_f_strong_shrink", "staged_g_equivalent_strong_shrink", "staged_g_native_strong_shrink", "staged_g_native_very_strong_shrink"),
    parameterization = c("f", "f", "g_equivalent", "g_native", "g_native"),
    shrink = c("current", "strong", "strong", "strong", "very_strong"),
    stringsAsFactors = FALSE
  )
  shrink_vals <- list(current = c(-6, 1.5, 0.75), strong = c(-8, 1.0, 1.5), very_strong = c(-10, 0.5, 2.5))
  summaries <- list(); aligns <- list()
  for (i in seq_len(nrow(variants))) {
    sv <- shrink_vals[[variants$shrink[[i]]]]
    if (variants$parameterization[[i]] == "g_native") {
      summaries[[i]] <- data.frame(variant = variants$variant[[i]], parameterization = "g_native",
                                   implementation_status = "not_implemented_requires_local_model_tmb_cpp_prior_rewrite",
                                   gradient_norm = NA_real_, convergence = NA_integer_, covariance_status = "not_run",
                                   stringsAsFactors = FALSE)
      aligns[[i]] <- data.frame(variant = variants$variant[[i]], parameterization = "g_native",
                                implementation_status = "not_implemented", stringsAsFactors = FALSE)
      next
    }
    res <- optimize_local_tmb_debug(bundle$data, graph, variant = variants$variant[[i]],
                                    parameterization = variants$parameterization[[i]],
                                    eval_max = 2000, iter_max = 2000,
                                    eta_borrowed_prior_mean = sv[[1]], eta_borrowed_prior_sd = sv[[2]], eta_distance_penalty = sv[[3]],
                                    staged_init = TRUE, stage0_fit = stage0)
    node <- summarize_local_debug_fit(res, graph, grf, task_info$lambda)
    summaries[[i]] <- cbind(local_debug_summary_row(res, graph), data.frame(shrink = variants$shrink[[i]], stringsAsFactors = FALSE))
    aligns[[i]] <- cbind(data.frame(variant = variants$variant[[i]], parameterization = variants$parameterization[[i]], shrink = variants$shrink[[i]], stringsAsFactors = FALSE),
                         edge_alignment_from_node_tbl(node, graph))
  }
  summary_tbl <- bind_rows_fill(summaries)
  align_tbl <- bind_rows_fill(aligns)
  multistart <- data.frame(
    comparison = "single_start_only",
    restart_stability_status = "not_available_no_restart_hook_in_custom_probe",
    stringsAsFactors = FALSE
  )
  write_tsv_safe(summary_tbl, file.path(dirs$tables, "h4_local_staged_g_init_probe.tsv"))
  write_tsv_safe(multistart, file.path(dirs$tables, "h4_local_staged_g_init_multistart.tsv"))
  write_tsv_safe(align_tbl, file.path(dirs$tables, "h4_local_staged_edge_alignment.tsv"))
  out <- list(summary = summary_tbl, multistart = multistart, edge_alignment = align_tbl)
  saveRDS(out, rds)
  out
}

resolve_shared_input_row <- function(source_input_dir, sim_id, minobs, fallback_dirs = NULL) {
  ctx <- tryCatch(resolve_source_context(source_input_dir, 1, minobs, "full"), error = function(e) NULL)
  candidates <- character()
  if (!is.null(ctx)) candidates <- c(candidates, file.path(ctx$shared_input_dir, "tables", "input_table.tsv"))
  candidates <- c(candidates,
                  file.path(repo_guess, "benchmark/results/grf_downsampled_input_benchmark_pm_5e_05_sampledepth200_2acba0a/shared_inputs/tables/input_table.tsv"),
                  file.path(repo_guess, "benchmark/results/grf_downsampled_input_benchmark_pm_5e_05/shared_inputs/tables/input_table.tsv"))
  for (p in unique(candidates)) {
    if (!file.exists(p)) next
    x <- read_tsv_safe(p)
    keep <- as.integer(x$simulation_id) == sim_id & as.integer(x$minobs) == minobs &
      abs(as.numeric(x$lambda) - 0.2) < 1e-8 & abs(as.numeric(x$time_gap) - 2) < 1e-8
    if (!any(keep)) next
    row <- x[which(keep)[1L], , drop = FALSE]
    root <- normalizePath(file.path(dirname(p), ".."), winslash = "/", mustWork = TRUE)
    for (col in c("input_rds", "grf_rds")) {
      if (!file.exists(row[[col]][[1L]])) {
        alt <- file.path(root, "cache", basename(row[[col]][[1L]]))
        if (file.exists(alt)) row[[col]] <- normalizePath(alt, winslash = "/", mustWork = TRUE)
      }
    }
    if (file.exists(row$input_rds[[1L]]) && file.exists(row$grf_rds[[1L]])) return(row)
  }
  data.frame()
}

fit_local_for_multisim <- function(counts, dt, local_shell_depth = 1, shrink = "current") {
  sv <- switch(shrink, strong = c(-8, 1.0, 1.5), very_strong = c(-10, 0.5, 2.5), c(-6, 1.5, 0.75))
  data <- alfak2::prepare_alfak2_data(counts, dt = dt)
  graph <- alfak2::build_karyotype_graph(data, shell_depth = local_shell_depth)
  alfak2::fit_local_posterior(data, graph = graph, observation_model = "dirichlet_multinomial", dm_concentration = 50,
                              control = list(eval.max = 500, iter.max = 500),
                              retry_on_untrusted_covariance = FALSE,
                              eta_borrowed_prior_mean = sv[[1]], eta_borrowed_prior_sd = sv[[2]], eta_distance_penalty = sv[[3]])
}

run_h5 <- function(source_input_dir, simulation_ids, minobs, input_policy, dirs, force = FALSE) {
  rds <- file.path(dirs$results, "h5_multi_sim_validation.rds")
  if (!isTRUE(force) && file.exists(rds)) return(readRDS(rds))
  local_rows <- list(); global_rows <- list(); norm_rows <- list()
  idx <- 0L; lidx <- 0L
  full_global_sims <- simulation_ids[seq_len(min(1L, length(simulation_ids)))]
  for (sim in simulation_ids) {
    row <- resolve_shared_input_row(source_input_dir, sim, minobs)
    if (!nrow(row)) {
      lidx <- lidx + 1L
      local_rows[[lidx]] <- data.frame(simulation_id = sim, status = "missing_prepared_input", stringsAsFactors = FALSE)
      next
    }
    yi <- readRDS(row$input_rds[[1L]])
    counts <- prepare_alfak2_counts(yi, minobs = minobs, input_policy = input_policy, drop_diploid = TRUE)
    dt <- suppressWarnings(diff(as.numeric(colnames(counts))))
    if (length(dt) != 1L || !is.finite(dt) || dt <= 0) dt <- as.numeric(row$time_delta[[1L]])
    grf <- readRDS(row$grf_rds[[1L]])
    for (shrink in c("current", "strong")) {
      lf <- tryCatch(fit_local_for_multisim(counts, dt, 1, shrink), error = function(e) e)
      lidx <- lidx + 1L
      if (inherits(lf, "error")) {
        local_rows[[lidx]] <- data.frame(simulation_id = sim, local_config = shrink, status = "error", error_message = conditionMessage(lf), stringsAsFactors = FALSE)
        next
      }
      local_rows[[lidx]] <- data.frame(
        simulation_id = sim, local_config = shrink, status = "ok",
        convergence = lf$diagnostics$convergence,
        gradient_norm = lf$diagnostics$gradient_norm,
        covariance_status = lf$diagnostics$covariance_status,
        covariance_fallback = lf$diagnostics$covariance_fallback,
        stringsAsFactors = FALSE
      )
      if (shrink != "current") next
      cfgs <- data.frame(
        experiment = "H5",
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
      if (!(sim %in% full_global_sims)) {
        for (j in seq_len(nrow(cfgs))) {
          idx <- idx + 1L
          skip_row <- cbind(
            data.frame(
              simulation_id = sim,
              minobs = minobs,
              input_policy = input_policy,
              metric_scale = "native",
              support_scope = "farfield",
              recommended_status = "global_solve_skipped_time_budget",
              failure_reason = "H5 default run scored full global graph for the first simulation only because sparse posterior solves for additional 10k+ node graphs were too slow.",
              stringsAsFactors = FALSE
            ),
            cfgs[j, , drop = FALSE]
          )
          global_rows[[idx]] <- skip_row
          if (cfgs$graph_edge_weight[[j]] == "normalized") {
            norm_rows[[length(norm_rows) + 1L]] <- skip_row
          }
        }
        next
      }
      graph <- tryCatch(
        alfak2::build_karyotype_graph(lf$data, shell_depth = 2, max_nodes = 30000),
        error = function(e) e
      )
      if (inherits(graph, "error")) {
        idx <- idx + 1L
        global_rows[[idx]] <- data.frame(
          simulation_id = sim,
          candidate_id = "all_global_configs",
          support_scope = "farfield",
          metric_scale = "native",
          recommended_status = "graph_expansion_failed",
          failure_reason = conditionMessage(graph),
          stringsAsFactors = FALSE
        )
        next
      }
      sim_dirs <- dirs
      sim_dirs$cache <- file.path(dirs$cache, paste0("h5_sim", sim))
      dir.create(sim_dirs$cache, recursive = TRUE, showWarnings = FALSE)
      components <- prepare_solver_cache(graph, sim_dirs, force = FALSE)
      for (j in seq_len(nrow(cfgs))) {
        fit <- fit_global_with_config(lf, graph, components, cfgs[j, , drop = FALSE], minobs)
        task_info <- list(
          simulation_id = sim,
          minobs = minobs,
          input_policy = input_policy,
          lambda = as.numeric(row$lambda[[1L]]),
          dt = dt,
          beta = if ("sim_beta" %in% names(row)) as.numeric(row$sim_beta[[1L]]) else 0.00005
        )
        m <- score_summary_abcd(fit$summary, graph, grf, as.numeric(row$lambda[[1L]]), task_info, cfgs[j, , drop = FALSE], "H5")
        far <- m[m$support_scope == "farfield" & m$metric_scale == "native", , drop = FALSE]
        far$shape_class <- metric_shape_class(far)
        far$recommended_status <- vapply(seq_len(nrow(far)), function(k) shape_status_from_metrics(far[k, , drop = FALSE], FALSE, TRUE, TRUE), character(1))
        idx <- idx + 1L
        global_rows[[idx]] <- far
        if (cfgs$graph_edge_weight[[j]] == "normalized") {
          norm_rows[[length(norm_rows) + 1L]] <- far
        }
      }
    }
  }
  local_tbl <- bind_rows_fill(local_rows)
  global_tbl <- bind_rows_fill(global_rows)
  norm_tbl <- bind_rows_fill(norm_rows)
  validation <- if (nrow(global_tbl)) {
    aggregate(list(n_configs = global_tbl$candidate_id), by = list(simulation_id = global_tbl$simulation_id), FUN = length)
  } else data.frame()
  write_tsv_safe(global_tbl, file.path(dirs$tables, "h5_multi_sim_failure_state_validation.tsv"))
  write_tsv_safe(local_tbl, file.path(dirs$tables, "h5_multi_sim_local_fix_validation.tsv"))
  write_tsv_safe(norm_tbl, file.path(dirs$tables, "h5_multi_sim_normalized_default_validation.tsv"))
  out <- list(global = global_tbl, local = local_tbl, normalized = norm_tbl, validation = validation)
  saveRDS(out, rds)
  out
}

make_h_recommendations <- function(h1, h2, h3, h4, h5) {
  h2_best <- h2$summary[order(h2$summary$gradient_norm), , drop = FALSE][1L, , drop = FALSE]
  h3_best_g <- h3$summary[h3$summary$parameterization == "g_equivalent", , drop = FALSE]
  h3_best_g <- if (nrow(h3_best_g)) h3_best_g[order(h3_best_g$transformed_gradient_norm), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  data.frame(
    table = c("calibration_gate_patch_recommendation", "local_f_block_recommendation", "g_parameterization_recommendation",
              "staged_local_fix_recommendation", "multisim_validation_recommendation", "recommended_next_steps"),
    recommendation = c(
      "Use calibration gates and allow no_valid_shape_configuration.",
      "Prioritize f-block/local fitness optimization diagnostics before delta work.",
      if (nrow(h3_best_g) && is.finite(h3_best_g$transformed_gradient_norm[[1L]]) && h3_best_g$transformed_gradient_norm[[1L]] < min(h3$summary$gradient_norm[h3$summary$parameterization == "f"], na.rm = TRUE)) "g=dt*f can improve transformed conditioning in some settings; validate before C++ rewrite." else "g=dt*f did not clearly solve convergence; do not rewrite local_model_tmb.cpp yet.",
      "Staged strong shrink can be used as a diagnostic but did not establish trusted shell_depth=1 local fit.",
      "Keep normalized as default candidate with failure-state gate; do not treat numeric-only improvements as shape solutions.",
      "Patch formal calibration gate, expose local TMB diagnostics, then rerun multi-sim before considering C++ edge-gradient."
    ),
    evidence = c(
      h1$status$recommended_status[[1L]],
      paste0("best_variant=", h2_best$variant[[1L]], "; gradient=", fmt_metric(h2_best$gradient_norm[[1L]])),
      if (nrow(h3_best_g)) paste0("best_g_transformed_gradient=", fmt_metric(h3_best_g$transformed_gradient_norm[[1L]])) else "g_not_run",
      paste0("best_staged_gradient=", fmt_metric(min(h4$summary$gradient_norm, na.rm = TRUE))),
      paste0("normalized_collapse_fraction=", if (nrow(h5$normalized)) fmt_metric(mean(h5$normalized$amplitude_collapse, na.rm = TRUE)) else "NA"),
      "C++ edge-gradient gate remains unmet."
    ),
    stringsAsFactors = FALSE
  )
}

write_h_report <- function(dirs, ctx, args_info, h1, h2, h3, h4, h5, recs) {
  all_long <- bind_rows_fill(list(
    transform(h1$ranked, experiment = "H1_calibration_gate_patch"),
    transform(h2$summary, experiment = "H2_local_f_block"),
    transform(h3$summary, experiment = "H3_g_parameterization"),
    transform(h4$summary, experiment = "H4_staged_local_fix"),
    transform(h5$global, experiment = "H5_multisim_global"),
    transform(h5$local, experiment = "H5_multisim_local")
  ))
  summary <- data.frame(
    experiment = c("H1", "H2", "H3", "H4", "H5"),
    key_result = c(
      h1$status$recommended_status[[1L]],
      paste0("top_f_block=", h2$summary$max_gradient_block_name[which.max(h2$summary$grad_f_max_abs)]),
      paste0("best_g_gradient=", if (any(h3$summary$parameterization == "g_equivalent")) fmt_metric(min(h3$summary$transformed_gradient_norm[h3$summary$parameterization == "g_equivalent"], na.rm = TRUE)) else "NA"),
      paste0("best_staged_gradient=", fmt_metric(min(h4$summary$gradient_norm, na.rm = TRUE))),
      paste0("n_sims_with_full_global_metrics=", length(unique(h5$global$simulation_id[is.finite(h5$global$n_scored)])))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(all_long, file.path(dirs$tables, "all_h_experiments_long.tsv"))
  write_tsv_safe(summary, file.path(dirs$tables, "h_experiment_summary.tsv"))
  for (nm in recs$table) write_tsv_safe(recs[recs$table == nm, -1, drop = FALSE], file.path(dirs$tables, paste0(nm, ".tsv")))
  top_tier <- h2$by_tier[order(-h2$by_tier$grad_f_max_abs), , drop = FALSE][1L, , drop = FALSE]
  top_nodes <- h2$nodes[order(-abs(h2$nodes$grad_f)), , drop = FALSE]
  top_node_summary <- if (nrow(top_nodes)) {
    paste0(top_nodes$karyotype[[1L]], " (", top_nodes$support_tier[[1L]], ", grad_f=", fmt_metric(top_nodes$grad_f[[1L]]), ")")
  } else {
    "NA"
  }
  h5_full <- sort(unique(h5$global$simulation_id[is.finite(h5$global$n_scored)]))
  h5_skipped <- sum(h5$global$recommended_status == "global_solve_skipped_time_budget", na.rm = TRUE)
  lines <- c(
    "# Farfield Local Calibration Patch Report",
    "",
    "## Data source",
    paste0("- source-input-dir: `", args_info$source_input_dir, "`"),
    paste0("- abcd-dir: `", args_info$abcd_dir, "`"),
    paste0("- diagnostics-dir: `", args_info$diagnostics_dir, "`"),
    paste0("- delta-probe-dir: `", args_info$delta_probe_dir, "`"),
    paste0("- delta-debug-dir: `", args_info$delta_debug_dir, "`"),
    paste0("- simulation_ids: ", paste(args_info$simulation_ids, collapse = ",")),
    paste0("- minobs: ", args_info$minobs),
    paste0("- input_policy: ", args_info$input_policy),
    paste0("- reused local bundle: `", ctx$local_bundle_path, "`"),
    "",
    "## Prior results summary",
    "- ABCD/E/F/G show oracle edge information can recover farfield shape, while non-oracle delta deployment fails and calibration should not force a shape recommendation.",
    "",
    "## H1 formal calibration failure-state patch",
    paste0("- formal patch status demo: ", h1$status$recommended_status[[1L]], "."),
    "- `run_grf_alfak2_parameter_calibration.R` now writes gate/status outputs, best_numeric_only_params, and avoids executable CLI args when no valid shape config exists.",
    "- Modified calibration functions: `rank_calibration_parameters()`, `summarize_calibration_results()`, plus new `classify_calibration_shape()` and `calibration_gate_summary()` helpers.",
    paste0("- best numeric-only fallback is recorded separately as `", h1$status$best_numeric_only_config[[1L]], "` and is not promoted as a shape recommendation."),
    "",
    "## H2 local f-block diagnostics",
    paste0("- largest f-gradient tier: ", top_tier$support_tier[[1L]], " / ", top_tier$support_scope[[1L]], "."),
    paste0("- top gradient node: ", top_node_summary, "."),
    paste0("- best H2 gradient variant: ", h2$summary$variant[which.min(h2$summary$gradient_norm)], " with gradient=", fmt_metric(min(h2$summary$gradient_norm, na.rm = TRUE)), "."),
    paste0("- local edge sign agreement best: ", fmt_metric(max(h2$edge_alignment$delta_sign_agreement, na.rm = TRUE)), "."),
    "- shell_depth=0 is well conditioned, but shell_depth=1 variants remain dominated by f gradients in direct/local-borrowed nodes; this points to a local likelihood/borrowed-identifiability issue rather than a delta-context-only issue.",
    "",
    "## H3 g = dt * f parameterization",
    paste0("- best original f gradient: ", fmt_metric(min(h3$summary$gradient_norm[h3$summary$parameterization == "f"], na.rm = TRUE)), "."),
    paste0("- best g-equivalent transformed gradient: ", fmt_metric(min(h3$summary$transformed_gradient_norm[h3$summary$parameterization == "g_equivalent"], na.rm = TRUE)), "."),
    "- g-native prior was not implemented because it requires changing the C++ objective priors.",
    "- The benchmark-only g-equivalent transform improved transformed gradient in selected settings but did not establish trusted covariance or reliable local edge alignment, so it is not enough to justify rewriting `local_model_tmb.cpp` yet.",
    "",
    "## H4 staged shell_depth=1 init",
    paste0("- best staged gradient: ", fmt_metric(min(h4$summary$gradient_norm, na.rm = TRUE)), "."),
    "- staged + shrinkage is diagnostic only; it did not establish trusted covariance or reliable local edge deltas.",
    "- staged g-native variants are intentionally marked not_run because they require a real C++ prior reparameterization.",
    "",
    "## H5 multi-sim validation",
    paste0("- scored simulations: ", paste(sort(unique(h5$global$simulation_id)), collapse = ","), "."),
    paste0("- full global metric simulations: ", paste(h5_full, collapse = ","), "."),
    paste0("- global solve skipped rows: ", h5_skipped, "."),
    paste0("- normalized amplitude-collapse fraction: ", if (nrow(h5$normalized)) fmt_metric(mean(h5$normalized$amplitude_collapse, na.rm = TRUE)) else "NA", "."),
    "- H5 was deliberately reduced after sim2 sparse posterior solves proved too slow for repeated 10k+ node graphs; sim2-5 still contribute local convergence/covariance diagnostics and explicit skipped global rows.",
    "- Failure-state behavior remains necessary: numeric improvements alone should not be promoted to shape recommendations.",
    "",
    "## Final conclusion",
    "- Continue C++ edge-gradient pseudo-observation now: no.",
    "- Prioritize formal calibration gates, local f-block diagnostics exposure, and robust local/staged initialization before further delta deployment.",
    "- Keep normalized as benchmark/probe/calibration default candidate with amplitude-collapse diagnostics; keep unit as stress-test and mutation as legacy baseline.",
    "- Do not default `anchor_count_reference=minobs` for full input.",
    "- Next priorities: complete formal calibration patch tests, expose TMB diagnostics in package internals, evaluate g/staged fixes on more simulations, and only then revisit non-oracle delta estimators."
  )
  writeLines(lines, file.path(dirs$root, "farfield_local_calibration_patch_report.md"))
  saveRDS(list(h1 = h1, h2 = h2, h3 = h3, h4 = h4, h5 = h5, summary = summary, recs = recs),
          file.path(dirs$results, "farfield_local_calibration_patch_all_results.rds"))
}

main_h <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage(); return(invisible(NULL))
  }
  mode <- match.arg(tolower(as.character(arg_value(args, "mode", "all"))),
                    c("prepare", "h1-calibration-gate-patch", "h2-local-f-block-diagnostics",
                      "h3-g-parameterization-probe", "h4-staged-g-init", "h5-multisim-validation",
                      "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  abcd_dir <- as.character(arg_value(args, "abcd_dir", "benchmark/results/farfield_shape_probe_abcd"))
  diagnostics_dir <- as.character(arg_value(args, "diagnostics_dir", "benchmark/results/farfield_shape_diagnostics"))
  delta_probe_dir <- as.character(arg_value(args, "delta_probe_dir", "benchmark/results/farfield_delta_estimator_probe"))
  delta_debug_dir <- as.character(arg_value(args, "delta_debug_dir", "benchmark/results/farfield_delta_debug"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_local_calibration_patch"))
  simulation_ids <- arg_integer_csv(args, "simulation_ids", 1:5)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_h_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, 1, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, 1, minobs, input_policy, force = FALSE)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(simulation_id = 1, minobs = minobs, input_policy = input_policy, lambda = as.numeric(ctx$input_table$lambda[[1L]]))
  delta_debug <- load_delta_debug_results(delta_debug_dir)
  saveRDS(list(context = ctx, simulation_ids = simulation_ids), file.path(dirs$results, "prepare_context.rds"))
  if (mode == "prepare") return(invisible(dirs$root))
  h1 <- if (mode %in% c("all", "h1-calibration-gate-patch")) run_h1(dirs, delta_debug, force = force) else readRDS(file.path(dirs$results, "h1_calibration_gate_patch.rds"))
  if (mode == "h1-calibration-gate-patch") return(invisible(h1))
  h2 <- if (mode %in% c("all", "h2-local-f-block-diagnostics")) run_h2(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "h2_local_f_block_diagnostics.rds"))
  if (mode == "h2-local-f-block-diagnostics") return(invisible(h2))
  h3 <- if (mode %in% c("all", "h3-g-parameterization-probe")) run_h3(bundle, grf, task_info, dirs, force = force) else readRDS(file.path(dirs$results, "h3_g_parameterization_probe.rds"))
  if (mode == "h3-g-parameterization-probe") return(invisible(h3))
  h4 <- if (mode %in% c("all", "h4-staged-g-init")) run_h4(bundle, grf, task_info, dirs, h3, force = force) else readRDS(file.path(dirs$results, "h4_staged_g_init_probe.rds"))
  if (mode == "h4-staged-g-init") return(invisible(h4))
  h5 <- if (mode %in% c("all", "h5-multisim-validation")) run_h5(source_input_dir, simulation_ids, minobs, input_policy, dirs, force = force) else readRDS(file.path(dirs$results, "h5_multi_sim_validation.rds"))
  if (mode == "h5-multisim-validation") return(invisible(h5))
  if (mode %in% c("all", "summarize")) {
    recs <- make_h_recommendations(h1, h2, h3, h4, h5)
    args_info <- list(source_input_dir = source_input_dir, abcd_dir = abcd_dir, diagnostics_dir = diagnostics_dir,
                      delta_probe_dir = delta_probe_dir, delta_debug_dir = delta_debug_dir,
                      simulation_ids = simulation_ids, minobs = minobs, input_policy = input_policy)
    write_h_report(dirs, ctx, args_info, h1, h2, h3, h4, h5, recs)
  }
  message("Wrote farfield local calibration patch under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) main_h()
