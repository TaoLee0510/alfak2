#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "benchmark/scr/run_farfield_shape_abcd.R"
script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
repo_guess <- normalizePath(file.path(dirname(script_file), "../.."), winslash = "/", mustWork = FALSE)
source(file.path(repo_guess, "benchmark", "scr", "probe_farfield_shape_configs.R"))

usage <- function() {
  cat(
    "Run far-field shape experiments A-D on cached prepared GRF inputs.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_farfield_shape_abcd.R --mode=all \\\n",
    "    --source-input-dir=benchmark/results/farfield_shape_probe_default \\\n",
    "    --output-dir=benchmark/results/farfield_shape_probe_abcd \\\n",
    "    --simulation-id=1 --minobs=5 --input-policy=full\n\n",
    "Modes:\n",
    "  prepare, experiment-a, experiment-b, experiment-c, experiment-d, summarize, all\n\n",
    "Options:\n",
    "  --quick=auto                 auto|true|false; auto keeps the full A grid but limits expensive B-D fanout if needed.\n",
    "  --cv-splits=20               direct-anchor holdout CV splits before quick-mode limiting.\n",
    "  --max-cv-candidates=40        maximum CV candidates before quick-mode limiting.\n",
    "  --force=false                recompute cached experiment fits.\n",
    sep = ""
  )
}

read_tsv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
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

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

path_under_benchmark <- function(path) {
  path <- as.character(path)
  if (!grepl("^/", path)) path <- file.path(repo_guess, path)
  norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  bench <- normalizePath(file.path(repo_guess, "benchmark"), winslash = "/", mustWork = TRUE)
  startsWith(norm, paste0(bench, "/")) || identical(norm, bench)
}

as_bool_or_auto <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("auto", "a")) return("auto")
  if (x %in% c("true", "t", "1", "yes", "y")) return("true")
  if (x %in% c("false", "f", "0", "no", "n")) return("false")
  stop("Invalid --quick value: ", x, call. = FALSE)
}

config_id <- function(prefix, graph_edge_weight, lambda_l, lambda_e, sigma_obs,
                      anchor_var_mode = "current", prior_mean_mode = "zero",
                      prior_mean_scale = 0, anchor_count_reference_mode = "none") {
  token <- function(x) gsub("-", "m", gsub("\\.", "p", format(x, scientific = TRUE, trim = TRUE)))
  paste(
    prefix,
    graph_edge_weight,
    paste0("ll", token(lambda_l)),
    paste0("le", token(lambda_e)),
    paste0("so", token(sigma_obs)),
    anchor_var_mode,
    prior_mean_mode,
    paste0("pm", token(prior_mean_scale)),
    anchor_count_reference_mode,
    sep = "__"
  )
}

make_baseline_config <- function(prefix = "baseline") {
  data.frame(
    experiment = prefix,
    candidate_id = "baseline_mutation_ll0p2_le1_so0p05",
    graph_edge_weight = "mutation",
    lambda_l = 0.2,
    lambda_e = 1,
    sigma_obs = 0.05,
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    solver = "matrix_mean",
    stringsAsFactors = FALSE
  )
}

make_probe_reference_config <- function(prefix = "baseline_probe_reference") {
  out <- make_baseline_config(prefix)
  out$candidate_id <- "baseline_probe_reference_mutation_ll0p2_le1_so0p02"
  out$sigma_obs <- 0.02
  out
}

make_experiment_a_grid <- function() {
  grid <- expand.grid(
    graph_edge_weight = c("normalized", "unit"),
    lambda_l = c(0.05, 0.2, 1, 5, 20, 100),
    lambda_e = c(0.005, 0.01, 0.05, 0.25),
    sigma_obs = c(0.02, 0.05, 0.1, 0.2),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$experiment <- "A"
  grid$anchor_var_mode <- "current"
  grid$prior_mean_mode <- "zero"
  grid$prior_mean_scale <- 0
  grid$anchor_count_reference_mode <- "none"
  grid$solver <- "matrix_mean"
  grid$candidate_id <- mapply(
    config_id,
    prefix = "A",
    graph_edge_weight = grid$graph_edge_weight,
    lambda_l = grid$lambda_l,
    lambda_e = grid$lambda_e,
    sigma_obs = grid$sigma_obs,
    anchor_var_mode = grid$anchor_var_mode,
    prior_mean_mode = grid$prior_mean_mode,
    prior_mean_scale = grid$prior_mean_scale,
    anchor_count_reference_mode = grid$anchor_count_reference_mode,
    USE.NAMES = FALSE
  )
  baseline <- make_baseline_config("A")
  probe_ref <- make_probe_reference_config("A_probe")
  out <- rbind(baseline, probe_ref, grid[, names(baseline), drop = FALSE])
  out[!duplicated(out$candidate_id), , drop = FALSE]
}

safe_num <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y)] <- default
  y
}

safe_cor2 <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

quantile_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  unname(stats::quantile(x, 0.95, names = FALSE) - stats::quantile(x, 0.05, names = FALSE))
}

pairwise_rank_loss <- function(pred, truth, margin = 0.01) {
  ok <- is.finite(pred) & is.finite(truth)
  pred <- pred[ok]
  truth <- truth[ok]
  if (length(pred) < 3L) return(NA_real_)
  cmb <- utils::combn(seq_along(pred), 2L)
  truth_sign <- sign(truth[cmb[1L, ]] - truth[cmb[2L, ]])
  keep <- truth_sign != 0
  if (!any(keep)) return(NA_real_)
  pred_margin <- truth_sign[keep] * (pred[cmb[1L, keep]] - pred[cmb[2L, keep]])
  mean(pmax(0, margin - pred_margin)^2)
}

shape_score_from_parts <- function(centered_rmse, spearman, estimate_sd_ratio, false_high_rate) {
  cr <- ifelse(is.finite(centered_rmse), centered_rmse, 10)
  sp <- ifelse(is.finite(spearman), spearman, -1)
  amp <- abs(log(pmax(ifelse(is.finite(estimate_sd_ratio), estimate_sd_ratio, 1e-4), 1e-4)))
  fh <- ifelse(is.finite(false_high_rate), false_high_rate, 1)
  cr - 0.25 * sp + 0.25 * amp + 0.10 * fh
}

compute_scope_metrics <- function(est, est_sd, truth, keep, n_nodes, scope_name) {
  keep <- keep & is.finite(est) & is.finite(truth)
  if (!any(keep)) {
    return(data.frame(
      support_scope = scope_name,
      n_nodes = as.integer(n_nodes),
      n_scored = 0L,
      mae = NA_real_,
      rmse = NA_real_,
      centered_mae = NA_real_,
      centered_rmse = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_,
      signed_bias = NA_real_,
      false_high_rate = NA_real_,
      sign_accuracy = NA_real_,
      estimate_sd = NA_real_,
      truth_sd = NA_real_,
      estimate_sd_ratio = NA_real_,
      estimate_range_ratio = NA_real_,
      estimate_iqr_ratio = NA_real_,
      mean_estimated_sd = NA_real_,
      shape_score = shape_score_from_parts(NA_real_, NA_real_, NA_real_, NA_real_),
      amplitude_collapse = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  x <- est[keep]
  y <- truth[keep]
  s <- est_sd[keep]
  err <- x - y
  xc <- x - mean(x)
  yc <- y - mean(y)
  centered_err <- xc - yc
  est_sd_val <- stats::sd(x)
  truth_sd_val <- stats::sd(y)
  sd_ratio <- if (is.finite(truth_sd_val) && truth_sd_val > 0) est_sd_val / truth_sd_val else NA_real_
  range_truth <- quantile_range(y)
  iqr_truth <- stats::IQR(y, na.rm = TRUE)
  range_ratio <- if (is.finite(range_truth) && range_truth > 0) quantile_range(x) / range_truth else NA_real_
  iqr_ratio <- if (is.finite(iqr_truth) && iqr_truth > 0) stats::IQR(x, na.rm = TRUE) / iqr_truth else NA_real_
  spearman <- safe_cor2(x, y, "spearman")
  centered_rmse <- sqrt(mean(centered_err^2))
  false_high_rate <- mean(xc > 0 & yc <= 0, na.rm = TRUE)
  data.frame(
    support_scope = scope_name,
    n_nodes = as.integer(n_nodes),
    n_scored = as.integer(length(x)),
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2)),
    centered_mae = mean(abs(centered_err)),
    centered_rmse = centered_rmse,
    pearson = safe_cor2(x, y, "pearson"),
    spearman = spearman,
    signed_bias = mean(err),
    false_high_rate = false_high_rate,
    sign_accuracy = mean(sign(xc) == sign(yc), na.rm = TRUE),
    estimate_sd = est_sd_val,
    truth_sd = truth_sd_val,
    estimate_sd_ratio = sd_ratio,
    estimate_range_ratio = range_ratio,
    estimate_iqr_ratio = iqr_ratio,
    mean_estimated_sd = if (any(is.finite(s))) mean(s, na.rm = TRUE) else NA_real_,
    shape_score = shape_score_from_parts(centered_rmse, spearman, sd_ratio, false_high_rate),
    amplitude_collapse = !is.finite(sd_ratio) || sd_ratio < 0.02,
    stringsAsFactors = FALSE
  )
}

score_summary_abcd <- function(summary, graph, grf, lambda, task_info, config, prior_mean_status = "none") {
  summary <- add_legacy_scale_if_possible(summary, graph, task_info$dt, task_info$beta)
  truth_map <- compute_grf_truth(summary$karyotype, grf$centroids, lambda)
  truth <- as.numeric(truth_map[as.character(summary$karyotype)])
  tier <- as.character(summary$support_tier)
  distance <- as.integer(summary$support_distance)
  scopes <- list(
    direct = tier == "directly_informed",
    local_borrowed = tier == "local_borrowed",
    weakly_supported = tier == "weakly_supported",
    farfield = distance >= 2L | tier %in% c("weakly_supported", "graph_borrowed", "prior_dominated"),
    all = rep(TRUE, nrow(summary))
  )
  for (st in sort(unique(tier))) scopes[[paste0("tier_", st)]] <- tier == st
  for (dd in sort(unique(distance[is.finite(distance)]))) scopes[[paste0("distance_", dd)]] <- distance == dd
  scales <- c(native = "fitness_mean")
  if ("fitness_mean_alfakR_scale" %in% names(summary)) scales <- c(scales, alfakR_scale = "fitness_mean_alfakR_scale")
  rows <- list()
  idx <- 0L
  for (scale_name in names(scales)) {
    est <- as.numeric(summary[[scales[[scale_name]]]])
    est_sd <- if (identical(scale_name, "alfakR_scale") && "fitness_sd_alfakR_scale" %in% names(summary)) {
      as.numeric(summary$fitness_sd_alfakR_scale)
    } else {
      as.numeric(summary$fitness_sd)
    }
    for (scope_name in names(scopes)) {
      keep0 <- scopes[[scope_name]]
      m <- compute_scope_metrics(est, est_sd, truth, keep0, sum(keep0, na.rm = TRUE), scope_name)
      idx <- idx + 1L
      rows[[idx]] <- cbind(
        data.frame(
          simulation_id = task_info$simulation_id,
          minobs = task_info$minobs,
          input_policy = task_info$input_policy,
          experiment = config$experiment,
          candidate_id = config$candidate_id,
          graph_edge_weight = config$graph_edge_weight,
          lambda_l = config$lambda_l,
          lambda_e = config$lambda_e,
          sigma_obs = config$sigma_obs,
          anchor_var_mode = config$anchor_var_mode,
          prior_mean_mode = config$prior_mean_mode,
          prior_mean_scale = config$prior_mean_scale,
          anchor_count_reference_mode = config$anchor_count_reference_mode,
          metric_scale = scale_name,
          prior_mean_status = prior_mean_status,
          stringsAsFactors = FALSE
        ),
        m
      )
    }
  }
  do.call(rbind, rows)
}

resolve_source_context <- function(source_input_dir, simulation_id, minobs, input_policy) {
  source_input_dir <- normalizePath(source_input_dir, winslash = "/", mustWork = TRUE)
  task_id <- paste0("sim", simulation_id, "_minobs", minobs, "_", input_policy)
  probe_rds <- file.path(source_input_dir, "farfield_shape_probe.rds")
  source_probe_dir <- if (file.exists(probe_rds)) source_input_dir else NULL
  shared_dir <- source_input_dir
  if (!is.null(source_probe_dir)) {
    probe <- readRDS(probe_rds)
    shared_dir <- normalizePath(probe$source_input_dir, winslash = "/", mustWork = TRUE)
  } else {
    probe <- NULL
  }
  local_bundle_path <- if (!is.null(source_probe_dir)) {
    file.path(source_probe_dir, "fits", task_id, "local_bundle.rds")
  } else {
    NA_character_
  }
  input_table_path <- file.path(shared_dir, "tables", "input_table.tsv")
  if (!file.exists(input_table_path)) stop("Missing prepared input table: ", input_table_path, call. = FALSE)
  input_tbl <- read_tsv_safe(input_table_path)
  input_tbl <- input_tbl[input_tbl$simulation_id == simulation_id & input_tbl$minobs == minobs, , drop = FALSE]
  if (!nrow(input_tbl)) stop("No prepared input row for simulation/minobs.", call. = FALSE)
  list(
    source_probe_dir = source_probe_dir,
    shared_input_dir = shared_dir,
    local_bundle_path = local_bundle_path,
    input_table = input_tbl[1L, , drop = FALSE],
    probe = probe
  )
}

prepare_abcd_bundle <- function(ctx, dirs, simulation_id, minobs, input_policy,
                                local_shell_depth = 1L, global_extra_shell = 1L,
                                max_nodes = 150000L, force = FALSE) {
  cache_path <- file.path(dirs$cache, paste0("local_bundle_sim", simulation_id, "_minobs", minobs, "_", input_policy, ".rds"))
  if (!isTRUE(force) && file.exists(cache_path)) return(readRDS(cache_path))
  if (length(ctx$local_bundle_path) && !is.na(ctx$local_bundle_path) && file.exists(ctx$local_bundle_path)) {
    bundle <- readRDS(ctx$local_bundle_path)
    saveRDS(bundle, cache_path)
    return(bundle)
  }
  message("No cached local bundle found; rebuilding one local fit.")
  task <- ctx$input_table
  yi <- readRDS(task$input_rds)
  counts <- prepare_counts(yi, minobs = minobs, input_policy = input_policy)
  selected_times <- suppressWarnings(as.numeric(colnames(counts)))
  dt <- if (length(selected_times) == 2L && all(is.finite(selected_times))) diff(selected_times) else NA_real_
  if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
  beta <- as.numeric(task$sim_pm)
  pkgload::load_all(repo_guess, quiet = TRUE)
  prep_input_depth <- get("prepare_counts_for_input_depth", envir = asNamespace("alfak2"))
  resolve_obs <- get("resolve_fit_observation_controls", envir = asNamespace("alfak2"))
  k_mat <- parse_karyotype_ids(rownames(counts))
  max_cn <- max(k_mat, na.rm = TRUE) + local_shell_depth + global_extra_shell
  obs <- resolve_obs("effective", "dirichlet_multinomial", 50)
  data <- prep_input_depth(
    counts,
    dt = dt,
    beta = beta,
    input_depth = "effective",
    effective_depth = NULL,
    effective_depth_mode = "min",
    effective_depth_rounding = "hash",
    effective_depth_seed = NULL
  )
  local_graph <- alfak2::build_karyotype_graph(
    data,
    transition_kernel = "exact",
    shell_depth = local_shell_depth,
    min_cn = 0L,
    max_cn = as.integer(max_cn),
    max_nodes = max_nodes
  )
  local <- alfak2::fit_local_posterior(
    data,
    local_graph,
    observation_model = obs$observation_model,
    dm_concentration = obs$dm_concentration,
    observation_weight_mode = "likelihood",
    control = list(eval.max = 500L, iter.max = 500L),
    retry_control = list(eval.max = 2000L, iter.max = 2000L)
  )
  global_graph <- alfak2::build_karyotype_graph(
    data,
    transition_kernel = "exact",
    shell_depth = local_shell_depth + global_extra_shell,
    min_cn = 0L,
    max_cn = as.integer(max_cn),
    max_nodes = max_nodes
  )
  bundle <- list(
    data = data,
    local = local,
    global_graph = global_graph,
    counts = counts,
    minobs = minobs,
    input_policy = input_policy,
    local_shell_depth = local_shell_depth,
    global_extra_shell = global_extra_shell
  )
  saveRDS(bundle, cache_path)
  bundle
}

row_keys <- function(mat) {
  apply(mat, 1L, paste, collapse = ".")
}

make_sparse_penalty <- function(n, rows, cols, vals) {
  Matrix::sparseMatrix(
    i = as.integer(rows),
    j = as.integer(cols),
    x = as.numeric(vals),
    dims = c(n, n),
    giveCsparse = TRUE
  )
}

add_penalty_terms_r <- function(ii, jj, xx, coeff_index, coeff_value, scale) {
  k <- length(coeff_index)
  if (!k || !is.finite(scale) || scale == 0) return(list(i = ii, j = jj, x = xx))
  ii <- c(ii, rep(as.integer(coeff_index), each = k))
  jj <- c(jj, rep(as.integer(coeff_index), times = k))
  xx <- c(xx, as.numeric(scale) * rep(as.numeric(coeff_value), each = k) * rep(as.numeric(coeff_value), times = k))
  list(i = ii, j = jj, x = xx)
}

build_curvature_penalty <- function(graph) {
  kary <- graph$karyotypes
  n <- nrow(kary)
  p <- ncol(kary)
  keys <- row_keys(kary)
  idx_env <- new.env(parent = emptyenv(), hash = TRUE, size = n * 2L)
  for (i in seq_len(n)) assign(keys[[i]], i, envir = idx_env)
  ii_list <- vector("list", n)
  jj_list <- vector("list", n)
  xx_list <- vector("list", n)
  for (i in seq_len(n)) {
    x <- kary[i, ]
    ii <- integer()
    jj <- integer()
    xx <- numeric()
    for (cc in seq_len(p)) {
      lo <- x
      hi <- x
      lo[[cc]] <- lo[[cc]] - 1L
      hi[[cc]] <- hi[[cc]] + 1L
      ilo <- get0(paste(lo, collapse = "."), envir = idx_env, inherits = FALSE, ifnotfound = NA_integer_)
      ihi <- get0(paste(hi, collapse = "."), envir = idx_env, inherits = FALSE, ifnotfound = NA_integer_)
      if (is.finite(ilo) && is.finite(ihi)) {
        tmp <- add_penalty_terms_r(ii, jj, xx, c(ilo, i, ihi), c(1, -2, 1), 1)
        ii <- tmp$i
        jj <- tmp$j
        xx <- tmp$x
      }
    }
    for (c1 in seq_len(p - 1L)) {
      for (c2 in seq.int(c1 + 1L, p)) {
        xc <- x
        xd <- x
        xcd <- x
        xc[[c1]] <- xc[[c1]] + 1L
        xd[[c2]] <- xd[[c2]] + 1L
        xcd[[c1]] <- xcd[[c1]] + 1L
        xcd[[c2]] <- xcd[[c2]] + 1L
        ic <- get0(paste(xc, collapse = "."), envir = idx_env, inherits = FALSE, ifnotfound = NA_integer_)
        id <- get0(paste(xd, collapse = "."), envir = idx_env, inherits = FALSE, ifnotfound = NA_integer_)
        icd <- get0(paste(xcd, collapse = "."), envir = idx_env, inherits = FALSE, ifnotfound = NA_integer_)
        if (is.finite(ic) && is.finite(id) && is.finite(icd)) {
          tmp <- add_penalty_terms_r(ii, jj, xx, c(i, ic, id, icd), c(1, -1, -1, 1), 0.5)
          ii <- tmp$i
          jj <- tmp$j
          xx <- tmp$x
        }
      }
    }
    if (length(ii)) {
      ii_list[[i]] <- ii
      jj_list[[i]] <- jj
      xx_list[[i]] <- xx
    }
  }
  ii <- unlist(ii_list, use.names = FALSE)
  jj <- unlist(jj_list, use.names = FALSE)
  xx <- unlist(xx_list, use.names = FALSE)
  make_sparse_penalty(n, ii, jj, xx)
}

build_edge_penalty <- function(graph, mode) {
  weight_fun <- get("graph_edge_weights", envir = asNamespace("alfak2"))
  w <- weight_fun(graph$edge_weight, mode)
  from <- as.integer(graph$edge_from)
  to <- as.integer(graph$edge_to)
  keep <- is.finite(w) & w > 0 & from != to
  from <- from[keep]
  to <- to[keep]
  w <- w[keep]
  make_sparse_penalty(
    length(graph$labels),
    c(from, to, from, to),
    c(from, to, to, from),
    c(w, w, -w, -w)
  )
}

prepare_solver_cache <- function(graph, dirs, force = FALSE) {
  cache_path <- file.path(dirs$cache, "matrix_solver_components.rds")
  if (!isTRUE(force) && file.exists(cache_path)) return(readRDS(cache_path))
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.", call. = FALSE)
  message("Building Matrix solver components for ", length(graph$labels), " global nodes.")
  started <- Sys.time()
  components <- list(
    labels = as.character(graph$labels),
    edge = list(
      mutation = build_edge_penalty(graph, "mutation"),
      normalized = build_edge_penalty(graph, "normalized"),
      unit = build_edge_penalty(graph, "unit")
    ),
    curvature = build_curvature_penalty(graph),
    eps = 1e-5,
    built_at = as.character(Sys.time())
  )
  components$elapsed_sec <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  saveRDS(components, cache_path)
  components
}

clone_local_with_anchor_var_mode <- function(local_fit, mode) {
  out <- local_fit
  mode <- as.character(mode)
  if (mode %in% c("current", "count_inflated")) return(out)
  value <- switch(
    mode,
    constant_0.05 = 0.05,
    constant_0.10 = 0.10,
    constant_0.20 = 0.20,
    stop("Unsupported anchor_var_mode: ", mode, call. = FALSE)
  )
  direct <- as.character(out$summary$support_tier) == "directly_informed" & is.finite(out$summary$fitness_mean)
  out$summary$fitness_sd[direct] <- value
  out$summary$fitness_sd_source[direct] <- paste0("constant_", format(value, trim = TRUE))
  out$summary$covariance_status[direct] <- "TMB_sdreport"
  if (!is.null(out$diagnostics)) {
    out$diagnostics$covariance_status <- "TMB_sdreport_constant_anchor_override"
    out$diagnostics$covariance_fallback <- FALSE
    out$diagnostics$fitness_sd_source <- paste0("constant_", format(value, trim = TRUE))
  }
  out
}

anchor_data_for_fit <- function(local_fit, graph, anchor_count_reference = NULL,
                                anchor_count_power = 1, anchor_exclude = character()) {
  anchor_match <- match(local_fit$summary$karyotype, as.character(graph$labels))
  tier_ok <- as.character(local_fit$summary$support_tier) == "directly_informed"
  if (length(anchor_exclude)) tier_ok <- tier_ok & !(as.character(local_fit$summary$karyotype) %in% as.character(anchor_exclude))
  anchor_count_all <- if ("effective_count_total" %in% names(local_fit$summary)) {
    as.numeric(local_fit$summary$effective_count_total)
  } else if ("count_total" %in% names(local_fit$summary)) {
    as.numeric(local_fit$summary$count_total)
  } else {
    rep(NA_real_, nrow(local_fit$summary))
  }
  count_ok <- is.finite(anchor_count_all) & anchor_count_all > 0
  keep <- which(!is.na(anchor_match) & is.finite(local_fit$summary$fitness_mean) & tier_ok & count_ok)
  if (!length(keep)) stop("No usable direct anchors.", call. = FALSE)
  cov_mult_fun <- get("anchor_covariance_multiplier", envir = asNamespace("alfak2"))
  count_mult_fun <- get("count_anchor_multiplier", envir = asNamespace("alfak2"))
  covariance_status <- as.character(local_fit$summary$covariance_status[keep])
  covariance_mult <- cov_mult_fun(covariance_status)
  count_mult <- count_mult_fun(anchor_count_all[keep], anchor_count_reference, anchor_count_power)
  var_base <- as.numeric(local_fit$summary$fitness_sd[keep])^2
  data.frame(
    node_id = as.integer(anchor_match[keep]),
    karyotype = as.character(local_fit$summary$karyotype[keep]),
    mean = as.numeric(local_fit$summary$fitness_mean[keep]),
    variance_base = var_base,
    variance = var_base * covariance_mult * count_mult,
    variance_multiplier = covariance_mult * count_mult,
    covariance_status = covariance_status,
    count_total = if ("count_total" %in% names(local_fit$summary)) as.numeric(local_fit$summary$count_total[keep]) else NA_real_,
    effective_count_total = anchor_count_all[keep],
    stringsAsFactors = FALSE
  )
}

support_tier_from_graph <- function(graph, posterior_sd) {
  tier <- as.character(graph$support_tier)
  tier[graph$support_distance > 2L] <- "graph_borrowed"
  if (any(is.finite(posterior_sd))) {
    prior_dominated <- posterior_sd > stats::quantile(posterior_sd, 0.9, na.rm = TRUE) &
      !(tier %in% c("directly_informed", "local_borrowed", "weakly_supported"))
    tier[prior_dominated] <- "prior_dominated"
  }
  tier
}

fit_global_with_config <- function(local_fit, graph, components, config, minobs,
                                   prior_mean = NULL, anchor_exclude = character()) {
  local_use <- clone_local_with_anchor_var_mode(local_fit, config$anchor_var_mode)
  anchor_ref <- if (identical(as.character(config$anchor_var_mode), "count_inflated") ||
                    identical(as.character(config$anchor_count_reference_mode), "minobs")) {
    as.numeric(minobs)
  } else {
    NULL
  }
  if (!is.null(prior_mean)) {
    idx <- match(as.character(local_use$summary$karyotype), as.character(graph$labels))
    ok <- !is.na(idx) & is.finite(local_use$summary$fitness_mean) & is.finite(prior_mean[idx])
    local_use$summary$fitness_mean[ok] <- local_use$summary$fitness_mean[ok] - prior_mean[idx[ok]]
  }
  anchors <- anchor_data_for_fit(local_use, graph, anchor_count_reference = anchor_ref,
                                 anchor_count_power = 1, anchor_exclude = anchor_exclude)
  n <- length(graph$labels)
  precision <- 1 / pmax(1e-10, anchors$variance + as.numeric(config$sigma_obs)^2)
  q_anchor <- Matrix::sparseMatrix(
    i = anchors$node_id,
    j = anchors$node_id,
    x = precision,
    dims = c(n, n),
    giveCsparse = TRUE
  )
  q <- as.numeric(config$lambda_l) * components$edge[[as.character(config$graph_edge_weight)]] +
    as.numeric(config$lambda_e) * components$curvature +
    q_anchor +
    Matrix::Diagonal(n, components$eps)
  q <- Matrix::forceSymmetric(q, uplo = "U")
  rhs <- numeric(n)
  rhs[anchors$node_id] <- rhs[anchors$node_id] + anchors$mean * precision
  started <- Sys.time()
  fit <- tryCatch({
    chol <- Matrix::Cholesky(q, LDL = TRUE, perm = TRUE)
    mean <- as.numeric(Matrix::solve(chol, rhs))
    approx_sd <- sqrt(1 / pmax(as.numeric(Matrix::diag(q)), 1e-10))
    if (!is.null(prior_mean)) {
      mean <- mean + prior_mean
    }
    summary <- data.frame(
      node_id = seq_along(graph$labels),
      karyotype = as.character(graph$labels),
      support_tier = support_tier_from_graph(graph, approx_sd),
      support_distance = as.integer(graph$support_distance),
      fitness_mean = mean,
      fitness_sd = approx_sd,
      fitness_sd_source = "diagonal_precision_approx",
      conf_low = mean - 1.959963984540054 * approx_sd,
      conf_high = mean + 1.959963984540054 * approx_sd,
      stringsAsFactors = FALSE
    )
    list(
      graph = graph,
      summary = summary,
      anchors = anchors,
      hyperparameters = list(
        lambda_l = as.numeric(config$lambda_l),
        lambda_e = as.numeric(config$lambda_e),
        sigma_obs = as.numeric(config$sigma_obs),
        cv_score = NA_real_,
        cv_status = "fixed_hyperparameters_matrix_mean",
        graph_edge_weight = as.character(config$graph_edge_weight),
        anchor_count_reference = anchor_ref,
        anchor_count_power = 1,
        anchor_excluded = length(anchor_exclude)
      ),
      diagnostics = list(
        factorization_status = "ok",
        solver = "Matrix::Cholesky_posterior_mean",
        posterior_sd_status = "diagonal_precision_approx",
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    )
  }, error = function(e) {
    list(
      graph = graph,
      summary = data.frame(),
      anchors = anchors,
      hyperparameters = list(
        lambda_l = as.numeric(config$lambda_l),
        lambda_e = as.numeric(config$lambda_e),
        sigma_obs = as.numeric(config$sigma_obs),
        cv_status = "error",
        graph_edge_weight = as.character(config$graph_edge_weight),
        anchor_count_reference = anchor_ref
      ),
      diagnostics = list(
        factorization_status = "error",
        solver = "Matrix::Cholesky_posterior_mean",
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    )
  })
  fit
}

fit_cached <- function(local_fit, graph, components, config, grf, task_info, dirs,
                       prior_mean = NULL, prior_mean_status = "none", force = FALSE) {
  fit_path <- file.path(dirs$fits, paste0(config$candidate_id, ".rds"))
  metric_path <- file.path(dirs$results, paste0(config$candidate_id, "_metrics.rds"))
  if (!isTRUE(force) && file.exists(fit_path) && file.exists(metric_path)) {
    return(list(fit = readRDS(fit_path), metrics = readRDS(metric_path), fit_path = fit_path))
  }
  fit <- fit_global_with_config(local_fit, graph, components, config, task_info$minobs, prior_mean = prior_mean)
  metrics <- if (nrow(fit$summary)) {
    score_summary_abcd(fit$summary, fit$graph, grf, task_info$lambda, task_info, config, prior_mean_status)
  } else {
    data.frame()
  }
  saveRDS(fit, fit_path)
  saveRDS(metrics, metric_path)
  list(fit = fit, metrics = metrics, fit_path = fit_path)
}

build_prior_mean_zero <- function(graph) {
  list(mean = rep(0, length(graph$labels)), status = "zero")
}

propagate_prior_by_context <- function(graph, seed_mean, delta_context, scale) {
  out <- rep(NA_real_, length(graph$labels))
  seed_ok <- is.finite(seed_mean)
  out[seed_ok] <- seed_mean[seed_ok]
  max_d <- max(graph$support_distance, na.rm = TRUE)
  from <- as.integer(graph$parent_from0) + 1L
  to <- as.integer(graph$parent_to0) + 1L
  ctx <- as.integer(graph$parent_context0) + 1L
  w <- as.numeric(graph$parent_weight)
  for (dd in seq_len(max_d)) {
    children <- which(graph$support_distance == dd & !is.finite(out))
    if (!length(children)) next
    for (child in children) {
      e <- which(to == child & is.finite(out[from]) & ctx >= 1L & ctx <= length(delta_context))
      if (!length(e)) next
      vals <- out[from[e]] + as.numeric(scale) * delta_context[ctx[e]]
      ww <- w[e]
      ww[!is.finite(ww) | ww <= 0] <- 1
      out[child] <- stats::weighted.mean(vals, ww, na.rm = TRUE)
    }
  }
  out[!is.finite(out)] <- 0
  out
}

direct_seed_mean <- function(local_fit, graph) {
  seed <- rep(NA_real_, length(graph$labels))
  direct <- as.character(local_fit$summary$support_tier) == "directly_informed" & is.finite(local_fit$summary$fitness_mean)
  idx <- match(as.character(local_fit$summary$karyotype[direct]), as.character(graph$labels))
  ok <- !is.na(idx)
  seed[idx[ok]] <- as.numeric(local_fit$summary$fitness_mean[direct][ok])
  seed
}

build_prior_mean_local_context_delta <- function(local_fit, graph, scale) {
  delta <- local_fit$parameter_mode$delta_context
  if (is.null(delta) || !length(delta) || !any(is.finite(delta))) {
    return(list(mean = rep(0, length(graph$labels)), status = "local_context_delta_missing_fallback_zero"))
  }
  seed <- direct_seed_mean(local_fit, graph)
  list(
    mean = propagate_prior_by_context(graph, seed, as.numeric(delta), scale),
    status = "local_context_delta"
  )
}

build_prior_mean_empirical_edge_delta <- function(local_fit, graph, scale) {
  local_graph <- local_fit$graph
  ls <- local_fit$summary
  f <- setNames(as.numeric(ls$fitness_mean), as.character(ls$karyotype))
  labels <- as.character(local_graph$labels)
  from <- as.integer(local_graph$parent_from0) + 1L
  to <- as.integer(local_graph$parent_to0) + 1L
  ctx <- as.integer(local_graph$parent_context0) + 1L
  ok <- from >= 1L & to >= 1L & from <= length(labels) & to <= length(labels) &
    is.finite(f[labels[from]]) & is.finite(f[labels[to]]) &
    ctx >= 1L
  n_ctx <- length(graph$context_label)
  delta_context <- rep(NA_real_, n_ctx)
  if (any(ok)) {
    delta <- as.numeric(f[labels[to[ok]]] - f[labels[from[ok]]])
    ctx_ok <- ctx[ok]
    for (cc in seq_len(n_ctx)) {
      vals <- delta[ctx_ok == cc]
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 2L) delta_context[[cc]] <- stats::median(vals)
    }
    global_med <- stats::median(delta[is.finite(delta)], na.rm = TRUE)
    if (!is.finite(global_med)) global_med <- 0
    delta_context[!is.finite(delta_context)] <- global_med
    status <- "empirical_edge_delta"
  } else {
    delta_context[] <- 0
    status <- "empirical_edge_delta_no_edges_fallback_zero"
  }
  seed <- direct_seed_mean(local_fit, graph)
  list(
    mean = propagate_prior_by_context(graph, seed, delta_context, scale),
    status = status,
    components = data.frame(
      context_index = seq_along(delta_context),
      context_label = as.character(graph$context_label),
      empirical_delta = delta_context,
      stringsAsFactors = FALSE
    )
  )
}

build_prior_mean <- function(mode, local_fit, graph, scale) {
  mode <- as.character(mode)
  if (identical(mode, "zero")) return(build_prior_mean_zero(graph))
  if (identical(mode, "local_context_delta")) return(build_prior_mean_local_context_delta(local_fit, graph, scale))
  if (identical(mode, "empirical_edge_delta")) return(build_prior_mean_empirical_edge_delta(local_fit, graph, scale))
  stop("Unsupported prior_mean_mode: ", mode, call. = FALSE)
}

select_farfield_native <- function(metrics) {
  metrics[metrics$support_scope == "farfield" & metrics$metric_scale == "native", , drop = FALSE]
}

dedupe_configs <- function(x) {
  key <- interaction(
    x$graph_edge_weight, x$lambda_l, x$lambda_e, x$sigma_obs,
    x$anchor_var_mode, x$prior_mean_mode, x$prior_mean_scale, x$anchor_count_reference_mode,
    drop = TRUE
  )
  x[!duplicated(key), , drop = FALSE]
}

bind_rows_fill <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && nrow(x), rows)
  if (!length(rows)) return(data.frame())
  nms <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(nms, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, nms, drop = FALSE]
  })
  do.call(rbind, rows)
}

run_experiment_a <- function(bundle, components, grf, task_info, dirs, force = FALSE) {
  configs <- make_experiment_a_grid()
  metrics <- list()
  manifest <- list()
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    message("[A] ", i, "/", nrow(configs), " ", cfg$candidate_id)
    started <- Sys.time()
    res <- fit_cached(bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs, force = force)
    metrics[[i]] <- res$metrics
    manifest[[i]] <- data.frame(
      experiment = "A",
      candidate_id = cfg$candidate_id,
      status = if (nrow(res$fit$summary)) "ok" else "error",
      fit_path = res$fit_path,
      elapsed_sec = res$fit$diagnostics$elapsed_sec %||% as.numeric(difftime(Sys.time(), started, units = "secs")),
      solver = "Matrix::Cholesky_posterior_mean",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, metrics)
  write_tsv_safe(out, file.path(dirs$tables, "experiment_A_global_scale_grid.tsv"))
  far <- select_farfield_native(out)
  top_score <- far[order(far$shape_score, far$centered_rmse), , drop = FALSE]
  top_score$ranking_mode <- "shape_score"
  top_spearman <- far[order(-far$spearman, far$centered_rmse), , drop = FALSE]
  top_spearman$ranking_mode <- "spearman"
  top_pearson <- far[order(-far$pearson, far$centered_rmse), , drop = FALSE]
  top_pearson$ranking_mode <- "pearson"
  top_amp <- far[order(abs(log(pmax(far$estimate_sd_ratio, 1e-4)))), , drop = FALSE]
  top_amp$ranking_mode <- "estimate_sd_ratio_closest_1"
  top <- rbind(head(top_score, 20L), head(top_spearman, 20L), head(top_pearson, 20L), head(top_amp, 20L))
  top$all_amplitude_collapse <- all(far$estimate_sd_ratio < 0.02, na.rm = TRUE)
  write_tsv_safe(top, file.path(dirs$tables, "experiment_A_top_configs.tsv"))
  saveRDS(list(configs = configs, metrics = out, manifest = do.call(rbind, manifest)),
          file.path(dirs$results, "experiment_A_global_scale_grid.rds"))
  out
}

configs_from_metrics <- function(metrics, experiment = "candidate") {
  cols <- c(
    "candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
    "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode"
  )
  out <- unique(metrics[, intersect(cols, names(metrics)), drop = FALSE])
  out$experiment <- experiment
  out$solver <- "matrix_mean"
  out[, c("experiment", "candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
          "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode", "solver"), drop = FALSE]
}

select_b_candidates <- function(a_metrics) {
  far <- select_farfield_native(a_metrics)
  eligible <- far[is.finite(far$estimate_sd_ratio) & far$estimate_sd_ratio >= 0.02, , drop = FALSE]
  pick <- if (nrow(eligible)) eligible[order(eligible$shape_score), , drop = FALSE] else far[order(far$shape_score), , drop = FALSE]
  cfg <- configs_from_metrics(head(pick, 10L), "B")
  best_by_weight <- do.call(rbind, lapply(split(far, far$graph_edge_weight), function(x) head(x[order(x$shape_score), , drop = FALSE], 1L)))
  cfg <- rbind(cfg, configs_from_metrics(best_by_weight, "B"), make_baseline_config("B"))
  cfg <- dedupe_configs(cfg)
  cfg
}

run_experiment_b <- function(bundle, components, grf, task_info, dirs, a_metrics, force = FALSE) {
  base_cfg <- select_b_candidates(a_metrics)
  modes <- c("current", "constant_0.05", "constant_0.10", "constant_0.20", "count_inflated")
  configs <- do.call(rbind, lapply(modes, function(mode) {
    x <- base_cfg
    x$experiment <- "B"
    x$anchor_var_mode <- mode
    x$anchor_count_reference_mode <- ifelse(mode == "count_inflated", "minobs", "none")
    x$candidate_id <- mapply(
      config_id,
      prefix = "B",
      graph_edge_weight = x$graph_edge_weight,
      lambda_l = x$lambda_l,
      lambda_e = x$lambda_e,
      sigma_obs = x$sigma_obs,
      anchor_var_mode = x$anchor_var_mode,
      prior_mean_mode = x$prior_mean_mode,
      prior_mean_scale = x$prior_mean_scale,
      anchor_count_reference_mode = x$anchor_count_reference_mode,
      USE.NAMES = FALSE
    )
    x
  }))
  configs <- dedupe_configs(configs)
  metrics <- list()
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    message("[B] ", i, "/", nrow(configs), " ", cfg$candidate_id)
    metrics[[i]] <- fit_cached(bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs, force = force)$metrics
  }
  out <- do.call(rbind, metrics)
  write_tsv_safe(out, file.path(dirs$tables, "experiment_B_anchor_variance.tsv"))
  far <- select_farfield_native(out)
  summary <- aggregate(
    far[, c("centered_rmse", "pearson", "spearman", "estimate_sd_ratio", "shape_score")],
    by = list(anchor_var_mode = far$anchor_var_mode),
    FUN = function(x) median(x, na.rm = TRUE)
  )
  summary$n_configs <- as.integer(table(far$anchor_var_mode)[summary$anchor_var_mode])
  write_tsv_safe(summary, file.path(dirs$tables, "experiment_B_anchor_variance_summary.tsv"))
  saveRDS(list(configs = configs, metrics = out, summary = summary),
          file.path(dirs$results, "experiment_B_anchor_variance.rds"))
  out
}

select_c_base_configs <- function(a_metrics, b_metrics) {
  far <- rbind(select_farfield_native(a_metrics), select_farfield_native(b_metrics))
  far <- far[order(far$shape_score, far$centered_rmse), , drop = FALSE]
  cfg <- configs_from_metrics(head(far, 10L), "C")
  forced <- data.frame(
    experiment = "C",
    candidate_id = c("forced_normalized_ll0p2_le0p01_so0p05", "forced_unit_ll0p2_le0p01_so0p05"),
    graph_edge_weight = c("normalized", "unit"),
    lambda_l = c(0.2, 0.2),
    lambda_e = c(0.01, 0.01),
    sigma_obs = c(0.05, 0.05),
    anchor_var_mode = "current",
    prior_mean_mode = "zero",
    prior_mean_scale = 0,
    anchor_count_reference_mode = "none",
    solver = "matrix_mean",
    stringsAsFactors = FALSE
  )
  cfg <- rbind(cfg, forced, make_baseline_config("C"))
  dedupe_configs(cfg)
}

run_experiment_c <- function(bundle, components, grf, task_info, dirs, a_metrics, b_metrics, force = FALSE) {
  base_cfg <- select_c_base_configs(a_metrics, b_metrics)
  priors <- rbind(
    data.frame(prior_mean_mode = "zero", prior_mean_scale = 0, stringsAsFactors = FALSE),
    expand.grid(
      prior_mean_mode = c("local_context_delta", "empirical_edge_delta"),
      prior_mean_scale = c(0.25, 0.5, 1.0),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  configs <- do.call(rbind, lapply(seq_len(nrow(priors)), function(i) {
    x <- base_cfg
    x$experiment <- "C"
    x$prior_mean_mode <- priors$prior_mean_mode[[i]]
    x$prior_mean_scale <- priors$prior_mean_scale[[i]]
    x$candidate_id <- mapply(
      config_id,
      prefix = "C",
      graph_edge_weight = x$graph_edge_weight,
      lambda_l = x$lambda_l,
      lambda_e = x$lambda_e,
      sigma_obs = x$sigma_obs,
      anchor_var_mode = x$anchor_var_mode,
      prior_mean_mode = x$prior_mean_mode,
      prior_mean_scale = x$prior_mean_scale,
      anchor_count_reference_mode = x$anchor_count_reference_mode,
      USE.NAMES = FALSE
    )
    x
  }))
  configs <- dedupe_configs(configs)
  metrics <- list()
  components_rows <- list()
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    message("[C] ", i, "/", nrow(configs), " ", cfg$candidate_id)
    pm <- build_prior_mean(cfg$prior_mean_mode, bundle$local, bundle$global_graph, cfg$prior_mean_scale)
    if (!is.null(pm$components)) {
      cc <- pm$components
      cc$candidate_id <- cfg$candidate_id
      cc$prior_mean_mode <- cfg$prior_mean_mode
      cc$prior_mean_scale <- cfg$prior_mean_scale
      components_rows[[length(components_rows) + 1L]] <- cc
    }
    metrics[[i]] <- fit_cached(
      bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs,
      prior_mean = pm$mean,
      prior_mean_status = pm$status,
      force = force
    )$metrics
  }
  out <- do.call(rbind, metrics)
  write_tsv_safe(out, file.path(dirs$tables, "experiment_C_prior_mean_edge_slope.tsv"))
  comp <- if (length(components_rows)) do.call(rbind, components_rows) else data.frame()
  write_tsv_safe(comp, file.path(dirs$tables, "experiment_C_prior_mean_components.tsv"))
  far <- select_farfield_native(out)
  top <- far[order(far$shape_score, far$centered_rmse), , drop = FALSE]
  write_tsv_safe(head(top, 40L), file.path(dirs$tables, "experiment_C_top_configs.tsv"))
  saveRDS(list(configs = configs, metrics = out, components = comp),
          file.path(dirs$results, "experiment_C_prior_mean_edge_slope.rds"))
  out
}

select_cv_candidates <- function(a_metrics, c_metrics, max_candidates = 40L) {
  far_a <- select_farfield_native(a_metrics)
  far_c <- select_farfield_native(c_metrics)
  a_top <- head(far_a[order(far_a$shape_score), , drop = FALSE], 20L)
  c_top <- head(far_c[order(far_c$shape_score), , drop = FALSE], 20L)
  cfg <- rbind(configs_from_metrics(a_top, "D"), configs_from_metrics(c_top, "D"), make_baseline_config("D"))
  cfg <- dedupe_configs(cfg)
  cfg <- cfg[seq_len(min(nrow(cfg), as.integer(max_candidates))), , drop = FALSE]
  cfg
}

cv_metric_row <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L) {
    return(data.frame(
      cv_mse = NA_real_,
      cv_centered_rmse = NA_real_,
      cv_spearman = NA_real_,
      cv_pairwise_rank_loss = NA_real_,
      cv_estimate_sd_ratio = NA_real_,
      cv_shape_score = shape_score_from_parts(NA_real_, NA_real_, NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ))
  }
  pred <- pred[ok]
  truth <- truth[ok]
  pc <- pred - mean(pred)
  tc <- truth - mean(truth)
  sd_ratio <- if (stats::sd(truth) > 0) stats::sd(pred) / stats::sd(truth) else NA_real_
  data.frame(
    cv_mse = mean((pred - truth)^2),
    cv_centered_rmse = sqrt(mean((pc - tc)^2)),
    cv_spearman = safe_cor2(pred, truth, "spearman"),
    cv_pairwise_rank_loss = pairwise_rank_loss(pred, truth),
    cv_estimate_sd_ratio = sd_ratio,
    cv_shape_score = shape_score_from_parts(sqrt(mean((pc - tc)^2)), safe_cor2(pred, truth, "spearman"), sd_ratio, mean(pc > 0 & tc <= 0)),
    stringsAsFactors = FALSE
  )
}

run_direct_anchor_holdout_cv <- function(bundle, components, configs, task_info, cv_splits = 20L, seed = 29011L) {
  local <- bundle$local
  graph <- bundle$global_graph
  direct <- local$summary[
    as.character(local$summary$support_tier) == "directly_informed" &
      is.finite(local$summary$fitness_mean),
    ,
    drop = FALSE
  ]
  labels <- intersect(as.character(direct$karyotype), as.character(graph$labels))
  if (length(labels) < 4L) cv_splits <- min(cv_splits, length(labels))
  set.seed(seed)
  split_list <- vector("list", cv_splits)
  for (s in seq_len(cv_splits)) {
    n_hold <- if (length(labels) <= 6L) 1L else max(2L, min(length(labels) - 3L, round(length(labels) * 0.25)))
    split_list[[s]] <- sample(labels, n_hold)
  }
  rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, , drop = FALSE]
    pm <- build_prior_mean(cfg$prior_mean_mode, local, graph, cfg$prior_mean_scale)
    truth_anchor <- setNames(as.numeric(local$summary$fitness_mean), as.character(local$summary$karyotype))
    if (!identical(cfg$prior_mean_mode, "zero")) {
      gidx <- match(names(truth_anchor), as.character(graph$labels))
      ok <- !is.na(gidx) & is.finite(pm$mean[gidx])
      truth_anchor[ok] <- truth_anchor[ok]
    }
    for (s in seq_along(split_list)) {
      hold <- split_list[[s]]
      fit <- fit_global_with_config(local, graph, components, cfg, task_info$minobs,
                                    prior_mean = pm$mean, anchor_exclude = hold)
      pred <- fit$summary$fitness_mean[match(hold, as.character(fit$summary$karyotype))]
      truth <- as.numeric(truth_anchor[hold])
      m <- cv_metric_row(pred, truth)
      idx <- idx + 1L
      rows[[idx]] <- cbind(
        data.frame(
          split_id = s,
          n_holdout = length(hold),
          holdout_labels = paste(hold, collapse = ","),
          candidate_id = cfg$candidate_id,
          graph_edge_weight = cfg$graph_edge_weight,
          lambda_l = cfg$lambda_l,
          lambda_e = cfg$lambda_e,
          sigma_obs = cfg$sigma_obs,
          anchor_var_mode = cfg$anchor_var_mode,
          prior_mean_mode = cfg$prior_mean_mode,
          prior_mean_scale = cfg$prior_mean_scale,
          anchor_count_reference_mode = cfg$anchor_count_reference_mode,
          prior_mean_status = pm$status,
          stringsAsFactors = FALSE
        ),
        m
      )
    }
  }
  do.call(rbind, rows)
}

aggregate_cv_candidates <- function(cv_splits) {
  key_cols <- c("candidate_id", "graph_edge_weight", "lambda_l", "lambda_e", "sigma_obs",
                "anchor_var_mode", "prior_mean_mode", "prior_mean_scale", "anchor_count_reference_mode")
  parts <- split(cv_splits, interaction(cv_splits[key_cols], drop = TRUE))
  rows <- lapply(parts, function(x) {
    data.frame(
      x[1L, key_cols, drop = FALSE],
      n_splits = nrow(x),
      cv_mse = median(x$cv_mse, na.rm = TRUE),
      cv_centered_rmse = median(x$cv_centered_rmse, na.rm = TRUE),
      cv_spearman = median(x$cv_spearman, na.rm = TRUE),
      cv_pairwise_rank_loss = median(x$cv_pairwise_rank_loss, na.rm = TRUE),
      cv_estimate_sd_ratio = median(x$cv_estimate_sd_ratio, na.rm = TRUE),
      cv_shape_score = median(x$cv_shape_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

select_by_cv_objectives <- function(cv_candidates) {
  objectives <- c("mse", "centered_rmse", "centered_rmse_amp", "pairwise_rank", "shape_score")
  rows <- list()
  for (obj in objectives) {
    x <- cv_candidates
    x$objective <- obj
    x$objective_value <- switch(
      obj,
      mse = x$cv_mse,
      centered_rmse = x$cv_centered_rmse,
      centered_rmse_amp = x$cv_centered_rmse + 0.25 * abs(log(pmax(x$cv_estimate_sd_ratio, 1e-4))),
      pairwise_rank = ifelse(is.finite(x$cv_pairwise_rank_loss), x$cv_pairwise_rank_loss, 10),
      shape_score = x$cv_shape_score
    )
    x <- x[order(x$objective_value, x$cv_centered_rmse), , drop = FALSE]
    rows[[obj]] <- x[1L, , drop = FALSE]
  }
  do.call(rbind, rows)
}

run_experiment_d <- function(bundle, components, grf, task_info, dirs, a_metrics, c_metrics,
                             cv_splits = 20L, max_cv_candidates = 40L, quick = "auto", force = FALSE) {
  configs <- select_cv_candidates(a_metrics, c_metrics, max_candidates = max_cv_candidates)
  if (quick %in% c("auto", "true")) {
    configs <- head(configs, min(nrow(configs), 20L))
    cv_splits <- min(as.integer(cv_splits), 8L)
  }
  message("[D] CV candidates=", nrow(configs), " splits=", cv_splits)
  cv_splits_tbl <- run_direct_anchor_holdout_cv(bundle, components, configs, task_info, cv_splits = cv_splits)
  write_tsv_safe(cv_splits_tbl, file.path(dirs$tables, "experiment_D_shape_cv_splits.tsv"))
  cv_candidates <- aggregate_cv_candidates(cv_splits_tbl)
  write_tsv_safe(cv_candidates, file.path(dirs$tables, "experiment_D_shape_cv_candidates.tsv"))
  selected <- select_by_cv_objectives(cv_candidates)
  write_tsv_safe(selected, file.path(dirs$tables, "experiment_D_selected_configs.tsv"))
  eval_metrics <- list()
  for (i in seq_len(nrow(selected))) {
    cfg <- selected[i, , drop = FALSE]
    cfg$experiment <- "D"
    cfg$solver <- "matrix_mean"
    cfg$candidate_id <- paste0("D_selected_", cfg$objective, "__", cfg$candidate_id)
    pm <- build_prior_mean(cfg$prior_mean_mode, bundle$local, bundle$global_graph, cfg$prior_mean_scale)
    fit <- fit_cached(bundle$local, bundle$global_graph, components, cfg, grf, task_info, dirs,
                      prior_mean = pm$mean, prior_mean_status = pm$status, force = force)
    mm <- fit$metrics
    mm$cv_objective <- cfg$objective
    mm$cv_objective_value <- cfg$objective_value
    eval_metrics[[i]] <- mm
  }
  eval_out <- do.call(rbind, eval_metrics)
  write_tsv_safe(eval_out, file.path(dirs$tables, "experiment_D_selected_farfield_evaluation.tsv"))
  saveRDS(list(configs = configs, cv_splits = cv_splits_tbl, cv_candidates = cv_candidates,
               selected = selected, evaluation = eval_out),
          file.path(dirs$results, "experiment_D_shape_cv.rds"))
  eval_out
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

make_summary_by_experiment <- function(all_metrics) {
  far <- all_metrics[all_metrics$support_scope == "farfield" & all_metrics$metric_scale == "native", , drop = FALSE]
  if (!nrow(far)) return(data.frame())
  do.call(rbind, lapply(split(far, far$experiment), function(x) {
    best <- x[order(x$shape_score, x$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
    data.frame(
      experiment = best$experiment,
      n_configs = length(unique(x$candidate_id)),
      best_candidate_id = best$candidate_id,
      best_graph_edge_weight = best$graph_edge_weight,
      best_lambda_l = best$lambda_l,
      best_lambda_e = best$lambda_e,
      best_sigma_obs = best$sigma_obs,
      best_anchor_var_mode = best$anchor_var_mode,
      best_prior_mean_mode = best$prior_mean_mode,
      best_prior_mean_scale = best$prior_mean_scale,
      best_centered_rmse = best$centered_rmse,
      best_pearson = best$pearson,
      best_spearman = best$spearman,
      best_estimate_sd_ratio = best$estimate_sd_ratio,
      best_shape_score = best$shape_score,
      any_noncollapsed = any(x$estimate_sd_ratio >= 0.02, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

make_recommended_configs <- function(all_metrics) {
  far <- all_metrics[all_metrics$support_scope == "farfield" & all_metrics$metric_scale == "native", , drop = FALSE]
  if (!nrow(far)) return(data.frame())
  best_shape <- far[order(far$shape_score, far$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
  best_rank <- far[order(-far$spearman, far$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
  best_amp <- far[order(abs(log(pmax(far$estimate_sd_ratio, 1e-4))), far$centered_rmse), , drop = FALSE][1L, , drop = FALSE]
  out <- rbind(best_shape, best_rank, best_amp)
  out$recommendation_role <- c("best_shape_score", "best_spearman", "best_amplitude")
  out
}

make_default_change_impact <- function(all_metrics) {
  far <- all_metrics[all_metrics$support_scope == "farfield" & all_metrics$metric_scale == "native", , drop = FALSE]
  baseline <- far[far$graph_edge_weight == "mutation" & far$lambda_l == 0.2 & far$lambda_e == 1 &
                    far$sigma_obs == 0.05 & far$prior_mean_mode == "zero" &
                    far$anchor_var_mode == "current", , drop = FALSE]
  norm <- far[far$graph_edge_weight == "normalized" & far$prior_mean_mode == "zero" &
                far$anchor_var_mode == "current", , drop = FALSE]
  norm_best <- if (nrow(norm)) norm[order(norm$shape_score), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  if (!nrow(baseline) || !nrow(norm_best)) return(data.frame())
  baseline <- baseline[1L, , drop = FALSE]
  data.frame(
    comparison = "baseline_mutation_to_best_normalized_zero_prior",
    baseline_candidate_id = baseline$candidate_id,
    normalized_candidate_id = norm_best$candidate_id,
    delta_centered_rmse = norm_best$centered_rmse - baseline$centered_rmse,
    delta_pearson = norm_best$pearson - baseline$pearson,
    delta_spearman = norm_best$spearman - baseline$spearman,
    delta_estimate_sd_ratio = norm_best$estimate_sd_ratio - baseline$estimate_sd_ratio,
    baseline_shape_score = baseline$shape_score,
    normalized_shape_score = norm_best$shape_score,
    stringsAsFactors = FALSE
  )
}

local_diagnostics_table <- function(bundle) {
  local <- bundle$local
  data.frame(
    convergence = local$diagnostics$convergence %||% NA,
    message = local$diagnostics$message %||% NA_character_,
    gradient_norm = local$diagnostics$gradient_norm %||% NA_real_,
    covariance_status = local$diagnostics$covariance_status %||% NA_character_,
    covariance_fallback = local$diagnostics$covariance_fallback %||% NA,
    fitness_sd_source = local$diagnostics$fitness_sd_source %||% NA_character_,
    n_local_nodes = nrow(local$summary),
    n_global_nodes = length(bundle$global_graph$labels),
    n_anchors = sum(as.character(local$summary$support_tier) == "directly_informed" & is.finite(local$summary$fitness_mean)),
    stringsAsFactors = FALSE
  )
}

fmt_metric <- function(x) {
  if (!is.finite(x)) "NA" else format(round(x, 4), nsmall = 4)
}

write_report <- function(dirs, ctx, bundle, all_metrics, summary_by_exp, recommended, default_impact, quick_mode) {
  far <- all_metrics[all_metrics$support_scope == "farfield" & all_metrics$metric_scale == "native", , drop = FALSE]
  direct <- all_metrics[all_metrics$support_scope == "direct" & all_metrics$metric_scale == "native", , drop = FALSE]
  baseline <- far[far$graph_edge_weight == "mutation" & far$lambda_l == 0.2 & far$lambda_e == 1 &
                    far$sigma_obs == 0.05 & far$anchor_var_mode == "current" &
                    far$prior_mean_mode == "zero", , drop = FALSE]
  probe_ref <- far[far$candidate_id == "baseline_probe_reference_mutation_ll0p2_le1_so0p02", , drop = FALSE]
  a <- far[far$experiment == "A", , drop = FALSE]
  b <- far[far$experiment == "B", , drop = FALSE]
  c <- far[far$experiment == "C", , drop = FALSE]
  d <- far[far$experiment == "D", , drop = FALSE]
  ld <- local_diagnostics_table(bundle)
  best <- function(x) if (nrow(x)) x[order(x$shape_score, x$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_a <- best(a)
  a_normunit <- a[a$graph_edge_weight %in% c("normalized", "unit"), , drop = FALSE]
  best_a_normunit <- if (nrow(a_normunit)) a_normunit[order(a_normunit$centered_rmse), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  best_b <- best(b)
  best_c <- best(c)
  best_d <- best(d)
  direct_c_best <- if (nrow(best_c) && nrow(direct)) {
    direct[direct$candidate_id == best_c$candidate_id, , drop = FALSE][1L, , drop = FALSE]
  } else data.frame()
  lines <- c(
    "# Farfield Shape ABCD Report",
    "",
    "## Data source",
    paste0("- source_input_dir: `", ctx$shared_input_dir, "`"),
    paste0("- source_probe_dir: `", ifelse(is.null(ctx$source_probe_dir), "none", ctx$source_probe_dir), "`"),
    "- condition: simulation_id=1, minobs=5, input_policy=full",
    paste0("- solver: Matrix posterior-mean solver in `benchmark/scr/run_farfield_shape_abcd.R`; posterior SD uses diagonal-precision approximation, so shape/rank conclusions focus on posterior mean metrics."),
    paste0("- quick_mode: ", quick_mode),
    "",
    "## Local fit diagnostics",
    paste0("- convergence: ", ld$convergence, "; gradient_norm: ", fmt_metric(ld$gradient_norm)),
    paste0("- covariance_status: ", ld$covariance_status, "; covariance_fallback: ", ld$covariance_fallback),
    paste0("- fitness_sd_source: ", ld$fitness_sd_source),
    paste0("- n_local_nodes: ", ld$n_local_nodes, "; n_global_nodes: ", ld$n_global_nodes, "; n_anchors: ", ld$n_anchors),
    "",
    "## Baseline reproduction",
    paste0("- requested baseline mutation/lambda_l=0.2/lambda_e=1/sigma_obs=0.05 farfield native: centered_rmse=", fmt_metric(baseline$centered_rmse[1]), ", pearson=", fmt_metric(baseline$pearson[1]), ", spearman=", fmt_metric(baseline$spearman[1]), ", estimate_sd_ratio=", fmt_metric(baseline$estimate_sd_ratio[1]), "."),
    paste0("- prior probe reference mutation/lambda_l=0.2/lambda_e=1/sigma_obs=0.02 farfield native: centered_rmse=", fmt_metric(probe_ref$centered_rmse[1]), ", pearson=", fmt_metric(probe_ref$pearson[1]), ", spearman=", fmt_metric(probe_ref$spearman[1]), ", estimate_sd_ratio=", fmt_metric(probe_ref$estimate_sd_ratio[1]), "."),
    "",
    "## Experiment A: global scale grid",
    paste0("- best farfield/native shape_score config: ", best_a$candidate_id[1], " with graph_edge_weight=", best_a$graph_edge_weight[1], ", lambda_l=", best_a$lambda_l[1], ", lambda_e=", best_a$lambda_e[1], ", sigma_obs=", best_a$sigma_obs[1], "."),
    paste0("- metrics: centered_rmse=", fmt_metric(best_a$centered_rmse[1]), ", pearson=", fmt_metric(best_a$pearson[1]), ", spearman=", fmt_metric(best_a$spearman[1]), ", estimate_sd_ratio=", fmt_metric(best_a$estimate_sd_ratio[1]), "."),
    paste0("- best normalized/unit centered_rmse config: ", best_a_normunit$candidate_id[1], " with centered_rmse=", fmt_metric(best_a_normunit$centered_rmse[1]), ", pearson=", fmt_metric(best_a_normunit$pearson[1]), ", spearman=", fmt_metric(best_a_normunit$spearman[1]), ", estimate_sd_ratio=", fmt_metric(best_a_normunit$estimate_sd_ratio[1]), "."),
    paste0("- amplitude collapse across A: ", all(a$estimate_sd_ratio < 0.02, na.rm = TRUE), "."),
    "",
    "## Experiment B: anchor variance isolation",
    paste0("- best anchor variance config: ", best_b$candidate_id[1], " with anchor_var_mode=", best_b$anchor_var_mode[1], "."),
    paste0("- metrics: centered_rmse=", fmt_metric(best_b$centered_rmse[1]), ", pearson=", fmt_metric(best_b$pearson[1]), ", spearman=", fmt_metric(best_b$spearman[1]), ", estimate_sd_ratio=", fmt_metric(best_b$estimate_sd_ratio[1]), "."),
    "- constant anchor variance is interpreted by comparing `experiment_B_anchor_variance_summary.tsv`; count_inflated is retained as a candidate rather than promoted to default.",
    "",
    "## Experiment C: prior mean / edge-slope prototype",
    paste0("- best prior mean config: ", best_c$candidate_id[1], " with prior_mean_mode=", best_c$prior_mean_mode[1], ", prior_mean_scale=", best_c$prior_mean_scale[1], "."),
    paste0("- farfield metrics: centered_rmse=", fmt_metric(best_c$centered_rmse[1]), ", pearson=", fmt_metric(best_c$pearson[1]), ", spearman=", fmt_metric(best_c$spearman[1]), ", estimate_sd_ratio=", fmt_metric(best_c$estimate_sd_ratio[1]), "."),
    paste0("- matching direct/native centered_rmse for that config: ", fmt_metric(direct_c_best$centered_rmse[1]), "."),
    "",
    "## Experiment D: shape-aware CV",
    paste0("- best selected objective by oracle farfield shape among CV selections: ", best_d$cv_objective[1], " using ", best_d$candidate_id[1], "."),
    paste0("- selected farfield metrics: centered_rmse=", fmt_metric(best_d$centered_rmse[1]), ", pearson=", fmt_metric(best_d$pearson[1]), ", spearman=", fmt_metric(best_d$spearman[1]), ", estimate_sd_ratio=", fmt_metric(best_d$estimate_sd_ratio[1]), "."),
    "- CV objectives are trained only on direct-anchor holdout labels; GRF truth is used afterward for farfield evaluation.",
    "",
    "## Recommendation",
    "- Current evidence supports changing benchmark/probe default `graph_edge_weight` from `mutation` to `normalized` for scale-sane GRF stress testing, but normalized is not a standalone fix: the farfield shape score with amplitude penalty still favors mutation/constant-anchor variants in this single simulation.",
    "- Unit edge weights remain useful as a synthetic stress-test but should not replace normalized as the benchmark default unless it consistently wins across more simulations.",
    "- Lowering `lambda_e` alone is not sufficient if the selected configuration still has low `estimate_sd_ratio`; amplitude diagnostics must be part of calibration ranking.",
    "- `anchor_count_reference=minobs` should stay a candidate, not a full-input default, because this probe does not show a reliable improvement.",
    "- The next core-code direction is an edge-gradient or prior-mean mechanism, but local convergence/covariance should also be fixed because the current local fit is `untrusted_nonconverged` and drives fallback anchor variance.",
    "",
    "## Output tables",
    "- `tables/all_experiments_long.tsv`",
    "- `tables/summary_by_experiment.tsv`",
    "- `tables/recommended_configs.tsv`",
    "- `tables/default_change_impact.tsv`"
  )
  writeLines(lines, file.path(dirs$root, "farfield_shape_abcd_report.md"))
}

make_dirs <- function(output_dir) {
  if (!grepl("^/", output_dir)) output_dir <- file.path(repo_guess, output_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!path_under_benchmark(output_dir)) {
    stop("Refusing to write benchmark outputs outside benchmark/: ", output_dir, call. = FALSE)
  }
  dirs <- list(
    root = output_dir,
    cache = file.path(output_dir, "cache"),
    fits = file.path(output_dir, "fits"),
    tables = file.path(output_dir, "tables"),
    results = file.path(output_dir, "results")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

main_abcd <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  mode <- match.arg(mode, c("prepare", "experiment-a", "experiment-b", "experiment-c", "experiment-d", "summarize", "all"))
  source_input_dir <- as.character(arg_value(args, "source_input_dir", "benchmark/results/farfield_shape_probe_default"))
  output_dir <- as.character(arg_value(args, "output_dir", "benchmark/results/farfield_shape_probe_abcd"))
  simulation_id <- arg_integer(args, "simulation_id", 1L)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  quick <- as_bool_or_auto(arg_value(args, "quick", "auto"))
  cv_splits <- arg_integer(args, "cv_splits", 20L)
  max_cv_candidates <- arg_integer(args, "max_cv_candidates", 40L)
  pkgload::load_all(repo_guess, quiet = TRUE)
  dirs <- make_dirs(output_dir)
  ctx <- resolve_source_context(source_input_dir, simulation_id, minobs, input_policy)
  bundle <- prepare_abcd_bundle(ctx, dirs, simulation_id, minobs, input_policy, force = force)
  grf <- readRDS(ctx$input_table$grf_rds[[1L]])
  task_info <- list(
    simulation_id = simulation_id,
    minobs = minobs,
    input_policy = input_policy,
    lambda = as.numeric(ctx$input_table$lambda[[1L]]),
    dt = as.numeric(ctx$input_table$time_delta[[1L]]),
    beta = as.numeric(ctx$input_table$sim_pm[[1L]])
  )
  write_tsv_safe(local_diagnostics_table(bundle), file.path(dirs$tables, "local_fit_diagnostics.tsv"))
  components <- prepare_solver_cache(bundle$global_graph, dirs, force = force)
  saveRDS(list(context = ctx, task_info = task_info, local_diagnostics = local_diagnostics_table(bundle)),
          file.path(dirs$results, "prepare_context.rds"))
  if (identical(mode, "prepare")) return(invisible(dirs$root))

  a_path <- file.path(dirs$tables, "experiment_A_global_scale_grid.tsv")
  b_path <- file.path(dirs$tables, "experiment_B_anchor_variance.tsv")
  c_path <- file.path(dirs$tables, "experiment_C_prior_mean_edge_slope.tsv")
  d_path <- file.path(dirs$tables, "experiment_D_selected_farfield_evaluation.tsv")

  a_metrics <- if (mode %in% c("all", "experiment-a")) {
    run_experiment_a(bundle, components, grf, task_info, dirs, force = force)
  } else read_tsv_safe(a_path)
  if (identical(mode, "experiment-a")) return(invisible(a_metrics))

  b_metrics <- if (mode %in% c("all", "experiment-b")) {
    if (!nrow(a_metrics)) stop("Experiment B requires experiment A metrics.", call. = FALSE)
    run_experiment_b(bundle, components, grf, task_info, dirs, a_metrics, force = force)
  } else read_tsv_safe(b_path)
  if (identical(mode, "experiment-b")) return(invisible(b_metrics))

  c_metrics <- if (mode %in% c("all", "experiment-c")) {
    if (!nrow(a_metrics) || !nrow(b_metrics)) stop("Experiment C requires A and B metrics.", call. = FALSE)
    run_experiment_c(bundle, components, grf, task_info, dirs, a_metrics, b_metrics, force = force)
  } else read_tsv_safe(c_path)
  if (identical(mode, "experiment-c")) return(invisible(c_metrics))

  d_metrics <- if (mode %in% c("all", "experiment-d")) {
    if (!nrow(a_metrics) || !nrow(c_metrics)) stop("Experiment D requires A and C metrics.", call. = FALSE)
    run_experiment_d(bundle, components, grf, task_info, dirs, a_metrics, c_metrics,
                     cv_splits = cv_splits, max_cv_candidates = max_cv_candidates,
                     quick = quick, force = force)
  } else read_tsv_safe(d_path)
  if (identical(mode, "experiment-d")) return(invisible(d_metrics))

  all_metrics <- bind_rows_fill(list(a_metrics, b_metrics, c_metrics, d_metrics))
  write_tsv_safe(all_metrics, file.path(dirs$tables, "all_experiments_long.tsv"))
  summary_by_exp <- make_summary_by_experiment(all_metrics)
  recommended <- make_recommended_configs(all_metrics)
  default_impact <- make_default_change_impact(all_metrics)
  write_tsv_safe(summary_by_exp, file.path(dirs$tables, "summary_by_experiment.tsv"))
  write_tsv_safe(recommended, file.path(dirs$tables, "recommended_configs.tsv"))
  write_tsv_safe(default_impact, file.path(dirs$tables, "default_change_impact.tsv"))
  saveRDS(list(context = ctx, task_info = task_info, metrics = all_metrics,
               summary_by_experiment = summary_by_exp, recommended = recommended,
               default_change_impact = default_impact),
          file.path(dirs$results, "farfield_shape_abcd_all_results.rds"))
  write_report(dirs, ctx, bundle, all_metrics, summary_by_exp, recommended, default_impact, quick)
  message("Wrote ABCD outputs under: ", dirs$root)
  invisible(dirs$root)
}

if (sys.nframe() == 0L) {
  main_abcd()
}
