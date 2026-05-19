#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Probe alfak2 far-field landscape-shape configurations on prepared GRF inputs.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/probe_farfield_shape_configs.R [options]\n\n",
    "Options:\n",
    "  --source-input-dir=<dir>     Prepared shared_inputs directory. Defaults to latest local/grf shared input.\n",
    "  --output-dir=<dir>           Output directory. Defaults to benchmark/results/farfield_shape_probe_<timestamp>.\n",
    "  --n-sim=1                    Number of simulation ids to run.\n",
    "  --minobs=5                   minobs row from input_table.tsv.\n",
    "  --input-policy=full          full|minobs_matched|soft_minobs.\n",
    "  --force=false                Refit even if cached files exist.\n",
    "  --local-shell-depth=1        Local shell depth.\n",
    "  --global-extra-shell=1       Global extra shell depth.\n",
    "  --max-nodes=150000           Graph max nodes.\n",
    "  --eval-max=500               nlminb eval.max for local fit.\n",
    "  --iter-max=500               nlminb iter.max for local fit.\n",
    "  --retry-max=2000             nlminb retry eval/iter max for local fit.\n",
    "  --candidate-set=default      default|minimal.\n",
    "\n",
    "Outputs:\n",
    "  tables/candidate_metrics.tsv\n",
    "  tables/farfield_ranking.tsv\n",
    "  tables/fit_manifest.tsv\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) stop("Unexpected positional argument: ", arg, call. = FALSE)
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      val <- "true"
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- val
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  if (!is.null(args[[name]])) args[[name]] else default
}

arg_integer <- function(args, name, default) {
  value <- suppressWarnings(as.integer(arg_value(args, name, default)))
  if (!is.finite(value)) stop("Invalid integer option --", gsub("_", "-", name), call. = FALSE)
  value
}

arg_numeric <- function(args, name, default) {
  value <- suppressWarnings(as.numeric(arg_value(args, name, default)))
  if (!is.finite(value)) stop("Invalid numeric option --", gsub("_", "-", name), call. = FALSE)
  value
}

arg_logical <- function(args, name, default = FALSE) {
  value <- tolower(as.character(arg_value(args, name, if (default) "true" else "false")))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Invalid logical option --", gsub("_", "-", name), call. = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

latest_shared_input_dir <- function(root = "benchmark/results") {
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[basename(dirs) == "shared_inputs" &
                 file.exists(file.path(dirs, "tables", "input_table.tsv"))]
  if (!length(dirs)) {
    stop("Could not find a shared_inputs directory under ", root, call. = FALSE)
  }
  dirs <- dirs[order(file.info(dirs)$mtime, decreasing = TRUE)]
  for (dir in dirs) {
    input_tbl <- try(
      utils::read.delim(file.path(dir, "tables", "input_table.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE),
      silent = TRUE
    )
    if (inherits(input_tbl, "try-error") || !nrow(input_tbl)) next
    has_paths <- all(c("input_rds", "grf_rds") %in% names(input_tbl))
    if (!has_paths) next
    ok <- file.exists(input_tbl$input_rds[1L]) && file.exists(input_tbl$grf_rds[1L])
    if (isTRUE(ok)) return(dir)
  }
  stop("Found shared_inputs directories, but none have readable input_rds/grf_rds paths.", call. = FALSE)
}

parse_karyotype_ids <- function(ids) {
  ids <- as.character(ids)
  pieces <- strsplit(ids, ".", fixed = TRUE)
  lens <- lengths(pieces)
  if (!length(ids) || length(unique(lens)) != 1L) {
    stop("Karyotype ids have inconsistent dimensions.", call. = FALSE)
  }
  mat <- matrix(as.integer(unlist(pieces, use.names = FALSE)), nrow = length(ids), byrow = TRUE)
  if (anyNA(mat)) stop("Karyotype ids contain non-integer fields.", call. = FALSE)
  rownames(mat) <- ids
  mat
}

compute_grf_truth <- function(karyotypes, centroids, lambda) {
  karyotypes <- unique(as.character(karyotypes))
  karyotypes <- karyotypes[nzchar(karyotypes)]
  if (!length(karyotypes)) return(stats::setNames(numeric(0), character(0)))
  k_mat <- parse_karyotype_ids(karyotypes)
  if (ncol(k_mat) != ncol(centroids)) {
    stop("Karyotype dimension does not match GRF centroids.", call. = FALSE)
  }
  out <- vapply(seq_len(nrow(k_mat)), function(i) {
    diffs <- sweep(centroids, 2L, as.numeric(k_mat[i, ]), FUN = "-")
    distances <- sqrt(rowSums(diffs^2))
    sum(sin(distances / lambda)) / (pi * sqrt(nrow(centroids)))
  }, numeric(1))
  stats::setNames(out, rownames(k_mat))
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

centered_rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  xc <- x[ok] - mean(x[ok])
  yc <- y[ok] - mean(y[ok])
  sqrt(mean((xc - yc)^2))
}

metric_row <- function(est, truth) {
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) {
    return(c(
      n = 0, mae = NA_real_, rmse = NA_real_, centered_rmse = NA_real_,
      pearson = NA_real_, spearman = NA_real_, estimate_sd = NA_real_,
      truth_sd = NA_real_, estimate_sd_ratio = NA_real_
    ))
  }
  est_ok <- est[ok]
  truth_ok <- truth[ok]
  truth_sd <- stats::sd(truth_ok)
  est_sd <- stats::sd(est_ok)
  c(
    n = sum(ok),
    mae = mean(abs(est_ok - truth_ok)),
    rmse = sqrt(mean((est_ok - truth_ok)^2)),
    centered_rmse = centered_rmse(est_ok, truth_ok),
    pearson = safe_cor(est_ok, truth_ok, "pearson"),
    spearman = safe_cor(est_ok, truth_ok, "spearman"),
    estimate_sd = est_sd,
    truth_sd = truth_sd,
    estimate_sd_ratio = if (is.finite(truth_sd) && truth_sd > 0) est_sd / truth_sd else NA_real_
  )
}

drop_diploid_counts <- function(counts) {
  if (!nrow(counts)) return(counts)
  k_dim <- length(strsplit(rownames(counts)[[1]], ".", fixed = TRUE)[[1]])
  diploid <- paste(rep(2L, k_dim), collapse = ".")
  counts[rownames(counts) != diploid, , drop = FALSE]
}

prepare_counts <- function(yi, minobs, input_policy) {
  counts <- as.matrix(yi$x)
  storage.mode(counts) <- "integer"
  counts <- drop_diploid_counts(counts)
  if (identical(input_policy, "minobs_matched")) {
    counts <- counts[rowSums(counts, na.rm = TRUE) >= minobs, , drop = FALSE]
  } else if (identical(input_policy, "soft_minobs")) {
    row_totals <- rowSums(counts, na.rm = TRUE)
    weights <- pmin(1, pmax(0, row_totals) / minobs)
    weights[!is.finite(weights)] <- 0
    weights <- cbind(t0 = weights, t1 = weights)
    rownames(weights) <- rownames(counts)
    attr(counts, "observation_weights") <- weights
    attr(counts, "soft_minobs") <- list(minobs = minobs, rule = "row_total_over_minobs")
  } else if (!identical(input_policy, "full")) {
    stop("Unsupported input policy: ", input_policy, call. = FALSE)
  }
  if (!nrow(counts)) stop("No counts remain after input-policy filtering.", call. = FALSE)
  counts
}

candidate_table <- function(candidate_set = "default") {
  if (identical(candidate_set, "minimal")) {
    return(data.frame(
      candidate_id = c("baseline_mutation", "normalized_low_curv", "normalized_low_curv_countref"),
      graph_edge_weight = c("mutation", "normalized", "normalized"),
      lambda_l = c(0.2, 1.0, 1.0),
      lambda_e = c(1.0, 0.01, 0.01),
      sigma_obs = c(0.02, 0.02, 0.02),
      anchor_count_reference = c(NA_real_, NA_real_, NA_real_),
      anchor_count_reference_mode = c("none", "none", "minobs"),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    candidate_id = c(
      "baseline_mutation",
      "normalized_same_smoothness",
      "normalized_low_curv_1e2",
      "normalized_low_curv_1e3",
      "unit_low_curv_1e2",
      "mutation_low_curv_1e4",
      "normalized_low_curv_countref"
    ),
    graph_edge_weight = c("mutation", "normalized", "normalized", "normalized", "unit", "mutation", "normalized"),
    lambda_l = c(0.2, 0.2, 1.0, 1.0, 1.0, 0.2, 1.0),
    lambda_e = c(1.0, 1.0, 0.01, 0.001, 0.01, 0.0001, 0.01),
    sigma_obs = c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
    anchor_count_reference = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    anchor_count_reference_mode = c("none", "none", "none", "none", "none", "none", "minobs"),
    stringsAsFactors = FALSE
  )
}

add_legacy_scale_if_possible <- function(summary, graph, dt, beta, n0 = 1e5, nb = 1e7) {
  add_scale <- get("add_alfakR_scale_to_summary", envir = asNamespace("alfak2"))
  weights <- as.numeric(as.character(summary$support_tier) == "directly_informed")
  if (!any(weights > 0 & is.finite(summary$fitness_mean))) return(summary)
  add_scale(summary, graph$karyotypes, weights, dt = dt, beta = beta, n0 = n0, nb = nb, correct_efflux = TRUE)
}

score_summary <- function(summary, graph, grf, lambda, task_info, candidate) {
  summary <- add_legacy_scale_if_possible(summary, graph, task_info$dt, task_info$beta)
  truth_map <- compute_grf_truth(summary$karyotype, grf$centroids, lambda)
  truth <- as.numeric(truth_map[as.character(summary$karyotype)])
  scopes <- list(
    direct = summary$support_tier == "directly_informed",
    local_borrowed = summary$support_tier == "local_borrowed",
    farfield = summary$support_tier %in% c("weakly_supported", "graph_borrowed", "prior_dominated"),
    whole = rep(TRUE, nrow(summary))
  )
  scales <- c(native = "fitness_mean")
  if ("fitness_mean_alfakR_scale" %in% names(summary)) {
    scales <- c(scales, alfakR_scale = "fitness_mean_alfakR_scale")
  }
  rows <- list()
  idx <- 0L
  for (scale_name in names(scales)) {
    est <- as.numeric(summary[[scales[[scale_name]]]])
    for (scope_name in names(scopes)) {
      keep <- scopes[[scope_name]]
      m <- metric_row(est[keep], truth[keep])
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        simulation_id = task_info$simulation_id,
        minobs = task_info$minobs,
        input_policy = task_info$input_policy,
        candidate_id = candidate$candidate_id,
        graph_edge_weight = candidate$graph_edge_weight,
        lambda_l = candidate$lambda_l,
        lambda_e = candidate$lambda_e,
        sigma_obs = candidate$sigma_obs,
        anchor_count_reference_mode = candidate$anchor_count_reference_mode,
        metric_scale = scale_name,
        support_scope = scope_name,
        n = as.integer(m[["n"]]),
        mae = m[["mae"]],
        rmse = m[["rmse"]],
        centered_rmse = m[["centered_rmse"]],
        pearson = m[["pearson"]],
        spearman = m[["spearman"]],
        estimate_sd = m[["estimate_sd"]],
        truth_sd = m[["truth_sd"]],
        estimate_sd_ratio = m[["estimate_sd_ratio"]],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }

  repo_dir <- normalizePath(".", mustWork = TRUE)
  source_input_dir <- arg_value(args, "source_input_dir", NULL)
  if (is.null(source_input_dir)) source_input_dir <- latest_shared_input_dir()
  source_input_dir <- normalizePath(source_input_dir, mustWork = TRUE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- arg_value(args, "output_dir", file.path("benchmark", "results", paste0("farfield_shape_probe_", timestamp)))
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  table_dir <- file.path(output_dir, "tables")
  fit_dir <- file.path(output_dir, "fits")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

  n_sim <- arg_integer(args, "n_sim", 1L)
  minobs <- arg_integer(args, "minobs", 5L)
  input_policy <- as.character(arg_value(args, "input_policy", "full"))
  force <- arg_logical(args, "force", FALSE)
  local_shell_depth <- arg_integer(args, "local_shell_depth", 1L)
  global_extra_shell <- arg_integer(args, "global_extra_shell", 1L)
  max_nodes <- arg_integer(args, "max_nodes", 150000L)
  eval_max <- arg_integer(args, "eval_max", 500L)
  iter_max <- arg_integer(args, "iter_max", 500L)
  retry_max <- arg_integer(args, "retry_max", 2000L)
  candidate_set <- as.character(arg_value(args, "candidate_set", "default"))

  pkgload::load_all(repo_dir, quiet = TRUE)
  prep_input_depth <- get("prepare_counts_for_input_depth", envir = asNamespace("alfak2"))
  resolve_obs <- get("resolve_fit_observation_controls", envir = asNamespace("alfak2"))

  input_path <- file.path(source_input_dir, "tables", "input_table.tsv")
  input_tbl <- utils::read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  input_tbl <- input_tbl[input_tbl$minobs == minobs, , drop = FALSE]
  sim_ids <- unique(input_tbl$simulation_id)
  sim_ids <- sim_ids[seq_len(min(length(sim_ids), n_sim))]
  input_tbl <- input_tbl[input_tbl$simulation_id %in% sim_ids, , drop = FALSE]
  if (!nrow(input_tbl)) stop("No matching input rows for minobs=", minobs, call. = FALSE)

  candidates <- candidate_table(candidate_set)
  write_tsv(candidates, file.path(table_dir, "candidate_table.tsv"))

  metric_rows <- list()
  manifest_rows <- list()
  metric_i <- 0L
  manifest_i <- 0L

  for (task_idx in seq_len(nrow(input_tbl))) {
    task <- input_tbl[task_idx, , drop = FALSE]
    message("[task ", task_idx, "/", nrow(input_tbl), "] simulation=", task$simulation_id,
            " minobs=", minobs, " input_policy=", input_policy)
    yi <- readRDS(task$input_rds)
    counts <- prepare_counts(yi, minobs = minobs, input_policy = input_policy)
    selected_times <- suppressWarnings(as.numeric(colnames(counts)))
    dt <- if (length(selected_times) == 2L && all(is.finite(selected_times))) diff(selected_times) else NA_real_
    if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
    beta <- as.numeric(task$sim_pm)
    k_mat <- parse_karyotype_ids(rownames(counts))
    max_cn <- max(k_mat, na.rm = TRUE) + local_shell_depth + global_extra_shell
    grf <- readRDS(task$grf_rds)

    task_id <- paste0("sim", task$simulation_id, "_minobs", minobs, "_", input_policy)
    task_fit_dir <- file.path(fit_dir, task_id)
    dir.create(task_fit_dir, recursive = TRUE, showWarnings = FALSE)
    local_bundle_path <- file.path(task_fit_dir, "local_bundle.rds")
    if (!force && file.exists(local_bundle_path)) {
      bundle <- readRDS(local_bundle_path)
      data <- bundle$data
      local <- bundle$local
      global_graph <- bundle$global_graph
    } else {
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
        control = list(eval.max = eval_max, iter.max = iter_max),
        retry_control = list(eval.max = retry_max, iter.max = retry_max)
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
      saveRDS(bundle, local_bundle_path)
    }

    task_info <- list(
      simulation_id = as.integer(task$simulation_id),
      minobs = minobs,
      input_policy = input_policy,
      dt = as.numeric(dt),
      beta = as.numeric(beta)
    )

    for (cand_idx in seq_len(nrow(candidates))) {
      cand <- candidates[cand_idx, , drop = FALSE]
      cand_path <- file.path(task_fit_dir, paste0(cand$candidate_id, ".rds"))
      started <- Sys.time()
      status <- "ok"
      error_message <- NA_character_
      global <- NULL
      if (!force && file.exists(cand_path)) {
        global <- readRDS(cand_path)
      } else {
        anchor_ref <- if (identical(as.character(cand$anchor_count_reference_mode), "minobs")) {
          as.numeric(minobs)
        } else if (is.finite(cand$anchor_count_reference)) {
          as.numeric(cand$anchor_count_reference)
        } else {
          NULL
        }
        global <- tryCatch(
          alfak2::fit_graph_posterior(
            local,
            global_graph,
            lambda_l_grid = as.numeric(cand$lambda_l),
            lambda_e_grid = as.numeric(cand$lambda_e),
            sigma_obs_grid = as.numeric(cand$sigma_obs),
            graph_edge_weight = as.character(cand$graph_edge_weight),
            anchor_support_tiers = "directly_informed",
            anchor_count_reference = anchor_ref,
            anchor_count_power = 1,
            anchor_min_effective_count = 0
          ),
          error = function(e) {
            status <<- "error"
            error_message <<- conditionMessage(e)
            NULL
          }
        )
        if (!is.null(global)) saveRDS(global, cand_path)
      }
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      manifest_i <- manifest_i + 1L
      manifest_rows[[manifest_i]] <- data.frame(
        simulation_id = as.integer(task$simulation_id),
        minobs = minobs,
        input_policy = input_policy,
        candidate_id = cand$candidate_id,
        status = status,
        error_message = error_message,
        elapsed_sec = elapsed,
        fit_path = cand_path,
        local_covariance_status = local$diagnostics$covariance_status,
        local_gradient_norm = local$diagnostics$gradient_norm,
        n_local_nodes = nrow(local$summary),
        n_global_nodes = length(global_graph$labels),
        stringsAsFactors = FALSE
      )
      if (!is.null(global)) {
        scored <- score_summary(
          global$summary,
          global$graph,
          grf,
          lambda = as.numeric(task$lambda),
          task_info = task_info,
          candidate = cand
        )
        scored$local_covariance_status <- local$diagnostics$covariance_status
        scored$global_cv_status <- global$hyperparameters$cv_status
        scored$global_cv_score <- global$hyperparameters$cv_score
        metric_i <- metric_i + 1L
        metric_rows[[metric_i]] <- scored
      }
      message("  ", cand$candidate_id, ": ", status, " (", round(elapsed, 1), " sec)")
    }
  }

  metrics <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame()
  manifest <- if (length(manifest_rows)) do.call(rbind, manifest_rows) else data.frame()
  write_tsv(metrics, file.path(table_dir, "candidate_metrics.tsv"))
  write_tsv(manifest, file.path(table_dir, "fit_manifest.tsv"))

  if (nrow(metrics)) {
    far <- metrics[metrics$support_scope == "farfield" & metrics$metric_scale == "native", , drop = FALSE]
    split_key <- interaction(far$candidate_id, far$graph_edge_weight, far$lambda_l, far$lambda_e,
                             far$sigma_obs, far$anchor_count_reference_mode, drop = TRUE)
    ranking <- do.call(rbind, lapply(split(far, split_key), function(x) {
      data.frame(
        candidate_id = x$candidate_id[[1L]],
        graph_edge_weight = x$graph_edge_weight[[1L]],
        lambda_l = x$lambda_l[[1L]],
        lambda_e = x$lambda_e[[1L]],
        sigma_obs = x$sigma_obs[[1L]],
        anchor_count_reference_mode = x$anchor_count_reference_mode[[1L]],
        n_tasks = nrow(x),
        farfield_centered_rmse_median = median(x$centered_rmse, na.rm = TRUE),
        farfield_pearson_median = median(x$pearson, na.rm = TRUE),
        farfield_spearman_median = median(x$spearman, na.rm = TRUE),
        farfield_estimate_sd_ratio_median = median(x$estimate_sd_ratio, na.rm = TRUE),
        rank_score = median(x$pearson + x$spearman - x$centered_rmse, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    ranking <- ranking[order(-ranking$rank_score, ranking$farfield_centered_rmse_median), , drop = FALSE]
    ranking$rank <- seq_len(nrow(ranking))
    write_tsv(ranking, file.path(table_dir, "farfield_ranking.tsv"))
  }

  saveRDS(
    list(
      source_input_dir = source_input_dir,
      output_dir = output_dir,
      input_table = input_tbl,
      candidates = candidates,
      metrics = metrics,
      manifest = manifest
    ),
    file.path(output_dir, "farfield_shape_probe.rds")
  )
  message("Wrote probe outputs under: ", output_dir)
  invisible(output_dir)
}

if (sys.nframe() == 0L) {
  main()
}
