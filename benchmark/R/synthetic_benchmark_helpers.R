source_benchmark_real_helpers <- function(repo_dir) {
  source(file.path(repo_dir, "benchmark", "R", "real_sample_helpers.R"))
}

source_reference_logistic_model <- function(repo_dir) {
  source(file.path(repo_dir, "benchmark", "R", "reference_logistic_model.R"))
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

top_k_overlap <- function(truth, pred, k = 10) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  truth <- truth[ok]
  pred <- pred[ok]
  k <- min(k, length(truth))
  top_truth <- order(truth, decreasing = TRUE)[seq_len(k)]
  top_pred <- order(pred, decreasing = TRUE)[seq_len(k)]
  length(intersect(top_truth, top_pred)) / k
}

truth_for_summary <- function(summary_df, landscape) {
  idx <- match(as.character(summary_df$karyotype), as.character(landscape$labels))
  out <- rep(NA_real_, nrow(summary_df))
  out[!is.na(idx)] <- landscape$fitness[idx[!is.na(idx)]]
  out
}

benchmark_token <- function(x) {
  x <- if (is.numeric(x)) format(x, scientific = TRUE, digits = 12) else as.character(x)
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

reference_standard_cache_dir <- function(repo_dir,
                                         subdir = "reference_standard_inputs") {
  path <- file.path(repo_dir, "benchmark", "data", subdir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

reference_standard_cache_file <- function(repo_dir,
                                          family,
                                          replicate,
                                          pm,
                                          n_chr,
                                          min_cn,
                                          max_cn,
                                          K,
                                          N0,
                                          step_size,
                                          census_size,
                                          landscape_seed,
                                          subdir = "reference_standard_inputs") {
  name <- paste(
    "standard",
    benchmark_token(family),
    paste0("rep", replicate),
    paste0("pm", benchmark_token(pm)),
    paste0("chr", n_chr),
    paste0("cn", min_cn, "-", max_cn),
    paste0("seed", landscape_seed),
    paste0("K", benchmark_token(K)),
    paste0("N0", benchmark_token(N0)),
    paste0("h", benchmark_token(step_size)),
    paste0("census", benchmark_token(census_size)),
    "itime1",
    sep = "__"
  )
  file.path(reference_standard_cache_dir(repo_dir, subdir), paste0(name, ".rds"))
}

load_or_build_reference_standard_data <- function(repo_dir,
                                                  family,
                                                  replicate,
                                                  pm,
                                                  n_chr = 6,
                                                  min_cn = 1,
                                                  max_cn = 8,
                                                  K = 1e9,
                                                  N0 = 100,
                                                  step_size = 0.2,
                                                  census_size = K,
                                                  landscape_seed = 1,
                                                  force_rebuild = FALSE,
                                                  cache_subdir = "reference_standard_inputs",
                                                  verbose = TRUE) {
  cache_file <- reference_standard_cache_file(
    repo_dir = repo_dir,
    family = family,
    replicate = replicate,
    pm = pm,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    K = K,
    N0 = N0,
    step_size = step_size,
    census_size = census_size,
    landscape_seed = landscape_seed,
    subdir = cache_subdir
  )

  if (file.exists(cache_file) && !force_rebuild) {
    if (isTRUE(verbose)) message("Loading cached reference standard data: ", cache_file)
    reference <- readRDS(cache_file)
    reference$cache <- list(file = cache_file, status = "loaded")
    return(reference)
  }

  if (isTRUE(verbose)) message("Building reference standard data: ", cache_file)
  landscape <- generate_reference_fitness_landscape(
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    family = family,
    seed = landscape_seed
  )
  reference <- simulate_logistic_missegregation_reference(
    landscape = landscape,
    pm = pm,
    K = K,
    N0 = N0,
    step_size = step_size,
    census_size = census_size
  )
  reference$benchmark <- list(
    family = family,
    replicate = replicate,
    landscape_seed = landscape_seed,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    pm = pm,
    K = K,
    N0 = N0,
    step_size = step_size,
    census_size = census_size
  )
  reference$cache <- list(file = cache_file, status = "built")
  saveRDS(reference, cache_file)
  reference
}

reference_summary_row <- function(reference, family, replicate, pm) {
  data.frame(
    family = family,
    replicate = replicate,
    pm = pm,
    n_reference_states = length(reference$landscape$labels),
    standard_input_rows = nrow(reference$counts),
    occupied_states = reference$sparsity$occupied_states,
    fraction_occupied = reference$sparsity$fraction_occupied,
    T = reference$time$T,
    dt = reference$time$dt,
    K = reference$K,
    N0 = reference$N0,
    census_size = reference$census_size,
    cache_status = reference$cache$status,
    cache_file = reference$cache$file,
    stringsAsFactors = FALSE
  )
}

metric_from_fit <- function(fit,
                            landscape,
                            family,
                            pm,
                            replicate,
                            downsample_fraction,
                            min_obs,
                            depth,
                            fraction_observed,
                            reference_occupied_fraction,
                            layer = c("global", "local"),
                            top_k = 10) {
  layer <- match.arg(layer)
  s <- alfak2::summarize_alfak2(fit, layer = layer)
  truth <- truth_for_summary(s, landscape)
  pred <- s$fitness_mean
  ok <- is.finite(truth) & is.finite(pred)
  err <- pred[ok] - truth[ok]
  coverage95 <- mean(truth[ok] >= s$conf_low[ok] & truth[ok] <= s$conf_high[ok], na.rm = TRUE)
  data.frame(
    family = family,
    pm = pm,
    replicate = replicate,
    downsample_fraction = downsample_fraction,
    min_obs = min_obs,
    depth = depth,
    layer = layer,
    sampling_fraction = downsample_fraction,
    fraction_observed = fraction_observed,
    reference_occupied_fraction = reference_occupied_fraction,
    n_reference_states = length(landscape$labels),
    n_predicted = nrow(s),
    n_scored = sum(ok),
    pearson = safe_cor(pred, truth, "pearson"),
    spearman = safe_cor(pred, truth, "spearman"),
    rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    mae = if (length(err)) mean(abs(err)) else NA_real_,
    coverage95 = coverage95,
    mean_width95 = mean(s$conf_high[ok] - s$conf_low[ok], na.rm = TRUE),
    topk_overlap = top_k_overlap(truth, pred, k = top_k),
    stringsAsFactors = FALSE
  )
}

estimate_elbow_grid <- function(x, y, min_side = 1L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L * min_side + 1L) {
    return(data.frame(breakpoint = NA_real_, rss = NA_real_, n = length(x), status = "too_few_points"))
  }
  candidates <- sort(unique(x))
  best <- data.frame(breakpoint = NA_real_, rss = Inf, n = length(x), status = "no_valid_breakpoint")
  for (b in candidates) {
    left <- which(x <= b)
    right <- which(x > b)
    if (length(left) < min_side || length(right) < min_side) next
    rss <- sum(stats::lm(y[left] ~ x[left])$residuals^2) +
      sum(stats::lm(y[right] ~ x[right])$residuals^2)
    if (rss < best$rss) {
      best <- data.frame(breakpoint = b, rss = rss, n = length(x), status = "ok")
    }
  }
  best
}

coerce_fraction_grid <- function(downsample_fraction_grid,
                                 depths,
                                 census_size) {
  if (!is.null(depths)) {
    return(sort(unique(as.numeric(depths) / as.numeric(census_size))))
  }
  downsample_fraction_grid <- sort(unique(as.numeric(downsample_fraction_grid)))
  if (any(!is.finite(downsample_fraction_grid)) ||
      any(downsample_fraction_grid <= 0) ||
      any(downsample_fraction_grid > 1)) {
    stop("`downsample_fraction_grid` must contain finite values in (0, 1].", call. = FALSE)
  }
  downsample_fraction_grid
}

run_synthetic_reference_task <- function(task,
                                         repo_dir,
                                         min_obs_grid,
                                         downsample_fraction_grid,
                                         n_chr,
                                         min_cn,
                                         max_cn,
                                         K,
                                         N0,
                                         step_size,
                                         census_size,
                                         dropout_prob,
                                         seed,
                                         fit_max_nodes,
                                         local_shell_depth,
                                         global_extra_shell,
                                         force_rebuild_reference,
                                         reference_cache_subdir,
                                         verbose) {
  family <- as.character(task$family)
  replicate <- as.integer(task$replicate)
  pm <- as.numeric(task$pm)
  landscape_seed <- as.integer(task$landscape_seed)

  reference <- load_or_build_reference_standard_data(
    repo_dir = repo_dir,
    family = family,
    replicate = replicate,
    pm = pm,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    K = K,
    N0 = N0,
    step_size = step_size,
    census_size = census_size,
    landscape_seed = landscape_seed,
    force_rebuild = force_rebuild_reference,
    cache_subdir = reference_cache_subdir,
    verbose = verbose
  )
  landscape <- reference$landscape
  reference_row <- reference_summary_row(reference, family, replicate, pm)
  rows <- list()
  fits <- list()
  row_idx <- 0L

  for (min_obs in min_obs_grid) {
    for (downsample_fraction in downsample_fraction_grid) {
      depth <- max(1L, as.integer(round(reference$census_size * downsample_fraction)))
      counts <- stratified_downsample_reference_population(
        reference,
        target_depth = depth,
        detection_threshold = min_obs,
        dropout_prob = dropout_prob,
        seed = seed + depth + 7919L * replicate + round(pm * 1e7) + 104729L * min_obs
      )
      fraction_observed <- nrow(counts) / length(landscape$labels)
      fit <- alfak2::fit_alfak2(
        counts,
        dt = reference$time$dt,
        beta = pm,
        min_cn = min_cn,
        max_cn = max_cn,
        local_shell_depth = local_shell_depth,
        global_extra_shell = global_extra_shell,
        max_nodes = fit_max_nodes,
        lambda_l_grid = c(1),
        lambda_e_grid = c(0.25),
        sigma_obs_grid = c(0.05),
        control = list(eval.max = 200, iter.max = 200)
      )
      fit_key <- paste(family, replicate, pm, min_obs, downsample_fraction, sep = "__")
      fits[[fit_key]] <- fit
      for (layer in c("local", "global")) {
        row_idx <- row_idx + 1L
        rows[[row_idx]] <- metric_from_fit(
          fit,
          landscape = landscape,
          family = family,
          pm = pm,
          replicate = replicate,
          downsample_fraction = downsample_fraction,
          min_obs = min_obs,
          depth = depth,
          fraction_observed = fraction_observed,
          reference_occupied_fraction = reference$sparsity$fraction_occupied,
          layer = layer
        )
      }
    }
  }

  list(
    metrics = do.call(rbind, rows),
    reference_summary = reference_row,
    fits = fits
  )
}

run_synthetic_sparsity_benchmark <- function(repo_dir,
                                             families = c("structured_epistatic", "rugged"),
                                             pm_grid = c(1e-4, 5e-4, 1e-3),
                                             min_obs_grid = c(1L),
                                             downsample_fraction_grid = c(2e-4, 5e-4, 1e-3),
                                             depths = NULL,
                                             replicates = 2,
                                             n_chr = 6,
                                             min_cn = 1,
                                             max_cn = 8,
                                             K = 1e9,
                                             N0 = 100,
                                             step_size = 0.2,
                                             census_size = K,
                                             dropout_prob = 0,
                                             seed = 1,
                                             fit_max_nodes = 15000,
                                             local_shell_depth = 0,
                                             global_extra_shell = 1,
                                             force_rebuild_reference = FALSE,
                                             reference_cache_subdir = "reference_standard_inputs",
                                             parallel_workers = 1L,
                                             verbose = TRUE) {
  dirs <- ensure_benchmark_dirs(repo_dir, "synthetic_sparsity")
  families <- unique(as.character(families))
  pm_grid <- unique(as.numeric(pm_grid))
  downsample_fraction_grid <- coerce_fraction_grid(downsample_fraction_grid, depths, census_size)
  min_obs_grid <- sort(unique(as.integer(min_obs_grid)))
  if (any(is.na(min_obs_grid)) || any(min_obs_grid < 1L)) {
    stop("`min_obs_grid` must contain positive integer thresholds.", call. = FALSE)
  }
  parallel_workers <- as.integer(parallel_workers)
  if (!is.finite(parallel_workers) || parallel_workers < 1L) {
    stop("`parallel_workers` must be a positive integer.", call. = FALSE)
  }

  tasks <- expand.grid(
    family = families,
    replicate = seq_len(replicates),
    pm = pm_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks$family_index <- match(tasks$family, families)
  tasks$landscape_seed <- seed + 1000L * tasks$family_index + tasks$replicate

  run_task <- function(i) {
    run_synthetic_reference_task(
      task = tasks[i, , drop = FALSE],
      repo_dir = repo_dir,
      min_obs_grid = min_obs_grid,
      downsample_fraction_grid = downsample_fraction_grid,
      n_chr = n_chr,
      min_cn = min_cn,
      max_cn = max_cn,
      K = K,
      N0 = N0,
      step_size = step_size,
      census_size = census_size,
      dropout_prob = dropout_prob,
      seed = seed,
      fit_max_nodes = fit_max_nodes,
      local_shell_depth = local_shell_depth,
      global_extra_shell = global_extra_shell,
      force_rebuild_reference = force_rebuild_reference,
      reference_cache_subdir = reference_cache_subdir,
      verbose = verbose
    )
  }

  if (parallel_workers > 1L && .Platform$OS.type != "windows") {
    workers <- min(parallel_workers, nrow(tasks))
    if (isTRUE(verbose)) message("Running synthetic benchmark with ", workers, " parallel workers.")
    task_results <- parallel::mclapply(seq_len(nrow(tasks)), run_task, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    if (parallel_workers > 1L && .Platform$OS.type == "windows") {
      warning("Parallel benchmark execution uses forked workers and is disabled on Windows; running serially.", call. = FALSE)
    }
    task_results <- lapply(seq_len(nrow(tasks)), run_task)
  }

  metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
  reference_summary <- do.call(rbind, lapply(task_results, `[[`, "reference_summary"))
  fits <- unlist(lapply(task_results, `[[`, "fits"), recursive = FALSE)
  elbow <- do.call(rbind, lapply(
    split(metrics, list(metrics$layer, metrics$pm, metrics$min_obs), drop = TRUE),
    function(df) {
      x_measure <- if (length(unique(df$fraction_observed[is.finite(df$fraction_observed)])) >= 3L) {
        "fraction_observed"
      } else {
        "downsample_fraction"
      }
      e <- estimate_elbow_grid(df[[x_measure]], df$spearman, min_side = 1L)
      data.frame(
        layer = unique(df$layer)[1],
        pm = unique(df$pm)[1],
        min_obs = unique(df$min_obs)[1],
        x_measure = x_measure,
        e,
        stringsAsFactors = FALSE
      )
    }
  ))
  utils::write.csv(reference_summary, file.path(dirs$tables, "synthetic_reference_input_summary.csv"), row.names = FALSE)
  utils::write.csv(metrics, file.path(dirs$tables, "synthetic_sparsity_metrics.csv"), row.names = FALSE)
  utils::write.csv(elbow, file.path(dirs$tables, "synthetic_sparsity_elbow.csv"), row.names = FALSE)
  saveRDS(
    list(metrics = metrics, elbow = elbow, reference_summary = reference_summary, fits = fits),
    file.path(dirs$results, "synthetic_sparsity_benchmark.rds")
  )
  list(metrics = metrics, elbow = elbow, reference_summary = reference_summary, fits = fits, dirs = dirs)
}

plot_sparsity_performance <- function(metrics,
                                      metric = "spearman",
                                      x = c("downsample_fraction", "fraction_observed")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  x <- match.arg(x)
  ggplot2::ggplot(metrics, ggplot2::aes(.data[[x]], .data[[metric]], color = .data$layer)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "loess", se = FALSE, formula = y ~ x) +
    ggplot2::facet_grid(min_obs + pm ~ family, labeller = ggplot2::label_both) +
    ggplot2::labs(x = x, y = metric, color = "Layer") +
    ggplot2::theme_bw()
}

plot_error_vs_depth <- function(metrics, metric = "rmse") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  ggplot2::ggplot(metrics, ggplot2::aes(.data$downsample_fraction, .data[[metric]], color = .data$layer)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group = interaction(.data$family, .data$replicate, .data$min_obs, .data$layer)), alpha = 0.35) +
    ggplot2::facet_grid(min_obs + pm ~ family, labeller = ggplot2::label_both) +
    ggplot2::labs(x = "Downsample fraction per timepoint", y = metric, color = "Layer") +
    ggplot2::theme_bw()
}
