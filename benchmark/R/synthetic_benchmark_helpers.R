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

benchmark_format_karyotypes <- function(karyotypes) {
  if (exists("format_karyotypes", mode = "function")) {
    return(as.character(format_karyotypes(karyotypes)))
  }
  if (requireNamespace("alfak2", quietly = TRUE) &&
      "format_karyotypes" %in% getNamespaceExports("alfak2")) {
    return(as.character(alfak2::format_karyotypes(karyotypes)))
  }
  stop("Benchmark karyotype formatting requires `format_karyotypes()`.", call. = FALSE)
}

benchmark_predict_landscape_fitness <- function(landscape, karyotypes) {
  if (exists("predict_landscape_fitness", mode = "function")) {
    return(as.numeric(predict_landscape_fitness(landscape, karyotypes)))
  }
  if (requireNamespace("alfak2", quietly = TRUE) &&
      "predict_landscape_fitness" %in% getNamespaceExports("alfak2")) {
    return(as.numeric(alfak2::predict_landscape_fitness(landscape, karyotypes)))
  }
  stop("Lazy benchmark landscapes require `predict_landscape_fitness()`.", call. = FALSE)
}

truth_for_summary <- function(summary_df, landscape) {
  labels <- as.character(summary_df$karyotype)
  out <- rep(NA_real_, nrow(summary_df))
  if (!is.null(landscape$labels) && !is.null(landscape$fitness)) {
    idx <- match(labels, as.character(landscape$labels))
    out[!is.na(idx)] <- landscape$fitness[idx[!is.na(idx)]]
    missing <- is.na(idx)
  } else {
    missing <- rep(TRUE, length(labels))
  }
  if (any(missing) && inherits(landscape, "alfak2_grf_landscape")) {
    out[missing] <- benchmark_predict_landscape_fitness(landscape, labels[missing])
  }
  out
}

benchmark_reference_state_count <- function(landscape) {
  if (!is.null(landscape$labels)) return(length(landscape$labels))
  if (!is.null(landscape$n_chr) && !is.null(landscape$min_cn) && !is.null(landscape$max_cn)) {
    return((as.numeric(landscape$max_cn) - as.numeric(landscape$min_cn) + 1)^as.numeric(landscape$n_chr))
  }
  NA_real_
}

materialize_benchmark_landscape <- function(landscape,
                                            n_chr,
                                            min_cn,
                                            max_cn,
                                            diploid_cn,
                                            diploid_fitness,
                                            max_materialized_states = 1e6) {
  required <- c("labels", "karyotypes", "fitness", "min_cn", "max_cn")
  if (is.list(landscape) && all(required %in% names(landscape))) {
    return(landscape)
  }
  states_per_chr <- as.numeric(max_cn) - as.numeric(min_cn) + 1
  n_states <- states_per_chr^as.numeric(n_chr)
  if (!is.finite(n_states) || n_states > max_materialized_states) {
    stop(
      "Synthetic logistic reference requires a materialized truth lattice; requested ",
      format(n_states, scientific = FALSE, big.mark = ","),
      " states, above `max_materialized_states = ",
      format(max_materialized_states, scientific = FALSE, big.mark = ","),
      "`.",
      call. = FALSE
    )
  }

  grid <- expand.grid(
    replicate(n_chr, seq.int(min_cn, max_cn), simplify = FALSE),
    KEEP.OUT.ATTRS = FALSE
  )
  karyotypes <- as.matrix(grid)
  storage.mode(karyotypes) <- "integer"
  labels <- benchmark_format_karyotypes(karyotypes)
  fitness <- benchmark_predict_landscape_fitness(landscape, karyotypes)
  diploid_index <- which(rowSums(abs(karyotypes - diploid_cn)) == 0)
  if (length(diploid_index) != 1L) {
    stop("Could not identify exactly one diploid state in materialized landscape.", call. = FALSE)
  }

  landscape$labels <- labels
  landscape$karyotypes <- karyotypes
  landscape$fitness <- fitness
  landscape$min_cn <- min_cn
  landscape$max_cn <- max_cn
  landscape$diploid_cn <- diploid_cn
  landscape$diploid_index <- diploid_index
  landscape$diploid_fitness <- diploid_fitness
  landscape$materialized_states <- nrow(karyotypes)
  class(landscape) <- unique(c("alfak2_materialized_landscape", class(landscape)))
  landscape
}

benchmark_token <- function(x) {
  x <- if (is.numeric(x)) format(x, scientific = TRUE, digits = 12) else as.character(x)
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

benchmark_or <- function(x, y) {
  if (is.null(x)) y else x
}

reference_standard_cache_dir <- function(repo_dir,
                                         subdir = "reference_standard_inputs") {
  path <- file.path(repo_dir, "benchmark", "data", subdir)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

reference_standard_cache_file <- function(repo_dir,
                                          replicate,
                                          pm,
                                          n_chr,
                                          min_cn,
                                          max_cn,
                                          diploid_cn,
                                          diploid_fitness,
                                          lower,
                                          upper,
                                          ell,
                                          use_full_range,
                                          K,
                                          N0,
                                          step_size,
                                          census_size,
                                          landscape_seed,
                                          subdir = "reference_standard_inputs") {
  name <- paste(
    "standard",
    "l1_gp",
    paste0("rep", replicate),
    paste0("pm", benchmark_token(pm)),
    paste0("chr", n_chr),
    paste0("cn", min_cn, "-", max_cn),
    paste0("dip", diploid_cn),
    paste0("dfit", benchmark_token(diploid_fitness)),
    paste0("range", benchmark_token(lower), "_", benchmark_token(upper)),
    paste0("ell", benchmark_token(ell)),
    paste0("ufr", benchmark_token(use_full_range)),
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

build_benchmark_truth_landscape <- function(n_chr,
                                            min_cn,
                                            max_cn,
                                            diploid_cn,
                                            diploid_fitness,
                                            lower,
                                            upper,
                                            ell,
                                            use_full_range,
                                            seed) {
  fn <- if (exists("simulate_l1_gp_landscape", mode = "function")) {
    simulate_l1_gp_landscape
  } else {
    alfak2::simulate_l1_gp_landscape
  }
  fn(
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    seed = seed,
    use_full_range = use_full_range,
    include_table = FALSE
  )
}

load_or_build_reference_standard_data <- function(repo_dir,
                                                  replicate,
                                                  pm,
                                                  n_chr = 6,
                                                  min_cn = 1,
                                                  max_cn = 8,
                                                  diploid_cn = 2,
                                                  diploid_fitness = 1,
                                                  lower = -5,
                                                  upper = 5,
                                                  ell = 2.5,
                                                  use_full_range = 0.95,
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
    replicate = replicate,
    pm = pm,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    use_full_range = use_full_range,
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
  landscape <- build_benchmark_truth_landscape(
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    use_full_range = use_full_range,
    seed = landscape_seed
  )
  landscape <- materialize_benchmark_landscape(
    landscape,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness
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
    landscape_model = "l1_gp",
    replicate = replicate,
    landscape_seed = landscape_seed,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    use_full_range = use_full_range,
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

reference_summary_row <- function(reference, replicate, pm) {
  cache <- benchmark_or(reference$cache, list(status = NA_character_, file = NA_character_))
  data.frame(
    landscape_model = "l1_gp",
    replicate = replicate,
    pm = pm,
    n_reference_states = benchmark_reference_state_count(reference$landscape),
    truth_family = benchmark_or(reference$landscape$family, "l1_gp"),
    covariance = benchmark_or(reference$landscape$covariance, NA_character_),
    standard_input_rows = nrow(reference$counts),
    occupied_states = reference$sparsity$occupied_states,
    fraction_occupied = reference$sparsity$fraction_occupied,
    T = reference$time$T,
    dt = reference$time$dt,
    K = reference$K,
    N0 = reference$N0,
    census_size = reference$census_size,
    cache_status = benchmark_or(cache$status, NA_character_),
    cache_file = benchmark_or(cache$file, NA_character_),
    stringsAsFactors = FALSE
  )
}

metric_from_fit <- function(fit,
                            landscape,
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
    landscape_model = "l1_gp",
    pm = pm,
    replicate = replicate,
    downsample_fraction = downsample_fraction,
    min_obs = min_obs,
    depth = depth,
    layer = layer,
    sampling_fraction = downsample_fraction,
    fraction_observed = fraction_observed,
    reference_occupied_fraction = reference_occupied_fraction,
    n_reference_states = benchmark_reference_state_count(landscape),
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
                                         diploid_cn,
                                         diploid_fitness,
                                         lower,
                                         upper,
                                         ell,
                                         use_full_range,
                                         K,
                                         N0,
                                         step_size,
                                         census_size,
                                         dropout_prob,
                                         seed,
                                         fit_max_nodes,
                                         local_shell_depth,
                                         global_extra_shell,
                                         input_depth,
                                         effective_depth,
                                         effective_depth_mode,
                                         effective_depth_rounding,
                                         effective_depth_seed,
                                         observation_model,
                                         dm_concentration,
                                         observation_weight_mode,
                                         force_rebuild_reference,
                                         reference_cache_subdir,
                                         verbose) {
  replicate <- as.integer(task$replicate)
  pm <- as.numeric(task$pm)
  landscape_seed <- as.integer(task$landscape_seed)

  reference <- load_or_build_reference_standard_data(
    repo_dir = repo_dir,
    replicate = replicate,
    pm = pm,
    n_chr = n_chr,
    min_cn = min_cn,
    max_cn = max_cn,
    diploid_cn = diploid_cn,
    diploid_fitness = diploid_fitness,
    lower = lower,
    upper = upper,
    ell = ell,
    use_full_range = use_full_range,
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
  reference_row <- reference_summary_row(reference, replicate, pm)
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
        input_depth = input_depth,
        effective_depth = effective_depth,
        effective_depth_mode = effective_depth_mode,
        effective_depth_rounding = effective_depth_rounding,
        effective_depth_seed = effective_depth_seed,
        observation_model = observation_model,
        dm_concentration = dm_concentration,
        observation_weight_mode = observation_weight_mode,
        control = list(eval.max = 200, iter.max = 200)
      )
      fit_key <- paste("l1_gp", replicate, pm, min_obs, downsample_fraction, sep = "__")
      fits[[fit_key]] <- fit
      for (layer in c("local", "global")) {
        row_idx <- row_idx + 1L
        rows[[row_idx]] <- metric_from_fit(
          fit,
          landscape = landscape,
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
                                             pm_grid = c(1e-4, 5e-4, 1e-3),
                                             min_obs_grid = c(1L),
                                             downsample_fraction_grid = c(2e-4, 5e-4, 1e-3),
                                             depths = NULL,
                                             replicates = 2,
                                             n_chr = 6,
                                             min_cn = 1,
                                             max_cn = 8,
                                             diploid_cn = 2,
                                             diploid_fitness = 1,
                                             lower = -5,
                                             upper = 5,
                                             ell = 2.5,
                                             use_full_range = 0.95,
                                             K = 1e9,
                                             N0 = 100,
                                             step_size = 0.2,
                                             census_size = K,
                                             dropout_prob = 0,
                                             seed = 1,
                                             fit_max_nodes = 15000,
                                             local_shell_depth = 0,
                                             global_extra_shell = 1,
                                             input_depth = "raw",
                                             effective_depth = NULL,
                                             effective_depth_mode = "min",
                                             effective_depth_rounding = "hash",
                                             effective_depth_seed = NULL,
                                             observation_model = NULL,
                                             dm_concentration = NULL,
                                             observation_weight_mode = "likelihood",
                                             force_rebuild_reference = FALSE,
                                             reference_cache_subdir = "reference_standard_inputs",
                                             parallel_workers = 1L,
                                             verbose = TRUE) {
  dirs <- ensure_benchmark_dirs(repo_dir, "synthetic_sparsity")
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
    replicate = seq_len(replicates),
    pm = pm_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks$landscape_seed <- seed + tasks$replicate

  run_task <- function(i) {
    run_synthetic_reference_task(
      task = tasks[i, , drop = FALSE],
      repo_dir = repo_dir,
      min_obs_grid = min_obs_grid,
      downsample_fraction_grid = downsample_fraction_grid,
      n_chr = n_chr,
      min_cn = min_cn,
      max_cn = max_cn,
      diploid_cn = diploid_cn,
      diploid_fitness = diploid_fitness,
      lower = lower,
      upper = upper,
      ell = ell,
      use_full_range = use_full_range,
      K = K,
      N0 = N0,
      step_size = step_size,
      census_size = census_size,
      dropout_prob = dropout_prob,
      seed = seed,
      fit_max_nodes = fit_max_nodes,
      local_shell_depth = local_shell_depth,
      global_extra_shell = global_extra_shell,
      input_depth = input_depth,
      effective_depth = effective_depth,
      effective_depth_mode = effective_depth_mode,
      effective_depth_rounding = effective_depth_rounding,
      effective_depth_seed = effective_depth_seed,
      observation_model = observation_model,
      dm_concentration = dm_concentration,
      observation_weight_mode = observation_weight_mode,
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
      stop("Forked parallelism is unavailable on Windows; set `parallel_workers = 1`.", call. = FALSE)
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
    ggplot2::facet_grid(min_obs ~ pm, labeller = ggplot2::label_both) +
    ggplot2::labs(x = x, y = metric, color = "Layer") +
    ggplot2::theme_bw()
}

plot_error_vs_depth <- function(metrics, metric = "rmse") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  ggplot2::ggplot(metrics, ggplot2::aes(.data$downsample_fraction, .data[[metric]], color = .data$layer)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group = interaction(.data$replicate, .data$min_obs, .data$pm, .data$layer)), alpha = 0.35) +
    ggplot2::facet_grid(min_obs ~ pm, labeller = ggplot2::label_both) +
    ggplot2::labs(x = "Downsample fraction per timepoint", y = metric, color = "Layer") +
    ggplot2::theme_bw()
}
