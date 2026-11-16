# The three multi-timepoint methods of fit_sparse_joint() compared on three
# simulated sets (README.md in this directory):
#   e1  the stage-2 passaged cultures of ../multi_timepoint (288 datasets),
#       linked, not copied;
#   e2  clonal sweeps: the same landscapes sampled from day 240 or 480, 120
#       days apart, where clones replace one another between time points;
#   e3  neutral: the stage-2 design, but from the first time point on every viable
#       karyotype divides at one rate, so every fitness difference a method
#       reports is a false positive.
# Inputs in results/inputs/<set>, fits and scores of a build in
# results/<build>/<set> (git-ignored). The simulator, landscapes (seed base
# 97000) and truth populations are those of ../multi_timepoint.

sm_dir <- if (exists("sm_script_dir")) sm_script_dir else getwd()
mt_script_dir <- normalizePath(file.path(sm_dir, "..", "multi_timepoint"))
source(file.path(mt_script_dir, "setup.R"))

sm_default_output <- file.path(si_repo, "validation", "sparse_joint", "series_methods", "results")
sm_sets <- c("e1", "e2", "e3")

# `build` NULL: the set's inputs (datasets stage); otherwise that build's fits.
sm_config <- function(set = sm_sets, build = NULL) {
  set <- match.arg(set)
  inputs <- if (set == "e1") file.path(mt_default_output, "series") else file.path(sm_default_output, "inputs", set)
  config <- mt_config(if (is.null(build)) inputs else file.path(sm_default_output, build, set), "series")
  config$set <- set
  config$build <- build
  config$input_dir <- inputs
  config$library_root <- mt_default_output
  if (set == "e2") {
    # Day 1440 is late enough that clones rarely replace one another within a
    # window; from day 240 or 480, 120 days apart, an intermediate population
    # lies mostly off the majors of the first and last time points in half or more of
    # the windows (truth populations, checked 2026-11-14).
    config$starts <- c(240, 480)
    config$spacings <- 120
    config$bottlenecks <- c(2e4, 2e3)
  }
  config$methods <- sm_methods()
  # Effect sizes for power (per day; true relative fitness of scored
  # karyotypes is typically 0.02-0.07 away from the reference), and the band
  # counted as no effect for false positives.
  config$power_effects <- c(0.02, 0.05, 0.1)
  config$null_band <- 0.005
  config
}

# Every method fits the same counts with the same settings as stage 2; only
# the multi-time-point method differs. `ends`: the first and last time point, a
# two-time-point fit.
sm_methods <- function() {
  # `options`: alfak2's internal options set during the fit (the variants
  # compared while choosing a method's details).
  series <- function(columns, ..., options = list()) {
    list(type = "series", birth_rate = 1, min_cn = 1L, max_cn = 8L, columns = columns, arguments = list(...),
         options = options)
  }
  methods <- list(
    zero = list(type = "zero"),
    ends = series("ends"),
    anchored = series("all", series = "anchored"),
    # Pairwise: a pair starts at balance over the history plus the
    # population's age at its first time point (default), or over the history
    # alone.
    pairwise = series("all", series = "pairwise"),
    pairwise_same = series("all", series = "pairwise", options = list(alfak2.pair_history = "same")),
    state_space = series("all", series = "state_space")
  )
  # Choosing the missegregation rate (validation/sparse_joint/beta_choice,
  # step 2): `ends` and `anchored` at the dataset's rate times each other
  # multiple of `beta = "auto"`, one fit per task; `--stage=compose` keeps, per
  # dataset, the candidate of lowest objective as `ends_auto` and
  # `anchored_auto`, which is what fit_sparse_joint(beta = "auto") returns.
  for (name in names(sm_beta_multiples)) {
    for (base in c("ends", "anchored")) {
      method <- methods[[base]]
      method$beta_multiple <- sm_beta_multiples[[name]]
      methods[[paste0(base, "_", name)]] <- method
    }
  }
  for (base in c("ends", "anchored")) {
    methods[[paste0(base, "_auto")]] <- list(type = "composed",
                                             from = c(base, paste0(base, "_", names(sm_beta_multiples))))
  }
  methods
}

sm_beta_multiples <- c(b1o16 = 1 / 16, b1o4 = 1 / 4, b4 = 4, b16 = 16, b64 = 64)

# Per dataset, the candidate fit of lowest objective of each composed method.
sm_stage_compose <- function(config) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  composed <- Filter(function(m) identical(m$type, "composed"), config$methods)
  for (name in names(composed)) {
    from <- composed[[name]]$from
    made <- missing <- 0L
    for (dataset_id in index$dataset_id) {
      paths <- vapply(from, function(m) si_fit_path(config, dataset_id, m), character(1))
      if (!all(file.exists(paths))) {
        missing <- missing + 1L
        next
      }
      fits <- lapply(paths, readRDS)
      ok <- vapply(fits, function(f) identical(f$status, "ok"), logical(1))
      objective <- vapply(fits, function(f) {
        o <- f$diagnostics$objective
        if (identical(f$status, "ok") && length(o) && is.finite(o)) as.numeric(o) else NA_real_
      }, numeric(1))
      multiple <- vapply(from, function(m) {
        b <- config$methods[[m]]$beta_multiple
        if (is.null(b)) 1 else b
      }, numeric(1))
      if (all(is.na(objective))) {
        result <- fits[[1L]]
      } else {
        best <- which.min(objective)
        result <- fits[[best]]
        result$beta_multiple <- multiple[[best]]
      }
      result$beta_profile <- data.frame(method = from, beta_multiple = multiple, objective = objective,
                                        status = vapply(fits, function(f) f$status, character(1)))
      result$seconds <- sum(vapply(fits, function(f) as.numeric(f$seconds), numeric(1)))
      result$method <- name
      path <- si_fit_path(config, dataset_id, name)
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, path)
      made <- made + 1L
    }
    cat(name, ": composed", made, "datasets;", missing, "with a candidate missing\n")
  }
  invisible()
}

# Links from a build's directory to the set's inputs.
sm_link <- function(config) {
  if (is.null(config$build)) return(invisible())
  dir.create(file.path(config$output_dir, "fits"), recursive = TRUE, showWarnings = FALSE)
  for (name in c("datasets", "nodes", "dataset_index.csv", "truth")) {
    path <- file.path(config$output_dir, name)
    target <- file.path(config$input_dir, name)
    if (!file.exists(target)) next
    if (!file.exists(path)) file.symlink(target, path)
    if (!file.exists(path)) stop("Link ", path, " -> ", target, " is broken.", call. = FALSE)
  }
  invisible()
}

# ---- Datasets (e2, e3) ------------------------------------------------------

sm_dataset_id <- function(set, landscape_id, start, bottleneck, spacing, depth, samples) {
  if (set == "e2") {
    sprintf("%s_t%04d_b%05d_s%03d_n%04d_k%d", landscape_id, as.integer(start), as.integer(bottleneck),
            as.integer(spacing), as.integer(depth), as.integer(samples))
  } else {
    mt_dataset_id(landscape_id, bottleneck, spacing, depth, samples)
  }
}

sm_window_plan <- function(config) {
  starts <- if (config$set == "e2") config$starts else config$t0
  plan <- list()
  for (i in seq_len(nrow(config$landscapes))) for (st in seq_along(starts)) {
    for (b in seq_along(config$bottlenecks)) for (s in seq_along(config$spacings)) {
      # Seeds apart from stage 2's (landscape_seed * 1000 + b * 10 + s).
      offset <- if (config$set == "e2") 500L + st * 100L else 700L
      plan[[length(plan) + 1L]] <- list(
        landscape = config$landscapes[i, ], start = starts[[st]], bottleneck = config$bottlenecks[[b]],
        spacing = config$spacings[[s]],
        seed = as.integer(config$landscapes$landscape_seed[[i]] * 1000L + offset + b * 10L + s)
      )
    }
  }
  plan
}

# stt_run_abm() with the study's copy of the simulator: `constant_fitness` > 0
# makes every viable karyotype divide at that rate.
sm_run_abm <- function(initial, centroids, wavelength, config, record_days, seed, beta, constant_fitness = 0) {
  abm <- config$abm
  record_steps <- as.integer(round(sort(unique(record_days)) / abm$delta_t))
  raw <- sm_abm_cpp(
    initial_karyotypes = matrix(as.integer(stt_parse_karyotypes(names(initial))), nrow = length(initial)),
    initial_counts = as.numeric(initial), centroids = centroids, wavelength = wavelength, beta = beta,
    delta_t = abm$delta_t, n_steps = max(record_steps), record_steps = record_steps,
    max_population = abm$max_pop, culling_survival = abm$culling_survival,
    min_cn = abm$min_cn, max_cn = abm$max_cn, seed = seed, constant_fitness = constant_fitness
  )
  lapply(stats::setNames(record_steps, as.character(record_steps * abm$delta_t)), function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
}

sm_load_abm <- function(config) {
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  if (!exists("sm_abm_cpp", mode = "function")) {
    cache <- file.path(config$output_dir, "cpp_cache")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    Rcpp::sourceCpp(file.path(sm_dir, "abm_constant.cpp"), cacheDir = cache, env = globalenv())
  }
  invisible(TRUE)
}

# The populations of one window: the truth population at the window's start,
# then passaged through `bottleneck` cells and recorded every `spacing` days.
# In e3 every viable karyotype divides at the start population's mean rate.
sm_make_window <- function(config, truth, item) {
  window_config <- config
  window_config$abm$culling_survival <- item$bottleneck / config$nb
  days <- item$start + item$spacing * (seq_len(config$max_samples) - 1L)
  start <- stt_snapshot(truth, item$start)
  constant <- 0
  if (config$set == "e3") {
    fitness <- stt_grf_fitness(names(start), truth$centroids, truth$wavelength)
    constant <- sum(start * pmax(fitness, 0)) / sum(start)
  }
  later <- sm_run_abm(start, truth$centroids, truth$wavelength, window_config,
                      record_days = days[-1L] - item$start, seed = item$seed, beta = truth$beta,
                      constant_fitness = constant)
  list(populations = c(list(start), unname(later)), constant_fitness = constant)
}

sm_stage_datasets <- function(config, workers) {
  if (config$set == "e1") stop("e1 links the stage-2 datasets; nothing to make.", call. = FALSE)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  truth_link <- file.path(config$output_dir, "truth")
  if (!file.exists(truth_link)) file.symlink(file.path(mt_default_output, "series", "truth"), truth_link)
  sm_load_abm(config)
  plan <- sm_window_plan(config)
  rows <- si_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
    ids <- character()
    window <- NULL
    for (d in seq_along(config$depths)) for (k in config$samples) {
      depth <- config$depths[[d]]
      id <- sm_dataset_id(config$set, truth$landscape_id, item$start, item$bottleneck, item$spacing, depth, k)
      ids <- c(ids, id)
      path <- stt_dataset_path(config, id)
      if (file.exists(path)) next
      if (is.null(window)) window <- sm_make_window(config, truth, item)
      populations <- window$populations
      counts <- mt_sample_columns(populations, depth, item$seed * 10L + d)
      counts <- counts[, seq_len(k), drop = FALSE]
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      counts <- counts[order(rowSums(counts), decreasing = TRUE), , drop = FALSE]
      times <- item$spacing * (seq_len(k) - 1L)
      colnames(counts) <- as.character(times)
      dataset <- list(
        dataset_id = id, set = config$set, landscape_id = truth$landscape_id, patient = 1L, split = truth$split,
        wavelength_class = truth$wavelength_class,
        scenario = sprintf("t%d_b%g_s%d_n%d_k%d", item$start, item$bottleneck, item$spacing, depth, k),
        t0 = item$start, times = times, interval = max(times), samples = k, spacing = item$spacing,
        bottleneck = item$bottleneck, depth = depth, seed = item$seed * 10L + d, beta = truth$beta,
        population_sizes = vapply(populations[seq_len(k)], sum, numeric(1)),
        counts = counts,
        truth = list(centroids = truth$centroids, wavelength = truth$wavelength, abm_delta_t = truth$abm_delta_t,
                     constant_fitness = if (window$constant_fitness > 0) window$constant_fitness,
                     frequencies = lapply(populations[seq_len(k)], function(p) p / sum(p)))
      )
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      temporary <- paste0(path, ".tmp", Sys.getpid())
      saveRDS(dataset, temporary)
      file.rename(temporary, path)
    }
    ids
  }, workers)
  ids <- unlist(rows)
  index <- do.call(rbind, lapply(ids, function(id) {
    dataset <- readRDS(stt_dataset_path(config, id))
    mt_dataset_nodes(config, dataset)
    counts <- dataset$counts
    data.frame(
      dataset_id = id, landscape_id = dataset$landscape_id, wavelength_class = dataset$wavelength_class,
      start = dataset$t0, bottleneck = dataset$bottleneck, spacing = dataset$spacing, depth = dataset$depth,
      samples = dataset$samples, n_karyotypes = nrow(counts),
      shared_all = sum(rowSums(counts > 0) == ncol(counts)), top_fraction = max(rowSums(counts)) / sum(counts),
      min_population = min(dataset$population_sizes), stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(index, si_index_path(config), row.names = FALSE)
  cat("Datasets:", nrow(index), " index:", si_index_path(config), "\n")
}

# How far each dataset's time points in between lie from the majors of the
# first and last time points (with the installed alfak2's own functions, defaults of
# fit_sparse_joint()): the gate of e2 is that clonal sweeps reproduce what the
# Salehi trajectories show.
sm_sweep_table <- function(config) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  rows <- lapply(index$dataset_id, function(id) {
    dataset <- readRDS(stt_dataset_path(config, id))
    counts <- dataset$counts
    k <- ncol(counts)
    coverage <- alfak2:::sparse_joint_domain_coverage(
      counts, alfak2:::sparse_joint_domain(counts, 0.01, 3, 5000, 1L, 8L, c(1L, k)), 0.01, 3)
    middle <- if (k > 2L) seq(2L, k - 1L) else integer(0)
    data.frame(dataset_id = id, samples = k,
               middle_no_major = sum(coverage$majors[middle] == 0),
               middle_no_domain = sum(coverage$cells_in_domain[middle] == 0),
               ends_no_domain = sum(coverage$cells_in_domain[c(1L, k)] == 0),
               min_middle_fraction = if (length(middle)) min(coverage$fraction[middle]) else NA_real_,
               stringsAsFactors = FALSE)
  })
  merge(index, do.call(rbind, rows), by = c("dataset_id", "samples"))
}

# ---- Fits ---------------------------------------------------------------------

sm_true_growth_rate <- function(labels, truth) {
  if (!is.null(truth$constant_fitness)) {
    return(rep(stt_effective_growth_rate(truth$constant_fitness, truth$abm_delta_t), length(labels)))
  }
  stt_true_growth_rate(labels, truth)
}

sm_fit_series <- function(config, dataset, nodes, method) {
  counts <- mt_counts(dataset, method$columns)
  times <- if (identical(method$columns, "ends")) dataset$times[c(1L, length(dataset$times))] else dataset$times
  if (length(method$options)) {
    old <- options(method$options)
    on.exit(options(old), add = TRUE)
  }
  warnings <- character()
  fit <- withCallingHandlers(
    do.call(alfak2::fit_sparse_joint, c(
      list(counts, times = times, beta = dataset$beta * if (is.null(method$beta_multiple)) 1 else method$beta_multiple,
           birth_rate = method$birth_rate, min_cn = method$min_cn,
           max_cn = method$max_cn),
      method$arguments)),
    warning = function(w) {
      warnings <<- c(warnings, paste(class(w)[[1L]], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  prediction <- stats::predict(fit, nodes$karyotype)
  d <- fit$diagnostics
  keep <- c("convergence", "objective", "finite_uncertainty", "n_states", "n_major", "dropped_karyotypes",
            "dropped_cells", "fitted_birth_rate", "history", "optimisation", "times", "intermediate_weight",
            "intermediate_dispersion", "dropped_samples", "series", "pair_history", "drift")
  list(status = "ok",
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                               shell = prediction$shell, stringsAsFactors = FALSE),
       diagnostics = d[intersect(keep, names(d))], warnings = warnings,
       history_profile = d$history_profile, domain_coverage = d$domain_coverage)
}

sm_run_method <- function(config, dataset, nodes, method) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    switch(method$type,
      zero = list(status = "ok", prediction = data.frame(karyotype = nodes$karyotype, fitness = 0,
                                                         fitness_sd = NA_real_, stringsAsFactors = FALSE)),
      series = sm_fit_series(config, dataset, nodes, method),
      composed = stop("Composed methods are made by --stage=compose.", call. = FALSE),
      stop("Unknown method type: ", method$type, call. = FALSE)
    ),
    error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e),
                             error_class = class(e)[[1L]])
  )
  result$seconds <- proc.time()[["elapsed"]] - started
  result$peak_mb <- sm_peak_mb()
  result
}

# Peak resident memory of this process so far (Linux; NA elsewhere).
sm_peak_mb <- function() {
  status <- "/proc/self/status"
  if (!file.exists(status)) return(NA_real_)
  line <- grep("^VmHWM:", readLines(status), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line)) / 1024
}

sm_stage_fit <- function(config, workers, method_names, datasets = NULL) {
  sm_link(config)
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  if (!is.null(datasets)) index <- index[grepl(datasets, index$dataset_id), , drop = FALSE]
  if (any(vapply(config$methods[method_names], function(m) m$type == "series", logical(1)))) {
    library_dir <- si_library(config$build, config$library_root)
    if (!dir.exists(file.path(library_dir, "alfak2"))) stop("No alfak2 build in ", library_dir, call. = FALSE)
    .libPaths(c(library_dir, .libPaths()))
    loadNamespace("alfak2")
    cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  }
  tasks <- list()
  for (dataset_id in index$dataset_id) for (method_name in method_names) {
    if (!file.exists(si_fit_path(config, dataset_id, method_name))) {
      tasks[[length(tasks) + 1L]] <- list(dataset_id = dataset_id, method_name = method_name)
    }
  }
  cat("Fits to run:", length(tasks), "\n")
  if (!length(tasks)) return(invisible())
  si_parallel(tasks, function(task) {
    dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
    nodes <- mt_dataset_nodes(config, dataset)
    set.seed(dataset$seed)
    result <- sm_run_method(config, dataset, nodes, config$methods[[task$method_name]])
    result$dataset_id <- task$dataset_id
    result$method <- task$method_name
    path <- si_fit_path(config, task$dataset_id, task$method_name)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(result, temporary)
    file.rename(temporary, path)
    TRUE
  }, workers)
  invisible()
}

# ---- Scores -------------------------------------------------------------------

# Power and false positives by distance from the observed karyotypes. A
# karyotype is called different from the reference when its 95% interval
# excludes zero. Power at effect e: share of karyotypes with |truth| >= e
# called different in the right direction; wrong_sign: called in the wrong
# direction. False positives: share called different among karyotypes with
# |truth| < null_band (every karyotype in e3).
sm_power <- function(score_nodes, effects, null_band) {
  x <- score_nodes
  called <- x$available & is.finite(x$pred_sd) & abs(x$pred_relative) > stats::qnorm(0.975) * x$pred_sd
  right <- called & sign(x$pred_relative) == sign(x$truth_relative)
  rows <- lapply(c("all", "d0", "d1", "d2"), function(stratum) {
    use <- if (stratum == "all") rep(TRUE, nrow(x)) else x$distance == as.integer(sub("d", "", stratum))
    if (!any(use)) return(NULL)
    null <- use & abs(x$truth_relative) < null_band
    out <- data.frame(stratum = stratum, n_nodes = sum(use), n_null = sum(null),
                      false_positive = if (any(null)) mean(called[null]) else NA_real_,
                      stringsAsFactors = FALSE)
    for (e in effects) {
      hit <- use & abs(x$truth_relative) >= e
      out[[sprintf("n_effect_%g", e)]] <- sum(hit)
      out[[sprintf("power_%g", e)]] <- if (any(hit)) mean(right[hit]) else NA_real_
      out[[sprintf("wrong_sign_%g", e)]] <- if (any(hit)) mean((called & !right)[hit]) else NA_real_
    }
    out
  })
  do.call(rbind, rows)
}

# Power at a controlled false-positive rate. The reported 95% intervals
# cannot exclude zero: their sd holds deficit x |relative fitness| (0.6, or
# 1.1 unobserved), so |estimate| / sd <= 1 / 0.6. Each method's own
# posterior z = estimate / model sd is used instead, with a threshold that
# report.R sets on e3 (no fitness differences) for a 5% false-positive rate.
# Here, per stratum: karyotypes with |truth| < null_band ("null") and with
# |truth| >= e, and how many exceed each threshold of `z_grid` (with the
# right sign for the latter; z = 0 gives the direction accuracy).
sm_z_grid <- c(seq(0, 6, by = 0.1), seq(6.5, 30, by = 0.5), 40, 60, 100)

sm_power_counts <- function(score_nodes, model_sd, effects, null_band, z_grid = sm_z_grid) {
  x <- score_nodes
  ok <- x$available & is.finite(model_sd) & model_sd > 0
  z <- ifelse(ok, abs(x$pred_relative) / model_sd, NA_real_)
  right <- ok & sign(x$pred_relative) == sign(x$truth_relative)
  rows <- list()
  for (stratum in c("all", "d0", "d1", "d2")) {
    use <- if (stratum == "all") rep(TRUE, nrow(x)) else x$distance == as.integer(sub("d", "", stratum))
    kinds <- c(list(null = use & abs(x$truth_relative) < null_band),
               stats::setNames(lapply(effects, function(e) use & abs(x$truth_relative) >= e),
                               sprintf("effect_%g", effects)))
    for (kind in names(kinds)) {
      member <- kinds[[kind]]
      if (!any(member)) next
      pass <- if (kind == "null") ok & member else right & member
      rows[[length(rows) + 1L]] <- data.frame(
        stratum = stratum, kind = kind, z = z_grid, n = sum(member),
        exceed = vapply(z_grid, function(t) sum(pass & z > t, na.rm = TRUE), numeric(1)),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

sm_stage_score <- function(config, workers) {
  sm_link(config)
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  methods <- intersect(names(config$methods), list.files(file.path(config$output_dir, "fits")))
  scored <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- mt_dataset_nodes(config, dataset)
    truth <- sm_true_growth_rate(nodes$karyotype, dataset$truth)
    rows <- lapply(methods, function(method) {
      path <- si_fit_path(config, dataset_id, method)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      key <- data.frame(dataset_id = dataset_id, method = method, status = fit$status, seconds = fit$seconds,
                        peak_mb = if (is.null(fit$peak_mb)) NA_real_ else fit$peak_mb, stringsAsFactors = FALSE)
      if (!identical(fit$status, "ok")) {
        key$detail <- fit$detail
        return(list(fit = key))
      }
      score <- stt_score_predictions(nodes, truth, fit$prediction, config$scoring$reference_min_count)
      d <- fit$diagnostics
      num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[[1L]]
      fit_row <- cbind(key, data.frame(
        convergence = num(d$convergence), history = num(d$history), n_states = num(d$n_states),
        n_major = num(d$n_major), dropped_samples = if (is.null(d$dropped_samples)) NA_real_ else NROW(d$dropped_samples),
        intermediate_weight = num(d$intermediate_weight), n_warnings = length(fit$warnings),
        stringsAsFactors = FALSE))
      model_sd <- if ("model_sd" %in% names(fit$prediction)) {
        fit$prediction$model_sd[match(score$nodes$karyotype, fit$prediction$karyotype)]
      } else rep(NA_real_, nrow(score$nodes))
      list(distance = cbind(key, score$metrics),
           power = cbind(key, sm_power(score$nodes, config$power_effects, config$null_band)),
           counts = cbind(key[c("dataset_id", "method")],
                          sm_power_counts(score$nodes, model_sd, config$power_effects, config$null_band)),
           fit = fit_row)
    })
    rows <- Filter(Negate(is.null), rows)
    lapply(c(distance = "distance", power = "power", counts = "counts", fit = "fit"),
           function(part) si_bind(lapply(rows, `[[`, part)))
  }, workers)
  saveRDS(merge(index, si_bind(lapply(scored, `[[`, "counts")), by = "dataset_id"),
          file.path(config$output_dir, "power_counts.rds"))
  for (part in c("distance", "power", "fit")) {
    table <- merge(index, si_bind(lapply(scored, `[[`, part)), by = "dataset_id")
    utils::write.csv(table, file.path(config$output_dir, sprintf("metrics_%s.csv", part)), row.names = FALSE)
  }
  cat("Scored", length(index$dataset_id), "datasets in", config$output_dir, "\n")
}
