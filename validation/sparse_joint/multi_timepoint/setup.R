# Stage 2 of the multi-timepoint work: one population sampled at several
# times, with passaging bottlenecks of different sizes between the time points
# (README.md in this directory).
#
# The data come from the development simulator of
# benchmark/sparse_two_timepoint, whose files are sourced unchanged, with
# this study's own landscapes (seed base 97000; held-out 98000, sealed),
# seeds and output directory. Stage code is that of ../short_interval.

mt_dir <- if (exists("mt_script_dir")) mt_script_dir else getwd()
si_script_dir <- normalizePath(file.path(mt_dir, "..", "short_interval"))
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_script_dir, "stages.R"))

mt_default_output <- file.path(si_repo, "validation", "sparse_joint", "multi_timepoint", "results")

mt_config <- function(output_dir = NULL, set = c("series", "heldout")) {
  set <- match.arg(set)
  if (is.null(output_dir)) output_dir <- file.path(mt_default_output, set)
  config <- stt_config(output_dir)
  classes <- data.frame(wavelength_class = c("w0.4", "w0.8", "w1.2", "w1.6"), wavelength = c(0.4, 0.8, 1.2, 1.6),
                        stringsAsFactors = FALSE)
  grid <- classes[rep(seq_len(nrow(classes)), each = 3L), ]
  index <- seq_len(nrow(grid))
  # Seed bases 71000 and 82000 belong to the development benchmark, 93000 to
  # 96000 to the short-interval and lattice-residual studies, 99000 to the
  # patient regrowth study.
  config$landscapes <- data.frame(
    landscape_id = sprintf(if (set == "series") "mt_L%02d" else "mh_L%02d", index), split = "mt",
    wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
    landscape_seed = (if (set == "series") 97000L else 98000L) + index, stringsAsFactors = FALSE
  )
  config$set <- set
  config$library_root <- mt_default_output
  config$patients_per_landscape <- 1L
  # The population grows on the development benchmark's settings (5e4 cells
  # at the start, thinned to 1% whenever it exceeds 2e6) until the first
  # time point at day 1440. From then on it is thinned to `bottleneck` cells,
  # nb / bottleneck growth per passage.
  config$abm$end_day <- 1440
  config$t0 <- 1440
  config$nb <- config$abm$max_pop
  config$bottlenecks <- c(2e4, 2e3, 2e2)
  # Days between time points: about one passage at the largest bottleneck (the
  # population divides about 0.9 times per day and grows 100-fold per
  # passage), and about ten.
  config$spacings <- c(6, 60)
  config$max_samples <- 8L
  config$samples <- c(4L, 8L)
  config$depths <- c(200, 1000)
  # Build cv3: the multi-timepoint model (mt1) with the convergence check of
  # stage C, whose fits are those of mt1 (the user stopped stage M and kept the
  # check, 2026-11-14).
  series <- function(columns, ..., build = "cv3") {
    list(type = "series", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L, columns = columns,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    # Every time point, and the first and last only (two timepoints: exactly
    # ff13cf4).
    all = series("all"),
    ends = series("ends"),
    # Every time point, the ones between the first and the last for the majors
    # only. `majors_any` (cv4): majors chosen from every time point; `majors1`
    # (cv5): majors and domain from the first and last time point, the time
    # points in between at full weight; `majors` (cv5): the same with their weight the
    # inverse of the majors' dispersion in a first fit. The pilot found the
    # fit on all time points worse than first+last, on karyotypes seen only in between
    # and on clones that became majors in between.
    majors_any = series("all", intermediate = "majors", build = "cv4"),
    majors1 = series("all", intermediate = "majors", intermediate_weight = 1, build = "cv5"),
    majors = series("all", intermediate = "majors", build = "cv5"),
    # alfakR dev2 f93b2d2 on every time point, as in the development benchmark
    # (minobs 20, 45 bootstraps), with the window's bottleneck as n0.
    alfakR = list(type = "alfakR", build = NULL, library = "alfakR_f93b2d2", minobs = 20L, nboot = 45L,
                  nn_prior = "empirical")
  )
  config
}

mt_dataset_id <- function(landscape_id, bottleneck, spacing, depth, samples) {
  sprintf("%s_b%05d_s%02d_n%04d_k%d", landscape_id, as.integer(bottleneck), as.integer(spacing), as.integer(depth),
          as.integer(samples))
}

# One ABM window per landscape, bottleneck and spacing, recorded at the
# first time point and every `spacing` days after it; each recorded population
# sampled at each depth. The four-time-point datasets are the first four
# columns of the eight-time-point ones (the same cells), so the two are paired.
mt_window_plan <- function(config) {
  plan <- list()
  for (i in seq_len(nrow(config$landscapes))) for (b in seq_along(config$bottlenecks)) {
    for (s in seq_along(config$spacings)) {
      plan[[length(plan) + 1L]] <- list(
        landscape = config$landscapes[i, ], bottleneck = config$bottlenecks[[b]], spacing = config$spacings[[s]],
        seed = as.integer(config$landscapes$landscape_seed[[i]] * 1000L + b * 10L + s)
      )
    }
  }
  plan
}

mt_make_window <- function(config, truth, item) {
  window_config <- config
  window_config$abm$culling_survival <- item$bottleneck / config$nb
  days <- config$t0 + item$spacing * (seq_len(config$max_samples) - 1L)
  start <- stt_snapshot(truth, config$t0)
  later <- stt_run_abm(start, truth$centroids, truth$wavelength, window_config,
                       record_days = days[-1L] - config$t0, seed = item$seed, beta = truth$beta)
  c(list(start), unname(later))
}

mt_sample_columns <- function(populations, depth, seed) {
  set.seed(seed)
  labels <- unique(unlist(lapply(populations, names)))
  counts <- vapply(populations, function(p) {
    out <- stats::setNames(numeric(length(labels)), labels)
    out[names(p)] <- stats::rmultinom(1L, depth, p / sum(p))[, 1L]
    out
  }, numeric(length(labels)))
  counts <- matrix(counts, nrow = length(labels), dimnames = list(labels, NULL))
  counts[rowSums(counts) > 0, , drop = FALSE]
}

mt_stage_datasets <- function(config, workers) {
  plan <- mt_window_plan(config)
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  rows <- si_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
    ids <- character()
    populations <- NULL
    for (d in seq_along(config$depths)) for (k in config$samples) {
      depth <- config$depths[[d]]
      id <- mt_dataset_id(truth$landscape_id, item$bottleneck, item$spacing, depth, k)
      ids <- c(ids, id)
      path <- stt_dataset_path(config, id)
      if (file.exists(path)) next
      if (is.null(populations)) populations <- mt_make_window(config, truth, item)
      counts <- mt_sample_columns(populations, depth, item$seed * 10L + d)
      counts <- counts[, seq_len(k), drop = FALSE]
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      counts <- counts[order(rowSums(counts), decreasing = TRUE), , drop = FALSE]
      times <- item$spacing * (seq_len(k) - 1L)
      colnames(counts) <- as.character(times)
      dataset <- list(
        dataset_id = id, landscape_id = truth$landscape_id, patient = 1L, split = truth$split,
        wavelength_class = truth$wavelength_class, scenario = sprintf("b%g_s%d_n%d_k%d", item$bottleneck,
                                                                      item$spacing, depth, k),
        t0 = config$t0, times = times, interval = max(times), samples = k, spacing = item$spacing,
        bottleneck = item$bottleneck, depth = depth, seed = item$seed * 10L + d, beta = truth$beta,
        population_sizes = vapply(populations[seq_len(k)], sum, numeric(1)),
        counts = counts,
        truth = list(centroids = truth$centroids, wavelength = truth$wavelength, abm_delta_t = truth$abm_delta_t,
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
      bottleneck = dataset$bottleneck, spacing = dataset$spacing, depth = dataset$depth, samples = dataset$samples,
      # Drift against sampling noise for a fitness difference, as estimated in
      # the plan: K^2 N_s / (5 n0), K passages one apart.
      drift_index = dataset$samples^2 * dataset$depth / (5 * dataset$bottleneck),
      n_karyotypes = nrow(counts), shared_all = sum(rowSums(counts > 0) == ncol(counts)),
      top_fraction = max(rowSums(counts)) / sum(counts),
      min_population = min(dataset$population_sizes), stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(index, si_index_path(config), row.names = FALSE)
  cat("Datasets:", nrow(index), " index:", si_index_path(config), "\n")
}

# Scored nodes from all time points: observed means seen at any time point, and
# the reference weights use the first time point and the sum of the later ones, so
# that every arm is scored on the same karyotypes against the same reference.
mt_dataset_nodes <- function(config, dataset) {
  path <- file.path(config$output_dir, "nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  counts <- dataset$counts
  collapsed <- cbind(counts[, 1L], rowSums(counts[, -1L, drop = FALSE]))
  rownames(collapsed) <- rownames(counts)
  scoring <- config$scoring
  nodes <- stt_scoring_nodes(collapsed, scoring$min_cn, scoring$max_cn, scoring$d1_cap, scoring$d2_cap, dataset$seed)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(nodes, temporary)
  file.rename(temporary, path)
  nodes
}
stt_dataset_nodes_si <- mt_dataset_nodes

mt_counts <- function(dataset, columns) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  if (identical(columns, "ends")) counts <- counts[, c(1L, ncol(counts)), drop = FALSE]
  counts[rowSums(counts) > 0, , drop = FALSE]
}

# Dispersion of the majors' counts about the fit's expectation at the later
# time points (the expectation at the first time point is free): mean of (y - E)^2 / E.
# About one under the model's Poisson sampling; passaging drift raises it.
mt_dispersion <- function(fit) {
  s <- fit$summary[fit$summary$major, , drop = FALSE]
  later <- grep("^count_t[1-9][0-9]*$", names(s), value = TRUE)
  y <- as.matrix(s[later])
  e <- as.matrix(s[sub("^count", "expected", later)])
  # Scale the expectation at each time point to its observed total over the majors,
  # so that the depth offsets do not count as dispersion.
  e <- sweep(e, 2L, colSums(y) / colSums(e), "*")
  mean((y - e)^2 / pmax(e, 1e-8))
}

mt_fit_series <- function(config, dataset, nodes, method) {
  counts <- mt_counts(dataset, method$columns)
  times <- if (identical(method$columns, "ends")) dataset$times[c(1L, length(dataset$times))] else dataset$times
  fit <- do.call(alfak2::fit_sparse_joint, c(
    list(counts, times = times, beta = dataset$beta, birth_rate = method$birth_rate, min_cn = method$min_cn,
         max_cn = method$max_cn),
    method$arguments
  ))
  prediction <- stats::predict(fit, nodes$karyotype)
  d <- fit$diagnostics
  list(status = "ok",
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                               shell = prediction$shell, stringsAsFactors = FALSE),
       majors = fit$summary$karyotype[fit$summary$major], scales = fit$trend$scales,
       dispersion = mt_dispersion(fit),
       diagnostics = d[intersect(c("convergence", "objective", "finite_uncertainty", "n_states", "n_major",
                                   "dropped_karyotypes", "dropped_cells", "fitted_birth_rate", "history",
                                   "optimisation", "times"), names(d))],
       history_profile = d$history_profile, domain_coverage = d$domain_coverage)
}

mt_fit_alfakR <- function(config, dataset, nodes, method) {
  work <- tempfile("alfakR_")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  counts <- mt_counts(dataset, "all")
  yi <- list(x = data.frame(counts, check.names = FALSE), dt = 1)
  suppressWarnings(suppressMessages(utils::capture.output(
    alfakR::alfak(yi = yi, outdir = work, passage_times = dataset$times, minobs = method$minobs,
                  nboot = method$nboot, n0 = dataset$bottleneck, nb = config$nb, pm = dataset$beta,
                  nn_prior = method$nn_prior, landscape_data_output = TRUE)
  )))
  landscape <- readRDS(file.path(work, "landscape.Rds"))
  index <- match(nodes$karyotype, as.character(landscape$k))
  fitness <- as.numeric(landscape$mean)[index]
  fitness_sd <- as.numeric(landscape$sd)[index]
  missing <- which(!is.finite(fitness))
  if (length(missing)) {
    kriging <- readRDS(file.path(work, "landscape_data.Rds"))$mean
    if (!is.null(kriging)) {
      fitness[missing] <- as.numeric(stats::predict(kriging, stt_parse_karyotypes(nodes$karyotype[missing])))
      fitness_sd[missing] <- NA_real_
    }
  }
  list(status = "ok", n_frequent = sum(landscape$fq),
       prediction = data.frame(karyotype = nodes$karyotype, fitness = fitness, fitness_sd = fitness_sd,
                               stringsAsFactors = FALSE))
}

si_run_method <- function(config, dataset, nodes, method) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    switch(method$type,
      zero = list(status = "ok", prediction = data.frame(karyotype = nodes$karyotype, fitness = 0,
                                                         fitness_sd = NA_real_, stringsAsFactors = FALSE)),
      series = mt_fit_series(config, dataset, nodes, method),
      alfakR = mt_fit_alfakR(config, dataset, nodes, method),
      stop("Unknown method type: ", method$type, call. = FALSE)
    ),
    error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e),
                             error_class = class(e)[[1L]])
  )
  result$seconds <- proc.time()[["elapsed"]] - started
  result
}

# alfakR runs from its own library; the alfak2 arms from an alfak2 build.
mt_stage_fit <- function(config, workers, method_names, datasets) {
  libraries <- unique(unlist(lapply(config$methods[method_names], `[[`, "library")))
  if (length(libraries)) {
    if (length(method_names) != 1L) stop("Run alfakR on its own.", call. = FALSE)
    .libPaths(c(file.path(config$library_root, "lib", libraries), .libPaths()))
    loadNamespace("alfakR")
    cat("alfakR from", dirname(find.package("alfakR")), "\n")
  }
  si_stage_fit(config, workers, method_names, datasets)
}

mt_stage_score <- function(config, workers) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  methods <- intersect(names(config$methods), list.files(file.path(config$output_dir, "fits")))
  scored <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- mt_dataset_nodes(config, dataset)
    truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
    rows <- lapply(methods, function(method) {
      path <- si_fit_path(config, dataset_id, method)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      key <- data.frame(dataset_id = dataset_id, method = method, status = fit$status, seconds = fit$seconds,
                        stringsAsFactors = FALSE)
      if (!identical(fit$status, "ok")) {
        key$detail <- fit$detail
        return(list(distance = NULL, fit = key))
      }
      score <- stt_score_predictions(nodes, truth, fit$prediction, config$scoring$reference_min_count)
      d <- fit$diagnostics
      fit_row <- cbind(key, data.frame(
        dispersion = if (is.null(fit$dispersion)) NA_real_ else fit$dispersion,
        convergence = if (is.null(d$convergence)) NA_real_ else d$convergence,
        history = if (is.null(d$history)) NA_real_ else d$history,
        n_states = if (is.null(d$n_states)) NA_real_ else d$n_states,
        n_major = if (is.null(d$n_major)) NA_real_ else d$n_major,
        n_frequent = if (is.null(fit$n_frequent)) NA_real_ else fit$n_frequent, stringsAsFactors = FALSE))
      list(distance = cbind(key, score$metrics), fit = fit_row)
    })
    rows <- Filter(Negate(is.null), rows)
    list(distance = si_bind(lapply(rows, `[[`, "distance")), fit = si_bind(lapply(rows, `[[`, "fit")))
  }, workers)
  distance <- merge(index, si_bind(lapply(scored, `[[`, "distance")), by = "dataset_id")
  fits <- merge(index, si_bind(lapply(scored, `[[`, "fit")), by = "dataset_id")
  utils::write.csv(distance, file.path(config$output_dir, "metrics_distance.csv"), row.names = FALSE)
  utils::write.csv(fits, file.path(config$output_dir, "metrics_fit.csv"), row.names = FALSE)
  cat("Scored", nrow(fits), "fits in", config$output_dir, "\n")
}
