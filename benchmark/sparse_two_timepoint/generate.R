# Ground-truth populations and two-timepoint count data.
#
# A "landscape" is one set of GRF centroids. A "patient" is one independent
# ABM trajectory on that landscape, started from the diploid founder with its
# own simulation seed, so patients of a landscape form a cohort that shares
# the fitness landscape but not the evolutionary history.

# Compile the simulator once, in the parent process, before forking workers.
stt_load_abm <- function(config, script_dir) {
  if (exists("stt_abm_cpp", mode = "function")) return(invisible(TRUE))
  cache <- file.path(config$output_dir, "cpp_cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(script_dir, "abm.cpp"), cacheDir = cache, env = globalenv())
  invisible(TRUE)
}

# Run the ABM from `initial` (named cell counts) and return the population at
# each requested day as a named count vector.
stt_run_abm <- function(initial, centroids, wavelength, config, record_days, seed,
                        beta = config$beta) {
  abm <- config$abm
  record_steps <- as.integer(round(sort(unique(record_days)) / abm$delta_t))
  raw <- stt_abm_cpp(
    initial_karyotypes = matrix(as.integer(stt_parse_karyotypes(names(initial))), nrow = length(initial)),
    initial_counts = as.numeric(initial),
    centroids = centroids, wavelength = wavelength, beta = beta,
    delta_t = abm$delta_t, n_steps = max(record_steps), record_steps = record_steps,
    max_population = abm$max_pop, culling_survival = abm$culling_survival,
    min_cn = abm$min_cn, max_cn = abm$max_cn, seed = seed
  )
  snapshots <- lapply(record_steps, function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
  names(snapshots) <- as.character(record_steps * abm$delta_t)
  snapshots
}

stt_truth_path <- function(config, landscape_id, patient) {
  file.path(config$output_dir, "truth", sprintf("%s_P%02d.rds", landscape_id, patient))
}

# One patient's trajectory. Every recorded snapshot is kept, so sampling
# windows can be chosen afterwards without resimulating.
stt_generate_truth <- function(config, landscape, patient, force = FALSE) {
  path <- stt_truth_path(config, landscape$landscape_id, patient)
  if (!force && file.exists(path)) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  founder <- rep(config$founder_cn, config$n_chr)
  set.seed(landscape$landscape_seed)
  centroids <- stt_generate_centroids(founder, config$n_waves, landscape$wavelength)
  abm <- config$abm
  initial <- stats::setNames(abm$pop_size, paste(founder, collapse = "."))
  snapshots <- stt_run_abm(initial, centroids, landscape$wavelength, config,
                           record_days = seq(0, abm$end_day, by = abm$record_days),
                           seed = landscape$landscape_seed * 100 + patient)
  truth <- list(
    landscape_id = landscape$landscape_id, patient = patient, split = landscape$split,
    wavelength_class = landscape$wavelength_class, wavelength = landscape$wavelength,
    centroids = centroids, beta = config$beta, abm_delta_t = abm$delta_t,
    snapshots = snapshots
  )
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(truth, temporary)
  file.rename(temporary, path)
  path
}

stt_snapshot <- function(truth, day) {
  population <- truth$snapshots[[as.character(day)]]
  if (is.null(population) || !length(population)) stop("No ABM snapshot at day ", day, ".", call. = FALSE)
  population
}

# Regrowth after a strong bottleneck: restart the ABM from `n_survivors` cells
# drawn from the day-`t0` population and run it for `interval` days.
stt_bottleneck_snapshot <- function(truth, config, t0, interval, n_survivors, seed) {
  start <- stt_snapshot(truth, t0)
  set.seed(seed)
  survivors <- stats::rmultinom(1L, n_survivors, start / sum(start))[, 1L]
  names(survivors) <- names(start)
  survivors <- survivors[survivors > 0]
  stt_run_abm(survivors, truth$centroids, truth$wavelength, config,
              record_days = interval, seed = seed, beta = truth$beta)[[1L]]
}

# Multinomial cell sampling at both timepoints. Rows are the karyotypes seen
# at least once; columns are named by day, with the first timepoint as day 0.
stt_sample_counts <- function(population_t0, population_t1, depth_t0, depth_t1, interval, seed) {
  set.seed(seed)
  labels <- union(names(population_t0), names(population_t1))
  draw <- function(population, depth) {
    out <- stats::setNames(numeric(length(labels)), labels)
    out[names(population)] <- stats::rmultinom(1L, depth, population / sum(population))[, 1L]
    out
  }
  counts <- cbind(draw(population_t0, depth_t0), draw(population_t1, depth_t1))
  colnames(counts) <- c("0", as.character(interval))
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  counts[order(rowSums(counts), decreasing = TRUE), , drop = FALSE]
}

stt_dataset_path <- function(config, dataset_id) {
  file.path(config$output_dir, "datasets", paste0(dataset_id, ".rds"))
}

# One two-timepoint dataset: the counts a method sees, and separately the
# population frequencies it is scored against.
stt_make_dataset <- function(config, truth, scenario, t0, seed, force = FALSE) {
  dataset_id <- sprintf("%s_P%02d_%s", truth$landscape_id, truth$patient, scenario$scenario)
  path <- stt_dataset_path(config, dataset_id)
  if (!force && file.exists(path)) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  population_t0 <- stt_snapshot(truth, t0)
  population_t1 <- if (isTRUE(scenario$bottleneck > 0)) {
    stt_bottleneck_snapshot(truth, config, t0, scenario$interval, scenario$bottleneck, seed + 1L)
  } else {
    stt_snapshot(truth, t0 + scenario$interval)
  }
  counts <- stt_sample_counts(population_t0, population_t1, scenario$depth_t0,
                              scenario$depth_t1, scenario$interval, seed)
  dataset <- list(
    dataset_id = dataset_id, landscape_id = truth$landscape_id, patient = truth$patient,
    split = truth$split, wavelength_class = truth$wavelength_class,
    scenario = scenario$scenario, t0 = t0, interval = scenario$interval,
    depth_t0 = scenario$depth_t0, depth_t1 = scenario$depth_t1,
    bottleneck = scenario$bottleneck, seed = seed, beta = truth$beta,
    counts = counts,
    truth = list(centroids = truth$centroids, wavelength = truth$wavelength,
                 abm_delta_t = truth$abm_delta_t,
                 frequency_t0 = population_t0 / sum(population_t0),
                 frequency_t1 = population_t1 / sum(population_t1))
  )
  saveRDS(dataset, path)
  path
}
