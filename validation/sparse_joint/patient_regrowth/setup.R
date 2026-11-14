# Patient regrowth: a recurrence that grows from a few surviving cells, and
# growth that starts some time after the first sample (README.md in this
# directory).
#
# The truth populations are those of the development benchmark's development
# landscapes (benchmark/results/sparse_two_timepoint/truth), linked and read
# in place. From each patient's population at day 1440 this study draws the
# survivors, regrows them with the benchmark's ABM and samples both
# timepoints with its own seeds (base 99000). Stage code, scored nodes and
# metrics are those of ../short_interval, sourced unchanged; fits go to this
# study's own results/.

pr_dir <- if (exists("pr_script_dir")) pr_script_dir else getwd()
si_script_dir <- normalizePath(file.path(pr_dir, "..", "short_interval"))
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_script_dir, "stages.R"))

pr_default_output <- file.path(si_repo, "validation", "sparse_joint", "patient_regrowth", "results")
pr_benchmark_output <- file.path(si_repo, "benchmark", "results", "sparse_two_timepoint")

# `build`: the alfak2 build of the multi-timepoint work the fits use
# (validation/sparse_joint/multi_timepoint/results/lib/alfak2_<build>). The
# study was first run on ff13cf4, in results/; the user moved it to the build
# with the convergence check (2026-11-13), whose fits go to results/<build>/,
# beside the shared datasets (pr_stage_link()).
pr_config <- function(output_dir = NULL, build = "ff13cf4") {
  if (is.null(output_dir)) {
    output_dir <- if (build == "ff13cf4") pr_default_output else file.path(pr_default_output, build)
  }
  config <- stt_config(output_dir)
  config$build <- build
  config$landscapes <- config$landscapes[config$landscapes$split == "dev", ]
  config$library_root <- file.path(si_repo, "validation", "sparse_joint", "multi_timepoint", "results")
  config$t0 <- 1440
  # What the model is told: 360 days between the samples, 500 and 700 cells,
  # as in the benchmark's typical and bottleneck scenarios.
  config$interval <- 360
  config$depth <- c(500, 700)
  # `survivors`: cells drawn from the day-1440 population that regrow (0: the
  # whole population goes on growing). `growing`: share of the 360 days during
  # which the population grows; before that the survivors neither divide nor
  # die. `cohort`: patients 2-8 of each landscape are generated as well.
  scenario <- function(name, survivors, growing, cohort = FALSE) {
    data.frame(scenario = name, survivors = survivors, growing = growing, cohort = cohort,
               stringsAsFactors = FALSE)
  }
  config$scenarios <- rbind(
    scenario("none", 0, 1, cohort = TRUE),
    scenario("s2000", 2000, 1),
    scenario("s0200", 200, 1, cohort = TRUE),
    scenario("s0020", 20, 1, cohort = TRUE),
    scenario("s0200_d50", 200, 0.5),
    scenario("s0200_d75", 200, 0.25),
    scenario("s0200_d875", 200, 0.125),
    scenario("s0020_d875", 20, 0.125)
  )
  config$cohort_size <- 8L
  sparse <- function(...) {
    list(type = "sparse_joint", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    # ff13cf4 at its defaults, told one division per day: the benchmark's
    # populations divide about 0.9 times per day while they grow.
    base = sparse(),
    # The birth rate averaged over the 360 days, growing share x 1. Run on the
    # dormant scenarios only (elsewhere it is `base`).
    equivalent = c(sparse(), list(rate = "growing")),
    # Rate known only as a range, point predictions at its lowest rate.
    rates = list(type = "rates", build = build, birth_rate_range = c(0.1, 1), min_cn = 1L, max_cn = 8L,
                 arguments = list()),
    # Eight patients of a landscape fitted together, scored on patient 1. Run
    # on the cohort scenarios only.
    cohort = list(type = "cohort", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                  arguments = list())
  )
  config
}

pr_dataset_id <- function(landscape_id, patient, scenario) {
  sprintf("%s_P%02d_%s", landscape_id, patient, scenario)
}

# Patient 1 of every landscape for every scenario (the scored datasets), and
# patients 2-8 for the cohort scenarios. One sampling seed per landscape and
# patient, shared by its scenarios: the first sample is then the same cells in
# every scenario, and the survivors of the same size are the same cells.
pr_dataset_plan <- function(config) {
  plan <- list()
  for (i in seq_len(nrow(config$landscapes))) {
    landscape <- config$landscapes[i, ]
    for (s in seq_len(nrow(config$scenarios))) {
      scenario <- config$scenarios[s, ]
      patients <- if (scenario$cohort) seq_len(config$cohort_size) else 1L
      for (patient in patients) {
        plan[[length(plan) + 1L]] <- list(
          landscape = landscape, patient = patient, scenario = scenario,
          seed = as.integer(99000000L + i * 1000L + patient * 10L)
        )
      }
    }
  }
  plan
}

# One dataset. The ABM is restarted from the survivors (or the whole day-1440
# population) with the same seed in every scenario, so a dormant scenario's
# second population is its immediate-regrowth scenario's population at the end
# of the shorter growth time.
pr_make_dataset <- function(config, truth, scenario, seed) {
  dataset_id <- pr_dataset_id(truth$landscape_id, truth$patient, scenario$scenario)
  path <- stt_dataset_path(config, dataset_id)
  if (file.exists(path)) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  population_t0 <- stt_snapshot(truth, config$t0)
  start <- population_t0
  if (scenario$survivors > 0) {
    set.seed(seed + 1L)
    drawn <- stats::rmultinom(1L, scenario$survivors, population_t0 / sum(population_t0))[, 1L]
    names(drawn) <- names(population_t0)
    start <- drawn[drawn > 0]
  }
  growth_days <- config$interval * scenario$growing
  population_t1 <- stt_run_abm(start, truth$centroids, truth$wavelength, config,
                               record_days = growth_days, seed = seed + 1L, beta = truth$beta)[[1L]]
  counts <- stt_sample_counts(population_t0, population_t1, config$depth[[1L]], config$depth[[2L]],
                              config$interval, seed)
  share <- start / sum(start)
  dataset <- list(
    dataset_id = dataset_id, landscape_id = truth$landscape_id, patient = truth$patient,
    split = truth$split, wavelength_class = truth$wavelength_class,
    scenario = scenario$scenario, t0 = config$t0, interval = config$interval,
    depth_t0 = config$depth[[1L]], depth_t1 = config$depth[[2L]],
    survivors = scenario$survivors, growing = scenario$growing, growth_days = growth_days,
    founder_karyotypes = length(start), founder_top_share = max(share),
    population_t1 = sum(population_t1),
    bottleneck = scenario$survivors, shift = 0, seed = seed, beta = truth$beta,
    counts = counts,
    truth = list(centroids = truth$centroids, wavelength = truth$wavelength,
                 abm_delta_t = truth$abm_delta_t,
                 frequency_t0 = population_t0 / sum(population_t0),
                 frequency_t1 = population_t1 / sum(population_t1))
  )
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(dataset, temporary)
  file.rename(temporary, path)
  path
}

# Link the benchmark's truth populations (read in place). For a build other
# than ff13cf4, link instead the datasets, scored nodes, index and `zero`
# fits of results/, which do not depend on the build, by relative links that
# hold on HPC too.
pr_stage_link <- function(config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!identical(config$build, "ff13cf4")) {
    dir.create(file.path(config$output_dir, "fits"), showWarnings = FALSE)
    links <- c(datasets = "../datasets", nodes = "../nodes", dataset_index.csv = "../dataset_index.csv",
               "fits/zero" = "../../fits/zero")
    for (name in names(links)) {
      path <- file.path(config$output_dir, name)
      if (!file.exists(path)) file.symlink(links[[name]], path)
      if (!file.exists(path)) stop("Link ", path, " -> ", links[[name]], " is broken.", call. = FALSE)
    }
    cat("Linked the shared inputs of", pr_default_output, "into", config$output_dir, "\n")
    return(invisible())
  }
  link <- file.path(config$output_dir, "truth")
  if (!file.exists(link)) file.symlink(file.path(pr_benchmark_output, "truth"), link)
  missing <- unlist(lapply(config$landscapes$landscape_id, function(id) {
    paths <- stt_truth_path(config, id, seq_len(config$cohort_size))
    paths[!file.exists(paths)]
  }))
  if (length(missing)) stop("Missing truth populations: ", paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
  cat("Truth populations linked from", file.path(pr_benchmark_output, "truth"), "\n")
}

# Datasets, scored nodes of patient 1, and the index of the scored datasets.
pr_stage_datasets <- function(config, workers) {
  plan <- pr_dataset_plan(config)
  stt_load_abm(config, file.path(si_repo, "benchmark", "sparse_two_timepoint"))
  paths <- si_parallel(plan, function(item) {
    truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, item$patient))
    path <- pr_make_dataset(config, truth, item$scenario, item$seed)
    if (item$patient == 1L) stt_dataset_nodes_si(config, readRDS(path))
    path
  }, workers)
  scored <- vapply(plan, function(item) item$patient == 1L, logical(1))
  index <- do.call(rbind, lapply(paths[scored], function(path) {
    dataset <- readRDS(path)
    data.frame(
      dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
      wavelength_class = dataset$wavelength_class, scenario = dataset$scenario,
      survivors = dataset$survivors, growing = dataset$growing,
      founder_karyotypes = dataset$founder_karyotypes, founder_top_share = dataset$founder_top_share,
      population_t1 = dataset$population_t1,
      t(stt_count_profile(dataset$counts)), stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(index, si_index_path(config), row.names = FALSE)
  cat("Datasets:", length(paths), "(", nrow(index), "scored ) index:", si_index_path(config), "\n")
}

# Fitted prediction in the columns the scoring reads.
pr_prediction_frame <- function(prediction) {
  data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
             fitness_sd = prediction$relative_sd,
             model_sd = if (is.null(prediction$relative_sd_model)) NA_real_ else prediction$relative_sd_model,
             shell = if (is.null(prediction$shell)) NA_real_ else prediction$shell, stringsAsFactors = FALSE)
}

pr_counts <- function(dataset) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  counts
}

pr_fit_sparse_joint <- function(config, dataset, nodes, method) {
  rate <- if (identical(method$rate, "growing")) method$birth_rate * dataset$growing else method$birth_rate
  fit <- do.call(alfak2::fit_sparse_joint, c(
    list(pr_counts(dataset), dt = dataset$interval, beta = dataset$beta, birth_rate = rate,
         min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments
  ))
  d <- fit$diagnostics
  # `optimisation`: how the scale optimisation was checked (builds from cv1).
  list(status = "ok", prediction = pr_prediction_frame(stats::predict(fit, nodes$karyotype)),
       majors = fit$summary$karyotype[fit$summary$major], scales = fit$trend$scales, birth_rate = rate,
       diagnostics = d[intersect(c("convergence", "objective", "finite_uncertainty", "n_states", "n_major",
                                   "dropped_karyotypes", "dropped_cells", "fitted_birth_rate", "history",
                                   "optimisation"), names(d))],
       history_profile = d$history_profile, domain_coverage = d$domain_coverage)
}

pr_fit_rates <- function(config, dataset, nodes, method) {
  fit <- do.call(alfak2::fit_sparse_joint_rates, c(
    list(pr_counts(dataset), dt = dataset$interval, beta = dataset$beta, min_cn = method$min_cn,
         max_cn = method$max_cn, birth_rate_range = method$birth_rate_range),
    method$arguments
  ))
  prediction <- stats::predict(fit, nodes$karyotype)
  lowest <- fit$fits[[1L]]
  list(status = "ok", prediction = pr_prediction_frame(prediction),
       majors = lowest$summary$karyotype[lowest$summary$major], rates = fit$rates,
       rank_agreement = attr(prediction, "rank_agreement"),
       diagnostics = list(convergence = max(vapply(fit$fits, function(f) f$diagnostics$convergence, numeric(1))),
                          optimisation = lapply(fit$fits, function(f) f$diagnostics$optimisation)))
}

pr_fit_cohort <- function(config, dataset, nodes, method) {
  others <- pr_dataset_id(dataset$landscape_id, setdiff(seq_len(config$cohort_size), dataset$patient),
                          dataset$scenario)
  members <- c(list(dataset), lapply(others, function(id) readRDS(stt_dataset_path(config, id))))
  samples <- lapply(members, pr_counts)
  names(samples) <- vapply(members, function(m) m$dataset_id, character(1))
  fit <- do.call(alfak2::fit_sparse_joint_cohort, c(
    list(samples, dt = vapply(members, function(m) m$interval, numeric(1)), beta = dataset$beta,
         birth_rate = method$birth_rate, min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments
  ))
  prediction <- stats::predict(fit, nodes$karyotype, sample = 1L)
  first <- fit$samples[[1L]]
  d <- fit$diagnostics
  list(status = "ok", prediction = pr_prediction_frame(prediction),
       majors = if (is.null(first$major)) NULL else first$karyotype[first$major],
       scales = fit$trend$scales,
       diagnostics = list(convergence = d$convergence, objective = d$objective, trend = d$trend,
                          trend_objective = d$trend_objective, n_states = d$n_states, n_nodes = d$n_nodes,
                          optimisation = d$optimisation, trend_convergence = d$trend_convergence))
}

# Replaces si_run_method() of the short-interval study (which also runs a
# cross-validation this study does not need).
si_run_method <- function(config, dataset, nodes, method) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    switch(method$type,
      zero = list(status = "ok", prediction = data.frame(karyotype = nodes$karyotype, fitness = 0,
                                                         fitness_sd = NA_real_, stringsAsFactors = FALSE)),
      sparse_joint = pr_fit_sparse_joint(config, dataset, nodes, method),
      rates = pr_fit_rates(config, dataset, nodes, method),
      cohort = pr_fit_cohort(config, dataset, nodes, method),
      stop("Unknown method type: ", method$type, call. = FALSE)
    ),
    error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e),
                             error_class = class(e)[[1L]])
  )
  result$seconds <- proc.time()[["elapsed"]] - started
  result
}

# Scores. The truth is per calendar day: the rate while growing times the
# growing share, which is what a fit told 360 days estimates. Skill, slope,
# rank correlation and coverage then compare across dormant shares; RMSE
# itself shrinks with the share and is compared within a scenario only.
pr_stage_score <- function(config, workers) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  methods <- intersect(names(config$methods), list.files(file.path(config$output_dir, "fits")))
  scored <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- stt_dataset_nodes_si(config, dataset)
    truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth) * dataset$growing
    rows <- lapply(methods, function(method) {
      path <- si_fit_path(config, dataset_id, method)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      key <- data.frame(dataset_id = dataset_id, method = method, status = fit$status,
                        seconds = fit$seconds, stringsAsFactors = FALSE)
      if (!identical(fit$status, "ok")) {
        key$detail <- fit$detail
        return(list(distance = NULL, fit = key))
      }
      score <- stt_score_predictions(nodes, truth, fit$prediction, config$scoring$reference_min_count)
      # Majors of the fit against the other observed karyotypes: where a
      # founder effect would show.
      table <- score$nodes
      observed <- table$distance == 0L
      major <- table$karyotype %in% fit$majors
      error <- table$pred_relative - table$truth_relative
      d <- fit$diagnostics
      fit_row <- cbind(key, data.frame(
        rmse_major = if (any(observed & major)) sqrt(mean(error[observed & major]^2)) else NA_real_,
        rmse_major_zero = if (any(observed & major)) sqrt(mean(table$truth_relative[observed & major]^2)) else NA_real_,
        n_major_scored = sum(observed & major),
        convergence = if (is.null(d$convergence)) NA_real_ else d$convergence,
        history = if (is.null(d$history)) NA_real_ else d$history,
        n_states = if (is.null(d$n_states)) NA_real_ else d$n_states,
        trend = if (is.null(d$trend)) NA_character_ else d$trend,
        rank_agreement = if (is.null(fit$rank_agreement)) NA_real_ else fit$rank_agreement,
        stringsAsFactors = FALSE))
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

# The model is invariant to rescaling time and rates together: a fit told
# 360 days at the average rate equals one told the growth time at the rate
# while growing, divided by the ratio. Checked on one dataset.
pr_stage_symmetry <- function(config, dataset_id) {
  .libPaths(c(si_library(config$build, config$library_root), .libPaths()))
  dataset <- readRDS(stt_dataset_path(config, dataset_id))
  nodes <- stt_dataset_nodes_si(config, dataset)
  fit_at <- function(dt, rate) {
    fit <- alfak2::fit_sparse_joint(pr_counts(dataset), dt = dt, beta = dataset$beta, birth_rate = rate,
                                    min_cn = 1L, max_cn = 8L)
    stats::predict(fit, nodes$karyotype)
  }
  told <- fit_at(dataset$interval, dataset$growing)
  growth <- fit_at(dataset$growth_days, 1)
  ratio <- dataset$interval / dataset$growth_days
  difference <- told$relative_fitness * ratio - growth$relative_fitness
  out <- c(max_abs_difference = max(abs(difference)), max_abs_relative = max(abs(growth$relative_fitness)),
           sd_ratio = stats::median(told$relative_sd * ratio / growth$relative_sd))
  print(out)
  saveRDS(list(dataset_id = dataset_id, summary = out, told = told, growth = growth),
          file.path(config$output_dir, "symmetry_check.rds"))
  invisible(out)
}
