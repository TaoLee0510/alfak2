# Method wrappers. Each takes a dataset and the scored nodes and returns
# list(status, prediction = data.frame(karyotype, fitness, fitness_sd),
# seconds, detail). `fitness` may be on any additive scale: scoring centres it.

stt_run_method <- function(method, dataset, nodes, work_dir) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    switch(method$type,
      zero = stt_method_zero(dataset, nodes),
      alfakR = stt_method_alfakR(dataset, nodes, method, work_dir),
      alfak2_two_layer = stt_method_alfak2_two_layer(dataset, nodes, method),
      alfak2_sparse_joint = stt_method_alfak2_sparse_joint(dataset, nodes, method),
      alfak2_sparse_joint_cohort = stt_method_alfak2_sparse_joint_cohort(dataset, nodes, method),
      alfak2_sparse_joint_rates = stt_method_alfak2_sparse_joint_rates(dataset, nodes, method),
      stop("Unknown method type: ", method$type, call. = FALSE)
    ),
    error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e))
  )
  result$seconds <- proc.time()[["elapsed"]] - started
  result
}

# The no-information answer: every karyotype is as fit as the reference.
stt_method_zero <- function(dataset, nodes) {
  list(status = "ok", detail = NA_character_,
       prediction = data.frame(karyotype = nodes$karyotype, fitness = 0, fitness_sd = NA_real_,
                               stringsAsFactors = FALSE))
}

# Original ALFA-K. Its landscape table covers the frequent karyotypes and two
# shells around them; any other scored node is predicted from the saved
# Kriging mean surface, without an uncertainty.
stt_method_alfakR <- function(dataset, nodes, method, work_dir) {
  out_dir <- file.path(work_dir, "alfakR")
  unlink(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  yi <- list(x = data.frame(dataset$counts, check.names = FALSE), dt = 1)
  suppressWarnings(suppressMessages(utils::capture.output(
    alfakR::alfak(
      yi = yi, outdir = out_dir, passage_times = c(0, dataset$interval),
      minobs = method$minobs, nboot = method$nboot, n0 = method$n0, nb = method$nb,
      pm = dataset$beta, nn_prior = method$nn_prior, landscape_data_output = TRUE
    )
  )))
  landscape <- readRDS(file.path(out_dir, "landscape.Rds"))
  index <- match(nodes$karyotype, as.character(landscape$k))
  fitness <- as.numeric(landscape$mean)[index]
  fitness_sd <- as.numeric(landscape$sd)[index]
  missing <- which(!is.finite(fitness))
  n_extrapolated <- 0L
  if (length(missing)) {
    kriging <- readRDS(file.path(out_dir, "landscape_data.Rds"))$mean
    if (!is.null(kriging)) {
      fitness[missing] <- as.numeric(stats::predict(kriging, stt_parse_karyotypes(nodes$karyotype[missing])))
      fitness_sd[missing] <- NA_real_
      n_extrapolated <- length(missing)
    }
  }
  list(status = "ok",
       detail = sprintf("frequent=%d kriging_only=%d", sum(landscape$fq), n_extrapolated),
       prediction = data.frame(karyotype = nodes$karyotype, fitness = fitness,
                               fitness_sd = fitness_sd, stringsAsFactors = FALSE))
}

# Current alfak2 two-layer pipeline, configured as in the June 22-chromosome
# benchmark and the real-sample workflow: local fit on observed karyotypes,
# graph Gaussian extrapolation with normalised edge weights. The global graph
# is the observed karyotypes plus the scored nodes.
#
# `method$graph_hyperparameters`, when given, fixes (lambda_l, lambda_e,
# sigma_obs). Without it the package default applies: leave-one-anchor-out
# selection over a 27-point grid, which refactorises the graph once per anchor
# and grid point and takes hours per dataset at 22 chromosomes.
stt_method_alfak2_two_layer <- function(dataset, nodes, method) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  karyotypes <- stt_parse_karyotypes(nodes$karyotype)
  targets <- nodes$karyotype[nodes$distance > 0L]
  arguments <- list(
    counts, dt = dataset$interval, beta = dataset$beta,
    transition_kernel = "exact",
    local_shell_depth = method$local_shell_depth,
    min_cn = 0, max_cn = as.integer(max(karyotypes)),
    max_nodes = method$max_nodes,
    graph_edge_weight = "normalized",
    extrapolation_method = "graph_gaussian_baseline",
    target_karyotypes = targets,
    observation_weight_mode = "fractional_count"
  )
  if (!is.null(method$graph_hyperparameters)) {
    fixed <- method$graph_hyperparameters
    arguments$lambda_l_grid <- fixed[["lambda_l"]]
    arguments$lambda_e_grid <- fixed[["lambda_e"]]
    arguments$sigma_obs_grid <- fixed[["sigma_obs"]]
  }
  fit <- suppressWarnings(do.call(alfak2::fit_alfak2, arguments))
  summary <- alfak2::summarize_alfak2(fit, layer = "global")
  index <- match(nodes$karyotype, as.character(summary$karyotype))
  list(status = "ok",
       detail = sprintf("local_covariance=%s nodes=%d",
                        fit$local$diagnostics$covariance_status, nrow(summary)),
       prediction = data.frame(karyotype = nodes$karyotype,
                               fitness = as.numeric(summary$fitness_mean)[index],
                               fitness_sd = as.numeric(summary$fitness_sd)[index],
                               stringsAsFactors = FALSE))
}

# Sparse joint count model. The baseline birth rate is an input; the benchmark
# supplies one value for all datasets and never the simulated truth.
stt_method_alfak2_sparse_joint <- function(dataset, nodes, method) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  # `method$arguments` holds any further arguments of the variant under test.
  fit <- do.call(alfak2::fit_sparse_joint, c(
    list(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
         min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments
  ))
  prediction <- stats::predict(fit, nodes$karyotype)
  diagnostics <- fit$diagnostics
  list(status = "ok",
       detail = sprintf("states=%d major=%d dropped=%d convergence=%d finite_sd=%s birth=%.2f",
                        diagnostics$n_states, diagnostics$n_major, diagnostics$dropped_karyotypes,
                        diagnostics$convergence, diagnostics$finite_uncertainty, diagnostics$fitted_birth_rate),
       # `model_sd` and `shell` are what the calibration is refitted from.
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                               shell = prediction$shell, stringsAsFactors = FALSE))
}

# Datasets fitted together with `dataset`, which comes first. "landscape":
# the other patients of its landscape. "other_landscapes": patient 1 of the
# following landscapes of the split, a cohort that does not share a landscape.
stt_cohort_members <- function(dataset, method) {
  others <- if (identical(method$members, "landscape")) {
    sprintf("%s_P%02d_%s", dataset$landscape_id, setdiff(seq_len(method$cohort_size), dataset$patient),
            dataset$scenario)
  } else {
    ids <- method$landscapes$landscape_id[method$landscapes$split == dataset$split]
    at <- match(dataset$landscape_id, ids)
    sprintf("%s_P01_%s", ids[(at + seq_len(method$cohort_size - 1L) - 1L) %% length(ids) + 1L], dataset$scenario)
  }
  paths <- file.path(method$dataset_dir, paste0(others, ".rds"))
  if (any(!file.exists(paths))) stop("Cohort datasets are missing: generate them with --cohort.", call. = FALSE)
  c(list(dataset), lapply(paths, readRDS))
}

# Cohort fit of the sparse joint model, predicted for the first member.
stt_method_alfak2_sparse_joint_cohort <- function(dataset, nodes, method) {
  members <- stt_cohort_members(dataset, method)
  samples <- lapply(members, function(member) {
    counts <- member$counts
    storage.mode(counts) <- "integer"
    counts
  })
  names(samples) <- vapply(members, function(member) member$dataset_id, character(1))
  fit <- do.call(alfak2::fit_sparse_joint_cohort, c(
    list(samples, dt = vapply(members, function(member) member$interval, numeric(1)), beta = dataset$beta,
         birth_rate = method$birth_rate, min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments
  ))
  prediction <- stats::predict(fit, nodes$karyotype, sample = 1L)
  diagnostics <- fit$diagnostics
  list(status = "ok",
       detail = sprintf("nodes=%d states=%d shared=%d convergence=%d finite_sd=%s deviation=%.4f trend=%s%s",
                        diagnostics$n_nodes, diagnostics$n_states, diagnostics$shared_nodes,
                        diagnostics$convergence, diagnostics$finite_uncertainty, fit$trend$scales[["deviation"]],
                        if (is.null(diagnostics$trend)) "shared" else diagnostics$trend,
                        if (is.null(diagnostics$trend_objective)) "" else
                          sprintf(" objective_shared=%.2f objective_sample=%.2f", diagnostics$trend_objective[["shared"]],
                                  diagnostics$trend_objective[["sample"]])),
       # `model_sd` and `shell` are what the calibration is refitted from.
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                               shell = prediction$shell, stringsAsFactors = FALSE))
}

# Sparse joint model with the baseline birth rate known only as a range: the
# prediction at the centre of the range, with an interval covering the range.
stt_method_alfak2_sparse_joint_rates <- function(dataset, nodes, method) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  fit <- alfak2::fit_sparse_joint_rates(counts, dt = dataset$interval, beta = dataset$beta,
                                        min_cn = method$min_cn, max_cn = method$max_cn,
                                        birth_rate_range = method$birth_rate_range)
  prediction <- stats::predict(fit, nodes$karyotype)
  list(status = "ok",
       detail = sprintf("rates=%s rank_agreement=%.3f", paste(signif(fit$rates, 3), collapse = "/"),
                        attr(prediction, "rank_agreement")),
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, stringsAsFactors = FALSE))
}
