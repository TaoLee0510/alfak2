#!/usr/bin/env Rscript
# Diagnostic information interventions on frozen v1; never a promotion test.
oracle_read_argument <- function(args, name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}

oracle_fingerprint <- function(x) {
  path <- tempfile("oracle_fingerprint_")
  on.exit(unlink(path))
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}

oracle_initial_groups <- function(dataset) {
  y <- matrix(0, length(dataset$target_labels), 2L)
  y[match(rownames(dataset$counts), dataset$target_labels), ] <- dataset$counts
  ifelse(y[, 1L] > 0, "t0_positive", ifelse(y[, 2L] > 0, "t1_only", "both_zero"))
}

oracle_point_diagnostics <- function(fit, dataset, arm, nodes, runtime) {
  index <- match(dataset$target_labels, fit$local$summary$karyotype)
  p0 <- fit$local$summary$pi0[index]
  mode_birth <- fit$local$parameter_mode$f[index]
  groups <- oracle_initial_groups(dataset)
  pi1 <- fit$local$summary$pi1[index]
  error <- nodes$pred_relative - nodes$truth_relative
  mass <- lapply(c("t0_positive", "t1_only", "both_zero"), function(group) {
    keep <- groups == group
    data.frame(dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
      family = dataset$family, sampling_replicate = dataset$sampling_replicate,
      arm = arm, initial_group = group, n_nodes = sum(keep),
      fitted_pi0_mass = sum(p0[keep]), truth_pi0_mass = sum(dataset$truth_pi0[keep]),
      fitted_pi1_mass = sum(pi1[keep]), truth_pi1_mass = sum(dataset$truth_pi1[keep]),
      relative_bias = if (any(keep)) mean(error[keep]) else NA_real_,
      relative_rmse = if (any(keep)) sqrt(mean(error[keep]^2)) else NA_real_)
  })
  list(mass = do.call(rbind, mass), baseline = data.frame(
    dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
    family = dataset$family, sampling_replicate = dataset$sampling_replicate, arm = arm,
    mode_pi0_weighted_birth = sum(p0 * mode_birth),
    true_initial_mean_birth = sum(dataset$truth_pi0 * dataset$truth),
    mode_predicted_log_growth = as.numeric(fit$local$posterior_predictive$log_total_growth),
    realized_log_growth = log(sum(dataset$latent_population_t1) / sum(dataset$latent_population_t0)),
    pi0_rmse = sqrt(mean((p0 - dataset$truth_pi0)^2)),
    pi1_rmse = sqrt(mean((pi1 - dataset$truth_pi1)^2)), runtime_seconds = runtime,
    stringsAsFactors = FALSE))
}

oracle_repeat_decomposition <- function(nodes) {
  # The first prespecified sampling replicate fixes a node's stratum for both
  # repeats. This avoids calling a changing support category a fixed estimand.
  reference <- nodes[nodes$arm == "count_only" & nodes$sampling_replicate == 1L,
                     c("landscape_id", "karyotype", "distance", "initial_group")]
  names(reference)[3:4] <- c("reference_distance", "reference_initial_group")
  paired <- merge(nodes, reference, by = c("landscape_id", "karyotype"), all.x = TRUE)
  paired$error <- paired$pred_relative - paired$truth_relative
  node_groups <- split(paired, interaction(paired$family, paired$landscape_id,
                                          paired$arm, paired$karyotype, drop = TRUE))
  moments <- lapply(node_groups, function(x) {
    e <- x$error[is.finite(x$error)]
    if (length(e) < 2L) return(NULL)
    data.frame(family = x$family[1L], landscape_id = x$landscape_id[1L], arm = x$arm[1L],
      karyotype = x$karyotype[1L], reference_distance = x$reference_distance[1L],
      reference_initial_group = x$reference_initial_group[1L], n_repeats = length(e),
      repeat_mean_error = mean(e), squared_repeat_mean_error = mean(e)^2,
      within_repeat_error_variance = mean((e - mean(e))^2), mse = mean(e^2))
  })
  moments <- do.call(rbind, moments)
  if (is.null(moments)) return(list(node_moments = data.frame(), distance = data.frame(), initial_group = data.frame()))
  summarize <- function(stratum) {
    formula <- stats::as.formula(paste(
      "cbind(repeat_mean_error, squared_repeat_mean_error, within_repeat_error_variance, mse) ~ family + arm +", stratum))
    ans <- stats::aggregate(formula, moments, mean)
    ans$rmse <- sqrt(ans$mse)
    ans$root_squared_repeat_mean_error <- sqrt(ans$squared_repeat_mean_error)
    ans$decomposition_residual <- ans$mse - ans$squared_repeat_mean_error - ans$within_repeat_error_variance
    ans
  }
  list(node_moments = moments, distance = summarize("reference_distance"),
       initial_group = summarize("reference_initial_group"))
}

oracle_main <- function(args = commandArgs(TRUE)) {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
  repo <- normalizePath(file.path(dirname(script), "..", ".."))
  run_root <- oracle_read_argument(args, "run-root", "/private/tmp/alfak2_iteration2_20261104_0518")
  data_dir <- oracle_read_argument(args, "data", file.path(run_root, "process-data"))
  output <- oracle_read_argument(args, "output", file.path(run_root, "oracle-results"))
  frozen <- oracle_read_argument(args, "v1-source", file.path(run_root, "v1-source"))
  dll <- oracle_read_argument(args, "dll", file.path(run_root, "current-lib/alfak2/libs/alfak2.so"))
  limit <- as.integer(oracle_read_argument(args, "limit", "2147483647"))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  options(alfak2.warning_log_dir = file.path(output, "warnings"))
  requireNamespace("Rcpp"); requireNamespace("TMB")
  binary <- dyn.load(normalizePath(dll))
  package <- new.env(parent = globalenv())
  routines <- getDLLRegisteredRoutines(binary)$`.Call`
  for (name in names(routines)) assign(name, routines[[name]], package)
  frozen_files <- sort(list.files(file.path(frozen, "R"), "[.]R$", full.names = TRUE))
  for (file in frozen_files) sys.source(file, package)
  # The only R departure from frozen v1 is the explicit oracle mechanism.
  sys.source(file.path(repo, "R/local-fit.R"), package)
  scoring_files <- file.path(repo, "validation/accuracy_upgrade", c("posterior_sampling.R", "metrics.R"))
  scoring <- new.env(parent = globalenv())
  for (file in scoring_files) sys.source(file, scoring)
  manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
  manifest <- head(manifest[manifest$mode == "finite_population", , drop = FALSE], limit)
  if (!nrow(manifest)) stop("No finite-population development datasets")
  hashes <- tools::md5sum(c(frozen_files, file.path(repo, "R/local-fit.R"), scoring_files, script, dll))
  saveRDS(list(source_hashes = hashes, data_hashes = tools::md5sum(manifest$dataset_file),
    session_info = utils::sessionInfo(), scope = "development_oracle_information_diagnostics_only",
    uncertainty = "uncalibrated_v1_joint_Laplace", population_log_growth_sd = 0.1),
    file.path(output, "oracle_provenance.rds"))
  all_nodes <- all_metrics <- all_status <- all_mass <- all_baseline <- list()
  for (i in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[i])
    if (!identical(dataset$split, "development")) stop("Oracle diagnostic only accepts development data")
    known_pi0 <- stats::setNames(dataset$truth_pi0, dataset$target_labels)
    population <- c(sum(dataset$latent_population_t0), sum(dataset$latent_population_t1))
    arms <- list(count_only = list(), known_pi0 = list(oracle_initial_probability = known_pi0),
      population_growth = list(population_sizes = population, population_log_growth_sd = 0.1),
      known_pi0_and_growth = list(oracle_initial_probability = known_pi0,
        population_sizes = population, population_log_growth_sd = 0.1))
    for (arm in names(arms)) {
      cache_file <- file.path(output, paste0("oracle_", dataset$dataset_id, "_", arm, ".rds"))
      key <- list(implementation = oracle_fingerprint(hashes),
                  dataset = unname(tools::md5sum(manifest$dataset_file[i])), arm = arm)
      result <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
      if (is.null(result) || !identical(result$cache_key, key)) {
        started <- proc.time()[["elapsed"]]
        call <- c(list(counts = dataset$counts, dt = dataset$config$dt, beta = dataset$config$beta,
          min_cn = dataset$config$min_cn, max_cn = dataset$config$max_cn,
          max_nodes = nrow(dataset$target_karyotypes), target_karyotypes = dataset$target_karyotypes,
          pipeline = "joint_count", forward_model = "division_process", inference_method = "laplace",
          sdreport_mode = "fitness_only", joint_control = list(baseline_sd = 0.5,
            main_sd = 0.07, curvature_sd = 0.04, interaction_sd = 0.04,
            residual_sd = 0.125, edge_sd = 0.25, collapse_gaussian = FALSE)), arms[[arm]])
        warnings <- character()
        fit <- tryCatch(withCallingHandlers(do.call(package$fit_alfak2, call), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        }), error = function(e) e)
        runtime <- proc.time()[["elapsed"]] - started
        if (inherits(fit, "error")) stop(sprintf("%s / %s failed: %s", dataset$dataset_id, arm, conditionMessage(fit)))
        if (!identical(fit$local$diagnostics$covariance_status, "TMB_sdreport"))
          stop(sprintf("Untrusted oracle covariance: %s / %s", dataset$dataset_id, arm))
        nodes <- scoring$extract_target_predictions(fit, dataset, arm)
        nodes$shared_distance <- nodes$distance
        nodes$distance <- as.numeric(dataset$native_target_distance)
        nodes$arm <- arm; nodes$sampling_replicate <- dataset$sampling_replicate
        nodes$initial_group <- oracle_initial_groups(dataset)
        native <- dataset; native$target_distance <- dataset$native_target_distance
        metrics <- scoring$score_dataset(nodes, native, runtime, "ok", fit$local$diagnostics$covariance_status)
        diagnostic <- oracle_point_diagnostics(fit, dataset, arm, nodes, runtime)
        result <- list(cache_key = key, nodes = nodes, metrics = metrics,
          mass = diagnostic$mass, baseline = diagnostic$baseline,
          status = data.frame(dataset_id = dataset$dataset_id, family = dataset$family,
            landscape_id = dataset$landscape_id, arm = arm, fit_status = "ok",
            covariance_status = fit$local$diagnostics$covariance_status,
            runtime_seconds = runtime, warnings = paste(unique(warnings), collapse = " | ")),
          diagnostics = fit$local$diagnostics)
        saveRDS(result, cache_file)
      }
      j <- length(all_nodes) + 1L
      all_nodes[[j]] <- result$nodes; all_metrics[[j]] <- result$metrics
      all_mass[[j]] <- result$mass; all_baseline[[j]] <- result$baseline
      all_status[[j]] <- result$status
      utils::write.csv(do.call(rbind, all_status), file.path(output, "oracle_fit_status.csv"), row.names = FALSE)
      cat(sprintf("oracle %d/%d %s %s %.2fs\n", i, nrow(manifest), dataset$dataset_id, arm,
                  result$status$runtime_seconds))
    }
  }
  nodes <- do.call(rbind, all_nodes); metrics <- do.call(rbind, all_metrics)
  mass <- do.call(rbind, all_mass); baseline <- do.call(rbind, all_baseline)
  utils::write.csv(nodes, file.path(output, "oracle_nodes.csv"), row.names = FALSE)
  utils::write.csv(metrics, file.path(output, "oracle_metrics.csv"), row.names = FALSE)
  utils::write.csv(mass, file.path(output, "oracle_initial_mass_groups.csv"), row.names = FALSE)
  utils::write.csv(baseline, file.path(output, "oracle_birth_and_growth.csv"), row.names = FALSE)
  family <- stats::aggregate(value ~ model + family + distance_group + metric, metrics, mean)
  utils::write.csv(family, file.path(output, "oracle_family_metrics.csv"), row.names = FALSE)
  decomposition <- oracle_repeat_decomposition(nodes)
  for (name in names(decomposition)) utils::write.csv(decomposition[[name]],
    file.path(output, paste0("oracle_repeat_", name, ".csv")), row.names = FALSE)
  saveRDS(manifest, file.path(output, "oracle_dataset_manifest.rds"))
  invisible(output)
}

if (sys.nframe() == 0L) oracle_main()
