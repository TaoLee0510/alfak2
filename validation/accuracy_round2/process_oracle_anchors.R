#!/usr/bin/env Rscript
# Oracle-anchor-only structure diagnostic. This is not a fit using count data,
# and is not a deployable competitor. Condition the frozen v1 Gaussian prior on
# exact log birth rates at actually observed nodes; no pseudo-observation SD.

process_condition_gaussian <- function(mean, covariance, observed, log_truth) {
  unknown <- setdiff(seq_along(mean), observed)
  out_mean <- as.numeric(mean)
  out_covariance <- matrix(0, length(mean), length(mean))
  out_mean[observed] <- log_truth[observed]
  if (length(unknown)) {
    Koo <- covariance[observed, observed, drop = FALSE]
    Kuo <- covariance[unknown, observed, drop = FALSE]
    regression <- t(solve(Koo, t(Kuo)))
    out_mean[unknown] <- mean[unknown] + as.numeric(regression %*% (log_truth[observed] - mean[observed]))
    out_covariance[unknown, unknown] <- covariance[unknown, unknown, drop = FALSE] -
      regression %*% covariance[observed, unknown, drop = FALSE]
  }
  out_covariance <- (out_covariance + t(out_covariance)) / 2
  stopifnot(max(abs(out_mean[observed] - log_truth[observed])) < 1e-12,
            min(eigen(out_covariance, symmetric = TRUE, only.values = TRUE)$values) > -1e-9)
  list(mean = out_mean, covariance = out_covariance)
}

process_pair_support <- function(prior, covariance, observed) {
  pair <- which(grepl('^interaction_', colnames(prior$design)))
  other <- setdiff(seq_len(ncol(prior$design)), pair)
  design_rank <- function(x) if (!ncol(x) || !nrow(x)) 0L else qr(x)$rank
  if (!length(pair)) return(list(pair_count = 0L, observed_pair_rank = 0L,
    incremental_pair_rank = 0L, nullspace_dimension = 0L,
    posterior_pair_sd_fraction_min = NA_real_, posterior_pair_sd_fraction_max = NA_real_,
    hidden_pair_log_variance = numeric(nrow(prior$design))))
  Xo <- prior$design[observed, pair, drop = FALSE]
  eigen <- eigen(crossprod(Xo), symmetric = TRUE)
  nonzero <- eigen$values > 1e-10 * max(1, max(eigen$values))
  nullspace <- eigen$vectors[, !nonzero, drop = FALSE]
  pair_variance <- prior$coef_sd[pair]^2
  cross <- sweep(Xo, 2L, pair_variance, '*')
  posterior <- diag(pair_variance, length(pair)) -
    crossprod(cross, solve(covariance[observed, observed, drop = FALSE], cross))
  fraction <- sqrt(pmax(0, diag(posterior)) / pair_variance)
  if (ncol(nullspace)) stopifnot(max(abs(crossprod(nullspace, posterior %*% nullspace) -
    diag(pair_variance[[1L]], ncol(nullspace)))) < 1e-8)
  hidden_variance <- if (ncol(nullspace))
    pair_variance[[1L]] * rowSums((prior$design[, pair, drop = FALSE] %*% nullspace)^2) else numeric(nrow(prior$design))
  list(pair_count = length(pair), observed_pair_rank = sum(nonzero),
    incremental_pair_rank = design_rank(prior$design[observed, , drop = FALSE]) -
      design_rank(prior$design[observed, other, drop = FALSE]),
    nullspace_dimension = sum(!nonzero),
    posterior_pair_sd_fraction_min = min(fraction), posterior_pair_sd_fraction_max = max(fraction),
    hidden_pair_log_variance = hidden_variance)
}

run_process_anchor_oracle <- function(data_dir, output_dir, frozen_source, repo,
                                      joint_overrides = list(),
                                      model_label = 'oracle_observed_truth_prior_only') {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prior_environment <- new.env(parent = globalenv())
  sys.source(file.path(frozen_source, 'R/data-validate.R'), prior_environment)
  # Equivalent label formatter avoids loading a native objective: this diagnostic
  # needs only the frozen prior constructor and exact Gaussian conditioning.
  prior_environment$format_karyotypes <- function(k) apply(k, 1L, paste, collapse = '.')
  sys.source(file.path(frozen_source, 'R/joint-landscape.R'), prior_environment)
  scoring <- new.env(parent = globalenv())
  sys.source(file.path(repo, 'validation/accuracy_upgrade/metrics.R'), scoring)
  manifest <- read.csv(file.path(data_dir, 'dataset_manifest.csv'), stringsAsFactors = FALSE)
  # ABM datasets only: each oracle knows fitness at this dataset's observed
  # nodes, not at states seen solely in the paired deterministic count sample.
  manifest <- manifest[manifest$mode == 'finite_population', , drop = FALSE]
  metrics <- native_metrics <- nodes <- status <- list()
  for (i in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[[i]])
    graph <- list(karyotypes = dataset$target_karyotypes, labels = dataset$target_labels,
                  min_cn = dataset$config$min_cn, max_cn = dataset$config$max_cn)
    prior_control <- utils::modifyList(dataset$config$joint_control, joint_overrides)
    prior <- do.call(prior_environment$build_joint_landscape_prior,
      c(list(graph = graph, forward_model = 'division_process'), prior_control))
    X <- sweep(prior$design, 2L, prior$coef_sd, '*')
    K <- tcrossprod(X) + solve(as.matrix(prior$precision))
    observed <- match(rownames(dataset$counts), graph$labels)
    support <- process_pair_support(prior, K, observed)
    third_order_feature <- apply(sweep(graph$karyotypes, 2L, prior$diagnostics$reference, '-'), 1L, prod)
    posterior <- process_condition_gaussian(prior$offset, K, observed, log(dataset$truth))
    mean_birth <- exp(posterior$mean + .5 * diag(posterior$covariance))
    covariance_birth <- tcrossprod(mean_birth) * expm1(posterior$covariance)
    dimnames(covariance_birth) <- dimnames(posterior$covariance) <- list(graph$labels, graph$labels)
    summary <- data.frame(karyotype = graph$labels, fitness_mean = mean_birth,
                           fitness_sd = sqrt(pmax(0, diag(covariance_birth))))
    synthetic_prediction <- list(global = list(summary = summary, fitness_covariance = covariance_birth),
      local = list(summary = summary, latent_fitness_mean = posterior$mean,
        latent_fitness_covariance = posterior$covariance,
        diagnostics = list(forward_model = 'division_process')))
    model <- model_label
    prediction <- scoring$extract_target_predictions(synthetic_prediction, dataset, model)
    # Exact anchors and the t0 reference are known by construction. Remove only
    # exp(log(x))/summation roundoff from their degenerate oracle intervals.
    prediction$pred_raw[observed] <- dataset$truth[observed]
    prediction$pred_relative[observed] <- dataset$truth_relative[observed]
    prediction$sd_raw[observed] <- prediction$sd_relative[observed] <- 0
    for (level in dataset$config$interval_levels) {
      suffix <- sprintf('%02d', round(100 * level))
      prediction[[paste0('relative_low_', suffix)]][observed] <- dataset$truth_relative[observed]
      prediction[[paste0('relative_high_', suffix)]][observed] <- dataset$truth_relative[observed]
    }
    prediction$is_oracle_anchor <- seq_along(graph$labels) %in% observed
    prediction$unobserved_pair_nullspace_log_variance <- support$hidden_pair_log_variance
    prediction$native_distance <- dataset$native_target_distance
    prediction$oracle_information <- 'true_log_birth_at_actual_observed_nodes_only'
    score <- scoring$score_dataset(prediction, dataset, 0,
      'diagnostic_only_exact_prior_conditioning', 'analytic_conditional_gaussian')
    native_dataset <- dataset
    native_dataset$target_distance <- dataset$native_target_distance
    native_prediction <- prediction
    native_prediction$distance <- dataset$native_target_distance
    native_score <- scoring$score_dataset(native_prediction, native_dataset, 0,
      'diagnostic_only_exact_prior_conditioning', 'analytic_conditional_gaussian')
    native_metrics[[i]] <- native_score
    metrics[[i]] <- score; nodes[[i]] <- prediction
    unknown <- setdiff(seq_along(graph$labels), observed)
    latent_error <- log(dataset$truth) - posterior$mean
    conditional_energy <- if (length(unknown)) as.numeric(crossprod(latent_error[unknown],
      solve(posterior$covariance[unknown, unknown, drop = FALSE], latent_error[unknown]))) / length(unknown) else NA_real_
    basis_coefficients <- as.numeric(qr.solve(prior$design, log(dataset$truth) - prior$offset))
    component_z <- abs(basis_coefficients / prior$coef_sd)
    maximum_z <- function(pattern) max(c(0, component_z[grepl(pattern, colnames(prior$design))]))
    status[[i]] <- data.frame(dataset_id = dataset$dataset_id, landscape_id = dataset$landscape_id,
      family = dataset$family, model = model, n_observed_anchors = length(observed),
      n_target = length(graph$labels), anchor_max_absolute_error = max(abs(mean_birth[observed] - dataset$truth[observed])),
      conditional_covariance_min_eigenvalue = min(eigen(posterior$covariance, symmetric = TRUE, only.values = TRUE)$values),
      conditional_mahalanobis_per_unknown = conditional_energy,
      truth_log_birth_basis_projection_rmse = sqrt(mean((log(dataset$truth) - prior$offset -
        as.numeric(prior$design %*% basis_coefficients))^2)),
      full_truth_OLS_main_max_prior_sd_units = maximum_z('^dosage_'),
      full_truth_OLS_curvature_max_prior_sd_units = maximum_z('^curvature_'),
      full_truth_OLS_interaction_max_prior_sd_units = maximum_z('^interaction_'),
      pair_columns = support$pair_count, observed_pair_rank = support$observed_pair_rank,
      incremental_pair_rank_given_other_basis = support$incremental_pair_rank,
      pair_nullspace_dimension = support$nullspace_dimension,
      posterior_pair_sd_fraction_min = support$posterior_pair_sd_fraction_min,
      posterior_pair_sd_fraction_max = support$posterior_pair_sd_fraction_max,
      third_order_nonzero_observed_nodes = sum(abs(third_order_feature[observed]) > 0),
      third_order_feature_unseen = all(abs(third_order_feature[observed]) == 0),
      diagnostic_only = TRUE)
    saveRDS(list(nodes = prediction, metrics = score, native_metrics = native_score, status = status[[i]],
                 latent_mean = posterior$mean, latent_covariance = posterior$covariance,
                 prior_covariance = K),
      file.path(output_dir, paste0(dataset$dataset_id, '_anchor_oracle.rds')))
  }
  write.csv(do.call(rbind, metrics), file.path(output_dir, 'oracle_anchor_metrics.csv'), row.names = FALSE)
  write.csv(do.call(rbind, native_metrics), file.path(output_dir, 'oracle_anchor_native_metrics.csv'), row.names = FALSE)
  write.csv(do.call(rbind, nodes), file.path(output_dir, 'oracle_anchor_nodes.csv'), row.names = FALSE)
  write.csv(do.call(rbind, status), file.path(output_dir, 'oracle_anchor_status.csv'), row.names = FALSE)
  saveRDS(list(diagnostic = 'oracle prior-only, no count likelihood or additional flow data',
    joint_overrides = joint_overrides, model_label = model_label,
    frozen_prior_source = frozen_source,
    source_md5 = tools::md5sum(c(file.path(frozen_source, 'R/joint-landscape.R'),
      file.path(repo, 'validation/accuracy_upgrade/metrics.R'),
      file.path(repo, 'validation/accuracy_round2/process_oracle_anchors.R'))),
    data_md5 = tools::md5sum(manifest$dataset_file)), file.path(output_dir, 'oracle_anchor_provenance.rds'))
  invisible(do.call(rbind, metrics))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  repo <- if (length(args) >= 4L) args[[4L]] else '/Users/4482173/Documents/GitHub/alfak2'
  data_dir <- if (length(args)) args[[1L]] else '/private/tmp/alfak2_iteration2_20261104_0518/process-data'
  output_dir <- if (length(args) >= 2L) args[[2L]] else '/private/tmp/alfak2_iteration2_20261104_0518/process-oracle'
  source <- if (length(args) >= 3L) args[[3L]] else '/private/tmp/alfak2_iteration2_20261104_0518/v1-source'
  result <- run_process_anchor_oracle(data_dir, output_dir, source, repo)
  cat(sprintf('Completed oracle-anchor-only structure diagnostics: %d metric rows\n', nrow(result)))
}
