#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
input <- round2_arg("matrices", file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "node_matrices.rds"))
manifest_file <- round2_arg("manifest", file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "test_dataset_manifest.csv"))
output <- round2_arg("output", file.path(script_dir, "results", "old_data_diagnostics"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(input)
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
stopifnot(identical(x$dataset_ids, manifest$dataset_id))
mapping <- c(original_baseline = "baseline", v1_joint = "division_v1", endpoint_ablation = "endpoint_v1",
              original_baseline_calibrated = "baseline_calibrated", v1_joint_calibrated = "division_v1_calibrated",
              endpoint_ablation_calibrated = "endpoint_v1_calibrated")
founder <- which(vapply(strsplit(x$target_labels, ".", fixed = TRUE), function(k) all(as.integer(k) == 2L), logical(1)))
stopifnot(length(founder) == 1L)
karyotypes <- do.call(rbind, lapply(strsplit(x$target_labels, ".", fixed = TRUE), as.integer))
distance_matrix <- as.matrix(stats::dist(karyotypes, method = "manhattan"))
edges <- which(distance_matrix == 1 & upper.tri(distance_matrix), arr.ind = TRUE)
safe_cor <- function(a, b, method) {
  valid <- is.finite(a) & is.finite(b)
  if (sum(valid) < 5L || stats::sd(a[valid]) < 1e-12 || stats::sd(b[valid]) < 1e-12) return(NA_real_)
  suppressWarnings(stats::cor(a[valid], b[valid], method = method))
}
case_rows <- decomposition <- slopes <- list()
for (model in names(mapping)) {
  z <- x$models[[mapping[[model]]]]
  error <- z$pred_relative - x$truth_relative
  calibrated <- grepl("_calibrated$", model)
  low <- if (calibrated) z$calibrated_low_95 else z$relative_low_95
  high <- if (calibrated) z$calibrated_high_95 else z$relative_high_95
  for (landscape in unique(manifest$landscape_id)) {
    columns <- which(manifest$landscape_id == landscape)
    reference_column <- columns[grepl("_A1_S1$", manifest$dataset_id[columns])]
    stopifnot(length(columns) == 4L, length(reference_column) == 1L)
    reference_distance <- x$reference_support_distance[, reference_column]
    family <- manifest$family[columns[1L]]
    for (column in columns) for (distance_view in c("actual_support", "fixed_A1S1_support")) {
      distance <- if (distance_view == "actual_support") x$reference_support_distance[, column] else reference_distance
      global_offset <- mean(error[, column])
      edge_distance <- pmax(distance[edges[, 1L]], distance[edges[, 2L]])
      for (d in c(-1L, sort(unique(distance)))) {
        ii <- if (d == -1L) seq_len(nrow(error)) else which(distance == d)
        e <- error[ii, column]; sd <- z$sd_relative[ii, column]
        centered_error <- e - global_offset
        shape_mse <- round2_mean(centered_error^2)
        cross_term <- 2 * global_offset * mean(centered_error)
        stopifnot(abs(mean(e^2) - shape_mse - global_offset^2 - cross_term) < 1e-12)
        edge_subset <- if (d == -1L) seq_len(nrow(edges)) else which(edge_distance == d)
        gradient_error <- error[edges[edge_subset, 2L], column] - error[edges[edge_subset, 1L], column]
        truth_sd <- if (length(ii) > 1L) stats::sd(x$truth_relative[ii, column]) else NA_real_
        pred_sd <- if (length(ii) > 1L) stats::sd(z$pred_relative[ii, column]) else NA_real_
        amplitude <- if (is.finite(truth_sd) && truth_sd > 0) pred_sd / truth_sd else NA_real_
        valid <- is.finite(e); interval <- valid & is.finite(low[ii, column]) & is.finite(high[ii, column])
        r <- -e / sd
        r[sd == 0 & abs(e) < 1e-12] <- 0
        lower_penalty <- 40 * pmax(low[ii, column] - x$truth_relative[ii, column], 0)
        upper_penalty <- 40 * pmax(x$truth_relative[ii, column] - high[ii, column], 0)
        width <- high[ii, column] - low[ii, column]
        case_rows[[length(case_rows) + 1L]] <- data.frame(model = model, family = family,
          landscape_id = landscape, dataset_id = manifest$dataset_id[column], abm_seed = manifest$abm_seed[column],
          sampling_seed = manifest$sampling_seed[column], distance_view = distance_view, distance = d,
          n_target = length(ii), prediction_fraction = mean(valid), interval_fraction = mean(interval),
          rmse = sqrt(round2_mean(e^2)), mse = round2_mean(e^2), mae = round2_mean(abs(e)), signed_bias = round2_mean(e),
          spatially_centered_rmse = sqrt(shape_mse), spatially_centered_mse = shape_mse,
          global_offset_squared = global_offset^2, global_offset_cross_term = cross_term,
          fixed_founder_rmse_diagnostic = sqrt(mean((e - error[founder, column])^2)),
          edge_gradient_rmse = sqrt(round2_mean(gradient_error^2)), n_edges = length(edge_subset),
          pearson = safe_cor(z$pred_relative[ii, column], x$truth_relative[ii, column], "pearson"),
          spearman = safe_cor(z$pred_relative[ii, column], x$truth_relative[ii, column], "spearman"),
          amplitude_ratio = amplitude, variance_ratio = amplitude^2, amplitude_fidelity_error = abs(amplitude - 1),
          amplitude_below_0p1 = if (is.finite(amplitude)) as.numeric(amplitude < 0.1) else NA_real_,
          absolute_error_q90 = round2_quantile(abs(e), 0.9), absolute_error_q99 = round2_quantile(abs(e), 0.99),
          coverage95 = if (any(interval)) mean(x$truth_relative[ii, column][interval] >= low[ii, column][interval] & x$truth_relative[ii, column][interval] <= high[ii, column][interval]) else NA_real_,
          lower_tail95 = if (any(interval)) mean(x$truth_relative[ii, column][interval] < low[ii, column][interval]) else NA_real_,
          upper_tail95 = if (any(interval)) mean(x$truth_relative[ii, column][interval] > high[ii, column][interval]) else NA_real_,
          width95 = round2_mean(width), width95_q90 = round2_quantile(width, 0.9),
          lower_penalty95 = round2_mean(lower_penalty), upper_penalty95 = round2_mean(upper_penalty),
          interval_score95 = round2_mean(width + lower_penalty + upper_penalty),
          standardized_q025 = round2_quantile(r, 0.025), standardized_q975 = round2_quantile(r, 0.975),
          standardized_abs_gt_1p96 = mean(abs(r) > 1.959963984540054, na.rm = TRUE))
      }
    }
    if (calibrated) next
    for (gauge_view in c("original_data_reference", "fixed_founder_reference_diagnostic")) {
      e <- error[, columns, drop = FALSE]
      if (gauge_view == "fixed_founder_reference_diagnostic") e <- sweep(e, 2L, e[founder, ], "-")
      grand <- rowMeans(e)
      abm <- unique(manifest$abm_seed[columns])
      group_columns <- lapply(abm, function(a) which(manifest$abm_seed[columns] == a))
      abm_mean <- vapply(group_columns, function(ii) rowMeans(e[, ii, drop = FALSE]), numeric(nrow(e)))
      within <- e
      for (a in seq_along(group_columns)) within[, group_columns[[a]]] <- sweep(e[, group_columns[[a]], drop = FALSE], 1L, abm_mean[, a], "-")
      between <- sweep(abm_mean, 1L, grand, "-")
      bias2 <- grand^2
      between_empirical <- rowMeans(between^2)
      within_empirical <- rowMeans(within^2)
      mse <- rowMeans(e^2)
      stopifnot(max(abs(mse - bias2 - between_empirical - within_empirical)) < 1e-12)
      ms_within <- rowSums(within^2) / 2
      ms_between <- 2 * rowSums(between^2)
      sigma_sampling <- ms_within
      sigma_abm <- (ms_between - ms_within) / 2
      bias2_corrected <- bias2 - sigma_abm / 2 - sigma_sampling / 4
      for (d in sort(unique(reference_distance))) {
        ii <- which(reference_distance == d)
        decomposition[[length(decomposition) + 1L]] <- data.frame(model = model, family = family,
          landscape_id = landscape, distance_view = "fixed_A1S1_support", distance = d, gauge_view = gauge_view,
          n_target = length(ii), n_abm = 2L, n_sampling_per_abm = 2L,
          empirical_mse = mean(mse[ii]), empirical_bias2_proxy = mean(bias2[ii]),
          empirical_between_abm = mean(between_empirical[ii]), empirical_within_sampling = mean(within_empirical[ii]),
          sampling_variance_anova = mean(sigma_sampling[ii]), abm_variance_anova_raw = mean(sigma_abm[ii]),
          bias2_anova_corrected_raw = mean(bias2_corrected[ii]),
          negative_abm_component_fraction = mean(sigma_abm[ii] < 0),
          negative_corrected_bias2_fraction = mean(bias2_corrected[ii] < 0))
      }
    }
    truth_mean <- rowMeans(x$truth_relative[, columns, drop = FALSE])
    pred_mean <- rowMeans(z$pred_relative[, columns, drop = FALSE])
    fit <- stats::lm(pred_mean ~ truth_mean)
    slopes[[length(slopes) + 1L]] <- data.frame(model = model, family = family, landscape_id = landscape,
      diagnostic_amplitude_slope = unname(stats::coef(fit)[2L]), diagnostic_intercept = unname(stats::coef(fit)[1L]),
      prediction_rms = sqrt(mean(pred_mean^2)), truth_rms = sqrt(mean(truth_mean^2)))
  }
}
case <- do.call(rbind, case_rows); decomposition <- do.call(rbind, decomposition); slopes <- do.call(rbind, slopes)
measure <- setdiff(names(case), c("model", "family", "landscape_id", "dataset_id", "abm_seed", "sampling_seed", "distance_view", "distance"))
land <- aggregate(case[, measure], case[, c("model", "family", "landscape_id", "distance_view", "distance")], round2_mean)
family <- aggregate(land[, measure], land[, c("model", "family", "distance_view", "distance")], round2_mean)
decomp_measure <- setdiff(names(decomposition), c("model", "family", "landscape_id", "distance_view", "distance", "gauge_view"))
decomp_family <- aggregate(decomposition[, decomp_measure], decomposition[, c("model", "family", "distance_view", "distance", "gauge_view")], round2_mean)
round2_write_csv(case, file.path(output, "case_metrics.csv"))
round2_write_csv(land, file.path(output, "landscape_metrics.csv"))
round2_write_csv(family, file.path(output, "family_distance_metrics.csv"))
round2_write_csv(decomposition, file.path(output, "landscape_bias_variance.csv"))
round2_write_csv(decomp_family, file.path(output, "family_bias_variance.csv"))
round2_write_csv(slopes, file.path(output, "landscape_amplitude_diagnostic.csv"))
round2_write_csv(data.frame(file = c(normalizePath(input), normalizePath(manifest_file)),
                           sha256 = unname(round2_sha256(c(input, manifest_file)))), file.path(output, "input_sha256.csv"))
round2_write_csv(data.frame(model = c("original_baseline", "v1_joint", "current_default_two_layer"),
                            status = c("old_exact_matrices_available", "old_exact_matrices_available", "prepared_not_run_no_surrogate")),
                 file.path(output, "control_status.csv"))
cat("Old-data attribution complete; no new fitting or sealed data generation.\n")
