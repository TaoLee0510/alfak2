#!/usr/bin/env Rscript
# Development diagnostic using already fitted, properly normalized component
# marginal likelihoods. No truth enters model weights.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "metrics.R"))
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: analyze_endpoint_bma.R SCALE_GRID_RESULTS DATA_DIRECTORY")
results_dir <- args[[1L]]; data_dir <- args[[2L]]
manifest <- utils::read.csv(file.path(data_dir, "dataset_manifest.csv"), stringsAsFactors = FALSE)
all_weights <- list()
for (include_two in c(FALSE, TRUE)) {
  scales <- c("0p05", "0p2", "0p5", "1", if (include_two) "2")
  model <- if (include_two) "endpoint_evidence_bma5_sensitivity" else "endpoint_evidence_bma4"
  metrics <- nodes <- statuses <- list()
  for (di in seq_len(nrow(manifest))) {
    dataset <- readRDS(manifest$dataset_file[di])
    files <- file.path(results_dir, paste0("joint_endpoint_s", scales, "_", dataset$dataset_id, ".rds"))
    component <- lapply(files, readRDS)
    objectives <- vapply(component, function(x) x$local_diagnostics$objective, numeric(1))
    if (any(!is.finite(objectives))) stop("Non-finite component evidence")
    weights <- exp(-(objectives - min(objectives)))
    weights <- weights / sum(weights)
    weight_row <- data.frame(model = model, dataset_id = dataset$dataset_id,
                             family = dataset$family, landscape_id = dataset$landscape_id,
                             component_scale = scales, objective = objectives, weight = weights)
    all_weights[[length(all_weights) + 1L]] <- weight_row
    z <- component[[1L]]$nodes
    stopifnot(all(vapply(component, function(x) identical(x$nodes$karyotype, z$karyotype), logical(1))))
    means <- vapply(component, function(x) x$nodes$pred_relative, numeric(nrow(z)))
    sds <- vapply(component, function(x) x$nodes$sd_relative, numeric(nrow(z)))
    if (any(!is.finite(means)) || any(!is.finite(sds)) || any(sds <= 0)) stop("Invalid component contrast moments")
    z$model <- model
    z$pred_relative <- as.numeric(means %*% weights)
    z$sd_relative <- sqrt(pmax(0, as.numeric((sds^2 + means^2) %*% weights) - z$pred_relative^2))
    z$pred_raw <- as.numeric(vapply(component, function(x) x$nodes$pred_raw, numeric(nrow(z))) %*% weights)
    raw_second <- vapply(component, function(x) x$nodes$sd_raw^2 + x$nodes$pred_raw^2, numeric(nrow(z)))
    z$sd_raw <- sqrt(pmax(0, as.numeric(raw_second %*% weights) - z$pred_raw^2))
    z$covariance_source <- "mixture_of_verified_component_contrast_marginals"
    z$interval_source <- "deterministic_Gaussian_mixture_contrast_quantiles"
    for (level in dataset$config$interval_levels) {
      suffix <- sprintf("%02d", round(100 * level))
      quantile <- function(i, probability) {
        lo <- min(means[i, ] - 12 * sds[i, ])
        hi <- max(means[i, ] + 12 * sds[i, ])
        stats::uniroot(function(q) sum(weights * stats::pnorm(q, means[i, ], sds[i, ])) - probability,
                       c(lo, hi), tol = 1e-9)$root
      }
      z[[paste0("relative_low_", suffix)]] <- vapply(seq_len(nrow(z)), quantile, numeric(1), probability = (1 - level) / 2)
      z[[paste0("relative_high_", suffix)]] <- vapply(seq_len(nrow(z)), quantile, numeric(1), probability = (1 + level) / 2)
    }
    runtime <- sum(vapply(component, function(x) x$status$runtime_seconds, numeric(1)))
    scored <- score_dataset(z, dataset, runtime, "ok", "component_Laplace_mixture")
    status <- data.frame(dataset_id = dataset$dataset_id, model = model, fit_status = "ok",
                          covariance_status = "component_Laplace_mixture", runtime_seconds = runtime,
                          warnings = "Endpoint-only development mixture; no new inference, no heldout count mixture available")
    metrics[[di]] <- scored; nodes[[di]] <- z; statuses[[di]] <- status
    saveRDS(list(metrics = scored, nodes = z, status = status, weights = weights,
                  component_scales = scales, component_objectives = objectives,
                  component_cache_keys = lapply(component, `[[`, "cache_key"),
                  weighting = "uniform model prior; exp(-Laplace marginal negative log likelihood); no truth"),
             file.path(results_dir, paste0(model, "_", dataset$dataset_id, ".rds")))
  }
  utils::write.csv(do.call(rbind, metrics), file.path(results_dir, paste0(model, "_metrics.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, nodes), file.path(results_dir, paste0(model, "_nodes.csv")), row.names = FALSE)
  utils::write.csv(do.call(rbind, statuses), file.path(results_dir, paste0(model, "_status.csv")), row.names = FALSE)
}
weights <- do.call(rbind, all_weights)
utils::write.csv(weights, file.path(results_dir, "evidence_bma_weights.csv"), row.names = FALSE)
groups <- split(weights, interaction(weights$model, weights$family, weights$component_scale, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) data.frame(
  model = x$model[1], family = x$family[1], component_scale = x$component_scale[1],
  mean_weight = mean(x$weight), min_weight = min(x$weight),
  q10_weight = unname(stats::quantile(x$weight, 0.1)), q90_weight = unname(stats::quantile(x$weight, 0.9)),
  max_weight = max(x$weight))))
utils::write.csv(summary, file.path(results_dir, "evidence_bma_weight_summary.csv"), row.names = FALSE)
artifact_dir <- file.path(script_dir, "results", "scale_grid_review")
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(weights, file.path(artifact_dir, "evidence_bma_weights.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(artifact_dir, "evidence_bma_weight_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
