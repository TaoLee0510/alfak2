#!/usr/bin/env Rscript
# Descriptive development attribution. Bootstrap units are six independent
# landscapes, stratified by the three fixed truth families; samples are nested.
args <- commandArgs(TRUE)
base <- if (length(args)) args[[1L]] else '/private/tmp/alfak2_iteration2_20261104_0518'
fit_dir <- file.path(base, 'process-results')
oracle_dir <- file.path(base, 'process-oracle')
output <- file.path(base, 'process-review')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
metrics <- read.csv(file.path(fit_dir, 'matched_metrics.csv'), stringsAsFactors = FALSE)
status <- read.csv(file.path(fit_dir, 'matched_status.csv'), stringsAsFactors = FALSE)
oracle <- read.csv(file.path(oracle_dir, 'oracle_anchor_metrics.csv'), stringsAsFactors = FALSE)
important <- c('relative_rmse', 'relative_bias', 'coverage_95', 'interval_score_95',
                'interval_width_95', 'edge_gradient_rmse')
metrics <- metrics[metrics$metric %in% important & grepl('^d[0-9]+$', metrics$distance_group), ]
landscape <- aggregate(value ~ mode + family + landscape_id + distance_group + metric, metrics, mean)
family <- aggregate(value ~ mode + family + distance_group + metric, landscape, mean)
overall <- aggregate(value ~ mode + distance_group + metric, landscape, mean)
write.csv(landscape, file.path(output, 'process_landscape_metrics.csv'), row.names = FALSE)
write.csv(family, file.path(output, 'process_family_metrics.csv'), row.names = FALSE)
write.csv(overall, file.path(output, 'process_overall_metrics.csv'), row.names = FALSE)

paired_development <- function(first, second, first_name, second_name,
                                keys = c('pair_id', 'family', 'landscape_id', 'distance_group', 'metric')) {
  first <- first[, c(keys, 'value')]; names(first)[names(first) == 'value'] <- 'first'
  second <- second[, c(keys, 'value')]; names(second)[names(second) == 'value'] <- 'second'
  pairs <- merge(first, second, by = keys, all = TRUE)
  pairs$delta <- pairs$second - pairs$first
  groups <- split(pairs, interaction(pairs$distance_group, pairs$metric, drop = TRUE))
  set.seed(921104L)
  summaries <- lapply(groups, function(z) {
    valid <- is.finite(z$first) & is.finite(z$second)
    finite <- z[valid, , drop = FALSE]
    if (!nrow(finite)) return(NULL)
    land <- aggregate(cbind(first, second, delta) ~ family + landscape_id, finite, mean)
    boot <- replicate(1000L, {
      sample_rows <- unlist(lapply(split(seq_len(nrow(land)), land$family), function(index)
        sample(index, length(index), replace = TRUE)), use.names = FALSE)
      mean(land$delta[sample_rows])
    })
    data.frame(distance_group = finite$distance_group[[1L]], metric = finite$metric[[1L]],
      first_arm = first_name, second_arm = second_name,
      first_mean = mean(land$first), second_mean = mean(land$second),
      second_minus_first = mean(land$delta),
      descriptive_delta_low = unname(quantile(boot, .025)),
      descriptive_delta_high = unname(quantile(boot, .975)),
      independent_landscapes = nrow(land), valid_dataset_pairs = nrow(finite),
      incomplete_dataset_pairs = sum(!valid),
      interpretation = 'development attribution, not an independent candidate test')
  })
  list(pairs = pairs, summary = do.call(rbind, summaries))
}
process <- paired_development(subset(metrics, mode == 'deterministic_mean'),
  subset(metrics, mode == 'finite_population'), 'deterministic_mean', 'finite_population')
write.csv(process$pairs, file.path(output, 'process_paired_dataset_deltas.csv'), row.names = FALSE)
write.csv(process$summary, file.path(output, 'process_paired_landscape_summary.csv'), row.names = FALSE)
oracle <- oracle[oracle$metric %in% important & grepl('^d[0-9]+$', oracle$distance_group), ]
oracle_landscape <- aggregate(value ~ family + landscape_id + distance_group + metric, oracle, mean)
oracle_family <- aggregate(value ~ family + distance_group + metric, oracle_landscape, mean)
write.csv(oracle_landscape, file.path(output, 'oracle_landscape_metrics.csv'), row.names = FALSE)
write.csv(oracle_family, file.path(output, 'oracle_family_metrics.csv'), row.names = FALSE)
oracle_comparison <- paired_development(subset(metrics, mode == 'finite_population'), oracle,
  'count_fit_with_fixed_eta_prior', 'oracle_observed_truth_prior_only',
  keys = c('dataset_id', 'family', 'landscape_id', 'distance_group', 'metric'))
write.csv(oracle_comparison$summary, file.path(output, 'oracle_vs_fitted_diagnostic.csv'), row.names = FALSE)

summary_file <- file.path(output, 'PROCESS_ATTRIBUTION.txt')
capture.output({
  cat('Process-noise and observed-fitness-oracle attribution: DEVELOPMENT ONLY\n\n')
  cat('Six independently generated landscapes, two per smooth/rough/interaction family.\n')
  cat('Each has one fixed initial integer population, one finite Gillespie path and\n')
  cat('two count samples. The deterministic arm uses exp(A*t) applied to that same\n')
  cat('initial population. Both arms use common sampling seeds and exactly identical\n')
  cat('t0 counts/reference. The entire 27-node three-chromosome CN1-3 domain is fixed.\n')
  cat('beta=.005, dt=1.5, N0=1200, n0=n1=150. Sampling repeats are not independent\n')
  cat('biological replicates. Distances use each matched pair\'s union observed set.\n\n')
  cat('Fits retain the frozen v1 landscape prior and division/count model, with a\n')
  cat('count-independent eta prior (mean=0, sd=5, distance penalty=0) for all nodes.\n')
  cat('This removes a count-dependent-prior difference between process arms; it is\n')
  cat('not v1 deployment performance and not simulation-based calibration.\n\n')
  cat(sprintf('Completed fits: %d/%d; trusted covariance: %d/%d.\n\n',
    sum(status$fit_status == 'ok'), nrow(status), sum(status$covariance_status == 'TMB_sdreport'), nrow(status)))
  print(overall[overall$metric %in% c('relative_rmse', 'coverage_95', 'interval_score_95'), ], row.names = FALSE)
  cat('\nPaired finite-population minus deterministic contrasts:\n')
  print(process$summary[process$summary$metric %in% c('relative_rmse', 'coverage_95', 'interval_score_95'),
    c('distance_group', 'metric', 'second_minus_first', 'descriptive_delta_low',
      'descriptive_delta_high', 'independent_landscapes', 'valid_dataset_pairs')], row.names = FALSE)
  cat('\nOracle structural diagnostic (NOT a deployable algorithm):\n')
  cat('Condition frozen Gaussian log-birth prior analytically on exact birth rates\n')
  cat('at actually observed ABM nodes. No count likelihood, pseudo-SD or flow data.\n')
  cat('This isolates what this prior predicts after the observed-fitness uncertainty\n')
  cat('is removed; it does not prove a statistical information impossibility.\n')
  print(oracle_family[oracle_family$metric %in% c('relative_rmse', 'coverage_95', 'interval_score_95'), ], row.names = FALSE)
  cat('\nAttribution and stopping implications:\n')
  cat('Undercoverage persists in the deterministic arm, so omitted finite-population\n')
  cat('noise cannot by itself explain the observed failures in these six landscapes.\n')
  cat('For the interaction family, far-node error and undercoverage remain even in\n')
  cat('the exact-anchor prior-only oracle. Frozen structure/shrinkage and uncertainty\n')
  cat('for unobserved interaction directions therefore merit a separate diagnostic.\n')
  cat('This is family-specific evidence of prior/structure mismatch and remaining\n')
  cat('uncertainty; it is not a universal impossibility theorem.\n')
  cat('Do not promote a process-noise extension or generic interval widening based\n')
  cat('on these results. At most two subsequent candidate changes should be selected\n')
  cat('after combining the independent pi0/growth oracle results, then tested on new\n')
  cat('held-out landscapes. No reserved round-two test data were generated here.\n')
}, file = summary_file)
cat(normalizePath(summary_file), '\n')
