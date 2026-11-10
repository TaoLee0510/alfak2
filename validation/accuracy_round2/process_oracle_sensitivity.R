#!/usr/bin/env Rscript
# Prespecified oracle structure sensitivity, never a trained candidate result.
script <- grep('^--file=', commandArgs(FALSE), value = TRUE)
directory <- dirname(normalizePath(sub('^--file=', '', script[[1L]])))
source(file.path(directory, 'process_oracle_anchors.R'))
args <- commandArgs(TRUE)
base <- if (length(args)) args[[1L]] else '/private/tmp/alfak2_iteration2_20261104_0518'
repo <- dirname(dirname(directory))
output <- file.path(base, 'process-oracle-sensitivity')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
variants <- list(pair_sd_004 = list(interaction_sd = .04),
                 pair_sd_020 = list(interaction_sd = .2),
                 pair_sd_050 = list(interaction_sd = .5),
                 no_pair = list(include_interactions = FALSE))
metrics <- native_metrics <- nodes <- statuses <- list()
for (variant in names(variants)) {
  folder <- file.path(output, variant)
  run_process_anchor_oracle(file.path(base, 'process-data'), folder,
    file.path(base, 'v1-source'), repo, joint_overrides = variants[[variant]],
    model_label = paste0('oracle_', variant))
  metrics[[variant]] <- read.csv(file.path(folder, 'oracle_anchor_metrics.csv'), stringsAsFactors = FALSE)
  metrics[[variant]]$variant <- variant
  native_metrics[[variant]] <- read.csv(file.path(folder, 'oracle_anchor_native_metrics.csv'), stringsAsFactors = FALSE)
  native_metrics[[variant]]$variant <- variant
  nodes[[variant]] <- read.csv(file.path(folder, 'oracle_anchor_nodes.csv'), stringsAsFactors = FALSE)
  nodes[[variant]]$variant <- variant
  statuses[[variant]] <- read.csv(file.path(folder, 'oracle_anchor_status.csv'), stringsAsFactors = FALSE)
  statuses[[variant]]$variant <- variant
}
metrics <- do.call(rbind, metrics); nodes <- do.call(rbind, nodes); status <- do.call(rbind, statuses)
native_metrics <- do.call(rbind, native_metrics)
keep <- metrics$metric %in% c('relative_rmse', 'coverage_95', 'interval_score_95', 'interval_width_95') &
  metrics$distance_group %in% c('d1', 'd2')
landscape <- aggregate(value ~ variant + family + landscape_id + distance_group + metric, metrics[keep, ], mean)
family <- aggregate(value ~ variant + family + distance_group + metric, landscape, mean)
overall <- aggregate(value ~ variant + distance_group + metric, landscape, mean)
write.csv(landscape, file.path(output, 'oracle_sensitivity_landscape.csv'), row.names = FALSE)
write.csv(family, file.path(output, 'oracle_sensitivity_family.csv'), row.names = FALSE)
write.csv(overall, file.path(output, 'oracle_sensitivity_overall.csv'), row.names = FALSE)
native_keep <- native_metrics$metric %in% c('relative_rmse', 'coverage_95', 'interval_score_95', 'interval_width_95') &
  grepl('^d[0-9]+$', native_metrics$distance_group)
native_landscape <- aggregate(value ~ variant + family + landscape_id + distance_group + metric,
  native_metrics[native_keep, ], mean)
native_family <- aggregate(value ~ variant + family + distance_group + metric, native_landscape, mean)
write.csv(native_landscape, file.path(output, 'oracle_sensitivity_native_landscape.csv'), row.names = FALSE)
write.csv(native_family, file.path(output, 'oracle_sensitivity_native_family.csv'), row.names = FALSE)
write.csv(status, file.path(output, 'oracle_observed_design_rank.csv'), row.names = FALSE)
write.csv(nodes, file.path(output, 'oracle_sensitivity_nodes.csv'), row.names = FALSE)
rank_summary <- aggregate(cbind(observed_pair_rank, incremental_pair_rank_given_other_basis,
  pair_nullspace_dimension, posterior_pair_sd_fraction_min, posterior_pair_sd_fraction_max) ~ variant + family,
  status, mean, na.action = na.pass)
write.csv(rank_summary, file.path(output, 'oracle_rank_summary.csv'), row.names = FALSE)
capture.output({
  cat('Oracle-only interaction-prior sensitivity: not deployable accuracy results.\n')
  cat('Same 12 development ABM datasets, exact observed true log-birth anchors.\n')
  cat('Only pair-interaction prior SD varies (.04/.2/.5), or interactions are omitted;\n')
  cat('main/curvature/residual remain frozen v1. No TMB fits or interval calibration.\n\n')
  cat('Main tables use the shared union-observed distance for paired-process\n')
  cat('comparability. The *_native_* tables additionally use distance to this ABM\n')
  cat('dataset\'s actual observed anchors; their d0 oracle error is exactly zero.\n\n')
  cat('Family-specific d1/d2 error, coverage and proper interval score:\n')
  print(family[family$metric != 'interval_width_95', ], row.names = FALSE)
  cat('\nRank and residual prior uncertainty for pair coefficients:\n')
  print(rank_summary, row.names = FALSE)
  cat('\nObserved pair-column rank is a design property, not full coefficient\n')
  cat('identifiability: main effects and the Gaussian residual remain confounded.\n')
  cat('Incremental rank projects out the other basis columns. A pair-design\n')
  cat('nullspace direction is exactly unconstrained by these oracle anchors.\n')
  cat('Posterior/prior pair SD ratios additionally account for residual uncertainty.\n')
  cat(sprintf('At v1 pair SD, %d/%d observed sets leave a pair-design nullspace;\n',
    sum(status$pair_nullspace_dimension[status$variant == 'pair_sd_004'] > 0),
    sum(status$variant == 'pair_sd_004')))
  cat(sprintf('%d/%d have no observed nonzero third-order product feature.\n',
    sum(status$third_order_feature_unseen[status$variant == 'pair_sd_004']),
    sum(status$variant == 'pair_sd_004')))
  cat('Pair-nullspace conditional variance equals its prior variance (checked\n')
  cat('analytically for every nonempty nullspace; tolerance 1e-8).\n')
  cat('Global pair widening improves interaction-family proper scores here but\n')
  cat('substantially worsens rough-family d2 point errors and proper scores.\n')
  cat('It cannot be treated as a uniform improvement or a replacement for a\n')
  cat('separate rough-residual explanation of observed local variation.\n')
  cat('These sensitivities diagnose the frozen prior and motivate a separately\n')
  cat('specified candidate; they are not independent validation or a method victory.\n')
}, file = file.path(output, 'ORACLE_STRUCTURE_SENSITIVITY.txt'))
cat(normalizePath(output), '\n')
