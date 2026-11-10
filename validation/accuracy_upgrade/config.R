# Frozen development pilot. Calibration/test landscapes are reserved, not run.
accuracy_upgrade_config <- function() {
  list(schema_version = 1L, n_chr = 4L, min_cn = 1L, max_cn = 3L,
       families = c("smooth", "rough", "interaction"),
       split_replicates = list(train = 1:2, calibration = 3:4, test = 5:8),
       run_split = "train", abm_replicates = 1:2, sampling_replicates = 1:2,
       beta = 0.005, dt = 1.5, time_step = 0.025,
       initial_population = 1200L, n0 = 150L, n1 = 150L,
       initial_neighbour_mass = 0.15,
       local_shell_depth = 1L, max_nodes = 81L,
       fit_arguments = list(lambda_l_grid = 1, lambda_e_grid = 0.1,
                            sigma_obs_grid = 0.05, graph_edge_weight = "unit",
                            control = list(eval.max = 500, iter.max = 500),
                            retry_control = list(eval.max = 1500, iter.max = 1500)),
       bootstrap_replicates = 1000L, bootstrap_seed = 911204L,
       interval_levels = c(0.5, 0.8, 0.95),
       node_weighting = "equal within each exact distance, equal across nonempty distances",
       gauge = "subtract the t0-count-weighted observed-node fitness mean",
       calibration_policy = "No truth-based calibration in development pilot; calibration/test seeds sealed")
}

accuracy_seed_manifest <- function(config) {
  rows <- list()
  for (family_id in seq_along(config$families)) {
    for (split in names(config$split_replicates)) {
      for (replicate in config$split_replicates[[split]]) {
        rows[[length(rows) + 1L]] <- data.frame(
          family = config$families[family_id], split = split,
          landscape_replicate = replicate,
          landscape_id = sprintf("%s_L%02d", config$families[family_id], replicate),
          landscape_seed = 110000L + 1000L * family_id + replicate,
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}
