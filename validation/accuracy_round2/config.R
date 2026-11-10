accuracy_round2_config <- function() {
  list(protocol_version = "round2-development-1", status = "development_only_new_splits_sealed",
    families = c("smooth", "rough", "interaction"),
    domains = list(
      cn4_1to3 = list(n_chr = 4L, min_cn = 1L, max_cn = 3L, n_target = 81L,
                      seed_namespace = 3100000L, primary_distance = "actual_observed_support", planned_distances = 0:3),
      cn3_1to4 = list(n_chr = 3L, min_cn = 1L, max_cn = 4L, n_target = 64L,
                      seed_namespace = 3200000L, primary_distance = "actual_observed_support", planned_distances = 0:5)),
    landscape_replicates = list(development = 1:8, calibration = 9:14, test_reserved = 15:34),
    planned_test_replicates = 15:26, max_test_landscapes_per_domain = 60L,
    calibration_landscapes_per_domain = 18L, provisional_test_landscapes_per_domain = 36L,
    abm_replicates = 1:2, sampling_replicates = 1:2,
    generator = "exact_gillespie_correlated_daughters",
    truth_versions = c(cn4_1to3 = "round1_linear_truth_unchanged", cn3_1to4 = "bounded_log_link_v1"),
    beta = 0.005, dt = 1.5, initial_population = 1200L,
    n0 = 150L, n1 = 150L, initial_neighbour_mass = 0.15,
    interval_levels = c(0.5, 0.8, 0.95),
    primary_effect = 0.05, family_noninferiority_margin = 0.05,
    primary_comparator = "v1_joint", additional_controls = c("current_default_two_layer", "zero_point_predictor"),
    historical_reference = "original_baseline",
    alpha = 0.05, target_power = 0.80, variance_inflation_sensitivity = 1.5,
    promotion_statement = "single_frozen_winner_intersection_union_one_sided_alpha_0p05",
    maximum_development_candidates = 2L, maximum_final_test_winners = 1L,
    multiple_cell_reporting = "simultaneous_or_Holm_separate_from_conjunctive_promotion",
    hard_metrics = c("rmse", "spatially_centered_rmse", "edge_gradient_rmse", "spearman", "amplitude_ratio", "variance_ratio", "interval_score95"),
    shape_primary_effect = 0.05, shape_family_noninferiority_margin = 0.05,
    shape_centering = "one spatial mean over the complete target domain; never recenter each shell",
    rank_amplitude_threshold_status = "must_freeze_on_development_before_any_new_test",
    simultaneous_reporting_budget = c(legacy_point_family = 16L, known_domain_with_shape_rank_amplitude = 76L,
                                       both_domains_with_shape_rank_amplitude = 190L),
    bootstrap_replicates = 5000L, power_seed = 610040L,
    gauge = "one t0-count-weighted observed reference per dataset, shared by every target",
    node_weighting = "equal within distance; equal landscapes within family; equal families",
    fit_slots = 2L, slot_directory = "/private/tmp/alfak2_iteration2_20261104_0518/fit_slots",
    freeze_policy = "No calibration generation before candidate freeze; no test generation before model/calibration/analysis freeze",
    family_usage = "generator and audit labels only; never passed as model inputs")
}

accuracy_round2_seed_manifest <- function(config) {
  rows <- list()
  for (domain_name in names(config$domains)) {
    domain <- config$domains[[domain_name]]
    for (fi in seq_along(config$families)) for (split in names(config$landscape_replicates)) {
      for (replicate in config$landscape_replicates[[split]]) {
        seed <- domain$seed_namespace + 1000L * fi + replicate
        rows[[length(rows) + 1L]] <- data.frame(domain = domain_name, family = config$families[fi],
          split = split, landscape_replicate = replicate,
          landscape_id = sprintf("%s_%s_L%02d", domain_name, config$families[fi], replicate),
          landscape_seed = seed, activated = split == "development" ||
            (split == "test_reserved" && replicate %in% config$planned_test_replicates),
          generation_allowed_now = split == "development", stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

accuracy_round2_model_arguments <- function(model) {
  switch(model,
    original_baseline = list(lambda_l_grid = 1, lambda_e_grid = 0.1, sigma_obs_grid = 0.05, graph_edge_weight = "unit"),
    zero_point_predictor = list(point_only = TRUE, relative_fitness = 0, interval_available = FALSE),
    current_default_two_layer = list(pipeline = "two_layer", extrapolation_control = list(compute_covariance = TRUE)),
    v1_joint = list(pipeline = "joint_count", forward_model = "division_process", inference_method = "laplace",
      sdreport_mode = "fitness_only", joint_control = list(baseline_sd = 0.5, main_sd = 0.07,
        curvature_sd = 0.04, interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25,
        interaction_support = "isotropic")),
    A = list(pipeline = "joint_count", forward_model = "division_process", inference_method = "laplace",
      sdreport_mode = "fitness_only", joint_control = list(baseline_sd = 0.5, main_sd = 0.07,
        curvature_sd = 0.04, interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25,
        interaction_support = "observed_nullspace", unidentified_interaction_sd = 0.2)),
    B = list(pipeline = "joint_count", forward_model = "division_process", inference_method = "laplace",
      sdreport_mode = "fitness_only", joint_control = list(baseline_sd = 0.5, main_sd = 0.07,
        curvature_sd = 0.04, interaction_sd = 0.04, residual_sd = 0.125, edge_sd = 0.25,
        interaction_support = "observed_nullspace", unidentified_interaction_sd = 0.2,
        prior_components = list(weak = list(residual_sd = 0.03, edge_sd = 0.06),
                                rough = list(residual_sd = 0.3, edge_sd = 0.6)), model_prior = c(0.75, 0.25))),
    stop("Unknown comparison model; candidate_v2 must be supplied as a frozen argument RDS"))
}
