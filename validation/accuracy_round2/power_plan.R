#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R")); source(file.path(script_dir, "config.R"))
config <- accuracy_round2_config()
input <- round2_arg("landscape-metrics", file.path(script_dir, "results", "old_data_diagnostics", "landscape_metrics.csv"))
output <- round2_arg("output", file.path(script_dir, "results", "power_planning"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
x <- utils::read.csv(input, stringsAsFactors = FALSE)
x <- x[x$distance_view == "actual_support" & x$distance >= 0, ]
reference <- x[x$model == "v1_joint", ]
proxy <- x[x$model == "endpoint_ablation", ]
risks <- c("rmse", "spatially_centered_rmse", "edge_gradient_rmse")
family_variance <- proxy_variance <- list()
for (risk in risks) for (family in config$families) for (distance in sort(unique(reference$distance))) {
  a <- reference[reference$family == family & reference$distance == distance, ]
  b <- proxy[proxy$family == family & proxy$distance == distance, ]
  b <- b[match(a$landscape_id, b$landscape_id), ]
  values <- a[[risk]]
  family_variance[[length(family_variance) + 1L]] <- data.frame(risk = risk, family = family, distance = distance,
    n_old_landscapes = length(values), reference_mean = mean(values), reference_sd = stats::sd(values))
  proxy_difference <- b[[risk]] - values
  proxy_variance[[length(proxy_variance) + 1L]] <- data.frame(risk = risk, family = family, distance = distance,
    comparator = "endpoint_ablation_minus_frozen_v1_proxy_only", n_old_landscapes = length(values),
    paired_sd = stats::sd(proxy_difference), risk_correlation = stats::cor(values, b[[risk]]),
    observed_difference_not_assumed_v2_effect = mean(proxy_difference))
}
family_variance <- do.call(rbind, family_variance); proxy_variance <- do.call(rbind, proxy_variance)
rows <- list()
budgets <- c(conjunction_only = 1L, point_family_16 = 16L, known_domain_shape_76 = 76L, both_domains_190 = 190L)
for (budget_name in names(budgets)) {
  critical <- stats::qnorm(1 - config$alpha / budgets[[budget_name]])
  for (rho in c(0, 0.5, 0.8, 0.95, 0.99)) for (inflation in c(1, config$variance_inflation_sensitivity)) {
    for (risk in risks) for (distance in sort(unique(reference$distance))) {
      f <- family_variance[family_variance$risk == risk & family_variance$distance == distance, ]
      # v2 risk variance is not yet observed. Equal marginal variance and rho
      # are explicit assumptions, not an estimate of candidate-v2 performance.
      paired_variance <- 2 * (1 - rho) * inflation * f$reference_sd^2
      target_effect <- config$primary_effect * mean(f$reference_mean)
      variance_coefficient <- sum(paired_variance) / 3
      needed <- 3 * ceiling((critical + stats::qnorm(config$target_power))^2 * variance_coefficient / target_effect^2 / 3)
      for (n_total in seq(18L, config$max_test_landscapes_per_domain, by = 3L)) {
        se <- sqrt(variance_coefficient / n_total)
        power <- if (se > 0) stats::pnorm(target_effect / se - critical) else 1
        rows[[length(rows) + 1L]] <- data.frame(gate = "family_balanced_superiority_5percent", risk = risk,
          family = "ALL", distance = distance, interval_budget = budget_name, simultaneous_gates = budgets[[budget_name]],
          assumed_v2_v1_correlation = rho, variance_multiplier = inflation,
          n_test_landscapes = n_total, n_per_family = n_total / 3,
          assumed_absolute_effect = target_effect, approximate_marginal_power = power,
          n_required_for_80percent_marginal = needed, required_exceeds_cap = needed > config$max_test_landscapes_per_domain)
      }
      for (fi in seq_len(nrow(f))) {
        family_margin <- config$family_noninferiority_margin * f$reference_mean[fi]
        # Assumes v2 and v1 have equal mean risk, hence no harm. The 5% NI
        # distance is to the upper margin, not a claim of a true 5% benefit.
        margin_multiplier <- 1 + config$family_noninferiority_margin
        family_variance_ni <- (1 + margin_multiplier^2 - 2 * margin_multiplier * rho) * inflation * f$reference_sd[fi]^2
        needed_ni <- 3 * ceiling((critical + stats::qnorm(config$target_power))^2 * family_variance_ni / family_margin^2)
        for (n_total in seq(18L, config$max_test_landscapes_per_domain, by = 3L)) {
          se <- sqrt(family_variance_ni / (n_total / 3))
          power <- if (se > 0) stats::pnorm(family_margin / se - critical) else 1
          rows[[length(rows) + 1L]] <- data.frame(gate = "family_noninferiority_5percent_assume_no_harm", risk = risk,
            family = f$family[fi], distance = distance, interval_budget = budget_name, simultaneous_gates = budgets[[budget_name]],
            assumed_v2_v1_correlation = rho, variance_multiplier = inflation,
            n_test_landscapes = n_total, n_per_family = n_total / 3,
            assumed_absolute_effect = family_margin, approximate_marginal_power = power,
            n_required_for_80percent_marginal = needed_ni, required_exceeds_cap = needed_ni > config$max_test_landscapes_per_domain)
        }
      }
    }
  }
}
planning <- do.call(rbind, rows)
round2_write_csv(planning, file.path(output, "power_curves.csv"))
round2_write_csv(family_variance, file.path(output, "old_v1_family_risk_variance.csv"))
round2_write_csv(proxy_variance, file.path(output, "endpoint_v1_pair_variance_proxy.csv"))
at_budget <- planning[planning$n_test_landscapes %in% c(36L, 60L) & planning$variance_multiplier == 1 &
                       planning$assumed_v2_v1_correlation %in% c(0.8, 0.95, 0.99) &
                       planning$interval_budget %in% c("conjunction_only", "known_domain_shape_76"), ]
summary <- aggregate(approximate_marginal_power ~ gate + risk + interval_budget + assumed_v2_v1_correlation + n_test_landscapes,
                     at_budget, min)
names(summary)[names(summary) == "approximate_marginal_power"] <- "minimum_marginal_planning_power_not_joint_power"
round2_write_csv(summary, file.path(output, "budget_summary.csv"))
decision <- list(status = "DEFER_TEST_N_NO_CANDIDATE_V2_OOF_PAIR_VARIANCE",
  primary_comparator = "frozen_v1", target_effect = 0.05, max_landscapes_per_domain = 60L,
  provisional_landscapes = 36L, no_test_generation_authorized = TRUE,
  limitations = c("Only four old independent landscapes per family", "Candidate-v2 loss covariance is unobserved",
    "Normal planning approximation; no power guarantee", "cn3_1to4 d4/d5 has no old-domain variance data",
    "IUT one-sided alpha=.05 is the primary conjunctive promotion rule; simultaneous budgets are reporting sensitivity only",
    "Actual candidate-v2 OOF loss vectors are required for family-stratified joint pass simulation",
    "Shape/rank/amplitude/uncertainty gates must be frozen before test", "No post-test expansion to chase significance"),
  source_sha256 = round2_sha256(c(input, file.path(script_dir, "power_plan.R"), file.path(script_dir, "config.R"))))
saveRDS(decision, file.path(output, "planning_decision.rds"))
cat(decision$status, "\n")
print(subset(summary, interval_budget == "conjunction_only" & assumed_v2_v1_correlation == 0.95), row.names = FALSE)
