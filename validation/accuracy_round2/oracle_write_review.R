#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
run_root <- if (length(args)) args[1L] else "/private/tmp/alfak2_iteration2_20261104_0518"
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
destination <- dirname(script)
data_dir <- file.path(run_root, "oracle-results")
reference_dir <- file.path(run_root, "oracle-reference")
artifact <- file.path(destination, "oracle_results")
dir.create(artifact, recursive = TRUE, showWarnings = FALSE)
read <- function(file) utils::read.csv(file.path(data_dir, file), stringsAsFactors = FALSE)
family <- read("oracle_family_metrics.csv")
mass <- read("oracle_initial_mass_groups.csv")
baseline <- read("oracle_birth_and_growth.csv")
status <- read("oracle_fit_status.csv")
decomposition <- read("oracle_repeat_distance.csv")
comparison <- utils::read.csv(file.path(reference_dir, "oracle_reference_comparison.csv"), stringsAsFactors = FALSE)
convergence <- utils::read.csv(file.path(reference_dir, "oracle_reference_convergence.csv"), stringsAsFactors = FALSE)
sampling <- utils::read.csv(file.path(reference_dir, "oracle_reference_sampling_check.csv"), stringsAsFactors = FALSE)
arms <- c("count_only", "known_pi0", "population_growth", "known_pi0_and_growth")
labels <- c(count_only = "Counts only", known_pi0 = "Exact initial frequencies",
            population_growth = "Population growth", known_pi0_and_growth = "Both information sources")
value <- function(arm, fam, distance, metric) {
  x <- family$value[family$model == arm & family$family == fam &
                      family$distance_group == paste0("d", distance) & family$metric == metric]
  if (length(x)) x[1L] else NA_real_
}
mass_summary <- stats::aggregate(cbind(fitted_pi0_mass, truth_pi0_mass) ~ arm + initial_group, mass, mean)
baseline_summary <- stats::aggregate(cbind(mode_pi0_weighted_birth, true_initial_mean_birth,
  mode_predicted_log_growth, realized_log_growth, pi0_rmse) ~ arm, baseline, mean)
lines <- c("# Round 2: initial-composition and growth-information oracles", "",
  "Exact initial-frequency information changes the inferred common birth baseline substantially in this development sample, but it does not recover the smooth or interaction d0 point errors. Population-growth information also shifts the baseline; its benefits and harms vary by family and distance. These interventions therefore do not isolate one sufficient explanation for the earlier near-field regressions.", "",
  "## Frozen information interventions", "",
  "The experiment uses six new three-chromosome landscapes (two smooth, two rough, two interaction), one exact Gillespie trajectory per landscape and two independently sampled count pairs per trajectory: 12 datasets and four paired information arms, 48 fits. Every arm uses exactly the same counts, full 27-node domain, v1 fitness prior and t0-count reference. Only the initial-frequency oracle and/or external total-population observation changes. The exact initial oracle fixes all probabilities, including true zeros; it supplies no pseudo-counts. Population growth uses the actual totals from the same trajectory and a prespecified log-growth SD of 0.1.", "",
  "All outputs are development diagnostics. Extra-information arms are not fair algorithmic improvements over count-only estimation. The v1 source tree is frozen; only the explicit diagnostic oracle mechanism is overlaid. No formal test is opened or promoted here, and no interval calibration is applied.", "",
  sprintf("All %d fits completed with trusted covariance. Single-fit runtimes ranged from %.2f to %.2f seconds.", nrow(status), min(status$runtime_seconds), max(status$runtime_seconds)), "",
  "## Common baseline and initial frequencies", "",
  "The birth column is the conditional mode averaged with the fitted initial frequencies; it is compared with the true initial-frequency-weighted birth rate. Log growth is the model's conditional-mode prediction compared with the realized finite-population log growth. These are distinct quantities.", "",
  "| Information | Mode-weighted birth | True weighted birth | Predicted log growth | Realized log growth | Initial-frequency RMSE |",
  "|---|---:|---:|---:|---:|---:|")
for (arm in arms) {
  x <- baseline_summary[baseline_summary$arm == arm, ]
  lines <- c(lines, sprintf("| %s | %.5f | %.5f | %.5f | %.5f | %.6f |", labels[arm],
    x$mode_pi0_weighted_birth, x$true_initial_mean_birth, x$mode_predicted_log_growth,
    x$realized_log_growth, x$pi0_rmse))
}
lines <- c(lines, "", "## Initial mass by observed appearance", "",
  "Groups use the unchanged count pair. A zero initial count is not asserted to mean biological absence: late-only rows retain nonzero true initial mass in this experiment. Means below include datasets where a group has no nodes, preserving total mass accounting.", "",
  "| Information | Count support | Estimated initial mass | True initial mass |",
  "|---|---|---:|---:|")
for (arm in arms) for (group in c("t0_positive", "t1_only", "both_zero")) {
  x <- mass_summary[mass_summary$arm == arm & mass_summary$initial_group == group, ]
  lines <- c(lines, sprintf("| %s | %s | %.7f | %.7f |", labels[arm], group,
                            x$fitted_pi0_mass, x$truth_pi0_mass))
}
lines <- c(lines, "", "## Fitness and raw interval results", "",
  "Distance is measured from the observed support in this finite-population dataset, identically for all four arms. The small development design contains d0–d2; it provides no d3 evidence. Scores use the full joint lognormal Laplace distribution before subtracting the original t0-count reference in each posterior draw.", "",
  "| Family | Distance | Information | Relative RMSE | 95% coverage | 95% interval score |",
  "|---|---|---|---:|---:|---:|")
for (fam in c("smooth", "rough", "interaction")) for (distance in 0:2) for (arm in arms)
  lines <- c(lines, sprintf("| %s | d%d | %s | %.5f | %.3f | %.5f |", fam, distance, labels[arm],
    value(arm, fam, distance, "relative_rmse"), value(arm, fam, distance, "coverage_95"),
    value(arm, fam, distance, "interval_score_95")))
lines <- c(lines, "", "## Repeat-error decomposition", "",
  "The CSV fixes each node's distance and appearance group using sampling replicate 1, then decomposes mean squared error across the two sampling repeats into squared repeat-mean error plus within-repeat error variance. Their numerical sum equals the reported MSE. Squared repeat-mean error is not an unbiased estimate of population bias squared; a shared ABM realization and finite sampling contribute to it. One biological path per landscape cannot separate systematic prior bias from shared process fluctuations. Native-distance metrics above and fixed-reference repeat strata are intentionally reported separately.", "",
  sprintf("Maximum algebraic residual in the repeat decomposition: %.3g.", max(abs(decomposition$decomposition_residual))), "",
  "## Independent small-graph posterior reference", "",
  "A separate fixed two-node graph uses a Gaussian log-birth prior and an initial-logit prior specified before the fixed count observations. The normalized coefficient-plus-GMRF prior is analytically reduced to its bivariate Gaussian log-birth distribution. At beta=0 the frequency and total-growth formulas are analytic. Gauss-Hermite importance quadrature integrates the original prior times that exact likelihood, including numerical integration of initial log odds when unknown; the broad Laplace-guided proposal does not replace the target density. Orders 40, 64, 96 and 160 check moment convergence, with a finer/wider log-odds integral at the higher orders. Independent prior importance sampling provides a second check with Monte Carlo standard errors.", "",
  "This is a fixed-data posterior comparison, not SBC or repeated-sampling coverage. Indicator-based Gauss-Hermite quantiles converge slowly and are not certified or used to make interval-coverage claims.", "",
  "| Information | Quantity | Reference mean | Laplace mean | Reference SD | Laplace / reference SD |",
  "|---|---|---:|---:|---:|---:|")
for (arm in arms) for (quantity in c("birth_1", "birth_2", "relative_contrast_1")) {
  x <- comparison[comparison$arm == arm & comparison$quantity == quantity, ]
  lines <- c(lines, sprintf("| %s | %s | %.6f | %.6f | %.6f | %.4f |", labels[arm], quantity,
    x$reference_mean, x$laplace_mean, x$reference_sd, x$sd_ratio))
}
fine <- convergence[convergence$coarse_order == 96L, ]
mean_z <- abs(sampling$sampling_mean - sampling$quadrature_mean) / sampling$sampling_mean_mcse
sd_z <- abs(sampling$sampling_sd - sampling$quadrature_sd) / sampling$sampling_sd_mcse
lines <- c(lines, "", sprintf("For quadrature order 96 versus 160, the maximum birth-mean difference is %.3g, covariance difference %.3g, and contrast-SD difference %.3g. The independent 500,000-prior-draw checks have ESS %.0f–%.0f; their maximum mean and SD discrepancies are %.2f and %.2f estimated Monte Carlo standard errors.",
  max(fine$max_mean_difference), max(fine$max_covariance_difference), max(fine$contrast_sd_difference),
  min(sampling$importance_ess), max(sampling$importance_ess), max(mean_z), max(sd_z)), "",
  "Without growth information, this example's Laplace/lognormal contrast variance is larger than the numerical posterior reference, not smaller. With growth information the contrast moments are much closer. Non-Gaussian baseline uncertainty is therefore a real approximation concern, but this example does not support attributing all observed undercoverage to Laplace underdispersion.", "",
  "## Limits and use in model design", "",
  "These six small landscapes do not reproduce the original four-chromosome family regressions exactly. They test information mechanisms while keeping each dataset and estimator fixed. Exact initial frequencies slightly help rough point estimation but do not repair smooth/interaction d0 point errors; growth information can tighten or shift intervals and sometimes worsens coverage. Candidate design should therefore consider both composition/baseline confounding and the fitness prior, with a separate process-noise assessment. No oracle gain is counted as an algorithmic promotion result.", "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
writeLines(lines, file.path(destination, "oracle_review.md"))
files <- c("oracle_family_metrics.csv", "oracle_initial_mass_groups.csv", "oracle_birth_and_growth.csv",
           "oracle_fit_status.csv", "oracle_repeat_distance.csv", "oracle_repeat_initial_group.csv")
file.copy(file.path(data_dir, files), artifact, overwrite = TRUE)
reference_files <- c("oracle_reference_comparison.csv", "oracle_reference_convergence.csv", "oracle_reference_sampling_check.csv")
file.copy(file.path(reference_dir, reference_files), artifact, overwrite = TRUE)
utils::write.csv(mass_summary, file.path(artifact, "oracle_initial_mass_summary.csv"), row.names = FALSE)
utils::write.csv(baseline_summary, file.path(artifact, "oracle_birth_growth_summary.csv"), row.names = FALSE)
cat(file.path(destination, "oracle_review.md"), "\n")
