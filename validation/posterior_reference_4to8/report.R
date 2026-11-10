#!/usr/bin/env Rscript
directory <- "validation/posterior_reference_4to8"
output <- Sys.getenv("ALFAK2_REFERENCE_OUTPUT", file.path(directory, "results"))
status <- read.csv(file.path(output, "case_status.csv"))
moments <- read.csv(file.path(output, "moment_comparison.csv"))
covariance <- read.csv(file.path(output, "covariance_comparison.csv"))
quantiles <- read.csv(file.path(output, "quantile_comparison.csv"))
checks <- read.csv(file.path(output, "independent_forward_checks.csv"))
tests <- read.csv(file.path(output, "validation/validation_tests.csv"))
stopifnot(nrow(status) == 16L, length(unique(moments$case_id)) == 16L,
  all(is.finite(moments$standardized_mean_difference)), all(is.finite(moments$sd_ratio)))
relative <- moments[moments$type == "relative", ]
range_text <- function(x, digits = 3) paste(formatC(range(x), digits = digits, format = "f"), collapse = "–")
rows <- lapply(status$case_id, function(id) {
  x <- relative[relative$case_id == id, ]; s <- status[status$case_id == id, ]
  C <- covariance[covariance$case_id == id & covariance$type == "relative", ]
  data.frame(case_id = id, n = s$n, profile = s$profile, beta = s$beta, growth = s$growth,
    reference_certified = s$passed, retained_per_chain = s$retained_per_chain,
    max_abs_standardized_mean_error = max(abs(x$standardized_mean_difference)),
    sd_ratio_min = min(x$sd_ratio), sd_ratio_max = max(x$sd_ratio),
    relative_covariance_error = C$laplace_relative_frobenius_error,
    reference_relative_covariance_mcse = C$reference_covariance_relative_mcse,
    any_material_mean = any(x$material_mean), any_material_sd = any(x$material_sd))
})
case_summary <- do.call(rbind, rows)
write.csv(case_summary, file.path(output, "case_comparison_summary.csv"), row.names = FALSE)
group <- do.call(rbind, lapply(c(4L, 8L), function(n) do.call(rbind, lapply(c(FALSE, TRUE), function(growth) {
  x <- relative[relative$n == n & relative$growth == growth, ]
  s <- case_summary[case_summary$n == n & case_summary$growth == growth, ]
  data.frame(n = n, growth = growth, certified_cases = sum(s$reference_certified), cases = nrow(s),
    sd_ratio_min = min(x$sd_ratio), sd_ratio_max = max(x$sd_ratio),
    max_abs_standardized_mean_error = max(abs(x$standardized_mean_difference)),
    cases_with_material_mean = sum(s$any_material_mean), cases_with_material_sd = sum(s$any_material_sd))
}))))
write.csv(group, file.path(output, "group_summary.csv"), row.names = FALSE)
lines <- c("# Independent posterior reference: 4 and 8 nodes", "",
  sprintf("%d/16 fixed-data references meet all preregistered Monte Carlo diagnostics. All 16 Laplace fits return trusted TMB covariance. This checkpoint compares algorithms under a fixed model; it does not certify biological fitness accuracy or repeated-sampling interval coverage.", sum(status$passed)), "",
  "## Scope and checks", "",
  "The full CN1–2 domains, counts, prior scales, beta values and information arms were fixed in [the protocol](PROTOCOL.md) before the primary comparisons. There are four independent MCMC chains per case. The exact aggregate Gaussian log-birth prior is analytically marginalized; the softmax mean-logit gauge is removed with the correct contrast covariance. Fitness priors use V1 scales; every initial-logit SD is explicitly 2. These are fixed-data mechanism probes, not an assessment of all deployment configurations.", "",
  sprintf("The independent sampler/target validation passes %d assertions with no failures, errors or warnings. It includes analytic Gaussian-prior moments, the earlier converged two-node quadrature, TMB pointwise likelihood agreement at 20 off-mode states, and explicit-versus-collapsed Gaussian-prior Laplace agreement on four graph/beta combinations.", sum(tests$passed)), "",
  sprintf("Across primary cases, maximum rank-normalized split Rhat is %.5f; minimum bulk/tail ESS is %.0f/%.0f. Retained draws per chain range from %d to %d. Maximum primary-mean/SD MCSE divided by reference SD is %.3f/%.3f; maximum relative covariance MCSE is %.3f. The case status table preserves failed attempts before any allowed length extension.", max(status$max_rhat), min(status$min_bulk_ess), min(status$min_tail_ess), min(status$retained_per_chain), max(status$retained_per_chain), max(status$max_mean_mcse_sd), max(status$max_sd_mcse_sd), max(status$max_covariance_relative_mcse)), "",
  "## Relative fitness comparison", "",
  "Every birth rate is contrasted with the same fixed t0-count-weighted reference. SD ratios are Laplace/MCMC. Mean errors are divided by the reference posterior SD. A material flag requires the predeclared numerical discrepancy threshold and a difference larger than four MCMC standard errors; it is not a population hypothesis test.", "",
  "| Nodes | Information | Certified cases | Relative SD ratio range | Largest absolute mean error / reference SD | Cases with material mean / SD error |",
  "|---:|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(group))) {
  x <- group[i, ]
  lines <- c(lines, sprintf("| %d | %s | %d/%d | %.3f–%.3f | %.3f | %d / %d |", x$n,
    if (x$growth) "Counts + growth" else "Counts only", x$certified_cases, x$cases,
    x$sd_ratio_min, x$sd_ratio_max, x$max_abs_standardized_mean_error,
    x$cases_with_material_mean, x$cases_with_material_sd))
}
lines <- c(lines, "", "The information-arm comparison probes whether constraining common growth reduces approximation error. Trusted Hessians alone do not diagnose that error. Growth adds information; any better agreement is a mechanism probe, not a fair algorithmic promotion over count-only data. These cases cannot attribute all far-field undercoverage to Laplace: a broad model posterior can still miss truths outside a restrictive model/prior.", "",
  "## All cases", "",
  "| Case | Reference certified | Draws/chain | Relative SD ratio | Max mean error / SD | Relative covariance error |",
  "|---|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(case_summary))) {
  x <- case_summary[i, ]
  lines <- c(lines, sprintf("| %s | %s | %d | %.3f–%.3f | %.3f | %.3f |", x$case_id,
    x$reference_certified, x$retained_per_chain, x$sd_ratio_min, x$sd_ratio_max,
    x$max_abs_standardized_mean_error, x$relative_covariance_error))
}
lines <- c(lines, "", "## Numerical forward checks", "",
  sprintf("An independent converged continuous-time Taylor solver checks 128 representative retained states per case (2,048 total). Maximum frequency and log-growth discrepancies are %.3g and %.3g; the maximum absolute log-likelihood discrepancy is %.3g. These are sampled posterior-state checks, not a global forward-error bound. Primary MCMC uses the same numerical likelihood as TMB; the continuous solver is a separate numerical-sensitivity check.", max(checks$maximum_pi1_difference), max(abs(checks$log_growth_difference)), max(abs(checks$exact_minus_split_log_likelihood))), "",
  "## Interpretation and next decision", "",
  "The tables separately compare latent log-rate moments, physical birth moments, relative/centered fitness and graph-edge contrasts. `moment_comparison.csv` decomposes total mean error into (i) using Laplace rather than reference log-rate moments and (ii) Gaussianizing the reference log-rate distribution before exponentiation. The terms sum exactly, but their magnitudes are descriptive rather than causal percentages. Complete covariance comparisons are retained; no shell-specific recentering is used.", "",
  "The scalar `posterior_squared_loss_ratio = 1 + standardized_mean_difference^2` is the conditional model-posterior quadratic-loss penalty of using the approximate mean. It is not an observed-truth RMSE ratio. Quantile tables report fixed-data reference posterior mass inside approximate intervals, not repeated-sampling coverage; both reference and iid approximate-quantile MCSE are explicit.", "",
  "The next computational candidate, if pursued, should integrate over common growth-scale uncertainty while preserving the full Gaussian prior and conditional shape uncertainty. Its accuracy must first be demonstrated against these references; changing an arbitrary centering constraint or merely rescaling intervals is not equivalent. Even an accurate posterior will not manufacture information about unvisited epistatic states. No production model, new accuracy calibration/test landscape, HPC job, commit or push was changed by this checkpoint.", "",
  "Full raw chains are kept in `/private/tmp/alfak2_reference_4to8_20261104`, with SHA256 and paths in `results/case_status.csv`. Reproducible cases, compact posterior summaries, diagnostics, covariance/quantile tables, source fingerprints and validation results are retained here.", "",
  "The original live runner emitted a trailing parser error after all 16 completed outputs were written: file-path options were edited while its already-parsed sampling loop was running. The saved results were subsequently verified case by case, including four-chain dimensions, retained lengths, convergence and raw-file SHA256 (`results/completion_integrity.txt`). Exact original computation sources were archived before those file-path edits. All target/sampler/case/fit/diagnostic function bodies remained unchanged, and a separate portable-path validation rerun passed the same 99 assertions. The current portable runner parses cleanly; the late process-exit error is recorded rather than hidden.", "",
  "The subsequent prior-preserving scalar-integration candidate and deployment-prior development pilot are reported separately in [the candidate review](../baseline_quadrature/REVIEW.md). Numerical accuracy and generating-truth accuracy remain separate decisions.", "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
writeLines(lines, Sys.getenv("ALFAK2_REFERENCE_REPORT", file.path(directory, "REVIEW.md")))
cat("Wrote REVIEW.md and comparison summaries.\n")
