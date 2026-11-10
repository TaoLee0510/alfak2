#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
root <- if (length(args)) args[1L] else "/private/tmp/alfak2_iteration2_20261104_0518"
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
destination <- dirname(script)
repo <- normalizePath(file.path(destination, "..", ".."))
input <- file.path(root, "oracle-fisher")
artifact <- file.path(destination, "oracle_results")
dir.create(artifact, recursive = TRUE, showWarnings = FALSE)
summary <- utils::read.csv(file.path(input, "oracle_fisher_summary.csv"), stringsAsFactors = FALSE)
checks <- utils::read.csv(file.path(input, "oracle_fisher_numerical_checks.csv"), stringsAsFactors = FALSE)
test_output <- capture.output(test_results <- testthat::test_file(
  file.path(repo, "tests/testthat/test-profile-fisher-math.R"), reporter = "summary"))
test_summary <- as.data.frame(test_results)
if (any(test_summary$failed > 0) || any(test_summary$error)) stop("Fisher mathematical regression failed")
columns <- intersect(c("file", "test", "nb", "failed", "skipped", "error", "warning", "passed"), names(test_summary))
utils::write.csv(test_summary[, columns, drop = FALSE], file.path(artifact, "oracle_fisher_math_tests.csv"), row.names = FALSE)
writeLines(test_output, file.path(artifact, "oracle_fisher_math_tests.txt"))
saveRDS(test_results, file.path(input, "oracle_fisher_math_tests.rds"))
value <- function(family, scenario, depth, column) {
  summary[summary$family == family & summary$scenario == scenario & summary$planned_depth == depth, column][1L]
}
labels <- c(unknown_pi0_profile = "Initial frequencies profiled",
            fixed_pi0_same_mode = "Initial frequencies fixed at same mode",
            fixed_true_pi0_same_birth = "True initial frequencies; same birth mode")
lines <- c("# Round 2: local expected Fisher information", "",
  "At the three examined v1 development modes, the 27-node count likelihood has algebraic fitness rank 26, yet the prior-relative amount of learning at 150 sampled cells per timepoint is only 0.54–0.71 effective directions after profiling initial frequencies. Scaling expected sampling information to 5,000 cells raises this diagnostic to 4.09–4.50 directions. These are local information calculations, not additional datasets, fitted methods or accuracy guarantees.", "",
  "## Expansion points and construction", "",
  "One prespecified finite-population dataset from the first landscape and first sampling replicate of each smooth, rough and interaction family supplies a frozen v1 conditional mode. Three 27-node count-only fits reproduce the required complete mode; no population-growth observation is added. All fitness priors retain v1 values. The independent forward calculation uses the exact continuous-time mean generator A=(2 P^T−I) diag(b), a converged positive matrix-exponential series, and central differences in g=log(b). Initial softmax derivatives and their propagation into final frequencies are analytic.", "",
  "For each timepoint t, standardize the probability Jacobian by sqrt(N_t) diag(1/sqrt(pi_t)). Stack the two initial-logit Jacobians into A_eta and the log-birth Jacobians into A_g (the initial-time birth block is zero). An SVD supplies an orthonormal basis Q_eta for the data-visible nuisance-column space. The projected design is R=A_g−Q_eta(Q_eta^T A_g), and the profiled expected multinomial Fisher matrix is F=R^T R. This PSD construction does not subtract nearly equal Fisher blocks or use an unstable inverse of the logit gauge.", "",
  "Let K be the complete frozen Gaussian prior covariance of g, including the feature coefficients and graph residual, and let K=L L^T. The prior-whitened matrix is (R L)^T(R L), with eigenvalues lambda_j. The reported effective learned directions are sum_j lambda_j/(1+lambda_j): a local Gaussian, prior-relative variance-learning diagnostic. It is not a literal count of globally identifiable parameters. Profiling eta omits its proper prior regularization, so this is not the full fitted posterior precision.", "",
  "## Information and depth planning", "",
  "Both timepoint depths are scaled together from 150 to 1,000 or 5,000, keeping the expansion point and probabilities fixed. Fisher information scales linearly with depth; the learned fractions saturate nonlinearly. Fixing pi0 at the same mode removes nuisance uncertainty at the identical expansion point and must increase information in PSD order. Replacing pi0 by its true value also changes the expansion point, and is reported separately rather than called a pure same-point information gain.", "",
  "| Family | Initial-frequency assumption | Learned dirs, 150 | 1,000 | 5,000 | Shape dirs, 150 | Shape dirs, 5,000 | Data rank |",
  "|---|---|---:|---:|---:|---:|---:|---:|")
for (family in c("smooth", "rough", "interaction")) for (scenario in c("unknown_pi0_profile", "fixed_pi0_same_mode"))
  lines <- c(lines, sprintf("| %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %d |", family, labels[scenario],
    value(family, scenario, 150, "effective_learned_directions"),
    value(family, scenario, 1000, "effective_learned_directions"),
    value(family, scenario, 5000, "effective_learned_directions"),
    value(family, scenario, 150, "shape_effective_directions"),
    value(family, scenario, 5000, "shape_effective_directions"),
    value(family, scenario, 150, "data_rank")))
lines <- c(lines, "", "| Family | True-pi0 sensitivity: learned dirs, 150 | 5,000 | Local data rank at 150 |",
  "|---|---:|---:|---:|")
for (family in c("smooth", "rough", "interaction"))
  lines <- c(lines, sprintf("| %s | %.3f | %.3f | %d |", family,
    value(family, "fixed_true_pi0_same_birth", 150, "effective_learned_directions"),
    value(family, "fixed_true_pi0_same_birth", 5000, "effective_learned_directions"),
    value(family, "fixed_true_pi0_same_birth", 150, "data_rank")))
lines <- c(lines, "", "## Physical null and shape directions", "",
  "With pi0 fixed, one final composition contains only n−1 independent probabilities for n free node birth rates. Consequently the data Fisher rank is at most 26 here; beta>0 does not remove this dimension bound. At beta=0, a common additive birth shift leaves frequencies unchanged, and its log-birth tangent is 1/b, not the all-ones vector. A common log-birth shift instead multiplies all rates and generally changes the evolutionary speed. Mutation can rotate the exact data-null direction into a mixture of baseline and shape effects.", "",
  "For the shape spectrum, map the common additive-birth direction into prior-whitened coordinates, a=L^(-1)(1/b)/||L^(-1)(1/b)||. Choose an orthonormal basis Q_shape perpendicular to a, then diagonalize Q_shape^T F_white Q_shape. This holds that particular prior-metric absolute coordinate fixed; it is not another fit or an unconditional separation of biological shape and absolute growth. Shape and absolute Rayleigh summaries are coupled and must not be added as independent information counts.", "",
  "| Family | Weakest direction: physical delta-birth alignment with a constant | Fixed-shape absolute Rayleigh at 150 |",
  "|---|---:|---:|")
for (family in c("smooth", "rough", "interaction"))
  lines <- c(lines, sprintf("| %s | %.3f | %.3f |", family,
    value(family, "unknown_pi0_profile", 150, "weakest_physical_absolute_alignment"),
    value(family, "unknown_pi0_profile", 150, "fixed_shape_absolute_rayleigh")))
lines <- c(lines, "", "The weakest physical directions therefore mix a common baseline with shape. A positive Rayleigh value when shape is held fixed does not establish absolute-rate identification when all shape directions are free. Full node loadings of the three weakest directions are retained in the CSV.", "",
  "## Numerical and mathematical checks", "",
  sprintf("Halving the central-difference step changes the log-birth Jacobian by at most %.3g and the profiled Fisher matrix by at most %.3g in relative Frobenius norm. Analytic initial/final logit Jacobians agree with finite differences to %.3g. The independent continuous-time mode frequencies agree with the fitted numerical solver to %.3g. Known-pi0 minus profiled information has minimum eigenvalue %.3g, consistent with PSD to numerical precision.",
    max(checks$logbirth_jacobian_epsilon_halving_max_change), max(checks$profiled_fisher_relative_change),
    max(checks$pi0_eta_jacobian_max_error, checks$pi1_eta_jacobian_max_error),
    max(checks$mode_cpp_vs_independent_exact_pi1), min(checks$known_minus_profile_min_eigenvalue)), "",
  sprintf("The independent mathematical regression file passes %d assertions: analytic beta-zero derivatives and additive-birth null, nuisance gauge/PSD ordering, mutation changing the null without breaking the composition rank bound, and linear information scaling with saturating learned fractions. Machine-readable results and the test transcript are saved beside the information tables.", sum(test_summary$passed)), "",
  "## Limits", "",
  "This is local expected information under the specified multinomial forward model at three prior-influenced modes. It omits finite-population process variance, does not use observed-Hessian curvature or posterior multimodality, and provides no global identifiability proof. Whitening and effective learned directions depend on K; increasing a prior variance can increase this relative diagnostic without increasing data information. Depth extrapolation does not simulate new data, refit modes, remove model bias or certify future RMSE/coverage. This diagnostic must not be substituted for independent candidate validation or used to relabel pair-feature support rank as statistical identification.", "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
writeLines(lines, file.path(destination, "oracle_fisher_review.md"))
files <- c("oracle_fisher_summary.csv", "oracle_fisher_spectra.csv", "oracle_fisher_numerical_checks.csv", "oracle_fisher_weak_directions.csv")
invisible(file.copy(file.path(input, files), artifact, overwrite = TRUE))
cat(file.path(destination, "oracle_fisher_review.md"), "\n")
