# Independent posterior reference: 4 and 8 nodes

16/16 fixed-data references meet all preregistered Monte Carlo diagnostics. All 16 Laplace fits return trusted TMB covariance. This checkpoint compares algorithms under a fixed model; it does not certify biological fitness accuracy or repeated-sampling interval coverage.

## Scope and checks

The full CN1–2 domains, counts, prior scales, beta values and information arms were fixed in [the protocol](PROTOCOL.md) before the primary comparisons. There are four independent MCMC chains per case. The exact aggregate Gaussian log-birth prior is analytically marginalized; the softmax mean-logit gauge is removed with the correct contrast covariance. Fitness priors use V1 scales; every initial-logit SD is explicitly 2. These are fixed-data mechanism probes, not an assessment of all deployment configurations.

The independent sampler/target validation passes 99 assertions with no failures, errors or warnings. It includes analytic Gaussian-prior moments, the earlier converged two-node quadrature, TMB pointwise likelihood agreement at 20 off-mode states, and explicit-versus-collapsed Gaussian-prior Laplace agreement on four graph/beta combinations.

Across primary cases, maximum rank-normalized split Rhat is 1.00431; minimum bulk/tail ESS is 1697/2331. Retained draws per chain range from 20000 to 160000. Maximum primary-mean/SD MCSE divided by reference SD is 0.025/0.024; maximum relative covariance MCSE is 0.041. The case status table preserves failed attempts before any allowed length extension.

## Relative fitness comparison

Every birth rate is contrasted with the same fixed t0-count-weighted reference. SD ratios are Laplace/MCMC. Mean errors are divided by the reference posterior SD. A material flag requires the predeclared numerical discrepancy threshold and a difference larger than four MCMC standard errors; it is not a population hypothesis test.

| Nodes | Information | Certified cases | Relative SD ratio range | Largest absolute mean error / reference SD | Cases with material mean / SD error |
|---:|---|---:|---:|---:|---:|
| 4 | Counts only | 4/4 | 1.417–1.817 | 0.996 | 4 / 4 |
| 4 | Counts + growth | 4/4 | 0.994–1.023 | 0.054 | 0 / 0 |
| 8 | Counts only | 4/4 | 1.562–3.284 | 1.379 | 4 / 4 |
| 8 | Counts + growth | 4/4 | 0.996–1.045 | 0.072 | 2 / 0 |

The information-arm comparison probes whether constraining common growth reduces approximation error. Trusted Hessians alone do not diagnose that error. Growth adds information; any better agreement is a mechanism probe, not a fair algorithmic promotion over count-only data. These cases cannot attribute all far-field undercoverage to Laplace: a broad model posterior can still miss truths outside a restrictive model/prior.

## All cases

| Case | Reference certified | Draws/chain | Relative SD ratio | Max mean error / SD | Relative covariance error |
|---|---|---:|---:|---:|---:|
| n4_balanced_beta0_counts | TRUE | 20000 | 1.506–1.564 | 0.487 | 1.399 |
| n4_balanced_beta0_growth | TRUE | 20000 | 0.994–1.020 | 0.016 | 0.032 |
| n4_balanced_beta0p005_counts | TRUE | 20000 | 1.516–1.616 | 0.506 | 1.547 |
| n4_balanced_beta0p005_growth | TRUE | 20000 | 1.003–1.023 | 0.026 | 0.034 |
| n4_sparse_beta0_counts | TRUE | 40000 | 1.610–1.817 | 0.996 | 2.111 |
| n4_sparse_beta0_growth | TRUE | 40000 | 1.008–1.021 | 0.044 | 0.029 |
| n4_sparse_beta0p005_counts | TRUE | 40000 | 1.417–1.590 | 0.604 | 1.412 |
| n4_sparse_beta0p005_growth | TRUE | 20000 | 1.005–1.022 | 0.054 | 0.038 |
| n8_balanced_beta0_counts | TRUE | 20000 | 1.607–1.732 | 0.515 | 1.835 |
| n8_balanced_beta0_growth | TRUE | 20000 | 0.996–1.045 | 0.038 | 0.049 |
| n8_balanced_beta0p005_counts | TRUE | 20000 | 1.562–1.706 | 0.484 | 1.732 |
| n8_balanced_beta0p005_growth | TRUE | 20000 | 1.010–1.036 | 0.054 | 0.067 |
| n8_sparse_beta0_counts | TRUE | 160000 | 1.984–3.284 | 1.379 | 7.323 |
| n8_sparse_beta0_growth | TRUE | 160000 | 1.013–1.027 | 0.065 | 0.052 |
| n8_sparse_beta0p005_counts | TRUE | 160000 | 1.596–2.313 | 0.764 | 3.227 |
| n8_sparse_beta0p005_growth | TRUE | 80000 | 1.005–1.027 | 0.072 | 0.034 |

## Numerical forward checks

An independent converged continuous-time Taylor solver checks 128 representative retained states per case (2,048 total). Maximum frequency and log-growth discrepancies are 8.97e-05 and 0.000121; the maximum absolute log-likelihood discrepancy is 0.00502. These are sampled posterior-state checks, not a global forward-error bound. Primary MCMC uses the same numerical likelihood as TMB; the continuous solver is a separate numerical-sensitivity check.

## Interpretation and next decision

The tables separately compare latent log-rate moments, physical birth moments, relative/centered fitness and graph-edge contrasts. `moment_comparison.csv` decomposes total mean error into (i) using Laplace rather than reference log-rate moments and (ii) Gaussianizing the reference log-rate distribution before exponentiation. The terms sum exactly, but their magnitudes are descriptive rather than causal percentages. Complete covariance comparisons are retained; no shell-specific recentering is used.

The scalar `posterior_squared_loss_ratio = 1 + standardized_mean_difference^2` is the conditional model-posterior quadratic-loss penalty of using the approximate mean. It is not an observed-truth RMSE ratio. Quantile tables report fixed-data reference posterior mass inside approximate intervals, not repeated-sampling coverage; both reference and iid approximate-quantile MCSE are explicit.

The next computational candidate, if pursued, should integrate over common growth-scale uncertainty while preserving the full Gaussian prior and conditional shape uncertainty. Its accuracy must first be demonstrated against these references; changing an arbitrary centering constraint or merely rescaling intervals is not equivalent. Even an accurate posterior will not manufacture information about unvisited epistatic states. No production model, new accuracy calibration/test landscape, HPC job, commit or push was changed by this checkpoint.

Full raw chains are kept in `/private/tmp/alfak2_reference_4to8_20261104`, with SHA256 and paths in `results/case_status.csv`. Reproducible cases, compact posterior summaries, diagnostics, covariance/quantile tables, source fingerprints and validation results are retained here.

The original live runner emitted a trailing parser error after all 16 completed outputs were written: file-path options were edited while its already-parsed sampling loop was running. The saved results were subsequently verified case by case, including four-chain dimensions, retained lengths, convergence and raw-file SHA256 (`results/completion_integrity.txt`). Exact original computation sources were archived before those file-path edits. All target/sampler/case/fit/diagnostic function bodies remained unchanged, and a separate portable-path validation rerun passed the same 99 assertions. The current portable runner parses cleanly; the late process-exit error is recorded rather than hidden.

The subsequent prior-preserving scalar-integration candidate and deployment-prior development pilot are reported separately in [the candidate review](../baseline_quadrature/REVIEW.md). Numerical accuracy and generating-truth accuracy remain separate decisions.

Generated: 2026-11-04 14:57:39 EST
