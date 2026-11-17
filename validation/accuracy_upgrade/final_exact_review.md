# alfak2: final independent exact-division validation

The frozen joint division candidate lowers average relative-fitness RMSE at every evaluated distance, but 4/4 distance-specific paired 95% intervals still include zero. Family-specific regressions and calibrated interval-score tradeoffs remain. These results do not establish uniformly better estimation for every karyotype or scenario.

The strict all-distance point-and-uncertainty promotion gate is not passed. The joint candidate remains an explicit opt-in; the default pipeline remains two_layer. No tuning or extra test repetitions were added after these results.

The test contains 12 independent landscapes and 48 datasets: three truth families, four reserved landscapes per family, two biological trajectories and two independent count samples per landscape. The entire 81-node four-chromosome CN1–3 domain is scored. Beta=0.005, interval=1.5, initial population=1,200 and count depths=150/150 are fixed.

## Point estimation: raw candidate versus raw baseline

All models use one t0-count-weighted observed-reference contrast. No per-distance centering, affine truth correction or calibrated change to a point prediction is used. Reported means first average biological/sampling repeats within each landscape, then weight landscapes equally. The confidence intervals resample those independent paired landscapes, with the prespecified 1,000-replicate bootstrap. They condition on the frozen candidate and calibration. One test dataset has no d3 target nodes; this is an absent stratum, not a missing prediction.

| Distance | Original RMSE | Joint division RMSE | RMSE reduction | Candidate minus original [95% landscape interval] | Paired landscapes/datasets |
|---|---:|---:|---:|---|---:|
| d0 | 0.15979 | 0.14735 | 7.8% | -0.01244 [-0.03426, 0.00406] | 12/48 |
| d1 | 0.20766 | 0.19535 | 5.9% | -0.01232 [-0.03081, 0.00174] | 12/48 |
| d2 | 0.25550 | 0.24412 | 4.5% | -0.01138 [-0.02853, 0.00159] | 12/48 |
| d3 | 0.29420 | 0.28174 | 4.2% | -0.01246 [-0.03179, 0.00059] | 12/47 |

## Uncertainty: original and candidate, before and after independent calibration

Both methods receive their own empirical calibration on the same separate six-landscape/24-dataset exact-simulation calibration split. Calibration changes interval endpoints only. It leaves every fitness prediction and SD unchanged. Proper interval scores penalize both excessive width and missed truth; lower is better. Generic calibration gains are not attributed to the new estimator.

| Distance | Model | 95% coverage | 95% interval score | Interval availability |
|---|---|---:|---:|---:|
| d0 | baseline | 0.622 | 2.78987 | 1.000 |
| d0 | division_v1 | 0.870 | 1.09178 | 1.000 |
| d0 | baseline_calibrated | 0.862 | 1.42458 | 1.000 |
| d0 | division_v1_calibrated | 0.959 | 1.03252 | 1.000 |
| d1 | baseline | 0.997 | 1.26803 | 1.000 |
| d1 | division_v1 | 0.843 | 1.18439 | 1.000 |
| d1 | baseline_calibrated | 0.868 | 1.64317 | 1.000 |
| d1 | division_v1_calibrated | 0.928 | 1.34307 | 1.000 |
| d2 | baseline | 0.997 | 1.47884 | 1.000 |
| d2 | division_v1 | 0.772 | 1.62747 | 1.000 |
| d2 | baseline_calibrated | 0.937 | 1.17902 | 1.000 |
| d2 | division_v1_calibrated | 0.962 | 1.23115 | 1.000 |
| d3 | baseline | 1.000 | 1.62914 | 1.000 |
| d3 | division_v1 | 0.723 | 2.30846 | 1.000 |
| d3 | baseline_calibrated | 0.935 | 1.41257 | 1.000 |
| d3 | division_v1_calibrated | 0.927 | 1.47378 | 1.000 |

The matched-calibration comparison is candidate_calibrated minus baseline_calibrated:

| Distance | Proper interval-score difference [95% landscape interval] | Paired landscapes |
|---|---|---:|
| d0 | -0.39206 [-1.21962, 0.20066] | 12 |
| d1 | -0.30010 [-1.09651, 0.19441] | 12 |
| d2 | 0.05212 [-0.05001, 0.15394] | 12 |
| d3 | 0.06121 [-0.09972, 0.21189] | 12 |

## Family-specific point errors

| Truth family | Distance | Original RMSE | Joint division RMSE | RMSE reduction |
|---|---|---:|---:|---:|
| smooth | d0 | 0.11413 | 0.11870 | -4.0% |
| smooth | d1 | 0.18697 | 0.18745 | -0.3% |
| smooth | d2 | 0.24880 | 0.24826 | 0.2% |
| smooth | d3 | 0.29161 | 0.28982 | 0.6% |
| rough | d0 | 0.31736 | 0.27153 | 14.4% |
| rough | d1 | 0.29279 | 0.25432 | 13.1% |
| rough | d2 | 0.27910 | 0.24379 | 12.6% |
| rough | d3 | 0.27054 | 0.23531 | 13.0% |
| interaction | d0 | 0.04789 | 0.05181 | -8.2% |
| interaction | d1 | 0.14323 | 0.14427 | -0.7% |
| interaction | d2 | 0.23861 | 0.24031 | -0.7% |
| interaction | d3 | 0.32045 | 0.32009 | 0.1% |

Remaining family/distance RMSE regressions: smooth d0, smooth d1, interaction d0, interaction d1, interaction d2.

## Count predictions and completion

Held-out multinomial count samples use a separate fixed seed and never enter fitting. These are per-cell log scores of the returned point probability forecasts, with a declared 1e-12 floor; they are not a fully integrated joint count predictive likelihood. Zero-probability hits and frequency RMSE are retained in the CSV.

| Model | t0 NLL/cell | t1 NLL/cell | Successful fits | Trusted local covariance |
|---|---:|---:|---:|---:|
| baseline | 0.74747 | 1.16925 | 48/48 | 7/48 |
| endpoint_v1 | 0.73601 | 1.16469 | 48/48 | 48/48 |
| division_v1 | 0.74091 | 1.15571 | 48/48 | 48/48 |

Failure and missing-interval cases remain in fit-status and availability denominators. Error averages and paired intervals use finite comparable predictions, so incomplete availability cannot be interpreted as successful universal accuracy.

## Freeze record and limits

Calibration frozen at: 2026-11-04 04:55:20 UTC ; the saved freeze record has test_generated=FALSE. The test landscape and sampling seeds were reserved before development. Candidate-v1 source/DLL and calibration-object checksums are retained in [the protocol artifact directory](results/protocol_v2/).

The initial tau-leap simulator has measurable finite-step growth bias. Formal calibration/test therefore use the independently implemented exact Gillespie process, chosen before test generation. The historical tau outcomes are retained as robustness evidence rather than overwritten. See [protocol v2](protocol_v2.md) and [simulator validation](results/simulator_validation.csv).

This result is confined to the specified simulation families, depth, beta, time interval and bounded target domain. It does not validate a complete 22-chromosome landscape or real biological measurements. Absolute birth-rate level remains prior-identified without external growth data. Calibration uses six independent landscapes and provides empirical intervals only, without a finite-sample or distribution-free coverage guarantee. Unknown distances/scenarios are explicitly unvalidated.

## Retained negative results and audit trail

- [Original pilot and Laplace ablation](pilot_review.md): repaired convergence did not automatically improve every fitness error.
- [Shared-scale and evidence-mixture screen](scale_grid_review.md): point/interval and truth-family tradeoffs remain visible.
- [Separated basis/residual screen](basis_residual_review.md) and [eight-dimensional Pareto table](results/basis_residual_pareto.csv).
- [Frozen tau robustness comparison](frozen_train_review.md): retained separately from the exact independent test.
- [Five-chromosome d0–d5 support-limit probe](farfield_review.md): remote states are not identified by one founder-only observation.
- [Machine-readable final results](results/final_exact/): paired raw and matched-calibration comparisons, all distance/family metrics, fit statuses, group completeness, the reserved test manifest and complete node-wise prediction/SD/interval matrices. Matrix checks verify calibrated and raw point predictions/SDs are identical.

Generated: 2026-11-04 00:07:14 EST
