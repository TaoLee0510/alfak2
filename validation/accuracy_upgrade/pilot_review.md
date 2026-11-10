# alfak2 accuracy development pilot

These are development results on six landscapes and 24 paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.

The finite target domain is all 81 four-chromosome karyotypes with copy numbers 1–3. Every method is scored on the same frozen target nodes and observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.

Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.

Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. Six landscapes are insufficient for a general claim of improvement at every distance.

## Completion

| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |
|---|---:|---:|---:|---:|
| ablation_laplace | 24 | 24 | 24 | 0.11 |
| baseline | 24 | 24 | 4 | 0.12 |

## Fitness and interval results

| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |
|---|---|---:|---:|---:|---:|---:|
| ablation_laplace | d0 | 0.1423 | 0.1171 | 1.000 | 0.8662 | 1.000 |
| ablation_laplace | d1 | 0.2012 | 0.1700 | 1.000 | 1.3856 | 1.000 |
| ablation_laplace | d2 | 0.2564 | 0.2209 | 1.000 | 1.5703 | 1.000 |
| ablation_laplace | d3 | 0.3013 | 0.2601 | 1.000 | 1.7043 | 1.000 |
| ablation_laplace | all_node_micro | 0.2354 | 0.1925 | 1.000 | 1.4089 | 1.000 |
| ablation_laplace | all_distance_macro | 0.2253 | 0.1920 | 1.000 | 1.3816 | 1.000 |
| baseline | d0 | 0.1379 | 0.1077 | 0.709 | 2.1510 | 1.000 |
| baseline | d1 | 0.1911 | 0.1551 | 0.998 | 1.2671 | 1.000 |
| baseline | d2 | 0.2521 | 0.2159 | 1.000 | 1.4803 | 1.000 |
| baseline | d3 | 0.2973 | 0.2605 | 1.000 | 1.6283 | 1.000 |
| baseline | all_node_micro | 0.2308 | 0.1839 | 0.950 | 1.5721 | 1.000 |
| baseline | all_distance_macro | 0.2196 | 0.1848 | 0.927 | 1.6317 | 1.000 |

## Independent held-out count predictions

Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.

| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |
|---|---:|---:|---:|---:|
| ablation_laplace | 0.7703 | 1.2195 | 0.00404 | 0.00453 |
| baseline | 0.7860 | 1.2226 | 0.00588 | 0.00623 |

## Fitness RMSE by truth family

| Model | Family | d0 | d1 | d2 | d3 | Equal-distance macro |
|---|---|---:|---:|---:|---:|---:|
| ablation_laplace | smooth | 0.1666 | 0.2462 | 0.3218 | 0.3726 | 0.2768 |
| ablation_laplace | rough | 0.1735 | 0.2021 | 0.1851 | 0.1531 | 0.1785 |
| ablation_laplace | interaction | 0.0867 | 0.1553 | 0.2623 | 0.3783 | 0.2207 |
| baseline | smooth | 0.1040 | 0.1762 | 0.2484 | 0.2940 | 0.2057 |
| baseline | rough | 0.2587 | 0.2581 | 0.2550 | 0.2321 | 0.2510 |
| baseline | interaction | 0.0510 | 0.1390 | 0.2528 | 0.3657 | 0.2021 |

## Paired candidate minus baseline differences

Negative RMSE and interval-score differences favor the candidate. Coverage differences have no universal favorable direction; assess distance from the nominal 0.95 and interval score jointly.

| Candidate | Distance | Metric | Mean difference | 95% landscape bootstrap interval | Landscapes |
|---|---|---|---:|---|---:|
| ablation_laplace | all_distance_macro | coverage_95 | 0.0732 | [0.0273, 0.1335] | 6 |
| ablation_laplace | d0 | coverage_95 | 0.2906 | [0.1024, 0.5124] | 6 |
| ablation_laplace | d1 | coverage_95 | 0.0022 | [0.0000, 0.0067] | 6 |
| ablation_laplace | d2 | coverage_95 | 0.0000 | [0.0000, 0.0000] | 6 |
| ablation_laplace | d3 | coverage_95 | 0.0000 | [0.0000, 0.0000] | 6 |
| ablation_laplace | all_distance_macro | interval_score_95 | -0.2501 | [-0.9235, 0.1580] | 6 |
| ablation_laplace | d0 | interval_score_95 | -1.2848 | [-3.9684, 0.3542] | 6 |
| ablation_laplace | d1 | interval_score_95 | 0.1184 | [0.1097, 0.1268] | 6 |
| ablation_laplace | d2 | interval_score_95 | 0.0900 | [0.0831, 0.0961] | 6 |
| ablation_laplace | d3 | interval_score_95 | 0.0760 | [0.0714, 0.0805] | 6 |
| ablation_laplace | all_distance_macro | relative_rmse | 0.0057 | [-0.0600, 0.0536] | 6 |
| ablation_laplace | d0 | relative_rmse | 0.0044 | [-0.0715, 0.0529] | 6 |
| ablation_laplace | d1 | relative_rmse | 0.0100 | [-0.0464, 0.0530] | 6 |
| ablation_laplace | d2 | relative_rmse | 0.0044 | [-0.0595, 0.0516] | 6 |
| ablation_laplace | d3 | relative_rmse | 0.0041 | [-0.0784, 0.0583] | 6 |

## Limits and next gate

The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.

Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.

Generated: 2026-11-03 23:29:15 EST
