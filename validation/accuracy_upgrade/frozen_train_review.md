# alfak2 accuracy development benchmark

These are development results on 6 landscapes and 24 paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.

The finite target domain contains 81 4-chromosome karyotypes. Every method is scored on the same frozen target nodes and reference observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.

Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.

Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. This small development sample is insufficient for a general claim of improvement at every distance.

## Completion

| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |
|---|---:|---:|---:|---:|
| baseline | 24 | 24 | 4 | 0.11 |
| division_v1 | 24 | 24 | 24 | 30.25 |
| endpoint_v1 | 24 | 24 | 24 | 0.09 |

## Fitness and interval results

| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |
|---|---|---:|---:|---:|---:|---:|
| baseline | d0 | 0.1379 | 0.1077 | 0.709 | 2.1510 | 1.000 |
| baseline | d1 | 0.1911 | 0.1551 | 0.998 | 1.2671 | 1.000 |
| baseline | d2 | 0.2521 | 0.2159 | 1.000 | 1.4803 | 1.000 |
| baseline | d3 | 0.2973 | 0.2605 | 1.000 | 1.6283 | 1.000 |
| baseline | all_node_micro | 0.2308 | 0.1839 | 0.950 | 1.5721 | 1.000 |
| baseline | all_distance_macro | 0.2196 | 0.1848 | 0.927 | 1.6317 | 1.000 |
| division_v1 | d0 | 0.1229 | 0.0957 | 0.907 | 0.6864 | 1.000 |
| division_v1 | d1 | 0.1822 | 0.1476 | 0.844 | 1.0330 | 1.000 |
| division_v1 | d2 | 0.2403 | 0.2050 | 0.765 | 1.6899 | 1.000 |
| division_v1 | d3 | 0.2880 | 0.2503 | 0.689 | 2.6563 | 1.000 |
| division_v1 | all_node_micro | 0.2196 | 0.1740 | 0.803 | 1.4387 | 1.000 |
| division_v1 | all_distance_macro | 0.2084 | 0.1746 | 0.801 | 1.5164 | 1.000 |
| endpoint_v1 | d0 | 0.1340 | 0.1065 | 0.912 | 0.7801 | 1.000 |
| endpoint_v1 | d1 | 0.1883 | 0.1557 | 0.872 | 0.9446 | 1.000 |
| endpoint_v1 | d2 | 0.2474 | 0.2123 | 0.796 | 1.4963 | 1.000 |
| endpoint_v1 | d3 | 0.2933 | 0.2530 | 0.698 | 2.6055 | 1.000 |
| endpoint_v1 | all_node_micro | 0.2261 | 0.1817 | 0.829 | 1.3158 | 1.000 |
| endpoint_v1 | all_distance_macro | 0.2158 | 0.1819 | 0.819 | 1.4566 | 1.000 |

## Independent held-out count predictions

Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.

| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |
|---|---:|---:|---:|---:|
| baseline | 0.7860 | 1.2226 | 0.00588 | 0.00623 |
| division_v1 | 0.7754 | 1.2024 | 0.00429 | 0.00411 |
| endpoint_v1 | 0.7739 | 1.2114 | 0.00503 | 0.00516 |

## Fitness RMSE by truth family

| Model | Family | d0 | d1 | d2 | d3 | Equal-distance macro |
|---|---|---:|---:|---:|---:|---:|
| baseline | smooth | 0.1040 | 0.1762 | 0.2484 | 0.2940 | 0.2057 |
| baseline | rough | 0.2587 | 0.2581 | 0.2550 | 0.2321 | 0.2510 |
| baseline | interaction | 0.0510 | 0.1390 | 0.2528 | 0.3657 | 0.2021 |
| division_v1 | smooth | 0.1101 | 0.1820 | 0.2542 | 0.2988 | 0.2113 |
| division_v1 | rough | 0.2042 | 0.2251 | 0.2137 | 0.1984 | 0.2104 |
| division_v1 | interaction | 0.0545 | 0.1394 | 0.2530 | 0.3667 | 0.2034 |
| endpoint_v1 | smooth | 0.1350 | 0.2108 | 0.2856 | 0.3336 | 0.2413 |
| endpoint_v1 | rough | 0.1980 | 0.2096 | 0.2007 | 0.1749 | 0.1958 |
| endpoint_v1 | interaction | 0.0692 | 0.1444 | 0.2559 | 0.3715 | 0.2102 |

## Paired candidate minus baseline differences

Negative RMSE and interval-score differences favor the candidate. Coverage differences have no universal favorable direction; assess distance from the nominal 0.95 and interval score jointly.

| Candidate | Distance | Metric | Mean difference | 95% landscape bootstrap interval | Landscapes |
|---|---|---|---:|---|---:|
| division_v1 | all_distance_macro | coverage_95 | -0.1259 | [-0.2007, -0.0206] | 6 |
| endpoint_v1 | all_distance_macro | coverage_95 | -0.1074 | [-0.1810, -0.0265] | 6 |
| division_v1 | d0 | coverage_95 | 0.1971 | [0.0758, 0.3668] | 6 |
| endpoint_v1 | d0 | coverage_95 | 0.2026 | [0.0955, 0.3196] | 6 |
| division_v1 | d1 | coverage_95 | -0.1538 | [-0.2211, -0.0895] | 6 |
| endpoint_v1 | d1 | coverage_95 | -0.1259 | [-0.1751, -0.0733] | 6 |
| division_v1 | d2 | coverage_95 | -0.2355 | [-0.3087, -0.1471] | 6 |
| endpoint_v1 | d2 | coverage_95 | -0.2044 | [-0.2772, -0.1362] | 6 |
| division_v1 | d3 | coverage_95 | -0.3112 | [-0.4634, -0.1623] | 6 |
| endpoint_v1 | d3 | coverage_95 | -0.3018 | [-0.4617, -0.1431] | 6 |
| division_v1 | all_distance_macro | interval_score_95 | -0.1153 | [-1.0165, 0.5753] | 6 |
| endpoint_v1 | all_distance_macro | interval_score_95 | -0.1751 | [-0.9691, 0.4475] | 6 |
| division_v1 | d0 | interval_score_95 | -1.4646 | [-3.8869, -0.0489] | 6 |
| endpoint_v1 | d0 | interval_score_95 | -1.3709 | [-3.7043, 0.0425] | 6 |
| division_v1 | d1 | interval_score_95 | -0.2341 | [-0.4122, -0.0488] | 6 |
| endpoint_v1 | d1 | interval_score_95 | -0.3226 | [-0.4737, -0.1638] | 6 |
| division_v1 | d2 | interval_score_95 | 0.2096 | [-0.1630, 0.5252] | 6 |
| endpoint_v1 | d2 | interval_score_95 | 0.0160 | [-0.3156, 0.3483] | 6 |
| division_v1 | d3 | interval_score_95 | 1.0280 | [-0.3113, 2.2301] | 6 |
| endpoint_v1 | d3 | interval_score_95 | 0.9771 | [-0.0894, 2.0083] | 6 |
| division_v1 | all_distance_macro | relative_rmse | -0.0112 | [-0.0403, 0.0052] | 6 |
| endpoint_v1 | all_distance_macro | relative_rmse | -0.0038 | [-0.0496, 0.0264] | 6 |
| division_v1 | d0 | relative_rmse | -0.0150 | [-0.0523, 0.0051] | 6 |
| endpoint_v1 | d0 | relative_rmse | -0.0039 | [-0.0541, 0.0267] | 6 |
| division_v1 | d1 | relative_rmse | -0.0090 | [-0.0328, 0.0049] | 6 |
| endpoint_v1 | d1 | relative_rmse | -0.0028 | [-0.0436, 0.0252] | 6 |
| division_v1 | d2 | relative_rmse | -0.0118 | [-0.0408, 0.0054] | 6 |
| endpoint_v1 | d2 | relative_rmse | -0.0047 | [-0.0487, 0.0261] | 6 |
| division_v1 | d3 | relative_rmse | -0.0093 | [-0.0344, 0.0054] | 6 |
| endpoint_v1 | d3 | relative_rmse | -0.0039 | [-0.0528, 0.0290] | 6 |

## Limits and next gate

The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.

Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.

Generated: 2026-11-03 23:50:02 EST
