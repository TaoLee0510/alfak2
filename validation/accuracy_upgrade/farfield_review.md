# alfak2 accuracy development benchmark

These are development results on 1 landscapes and 1 paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.

The finite target domain contains 243 5-chromosome karyotypes. Every method is scored on the same frozen target nodes and reference observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.

Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.

Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. This small development sample is insufficient for a general claim of improvement at every distance.

## Completion

| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |
|---|---:|---:|---:|---:|
| baseline | 1 | 1 | 1 | 0.41 |

## Fitness and interval results

| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |
|---|---|---:|---:|---:|---:|---:|
| baseline | d0 | 0.0000 | 0.0000 | 1.000 | 0.0000 | 1.000 |
| baseline | d1 | 0.0268 | 0.0120 | 1.000 | 97.4412 | 1.000 |
| baseline | d2 | 0.1020 | 0.0600 | 1.000 | 137.7891 | 1.000 |
| baseline | d3 | 0.1764 | 0.1390 | 1.000 | 168.7598 | 1.000 |
| baseline | d4 | 0.2546 | 0.2263 | 1.000 | 194.8835 | 1.000 |
| baseline | d5 | 0.3371 | 0.2814 | 1.000 | 217.9169 | 1.000 |
| baseline | all_node_micro | 0.2197 | 0.1677 | 1.000 | 175.1061 | 1.000 |
| baseline | all_distance_macro | 0.1495 | 0.1198 | 1.000 | 136.1317 | 1.000 |

## Independent held-out count predictions

Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.

| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |
|---|---:|---:|---:|---:|
| baseline | 0.0005 | 0.0005 | 0.00003 | 0.00003 |

## Fitness RMSE by truth family

| Model | Family | d0 | d1 | d2 | d3 | d4 | d5 | Equal-distance macro |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| baseline | interaction | 0.0000 | 0.0268 | 0.1020 | 0.1764 | 0.2546 | 0.3371 | 0.1495 |

## Limits and next gate

The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.

Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.

Generated: 2026-11-04 00:00:59 EST
