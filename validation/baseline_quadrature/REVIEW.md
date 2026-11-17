# Growth-scale integration: numerical and truth-development checks

The numerical candidate preserves the original Gaussian prior and integrates its existing independent intercept coefficient. It uses no additional growth observation in count-only fits. This separates a computational correction from changing the biological model or claiming new information.

## Small-graph numerical reference

All 16/16 cases pass quadrature refinement with every positive-weight inner node trusted. Against the independently checked MCMC references, relative-fitness SD ratios range 0.923–1.035, versus 0.994–3.284 for the original Gaussian/Laplace approximation. The largest standardized mean error falls from 1.379 to 0.122. Two nonzero-beta sparse count-only cases retain material mean discrepancies under the predeclared criterion; the approximation is improved but not exact.

These are fixed-model computational comparisons using an explicitly fixed SD-2 initial-logit prior. They are not generating-truth RMSE, repeated-sampling coverage, SBC or a blanket deployment guarantee. [Reference report](../posterior_reference_4to8/REVIEW.md); the complete per-case numerical tables are retained in `results/`.

## Deployment-prior truth pilot

Twelve unchanged, revealed development samples represent six finite-population landscapes, two per family. Both estimators retain deployment eta priors, the same V1 fitness prior, complete 27-node domain and native support distances. Every V1 mean/covariance is replayed and every conditional eta-prior vector is checked against the frozen fit. Both models receive 40,000 actual joint-distribution draws for reference-contrast intervals, with common random numbers. No calibration is learned here.

All 12/12 integration fits pass numerical refinement. Median integration time is 25.60s per dataset with at most two fit processes. Input hashes, exact consumed datasets/frozen fits and node-level outputs are archived in the two pilot-part directories.

Risks average sampling repeats within landscape, landscapes within family and families equally. Distance -1 covers the whole target domain. Shape removes one spatial error mean over the whole domain, never an individual shell.

| Distance | Model | RMSE | Shape RMSE | Edge RMSE | Raw IS95 | Coverage95 |
|---:|---|---:|---:|---:|---:|---:|
| -1 | baseline_quadrature | 0.19957 | 0.17298 | 0.20770 | 1.04444 | 0.926 |
| -1 | v1_joint | 0.20426 | 0.17320 | 0.20828 | 1.27001 | 0.886 |
| 0 | baseline_quadrature | 0.14031 | 0.12337 | 0.14889 | 0.82814 | 0.936 |
| 0 | v1_joint | 0.14427 | 0.12388 | 0.15161 | 1.08972 | 0.894 |
| 1 | baseline_quadrature | 0.20754 | 0.17926 | 0.21371 | 1.01265 | 0.940 |
| 1 | v1_joint | 0.21284 | 0.17908 | 0.21355 | 1.24317 | 0.890 |
| 2 | baseline_quadrature | 0.20337 | 0.18623 | 0.20447 | 1.36254 | 0.900 |
| 2 | v1_joint | 0.20655 | 0.18607 | 0.20430 | 1.44429 | 0.883 |

## Descriptive family/distance safeguards

Ratios above 1.05 flag observed harm relative to frozen V1. These six development landscapes cannot establish formal noninferiority; all cells remain visible.

| Family | Distance | RMSE ratio | Shape ratio | Edge ratio | Raw IS95 ratio |
|---|---:|---:|---:|---:|---:|
| interaction | -1 | 0.9960 | 1.0011 | 1.0008 | 0.9607 |
| rough | -1 | 0.9747 | 0.9961 | 0.9938 | 0.6547 |
| smooth | -1 | 0.9536 | 0.9982 | 0.9987 | 0.9513 |
| interaction | 0 | 1.0019 | 1.0136 | 1.0188 | 1.0451 |
| rough | 0 | 0.9701 | 0.9861 | 0.9812 | 0.5422 |
| smooth | 0 | 0.9478 | 0.9940 | 0.9512 | 1.0668 |
| interaction | 1 | 0.9875 | 1.0003 | 1.0005 | 0.9206 |
| rough | 1 | 0.9781 | 1.0022 | 0.9978 | 0.6701 |
| smooth | 1 | 0.9536 | 1.0006 | 1.0100 | 0.9840 |
| interaction | 2 | 1.0216 | 1.0038 | 1.0012 | 0.9475 |
| rough | 2 | 0.9533 | 0.9986 | 1.0005 | 1.1235 |
| smooth | 2 | 0.9533 | 0.9970 | 1.0005 | 0.8109 |

## Decision boundaries

The stage remains closed: no all-distance shape/uncertainty winner. Distances 1 and 2 have slightly higher centered shape/edge error. Raw interval-score ratios exceed 1.05 in smooth-d0 (1.0668) and rough-d2 (1.1235). The common-reference offset reduction explains 97.08% of the whole-domain MSE gain; the shape contribution is positive but small. This is a development stopping decision, not a formal independent-test rejection.

A more faithful posterior computation does not automatically improve estimation of truths outside the model's structural prior. The pilot is used only to decide whether an all-distance shape/uncertainty claim is plausible enough for further development. No production default is promoted, no new accuracy calibration/test is opened, and no final test N/power is certified. Numerical correction alone cannot establish the original all-distance accuracy objective.

The most important remaining distinctions are (i) residual inner-Laplace error for rare initial states, (ii) biological-prior mismatch and (iii) information absent from sparse two-timepoint observations. These must remain separate in any subsequent method or observation-design change.

Generated: 2026-11-04 14:57:40 EST
