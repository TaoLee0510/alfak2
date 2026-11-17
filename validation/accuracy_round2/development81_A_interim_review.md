# Round 2 finite-population candidate development comparison

Development screening only: 4 chromosomes, all 81 target nodes, 18 independent landscapes. Scores use each count sample's native observed-support distances, not common process-pair distances. Every old calibration/revealed-test input is explicitly relabeled development. No new calibration/test landscape is opened.

V1 is explicitly isotropic with the frozen small prior scales. A changes the interaction prior to its observed-design nullspace with unidentified SD 0.2. B adds the predeclared weak/rough residual components with prior weights 0.75/0.25. All use deployment eta priors, no population-growth/pi0 oracle and no generating-family model input. The current default two-layer API and an interval-free zero point predictor are separate controls.

Intervals for B sample the actual component mixture and then apply the shared reference contrast. Family labels are evaluation metadata only. Means first average count samples within a landscape, then landscapes within a family, then families equally.

## Completion

| Model | Returned | Successful | Trusted covariance | Median seconds |
|---|---:|---:|---:|---:|
| v1_joint | 18 | 18 | 18 | 0.000 |
| A | 18 | 18 | 18 | 18.831 |
| zero_point_predictor | 18 | 18 | 0 | 0.000 |

## Native-distance risk and shape

Spatial centering removes one error mean over the entire target domain before stratification. It never recenters individual shells. All unordered one-copy graph edges enter the gradient metric, assigned to the largest endpoint support distance. Distance -1 denotes the complete target domain.

| Distance | Model | RMSE | Centered shape RMSE | Edge gradient RMSE | Spearman | Amplitude ratio | 95% coverage | 95% interval score | Interval fraction |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | v1_joint | 0.22981 | 0.18073 | 0.18991 | 0.071 | 0.047 | 0.766 | 1.75373 | 1.000 |
| -1 | A | 0.22981 | 0.18073 | 0.18991 | 0.071 | 0.048 | 0.817 | 1.48379 | 1.000 |
| -1 | zero_point_predictor | 0.23580 | 0.18124 | 0.19081 | NA | 0.000 | NA | NA | 0.000 |
| 0 | v1_joint | 0.15000 | 0.12744 | 0.16531 | 0.125 | 0.199 | 0.848 | 1.26220 | 1.000 |
| 0 | A | 0.15000 | 0.12744 | 0.16530 | 0.128 | 0.199 | 0.853 | 1.26583 | 1.000 |
| 0 | zero_point_predictor | 0.15701 | 0.12943 | 0.17117 | NA | 0.000 | NA | NA | 0.000 |
| 1 | v1_joint | 0.20048 | 0.15174 | 0.16642 | 0.027 | 0.032 | 0.802 | 1.53515 | 1.000 |
| 1 | A | 0.20048 | 0.15174 | 0.16642 | 0.027 | 0.032 | 0.822 | 1.43528 | 1.000 |
| 1 | zero_point_predictor | 0.20665 | 0.15150 | 0.16687 | NA | 0.000 | NA | NA | 0.000 |
| 2 | v1_joint | 0.25362 | 0.20018 | 0.19690 | 0.074 | 0.024 | 0.730 | 1.98355 | 1.000 |
| 2 | A | 0.25361 | 0.20018 | 0.19690 | 0.074 | 0.024 | 0.802 | 1.60785 | 1.000 |
| 2 | zero_point_predictor | 0.25903 | 0.20040 | 0.19692 | NA | 0.000 | NA | NA | 0.000 |
| 3 | v1_joint | 0.27756 | 0.23693 | 0.21978 | 0.071 | 0.024 | 0.671 | 2.63899 | 1.000 |
| 3 | A | 0.27756 | 0.23694 | 0.21978 | 0.072 | 0.024 | 0.777 | 1.88498 | 1.000 |
| 3 | zero_point_predictor | 0.28210 | 0.23757 | 0.21991 | NA | 0.000 | NA | NA | 0.000 |

## Candidate versus V1 family safety screen

These are descriptive ratios, not a formal noninferiority test. Ratios above 1.05 flag the predeclared 5% family safety margin. All cells remain visible; this pilot does not certify the final test sample size.

| Family | Distance | Candidate | RMSE ratio | Shape RMSE ratio | Edge RMSE ratio | Proper interval-score ratio |
|---|---:|---|---:|---:|---:|---:|
| interaction | -1 | A | 1.0000 | 1.0000 | 1.0000 | 0.8600 |
| interaction | 0 | A | 1.0002 | 1.0002 | 1.0001 | 0.9986 |
| interaction | 1 | A | 1.0000 | 1.0000 | 1.0000 | 1.0185 |
| interaction | 2 | A | 1.0000 | 1.0000 | 1.0000 | 0.8587 |
| interaction | 3 | A | 1.0000 | 1.0000 | 1.0000 | 0.6988 |
| rough | -1 | A | 1.0000 | 1.0000 | 1.0000 | 0.8960 |
| rough | 0 | A | 0.9999 | 1.0000 | 0.9999 | 1.0010 |
| rough | 1 | A | 1.0000 | 1.0000 | 1.0000 | 0.9476 |
| rough | 2 | A | 1.0000 | 1.0000 | 1.0000 | 0.8371 |
| rough | 3 | A | 1.0000 | 1.0000 | 1.0000 | 0.8320 |
| smooth | -1 | A | 1.0000 | 1.0000 | 1.0000 | 0.7787 |
| smooth | 0 | A | 1.0001 | 1.0000 | 1.0001 | 1.0162 |
| smooth | 1 | A | 1.0000 | 1.0000 | 1.0000 | 0.8669 |
| smooth | 2 | A | 1.0000 | 1.0000 | 1.0000 | 0.7572 |
| smooth | 3 | A | 1.0000 | 1.0000 | 1.0000 | 0.6633 |

This is a development comparison, not an independent test or SBC study. Reused V1 results retain their frozen source-cache SHA and unchanged counts/reference; zero elapsed time denotes reuse, not a newly timed fit. The current default is a separate actual API control, not an endpoint-ablation surrogate. No prior weight is adjusted from a sealed test. A single winner, if any, must be frozen after development and evaluated once on new calibration/test landscapes. Calibration hyperparameter heldout scores used to select gamma/kappa are not by themselves nested OOF loss vectors for joint power certification.

Generated: 2026-11-04 01:05:18 EST
