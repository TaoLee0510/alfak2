# Round 2 finite-population candidate pilot

Development screening only: three chromosomes, CN1–3, all 27 target nodes. Scores use each count sample's native observed-support distances, not the common D/G process-pair distances. The pilot has only two landscapes per truth family and two count samples per finite-population trajectory. No new calibration/test landscape is opened.

V1 is explicitly isotropic with the frozen small prior scales. A changes the interaction prior to its observed-design nullspace with unidentified SD 0.2. B adds the predeclared weak/rough residual components with prior weights 0.75/0.25. All use deployment eta priors, no population-growth/pi0 oracle and no generating-family model input. The current default two-layer API and an interval-free zero point predictor are separate controls.

Intervals for B sample the actual component mixture and then apply the shared reference contrast. Family labels are evaluation metadata only. Means first average count samples within a landscape, then landscapes within a family, then families equally.

## Completion

| Model | Returned | Successful | Trusted covariance | Median seconds |
|---|---:|---:|---:|---:|
| v1_joint | 12 | 12 | 12 | 0.952 |
| A | 12 | 12 | 12 | 0.953 |
| B | 12 | 12 | 12 | 1.931 |
| current_default_two_layer | 12 | 12 | 12 | 0.044 |
| zero_point_predictor | 12 | 12 | 0 | 0.000 |

## Native-distance risk and shape

Spatial centering removes one error mean over the entire target domain before stratification. It never recenters individual shells. All unordered one-copy graph edges enter the gradient metric, assigned to the largest endpoint support distance. Distance -1 denotes the complete target domain.

| Distance | Model | RMSE | Centered shape RMSE | Edge gradient RMSE | Spearman | Amplitude ratio | 95% coverage | 95% interval score | Interval fraction |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | v1_joint | 0.20426 | 0.17320 | 0.20828 | 0.049 | 0.087 | 0.886 | 1.28536 | 1.000 |
| -1 | A | 0.20424 | 0.17320 | 0.20826 | 0.050 | 0.087 | 0.901 | 1.24210 | 1.000 |
| -1 | B | 0.19991 | 0.17269 | 0.20595 | 0.052 | 0.114 | 0.944 | 1.12117 | 1.000 |
| -1 | current_default_two_layer | 0.26809 | 0.21205 | 0.23229 | 0.038 | 0.588 | 0.994 | 2.28438 | 1.000 |
| -1 | zero_point_predictor | 0.20899 | 0.17391 | 0.21050 | NA | 0.000 | NA | NA | 0.000 |
| 0 | v1_joint | 0.14427 | 0.12388 | 0.15161 | 0.161 | 0.251 | 0.894 | 1.10949 | 1.000 |
| 0 | A | 0.14425 | 0.12386 | 0.15158 | 0.161 | 0.252 | 0.894 | 1.10839 | 1.000 |
| 0 | B | 0.13975 | 0.12265 | 0.14604 | 0.179 | 0.338 | 0.937 | 0.83704 | 1.000 |
| 0 | current_default_two_layer | 0.18717 | 0.14248 | 0.19012 | 0.108 | 0.785 | 0.985 | 0.97981 | 1.000 |
| 0 | zero_point_predictor | 0.15080 | 0.12794 | 0.15940 | NA | 0.000 | NA | NA | 0.000 |
| 1 | v1_joint | 0.21284 | 0.17908 | 0.21355 | -0.097 | 0.043 | 0.896 | 1.25247 | 1.000 |
| 1 | A | 0.21283 | 0.17908 | 0.21353 | -0.097 | 0.043 | 0.909 | 1.17336 | 1.000 |
| 1 | B | 0.20824 | 0.17927 | 0.21107 | -0.081 | 0.045 | 0.966 | 1.09930 | 1.000 |
| 1 | current_default_two_layer | 0.28948 | 0.22576 | 0.24014 | -0.036 | 0.560 | 1.000 | 2.58268 | 1.000 |
| 1 | zero_point_predictor | 0.21476 | 0.17770 | 0.21633 | NA | 0.000 | NA | NA | 0.000 |
| 2 | v1_joint | 0.20655 | 0.18607 | 0.20430 | 0.197 | 0.053 | 0.869 | 1.46202 | 1.000 |
| 2 | A | 0.20651 | 0.18605 | 0.20430 | 0.197 | 0.054 | 0.914 | 1.43784 | 1.000 |
| 2 | B | 0.20470 | 0.18628 | 0.20436 | 0.197 | 0.057 | 0.935 | 1.54753 | 1.000 |
| 2 | current_default_two_layer | 0.27150 | 0.23488 | 0.22173 | 0.112 | 1.316 | 1.000 | 3.66488 | 1.000 |
| 2 | zero_point_predictor | 0.21628 | 0.18667 | 0.20369 | NA | 0.000 | NA | NA | 0.000 |

## Candidate versus V1 family safety screen

These are descriptive ratios, not a formal noninferiority test. Ratios above 1.05 flag the predeclared 5% family safety margin. All cells remain visible; this pilot does not certify the final test sample size.

| Family | Distance | Candidate | RMSE ratio | Shape RMSE ratio | Edge RMSE ratio | Proper interval-score ratio |
|---|---:|---|---:|---:|---:|---:|
| interaction | -1 | A | 1.0000 | 1.0000 | 1.0000 | 0.9117 |
| interaction | 0 | A | 1.0001 | 1.0001 | 1.0001 | 1.0047 |
| interaction | 1 | A | 1.0000 | 1.0000 | 1.0000 | 0.9041 |
| interaction | 2 | A | 1.0000 | 1.0000 | 1.0000 | 0.8804 |
| rough | -1 | A | 0.9999 | 1.0000 | 0.9999 | 0.9773 |
| rough | 0 | A | 0.9997 | 0.9997 | 0.9997 | 0.9961 |
| rough | 1 | A | 1.0000 | 1.0001 | 0.9998 | 0.8980 |
| rough | 2 | A | 0.9998 | 0.9997 | 1.0000 | 1.3585 |
| smooth | -1 | A | 0.9997 | 0.9999 | 0.9998 | 1.0440 |
| smooth | 0 | A | 1.0000 | 0.9998 | 0.9999 | 1.0000 |
| smooth | 1 | A | 0.9997 | 0.9999 | 0.9998 | 1.0933 |
| smooth | 2 | A | 0.9996 | 0.9999 | 0.9997 | 0.9674 |
| interaction | -1 | B | 0.9970 | 1.0011 | 1.0008 | 0.8666 |
| interaction | 0 | B | 1.0164 | 1.0353 | 1.0363 | 0.9214 |
| interaction | 1 | B | 0.9921 | 1.0008 | 1.0015 | 0.8404 |
| interaction | 2 | B | 1.0140 | 0.9996 | 0.9999 | 0.8766 |
| rough | -1 | B | 0.9480 | 0.9896 | 0.9743 | 0.7636 |
| rough | 0 | B | 0.9306 | 0.9543 | 0.9237 | 0.5690 |
| rough | 1 | B | 0.9539 | 1.0022 | 0.9733 | 0.7745 |
| rough | 2 | B | 0.9361 | 1.0060 | 1.0001 | 1.7841 |
| smooth | -1 | B | 0.9958 | 1.0004 | 1.0029 | 1.1382 |
| smooth | 0 | B | 1.0072 | 1.0035 | 1.0078 | 1.1590 |
| smooth | 1 | B | 0.9939 | 0.9999 | 1.0019 | 1.2042 |
| smooth | 2 | B | 0.9952 | 0.9997 | 1.0020 | 0.9883 |

The pilot cannot demonstrate improvement in the new 64-node d4/d5 domain or the 81-node domain. It is not an SBC study and does not include full posterior integration of count forecasts. No model/prior weight is adjusted from a sealed test. A single winner, if any, must be frozen after development and evaluated once on new calibration/test landscapes.

Generated: 2026-11-04 00:50:53 EST
