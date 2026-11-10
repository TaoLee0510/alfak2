# Round 2 finite-population candidate development comparison

Development screening only: 4 chromosomes, all 81 target nodes, 18 independent landscapes. Scores use each count sample's native observed-support distances, not common process-pair distances. Every old calibration/revealed-test input is explicitly relabeled development. No new calibration/test landscape is opened.

V1 is explicitly isotropic with the frozen small prior scales. A changes the interaction prior to its observed-design nullspace with unidentified SD 0.2. B adds the predeclared weak/rough residual components with prior weights 0.75/0.25. All use deployment eta priors, no population-growth/pi0 oracle and no generating-family model input. The current default two-layer API and an interval-free zero point predictor are separate controls.

Intervals for B sample the actual component mixture and then apply the shared reference contrast. Family labels are evaluation metadata only. Means first average count samples within a landscape, then landscapes within a family, then families equally.

## Completion

| Model | Returned | Successful | Trusted covariance | Median seconds |
|---|---:|---:|---:|---:|
| v1_joint | 18 | 18 | 18 | 0.000 |
| A | 18 | 18 | 18 | 18.831 |
| B | 18 | 18 | 18 | 36.242 |
| current_default_two_layer | 18 | 18 | 18 | 0.164 |
| zero_point_predictor | 18 | 18 | 0 | 0.000 |

## Native-distance risk and shape

Spatial centering removes one error mean over the entire target domain before stratification. It never recenters individual shells. All unordered one-copy graph edges enter the gradient metric, assigned to the largest endpoint support distance. Distance -1 denotes the complete target domain.

| Distance | Model | RMSE | Centered shape RMSE | Edge gradient RMSE | Spearman | Amplitude ratio | 95% coverage | 95% interval score | Interval fraction |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | v1_joint | 0.22981 | 0.18073 | 0.18991 | 0.071 | 0.047 | 0.766 | 1.75373 | 1.000 |
| -1 | A | 0.22981 | 0.18073 | 0.18991 | 0.071 | 0.048 | 0.817 | 1.48379 | 1.000 |
| -1 | B | 0.22515 | 0.18088 | 0.19028 | 0.072 | 0.075 | 0.924 | 1.17725 | 1.000 |
| -1 | current_default_two_layer | 0.26649 | 0.20739 | 0.20282 | 0.050 | 0.535 | 0.999 | 2.35109 | 1.000 |
| -1 | zero_point_predictor | 0.23580 | 0.18124 | 0.19081 | NA | 0.000 | NA | NA | 0.000 |
| 0 | v1_joint | 0.15000 | 0.12744 | 0.16531 | 0.125 | 0.199 | 0.848 | 1.26220 | 1.000 |
| 0 | A | 0.15000 | 0.12744 | 0.16530 | 0.128 | 0.199 | 0.853 | 1.26583 | 1.000 |
| 0 | B | 0.14947 | 0.12886 | 0.16379 | 0.124 | 0.330 | 0.949 | 0.79086 | 1.000 |
| 0 | current_default_two_layer | 0.18001 | 0.14788 | 0.19021 | 0.128 | 0.825 | 0.994 | 1.03099 | 1.000 |
| 0 | zero_point_predictor | 0.15701 | 0.12943 | 0.17117 | NA | 0.000 | NA | NA | 0.000 |
| 1 | v1_joint | 0.20048 | 0.15174 | 0.16642 | 0.027 | 0.032 | 0.802 | 1.53515 | 1.000 |
| 1 | A | 0.20048 | 0.15174 | 0.16642 | 0.027 | 0.032 | 0.822 | 1.43528 | 1.000 |
| 1 | B | 0.19574 | 0.15198 | 0.16808 | 0.030 | 0.036 | 0.936 | 1.04268 | 1.000 |
| 1 | current_default_two_layer | 0.23914 | 0.17985 | 0.18274 | 0.006 | 0.573 | 1.000 | 2.06884 | 1.000 |
| 1 | zero_point_predictor | 0.20665 | 0.15150 | 0.16687 | NA | 0.000 | NA | NA | 0.000 |
| 2 | v1_joint | 0.25362 | 0.20018 | 0.19690 | 0.074 | 0.024 | 0.730 | 1.98355 | 1.000 |
| 2 | A | 0.25361 | 0.20018 | 0.19690 | 0.074 | 0.024 | 0.802 | 1.60785 | 1.000 |
| 2 | B | 0.24850 | 0.20023 | 0.19707 | 0.076 | 0.025 | 0.912 | 1.32183 | 1.000 |
| 2 | current_default_two_layer | 0.29178 | 0.22732 | 0.20790 | 0.068 | 0.533 | 1.000 | 2.81676 | 1.000 |
| 2 | zero_point_predictor | 0.25903 | 0.20040 | 0.19692 | NA | 0.000 | NA | NA | 0.000 |
| 3 | v1_joint | 0.27756 | 0.23693 | 0.21978 | 0.071 | 0.024 | 0.671 | 2.63899 | 1.000 |
| 3 | A | 0.27756 | 0.23694 | 0.21978 | 0.072 | 0.024 | 0.777 | 1.88498 | 1.000 |
| 3 | B | 0.27142 | 0.23695 | 0.21987 | 0.081 | 0.025 | 0.873 | 1.72088 | 1.000 |
| 3 | current_default_two_layer | 0.32442 | 0.26353 | 0.22918 | 0.053 | 0.647 | 1.000 | 3.47082 | 1.000 |
| 3 | zero_point_predictor | 0.28210 | 0.23757 | 0.21991 | NA | 0.000 | NA | NA | 0.000 |

Rank/correlation and amplitude/variance ratios average only their defined values; a small stratum or constant predictor can make them undefined while all point predictions remain available. Exact dataset/landscape denominators are retained in `descriptive_metric_denominators.csv`.
At native d3, Spearman is defined for 15/18 datasets and amplitude for 16/18 in each fitted model. The zero predictor has undefined rank at every distance. RMSE, centered shape RMSE, edge RMSE and proper interval score in the four fitted models each retain all 18/18 datasets at d0–d3.

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
| interaction | -1 | B | 0.9967 | 1.0007 | 1.0020 | 0.8901 |
| interaction | 0 | B | 1.0275 | 1.0524 | 0.9809 | 1.3480 |
| interaction | 1 | B | 0.9915 | 0.9996 | 1.0129 | 1.1782 |
| interaction | 2 | B | 0.9973 | 0.9999 | 1.0002 | 0.8168 |
| interaction | 3 | B | 0.9964 | 1.0001 | 1.0000 | 0.6993 |
| rough | -1 | B | 0.9757 | 0.9968 | 0.9950 | 0.4955 |
| rough | 0 | B | 0.9679 | 0.9818 | 0.9682 | 0.3221 |
| rough | 1 | B | 0.9773 | 1.0012 | 0.9947 | 0.4368 |
| rough | 2 | B | 0.9779 | 1.0001 | 1.0002 | 0.5516 |
| rough | 3 | B | 0.9674 | 1.0002 | 1.0003 | 0.7221 |
| smooth | -1 | B | 0.9689 | 1.0053 | 1.0190 | 0.7276 |
| smooth | 0 | B | 1.0388 | 1.0270 | 1.0570 | 1.4789 |
| smooth | 1 | B | 0.9639 | 1.0040 | 1.0422 | 0.8458 |
| smooth | 2 | B | 0.9652 | 1.0009 | 1.0039 | 0.6662 |
| smooth | 3 | B | 0.9647 | 0.9998 | 1.0014 | 0.5709 |

This is a development comparison, not an independent test or SBC study. Reused V1 results retain their frozen source-cache SHA and unchanged counts/reference; zero elapsed time denotes reuse, not a newly timed fit. The current default is a separate actual API control, not an endpoint-ablation surrogate. No prior weight is adjusted from a sealed test. A single winner, if any, must be frozen after development and evaluated once on new calibration/test landscapes. Calibration hyperparameter heldout scores used to select gamma/kappa are not by themselves nested OOF loss vectors for joint power certification.

Generated: 2026-11-04 01:33:01 EST
