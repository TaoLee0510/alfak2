# alfak2 accuracy development benchmark

These are development results on 6 landscapes and 24 paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.

The finite target domain contains 81 4-chromosome karyotypes. Every method is scored on the same frozen target nodes and observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.

Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.

Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. This small development sample is insufficient for a general claim of improvement at every distance.

## Completion

| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |
|---|---:|---:|---:|---:|
| baseline | 24 | 24 | 4 | 0.12 |
| endpoint_evidence_bma4 | 24 | 24 | 24 | 0.35 |
| endpoint_evidence_bma5_sensitivity | 24 | 24 | 24 | 0.44 |
| joint_endpoint_s0p05 | 24 | 24 | 24 | 0.09 |
| joint_endpoint_s0p2 | 24 | 24 | 24 | 0.09 |
| joint_endpoint_s0p5 | 24 | 24 | 24 | 0.09 |
| joint_endpoint_s1 | 24 | 24 | 24 | 0.09 |
| joint_endpoint_s2 | 24 | 24 | 24 | 0.09 |

## Fitness and interval results

| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |
|---|---|---:|---:|---:|---:|---:|
| baseline | d0 | 0.1379 | 0.1077 | 0.709 | 2.1510 | 1.000 |
| baseline | d1 | 0.1911 | 0.1551 | 0.998 | 1.2671 | 1.000 |
| baseline | d2 | 0.2521 | 0.2159 | 1.000 | 1.4803 | 1.000 |
| baseline | d3 | 0.2973 | 0.2605 | 1.000 | 1.6283 | 1.000 |
| baseline | all_node_micro | 0.2308 | 0.1839 | 0.950 | 1.5721 | 1.000 |
| baseline | all_distance_macro | 0.2196 | 0.1848 | 0.927 | 1.6317 | 1.000 |
| endpoint_evidence_bma4 | d0 | 0.1372 | 0.1090 | 0.988 | 0.7878 | 1.000 |
| endpoint_evidence_bma4 | d1 | 0.2000 | 0.1666 | 0.990 | 1.0802 | 1.000 |
| endpoint_evidence_bma4 | d2 | 0.2495 | 0.2151 | 0.979 | 1.2769 | 1.000 |
| endpoint_evidence_bma4 | d3 | 0.2901 | 0.2484 | 0.943 | 1.6569 | 1.000 |
| endpoint_evidence_bma4 | all_node_micro | 0.2325 | 0.1876 | 0.978 | 1.1772 | 1.000 |
| endpoint_evidence_bma4 | all_distance_macro | 0.2192 | 0.1847 | 0.975 | 1.2005 | 1.000 |
| endpoint_evidence_bma5_sensitivity | d0 | 0.1729 | 0.1411 | 0.988 | 1.0159 | 1.000 |
| endpoint_evidence_bma5_sensitivity | d1 | 0.2319 | 0.1954 | 0.995 | 1.4601 | 1.000 |
| endpoint_evidence_bma5_sensitivity | d2 | 0.2764 | 0.2377 | 0.993 | 1.7053 | 1.000 |
| endpoint_evidence_bma5_sensitivity | d3 | 0.2989 | 0.2581 | 0.962 | 2.0314 | 1.000 |
| endpoint_evidence_bma5_sensitivity | all_node_micro | 0.2617 | 0.2138 | 0.988 | 1.5396 | 1.000 |
| endpoint_evidence_bma5_sensitivity | all_distance_macro | 0.2450 | 0.2081 | 0.985 | 1.5532 | 1.000 |
| joint_endpoint_s0p05 | d0 | 0.1376 | 0.1075 | 0.376 | 3.4528 | 1.000 |
| joint_endpoint_s0p05 | d1 | 0.1909 | 0.1549 | 0.206 | 4.8434 | 1.000 |
| joint_endpoint_s0p05 | d2 | 0.2519 | 0.2157 | 0.115 | 6.9054 | 1.000 |
| joint_endpoint_s0p05 | d3 | 0.2971 | 0.2603 | 0.104 | 8.5171 | 1.000 |
| joint_endpoint_s0p05 | all_node_micro | 0.2306 | 0.1837 | 0.182 | 5.8698 | 1.000 |
| joint_endpoint_s0p05 | all_distance_macro | 0.2194 | 0.1846 | 0.200 | 5.9297 | 1.000 |
| joint_endpoint_s0p2 | d0 | 0.1345 | 0.1055 | 0.714 | 1.8237 | 1.000 |
| joint_endpoint_s0p2 | d1 | 0.1889 | 0.1531 | 0.580 | 2.1226 | 1.000 |
| joint_endpoint_s0p2 | d2 | 0.2494 | 0.2134 | 0.466 | 2.9605 | 1.000 |
| joint_endpoint_s0p2 | d3 | 0.2957 | 0.2577 | 0.483 | 4.1847 | 1.000 |
| joint_endpoint_s0p2 | all_node_micro | 0.2282 | 0.1815 | 0.540 | 2.6687 | 1.000 |
| joint_endpoint_s0p2 | all_distance_macro | 0.2171 | 0.1824 | 0.561 | 2.7729 | 1.000 |
| joint_endpoint_s0p5 | d0 | 0.1359 | 0.1085 | 0.958 | 0.7654 | 1.000 |
| joint_endpoint_s0p5 | d1 | 0.1921 | 0.1599 | 0.974 | 0.8664 | 1.000 |
| joint_endpoint_s0p5 | d2 | 0.2471 | 0.2121 | 0.953 | 1.0567 | 1.000 |
| joint_endpoint_s0p5 | d3 | 0.2947 | 0.2529 | 0.878 | 1.5109 | 1.000 |
| joint_endpoint_s0p5 | all_node_micro | 0.2278 | 0.1836 | 0.951 | 1.0089 | 1.000 |
| joint_endpoint_s0p5 | all_distance_macro | 0.2174 | 0.1833 | 0.941 | 1.0498 | 1.000 |
| joint_endpoint_s1 | d0 | 0.2002 | 0.1660 | 0.994 | 1.0813 | 1.000 |
| joint_endpoint_s1 | d1 | 0.2488 | 0.2135 | 0.998 | 1.5339 | 1.000 |
| joint_endpoint_s1 | d2 | 0.2924 | 0.2489 | 0.998 | 1.7658 | 1.000 |
| joint_endpoint_s1 | d3 | 0.3308 | 0.2849 | 0.997 | 1.9424 | 1.000 |
| joint_endpoint_s1 | all_node_micro | 0.2766 | 0.2293 | 0.997 | 1.5934 | 1.000 |
| joint_endpoint_s1 | all_distance_macro | 0.2681 | 0.2283 | 0.997 | 1.5809 | 1.000 |
| joint_endpoint_s2 | d0 | 0.3644 | 0.2997 | 0.991 | 1.7151 | 1.000 |
| joint_endpoint_s2 | d1 | 0.3762 | 0.3229 | 1.000 | 2.7882 | 1.000 |
| joint_endpoint_s2 | d2 | 0.4118 | 0.3481 | 1.000 | 3.2442 | 1.000 |
| joint_endpoint_s2 | d3 | 0.4142 | 0.3557 | 1.000 | 3.6028 | 1.000 |
| joint_endpoint_s2 | all_node_micro | 0.4011 | 0.3363 | 0.998 | 2.8784 | 1.000 |
| joint_endpoint_s2 | all_distance_macro | 0.3917 | 0.3316 | 0.998 | 2.8376 | 1.000 |

## Independent held-out count predictions

Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.

| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |
|---|---:|---:|---:|---:|
| baseline | 0.7860 | 1.2226 | 0.00588 | 0.00623 |
| endpoint_evidence_bma4 | NA | NA | NA | NA |
| endpoint_evidence_bma5_sensitivity | NA | NA | NA | NA |
| joint_endpoint_s0p05 | 0.7772 | 1.2135 | 0.00573 | 0.00597 |
| joint_endpoint_s0p2 | 0.7765 | 1.2130 | 0.00558 | 0.00579 |
| joint_endpoint_s0p5 | 0.7736 | 1.2115 | 0.00499 | 0.00511 |
| joint_endpoint_s1 | 0.7702 | 1.2106 | 0.00424 | 0.00451 |
| joint_endpoint_s2 | 0.7696 | 1.2124 | 0.00400 | 0.00464 |

## Fitness RMSE by truth family

| Model | Family | d0 | d1 | d2 | d3 | Equal-distance macro |
|---|---|---:|---:|---:|---:|---:|
| baseline | smooth | 0.1040 | 0.1762 | 0.2484 | 0.2940 | 0.2057 |
| baseline | rough | 0.2587 | 0.2581 | 0.2550 | 0.2321 | 0.2510 |
| baseline | interaction | 0.0510 | 0.1390 | 0.2528 | 0.3657 | 0.2021 |
| endpoint_evidence_bma4 | smooth | 0.1520 | 0.2274 | 0.3004 | 0.3377 | 0.2544 |
| endpoint_evidence_bma4 | rough | 0.1882 | 0.2260 | 0.1910 | 0.1627 | 0.1920 |
| endpoint_evidence_bma4 | interaction | 0.0715 | 0.1465 | 0.2570 | 0.3698 | 0.2112 |
| endpoint_evidence_bma5_sensitivity | smooth | 0.1922 | 0.2632 | 0.3336 | 0.3580 | 0.2867 |
| endpoint_evidence_bma5_sensitivity | rough | 0.2384 | 0.2793 | 0.2346 | 0.1668 | 0.2298 |
| endpoint_evidence_bma5_sensitivity | interaction | 0.0881 | 0.1532 | 0.2611 | 0.3719 | 0.2186 |
| joint_endpoint_s0p05 | smooth | 0.1043 | 0.1766 | 0.2488 | 0.2944 | 0.2060 |
| joint_endpoint_s0p05 | rough | 0.2576 | 0.2572 | 0.2539 | 0.2312 | 0.2500 |
| joint_endpoint_s0p05 | interaction | 0.0510 | 0.1390 | 0.2528 | 0.3657 | 0.2021 |
| joint_endpoint_s0p2 | smooth | 0.1085 | 0.1828 | 0.2557 | 0.3008 | 0.2119 |
| joint_endpoint_s0p2 | rough | 0.2429 | 0.2448 | 0.2397 | 0.2202 | 0.2369 |
| joint_endpoint_s0p2 | interaction | 0.0521 | 0.1392 | 0.2528 | 0.3662 | 0.2026 |
| joint_endpoint_s0p5 | smooth | 0.1403 | 0.2165 | 0.2903 | 0.3304 | 0.2444 |
| joint_endpoint_s0p5 | rough | 0.1924 | 0.2115 | 0.1926 | 0.1820 | 0.1946 |
| joint_endpoint_s0p5 | interaction | 0.0750 | 0.1483 | 0.2583 | 0.3717 | 0.2133 |
| joint_endpoint_s1 | smooth | 0.2400 | 0.3019 | 0.3719 | 0.3935 | 0.3268 |
| joint_endpoint_s1 | rough | 0.2024 | 0.2497 | 0.2137 | 0.2020 | 0.2170 |
| joint_endpoint_s1 | interaction | 0.1581 | 0.1949 | 0.2916 | 0.3969 | 0.2604 |
| joint_endpoint_s2 | smooth | 0.4297 | 0.4484 | 0.5097 | 0.4945 | 0.4706 |
| joint_endpoint_s2 | rough | 0.3511 | 0.3833 | 0.3477 | 0.2847 | 0.3417 |
| joint_endpoint_s2 | interaction | 0.3124 | 0.2969 | 0.3782 | 0.4634 | 0.3627 |

## Paired candidate minus baseline differences

Negative RMSE and interval-score differences favor the candidate. Coverage differences have no universal favorable direction; assess distance from the nominal 0.95 and interval score jointly.

| Candidate | Distance | Metric | Mean difference | 95% landscape bootstrap interval | Landscapes |
|---|---|---|---:|---|---:|
| endpoint_evidence_bma4 | all_distance_macro | coverage_95 | 0.0483 | [-0.0168, 0.1205] | 6 |
| endpoint_evidence_bma5_sensitivity | all_distance_macro | coverage_95 | 0.0578 | [0.0003, 0.1242] | 6 |
| joint_endpoint_s0p05 | all_distance_macro | coverage_95 | -0.7266 | [-0.7504, -0.7034] | 6 |
| joint_endpoint_s0p2 | all_distance_macro | coverage_95 | -0.3661 | [-0.4143, -0.3059] | 6 |
| joint_endpoint_s0p5 | all_distance_macro | coverage_95 | 0.0139 | [-0.0405, 0.0767] | 6 |
| joint_endpoint_s1 | all_distance_macro | coverage_95 | 0.0702 | [0.0249, 0.1343] | 6 |
| joint_endpoint_s2 | all_distance_macro | coverage_95 | 0.0710 | [0.0236, 0.1298] | 6 |
| endpoint_evidence_bma4 | d0 | coverage_95 | 0.2782 | [0.0980, 0.5044] | 6 |
| endpoint_evidence_bma5_sensitivity | d0 | coverage_95 | 0.2782 | [0.0979, 0.5044] | 6 |
| joint_endpoint_s0p05 | d0 | coverage_95 | -0.3338 | [-0.4550, -0.2116] | 6 |
| joint_endpoint_s0p2 | d0 | coverage_95 | 0.0044 | [-0.0252, 0.0368] | 6 |
| joint_endpoint_s0p5 | d0 | coverage_95 | 0.2486 | [0.1064, 0.4341] | 6 |
| joint_endpoint_s1 | d0 | coverage_95 | 0.2848 | [0.1008, 0.5310] | 6 |
| joint_endpoint_s2 | d0 | coverage_95 | 0.2816 | [0.0944, 0.5280] | 6 |
| endpoint_evidence_bma4 | d1 | coverage_95 | -0.0074 | [-0.0139, -0.0003] | 6 |
| endpoint_evidence_bma5_sensitivity | d1 | coverage_95 | -0.0027 | [-0.0089, 0.0044] | 6 |
| joint_endpoint_s0p05 | d1 | coverage_95 | -0.7914 | [-0.8514, -0.7315] | 6 |
| joint_endpoint_s0p2 | d1 | coverage_95 | -0.4175 | [-0.5120, -0.3297] | 6 |
| joint_endpoint_s0p5 | d1 | coverage_95 | -0.0237 | [-0.0382, -0.0095] | 6 |
| joint_endpoint_s1 | d1 | coverage_95 | -0.0000 | [-0.0034, 0.0033] | 6 |
| joint_endpoint_s2 | d1 | coverage_95 | 0.0022 | [0.0000, 0.0067] | 6 |
| endpoint_evidence_bma4 | d2 | coverage_95 | -0.0211 | [-0.0342, -0.0085] | 6 |
| endpoint_evidence_bma5_sensitivity | d2 | coverage_95 | -0.0066 | [-0.0123, -0.0013] | 6 |
| joint_endpoint_s0p05 | d2 | coverage_95 | -0.8852 | [-0.9205, -0.8502] | 6 |
| joint_endpoint_s0p2 | d2 | coverage_95 | -0.5340 | [-0.6174, -0.4402] | 6 |
| joint_endpoint_s0p5 | d2 | coverage_95 | -0.0470 | [-0.0736, -0.0178] | 6 |
| joint_endpoint_s1 | d2 | coverage_95 | -0.0015 | [-0.0046, 0.0000] | 6 |
| joint_endpoint_s2 | d2 | coverage_95 | 0.0000 | [0.0000, 0.0000] | 6 |
| endpoint_evidence_bma4 | d3 | coverage_95 | -0.0566 | [-0.1133, -0.0060] | 6 |
| endpoint_evidence_bma5_sensitivity | d3 | coverage_95 | -0.0378 | [-0.0887, 0.0000] | 6 |
| joint_endpoint_s0p05 | d3 | coverage_95 | -0.8957 | [-0.9530, -0.8382] | 6 |
| joint_endpoint_s0p2 | d3 | coverage_95 | -0.5171 | [-0.6639, -0.3495] | 6 |
| joint_endpoint_s0p5 | d3 | coverage_95 | -0.1221 | [-0.1849, -0.0604] | 6 |
| joint_endpoint_s1 | d3 | coverage_95 | -0.0026 | [-0.0078, 0.0000] | 6 |
| joint_endpoint_s2 | d3 | coverage_95 | 0.0000 | [0.0000, 0.0000] | 6 |
| endpoint_evidence_bma4 | all_distance_macro | interval_score_95 | -0.4312 | [-1.0342, -0.0145] | 6 |
| endpoint_evidence_bma5_sensitivity | all_distance_macro | interval_score_95 | -0.0785 | [-0.4921, 0.2528] | 6 |
| joint_endpoint_s0p05 | all_distance_macro | interval_score_95 | 4.2980 | [3.4059, 5.4859] | 6 |
| joint_endpoint_s0p2 | all_distance_macro | interval_score_95 | 1.1412 | [0.6351, 1.5931] | 6 |
| joint_endpoint_s0p5 | all_distance_macro | interval_score_95 | -0.5819 | [-1.2899, -0.0775] | 6 |
| joint_endpoint_s1 | all_distance_macro | interval_score_95 | -0.0508 | [-0.7474, 0.3545] | 6 |
| joint_endpoint_s2 | all_distance_macro | interval_score_95 | 1.2059 | [0.5138, 1.6184] | 6 |
| endpoint_evidence_bma4 | d0 | interval_score_95 | -1.3632 | [-4.0024, 0.2198] | 6 |
| endpoint_evidence_bma5_sensitivity | d0 | interval_score_95 | -1.1351 | [-3.6932, 0.4132] | 6 |
| joint_endpoint_s0p05 | d0 | interval_score_95 | 1.3018 | [0.6593, 1.9514] | 6 |
| joint_endpoint_s0p2 | d0 | interval_score_95 | -0.3273 | [-0.8414, -0.0328] | 6 |
| joint_endpoint_s0p5 | d0 | interval_score_95 | -1.3856 | [-3.8285, 0.1075] | 6 |
| joint_endpoint_s1 | d0 | interval_score_95 | -1.0697 | [-3.8525, 0.5812] | 6 |
| joint_endpoint_s2 | d0 | interval_score_95 | -0.4360 | [-3.2321, 1.2023] | 6 |
| endpoint_evidence_bma4 | d1 | interval_score_95 | -0.1870 | [-0.2878, -0.0459] | 6 |
| endpoint_evidence_bma5_sensitivity | d1 | interval_score_95 | 0.1930 | [-0.0446, 0.5323] | 6 |
| joint_endpoint_s0p05 | d1 | interval_score_95 | 3.5763 | [2.1561, 5.5019] | 6 |
| joint_endpoint_s0p2 | d1 | interval_score_95 | 0.8555 | [0.0389, 2.0484] | 6 |
| joint_endpoint_s0p5 | d1 | interval_score_95 | -0.4008 | [-0.4373, -0.3623] | 6 |
| joint_endpoint_s1 | d1 | interval_score_95 | 0.2667 | [0.2572, 0.2766] | 6 |
| joint_endpoint_s2 | d1 | interval_score_95 | 1.5210 | [1.5032, 1.5402] | 6 |
| endpoint_evidence_bma4 | d2 | interval_score_95 | -0.2034 | [-0.3185, -0.0575] | 6 |
| endpoint_evidence_bma5_sensitivity | d2 | interval_score_95 | 0.2249 | [-0.0547, 0.6246] | 6 |
| joint_endpoint_s0p05 | d2 | interval_score_95 | 5.4251 | [4.1653, 6.6288] | 6 |
| joint_endpoint_s0p2 | d2 | interval_score_95 | 1.4802 | [0.7317, 2.0891] | 6 |
| joint_endpoint_s0p5 | d2 | interval_score_95 | -0.4236 | [-0.5012, -0.3300] | 6 |
| joint_endpoint_s1 | d2 | interval_score_95 | 0.2854 | [0.2754, 0.2960] | 6 |
| joint_endpoint_s2 | d2 | interval_score_95 | 1.7639 | [1.7489, 1.7817] | 6 |
| endpoint_evidence_bma4 | d3 | interval_score_95 | 0.0285 | [-0.3212, 0.4368] | 6 |
| endpoint_evidence_bma5_sensitivity | d3 | interval_score_95 | 0.4031 | [0.0869, 0.7899] | 6 |
| joint_endpoint_s0p05 | d3 | interval_score_95 | 6.8887 | [4.5236, 8.8019] | 6 |
| joint_endpoint_s0p2 | d3 | interval_score_95 | 2.5564 | [0.9918, 3.8971] | 6 |
| joint_endpoint_s0p5 | d3 | interval_score_95 | -0.1174 | [-0.5598, 0.4586] | 6 |
| joint_endpoint_s1 | d3 | interval_score_95 | 0.3141 | [0.3030, 0.3239] | 6 |
| joint_endpoint_s2 | d3 | interval_score_95 | 1.9745 | [1.9513, 1.9907] | 6 |
| endpoint_evidence_bma4 | all_distance_macro | relative_rmse | -0.0004 | [-0.0517, 0.0379] | 6 |
| endpoint_evidence_bma5_sensitivity | all_distance_macro | relative_rmse | 0.0254 | [-0.0106, 0.0628] | 6 |
| joint_endpoint_s0p05 | all_distance_macro | relative_rmse | -0.0002 | [-0.0010, 0.0003] | 6 |
| joint_endpoint_s0p2 | all_distance_macro | relative_rmse | -0.0025 | [-0.0133, 0.0044] | 6 |
| joint_endpoint_s0p5 | all_distance_macro | relative_rmse | -0.0021 | [-0.0509, 0.0295] | 6 |
| joint_endpoint_s1 | all_distance_macro | relative_rmse | 0.0485 | [-0.0208, 0.1013] | 6 |
| joint_endpoint_s2 | all_distance_macro | relative_rmse | 0.1721 | [0.1122, 0.2311] | 6 |
| endpoint_evidence_bma4 | d0 | relative_rmse | -0.0007 | [-0.0605, 0.0395] | 6 |
| endpoint_evidence_bma5_sensitivity | d0 | relative_rmse | 0.0350 | [-0.0035, 0.0746] | 6 |
| joint_endpoint_s0p05 | d0 | relative_rmse | -0.0003 | [-0.0011, 0.0002] | 6 |
| joint_endpoint_s0p2 | d0 | relative_rmse | -0.0034 | [-0.0150, 0.0033] | 6 |
| joint_endpoint_s0p5 | d0 | relative_rmse | -0.0020 | [-0.0575, 0.0320] | 6 |
| joint_endpoint_s1 | d0 | relative_rmse | 0.0623 | [-0.0243, 0.1264] | 6 |
| joint_endpoint_s2 | d0 | relative_rmse | 0.2265 | [0.1452, 0.2994] | 6 |
| endpoint_evidence_bma4 | d1 | relative_rmse | 0.0089 | [-0.0264, 0.0387] | 6 |
| endpoint_evidence_bma5_sensitivity | d1 | relative_rmse | 0.0408 | [0.0147, 0.0721] | 6 |
| joint_endpoint_s0p05 | d1 | relative_rmse | -0.0002 | [-0.0009, 0.0003] | 6 |
| joint_endpoint_s0p2 | d1 | relative_rmse | -0.0022 | [-0.0129, 0.0047] | 6 |
| joint_endpoint_s0p5 | d1 | relative_rmse | 0.0010 | [-0.0413, 0.0304] | 6 |
| joint_endpoint_s1 | d1 | relative_rmse | 0.0577 | [0.0076, 0.1033] | 6 |
| joint_endpoint_s2 | d1 | relative_rmse | 0.1851 | [0.1398, 0.2431] | 6 |
| endpoint_evidence_bma4 | d2 | relative_rmse | -0.0026 | [-0.0569, 0.0383] | 6 |
| endpoint_evidence_bma5_sensitivity | d2 | relative_rmse | 0.0244 | [-0.0108, 0.0663] | 6 |
| joint_endpoint_s0p05 | d2 | relative_rmse | -0.0002 | [-0.0011, 0.0003] | 6 |
| joint_endpoint_s0p2 | d2 | relative_rmse | -0.0027 | [-0.0130, 0.0049] | 6 |
| joint_endpoint_s0p5 | d2 | relative_rmse | -0.0050 | [-0.0565, 0.0301] | 6 |
| joint_endpoint_s1 | d2 | relative_rmse | 0.0404 | [-0.0285, 0.0971] | 6 |
| joint_endpoint_s2 | d2 | relative_rmse | 0.1598 | [0.1032, 0.2217] | 6 |
| endpoint_evidence_bma4 | d3 | relative_rmse | -0.0072 | [-0.0709, 0.0332] | 6 |
| endpoint_evidence_bma5_sensitivity | d3 | relative_rmse | 0.0016 | [-0.0572, 0.0510] | 6 |
| joint_endpoint_s0p05 | d3 | relative_rmse | -0.0001 | [-0.0008, 0.0003] | 6 |
| joint_endpoint_s0p2 | d3 | relative_rmse | -0.0015 | [-0.0111, 0.0049] | 6 |
| joint_endpoint_s0p5 | d3 | relative_rmse | -0.0026 | [-0.0493, 0.0285] | 6 |
| joint_endpoint_s1 | d3 | relative_rmse | 0.0336 | [-0.0402, 0.0848] | 6 |
| joint_endpoint_s2 | d3 | relative_rmse | 0.1169 | [0.0340, 0.1876] | 6 |

## Limits and next gate

The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.

Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.

Generated: 2026-11-03 23:36:25 EST
