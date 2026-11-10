# alfak2 accuracy development benchmark

These are development results on 6 landscapes and 24 paired datasets. Calibration and test landscape seeds are reserved and unopened. No truth-based affine rescaling, tuning, or interval calibration is applied.

The finite target domain contains 81 4-chromosome karyotypes. Every method is scored on the same frozen target nodes and reference observed-set Manhattan distances. Each node uses one common t0-count-weighted observed-reference contrast. Exact-distance strata are weighted equally in the macro summary, preventing distant node counts from dominating.

Uncertainty intervals use the full covariance of that same contrast when available. Joint division-process intervals use 8,000 seeded draws from the correlated latent Gaussian Laplace approximation, transformed to birth rates and then to the fixed reference contrast. Gaussian graph intervals use the exact linear-contrast variance. The frozen baseline covariance is recovered from its unchanged Gaussian precision and accepted only when both means and diagonal variances match within 1e-7. Missing covariance produces missing interval scores. Raw fitness RMSE is retained for diagnostics; birth-rate and relative-fitness raw scales are not directly comparable.

Independent units are landscapes, not nodes or method runs. Paired 95% percentile bootstrap intervals resample landscape-level mean differences. This small development sample is insufficient for a general claim of improvement at every distance.

## Completion

| Model | Datasets | Successful fits | Trusted local covariance | Median seconds |
|---|---:|---:|---:|---:|
| baseline | 24 | 24 | 4 | 0.12 |
| joint_endpoint_b0p1_r0p2 | 24 | 24 | 24 | 0.23 |
| joint_endpoint_b0p1_r0p5 | 24 | 24 | 24 | 0.20 |
| joint_endpoint_b0p1_r1 | 24 | 24 | 24 | 0.18 |
| joint_endpoint_b0p2_r0p2 | 24 | 24 | 24 | 0.23 |
| joint_endpoint_b0p2_r0p5 | 24 | 24 | 24 | 0.20 |
| joint_endpoint_b0p2_r1 | 24 | 24 | 24 | 0.19 |
| joint_endpoint_b0p35_r0p2 | 24 | 24 | 24 | 0.22 |
| joint_endpoint_b0p35_r0p5 | 24 | 24 | 24 | 0.20 |
| joint_endpoint_b0p35_r1 | 24 | 24 | 24 | 0.18 |
| joint_endpoint_b0p5_r0p2 | 24 | 24 | 24 | 0.23 |
| joint_endpoint_b0p5_r0p5 | 24 | 24 | 24 | 0.20 |
| joint_endpoint_b0p5_r1 | 24 | 24 | 24 | 0.18 |

## Fitness and interval results

| Model | Distance | Relative RMSE | Relative MAE | 95% coverage | 95% interval score | Prediction coverage |
|---|---|---:|---:|---:|---:|---:|
| baseline | d0 | 0.1379 | 0.1077 | 0.709 | 2.1510 | 1.000 |
| baseline | d1 | 0.1911 | 0.1551 | 0.998 | 1.2671 | 1.000 |
| baseline | d2 | 0.2521 | 0.2159 | 1.000 | 1.4803 | 1.000 |
| baseline | d3 | 0.2973 | 0.2605 | 1.000 | 1.6283 | 1.000 |
| baseline | all_node_micro | 0.2308 | 0.1839 | 0.950 | 1.5721 | 1.000 |
| baseline | all_distance_macro | 0.2196 | 0.1848 | 0.927 | 1.6317 | 1.000 |
| joint_endpoint_b0p1_r0p2 | d0 | 0.1350 | 0.1059 | 0.690 | 2.0747 | 1.000 |
| joint_endpoint_b0p1_r0p2 | d1 | 0.1891 | 0.1532 | 0.505 | 2.7927 | 1.000 |
| joint_endpoint_b0p1_r0p2 | d2 | 0.2501 | 0.2139 | 0.321 | 4.1936 | 1.000 |
| joint_endpoint_b0p1_r0p2 | d3 | 0.2958 | 0.2580 | 0.246 | 5.6181 | 1.000 |
| joint_endpoint_b0p1_r0p2 | all_node_micro | 0.2287 | 0.1818 | 0.433 | 3.5719 | 1.000 |
| joint_endpoint_b0p1_r0p2 | all_distance_macro | 0.2175 | 0.1828 | 0.441 | 3.6698 | 1.000 |
| joint_endpoint_b0p1_r0p5 | d0 | 0.1340 | 0.1064 | 0.903 | 0.8053 | 1.000 |
| joint_endpoint_b0p1_r0p5 | d1 | 0.1881 | 0.1554 | 0.840 | 1.0229 | 1.000 |
| joint_endpoint_b0p1_r0p5 | d2 | 0.2478 | 0.2124 | 0.737 | 1.7230 | 1.000 |
| joint_endpoint_b0p1_r0p5 | d3 | 0.2933 | 0.2531 | 0.645 | 3.0037 | 1.000 |
| joint_endpoint_b0p1_r0p5 | all_node_micro | 0.2262 | 0.1817 | 0.787 | 1.4737 | 1.000 |
| joint_endpoint_b0p1_r0p5 | all_distance_macro | 0.2158 | 0.1818 | 0.781 | 1.6387 | 1.000 |
| joint_endpoint_b0p1_r1 | d0 | 0.1701 | 0.1420 | 0.988 | 0.8763 | 1.000 |
| joint_endpoint_b0p1_r1 | d1 | 0.2153 | 0.1836 | 0.983 | 1.0782 | 1.000 |
| joint_endpoint_b0p1_r1 | d2 | 0.2696 | 0.2322 | 0.944 | 1.2455 | 1.000 |
| joint_endpoint_b0p1_r1 | d3 | 0.3168 | 0.2766 | 0.835 | 1.7880 | 1.000 |
| joint_endpoint_b0p1_r1 | all_node_micro | 0.2493 | 0.2065 | 0.956 | 1.1777 | 1.000 |
| joint_endpoint_b0p1_r1 | all_distance_macro | 0.2430 | 0.2086 | 0.938 | 1.2470 | 1.000 |
| joint_endpoint_b0p2_r0p2 | d0 | 0.1345 | 0.1055 | 0.714 | 1.8237 | 1.000 |
| joint_endpoint_b0p2_r0p2 | d1 | 0.1889 | 0.1531 | 0.580 | 2.1226 | 1.000 |
| joint_endpoint_b0p2_r0p2 | d2 | 0.2494 | 0.2134 | 0.466 | 2.9605 | 1.000 |
| joint_endpoint_b0p2_r0p2 | d3 | 0.2957 | 0.2577 | 0.483 | 4.1847 | 1.000 |
| joint_endpoint_b0p2_r0p2 | all_node_micro | 0.2282 | 0.1815 | 0.540 | 2.6687 | 1.000 |
| joint_endpoint_b0p2_r0p2 | all_distance_macro | 0.2171 | 0.1824 | 0.561 | 2.7729 | 1.000 |
| joint_endpoint_b0p2_r0p5 | d0 | 0.1340 | 0.1065 | 0.912 | 0.7801 | 1.000 |
| joint_endpoint_b0p2_r0p5 | d1 | 0.1883 | 0.1557 | 0.872 | 0.9446 | 1.000 |
| joint_endpoint_b0p2_r0p5 | d2 | 0.2474 | 0.2123 | 0.796 | 1.4963 | 1.000 |
| joint_endpoint_b0p2_r0p5 | d3 | 0.2933 | 0.2530 | 0.698 | 2.6055 | 1.000 |
| joint_endpoint_b0p2_r0p5 | all_node_micro | 0.2261 | 0.1817 | 0.829 | 1.3158 | 1.000 |
| joint_endpoint_b0p2_r0p5 | all_distance_macro | 0.2158 | 0.1819 | 0.819 | 1.4566 | 1.000 |
| joint_endpoint_b0p2_r1 | d0 | 0.1709 | 0.1426 | 0.988 | 0.8832 | 1.000 |
| joint_endpoint_b0p2_r1 | d1 | 0.2159 | 0.1842 | 0.983 | 1.0921 | 1.000 |
| joint_endpoint_b0p2_r1 | d2 | 0.2696 | 0.2322 | 0.950 | 1.2410 | 1.000 |
| joint_endpoint_b0p2_r1 | d3 | 0.3165 | 0.2764 | 0.865 | 1.7067 | 1.000 |
| joint_endpoint_b0p2_r1 | all_node_micro | 0.2496 | 0.2068 | 0.960 | 1.1763 | 1.000 |
| joint_endpoint_b0p2_r1 | all_distance_macro | 0.2432 | 0.2089 | 0.946 | 1.2307 | 1.000 |
| joint_endpoint_b0p35_r0p2 | d0 | 0.1334 | 0.1048 | 0.796 | 1.3736 | 1.000 |
| joint_endpoint_b0p35_r0p2 | d1 | 0.1888 | 0.1535 | 0.781 | 1.2708 | 1.000 |
| joint_endpoint_b0p35_r0p2 | d2 | 0.2481 | 0.2125 | 0.751 | 1.5915 | 1.000 |
| joint_endpoint_b0p35_r0p2 | d3 | 0.2958 | 0.2569 | 0.698 | 2.5283 | 1.000 |
| joint_endpoint_b0p35_r0p2 | all_node_micro | 0.2275 | 0.1811 | 0.760 | 1.5989 | 1.000 |
| joint_endpoint_b0p35_r0p2 | all_distance_macro | 0.2165 | 0.1819 | 0.757 | 1.6910 | 1.000 |
| joint_endpoint_b0p35_r0p5 | d0 | 0.1345 | 0.1070 | 0.937 | 0.7589 | 1.000 |
| joint_endpoint_b0p35_r0p5 | d1 | 0.1894 | 0.1569 | 0.932 | 0.8471 | 1.000 |
| joint_endpoint_b0p35_r0p5 | d2 | 0.2468 | 0.2120 | 0.877 | 1.1735 | 1.000 |
| joint_endpoint_b0p35_r0p5 | d3 | 0.2936 | 0.2529 | 0.782 | 1.9317 | 1.000 |
| joint_endpoint_b0p35_r0p5 | all_node_micro | 0.2263 | 0.1822 | 0.895 | 1.0907 | 1.000 |
| joint_endpoint_b0p35_r0p5 | all_distance_macro | 0.2161 | 0.1822 | 0.882 | 1.1778 | 1.000 |
| joint_endpoint_b0p35_r1 | d0 | 0.1730 | 0.1443 | 0.991 | 0.9023 | 1.000 |
| joint_endpoint_b0p35_r1 | d1 | 0.2176 | 0.1860 | 0.989 | 1.1312 | 1.000 |
| joint_endpoint_b0p35_r1 | d2 | 0.2698 | 0.2325 | 0.969 | 1.2529 | 1.000 |
| joint_endpoint_b0p35_r1 | d3 | 0.3162 | 0.2762 | 0.905 | 1.5952 | 1.000 |
| joint_endpoint_b0p35_r1 | all_node_micro | 0.2507 | 0.2079 | 0.973 | 1.1878 | 1.000 |
| joint_endpoint_b0p35_r1 | all_distance_macro | 0.2442 | 0.2097 | 0.964 | 1.2204 | 1.000 |
| joint_endpoint_b0p5_r0p2 | d0 | 0.1330 | 0.1045 | 0.875 | 1.0885 | 1.000 |
| joint_endpoint_b0p5_r0p2 | d1 | 0.1899 | 0.1550 | 0.910 | 0.9412 | 1.000 |
| joint_endpoint_b0p5_r0p2 | d2 | 0.2474 | 0.2114 | 0.896 | 1.0914 | 1.000 |
| joint_endpoint_b0p5_r0p2 | d3 | 0.2967 | 0.2565 | 0.830 | 1.6936 | 1.000 |
| joint_endpoint_b0p5_r0p2 | all_node_micro | 0.2277 | 0.1811 | 0.885 | 1.1490 | 1.000 |
| joint_endpoint_b0p5_r0p2 | all_distance_macro | 0.2167 | 0.1819 | 0.878 | 1.2037 | 1.000 |
| joint_endpoint_b0p5_r0p5 | d0 | 0.1359 | 0.1085 | 0.958 | 0.7654 | 1.000 |
| joint_endpoint_b0p5_r0p5 | d1 | 0.1921 | 0.1599 | 0.974 | 0.8664 | 1.000 |
| joint_endpoint_b0p5_r0p5 | d2 | 0.2471 | 0.2121 | 0.953 | 1.0567 | 1.000 |
| joint_endpoint_b0p5_r0p5 | d3 | 0.2947 | 0.2529 | 0.878 | 1.5109 | 1.000 |
| joint_endpoint_b0p5_r0p5 | all_node_micro | 0.2278 | 0.1836 | 0.951 | 1.0089 | 1.000 |
| joint_endpoint_b0p5_r0p5 | all_distance_macro | 0.2174 | 0.1833 | 0.941 | 1.0498 | 1.000 |
| joint_endpoint_b0p5_r1 | d0 | 0.1766 | 0.1472 | 0.991 | 0.9330 | 1.000 |
| joint_endpoint_b0p5_r1 | d1 | 0.2210 | 0.1891 | 0.994 | 1.1997 | 1.000 |
| joint_endpoint_b0p5_r1 | d2 | 0.2711 | 0.2336 | 0.985 | 1.3178 | 1.000 |
| joint_endpoint_b0p5_r1 | d3 | 0.3165 | 0.2761 | 0.948 | 1.5202 | 1.000 |
| joint_endpoint_b0p5_r1 | all_node_micro | 0.2530 | 0.2101 | 0.986 | 1.2378 | 1.000 |
| joint_endpoint_b0p5_r1 | all_distance_macro | 0.2463 | 0.2115 | 0.980 | 1.2427 | 1.000 |

## Independent held-out count predictions

Held-out counts use each frozen latent population and an independent sampling seed (training sampling seed + 50,000,000); they are never used to fit the model. Scores are negative log likelihood per cell. Zero predictions are floored at 1e-12 and renormalized; positive held-out counts assigned zero probability are separately recorded in the CSV. Conditional expected cross entropy and full-target frequency RMSE are also retained.

| Model | t0 NLL/cell | t1 NLL/cell | t0 frequency RMSE | t1 frequency RMSE |
|---|---:|---:|---:|---:|
| baseline | 0.7860 | 1.2226 | 0.00588 | 0.00623 |
| joint_endpoint_b0p1_r0p2 | 0.7766 | 1.2130 | 0.00559 | 0.00581 |
| joint_endpoint_b0p1_r0p5 | 0.7740 | 1.2114 | 0.00504 | 0.00516 |
| joint_endpoint_b0p1_r1 | 0.7707 | 1.2101 | 0.00427 | 0.00452 |
| joint_endpoint_b0p2_r0p2 | 0.7765 | 1.2130 | 0.00558 | 0.00579 |
| joint_endpoint_b0p2_r0p5 | 0.7739 | 1.2114 | 0.00503 | 0.00516 |
| joint_endpoint_b0p2_r1 | 0.7707 | 1.2101 | 0.00427 | 0.00451 |
| joint_endpoint_b0p35_r0p2 | 0.7762 | 1.2130 | 0.00555 | 0.00576 |
| joint_endpoint_b0p35_r0p5 | 0.7738 | 1.2115 | 0.00502 | 0.00514 |
| joint_endpoint_b0p35_r1 | 0.7707 | 1.2102 | 0.00426 | 0.00451 |
| joint_endpoint_b0p5_r0p2 | 0.7759 | 1.2129 | 0.00551 | 0.00571 |
| joint_endpoint_b0p5_r0p5 | 0.7736 | 1.2115 | 0.00499 | 0.00511 |
| joint_endpoint_b0p5_r1 | 0.7706 | 1.2102 | 0.00426 | 0.00451 |

## Fitness RMSE by truth family

| Model | Family | d0 | d1 | d2 | d3 | Equal-distance macro |
|---|---|---:|---:|---:|---:|---:|
| baseline | smooth | 0.1040 | 0.1762 | 0.2484 | 0.2940 | 0.2057 |
| baseline | rough | 0.2587 | 0.2581 | 0.2550 | 0.2321 | 0.2510 |
| baseline | interaction | 0.0510 | 0.1390 | 0.2528 | 0.3657 | 0.2021 |
| joint_endpoint_b0p1_r0p2 | smooth | 0.1083 | 0.1823 | 0.2552 | 0.3014 | 0.2118 |
| joint_endpoint_b0p1_r0p2 | rough | 0.2448 | 0.2461 | 0.2421 | 0.2196 | 0.2382 |
| joint_endpoint_b0p1_r0p2 | interaction | 0.0519 | 0.1391 | 0.2529 | 0.3664 | 0.2026 |
| joint_endpoint_b0p1_r0p5 | smooth | 0.1344 | 0.2103 | 0.2852 | 0.3344 | 0.2411 |
| joint_endpoint_b0p1_r0p5 | rough | 0.1991 | 0.2098 | 0.2024 | 0.1740 | 0.1963 |
| joint_endpoint_b0p1_r0p5 | interaction | 0.0685 | 0.1441 | 0.2558 | 0.3717 | 0.2100 |
| joint_endpoint_b0p1_r1 | smooth | 0.2062 | 0.2720 | 0.3482 | 0.4011 | 0.3069 |
| joint_endpoint_b0p1_r1 | rough | 0.1763 | 0.2062 | 0.1901 | 0.1622 | 0.1837 |
| joint_endpoint_b0p1_r1 | interaction | 0.1279 | 0.1677 | 0.2705 | 0.3871 | 0.2383 |
| joint_endpoint_b0p2_r0p2 | smooth | 0.1085 | 0.1828 | 0.2557 | 0.3008 | 0.2119 |
| joint_endpoint_b0p2_r0p2 | rough | 0.2429 | 0.2448 | 0.2397 | 0.2202 | 0.2369 |
| joint_endpoint_b0p2_r0p2 | interaction | 0.0521 | 0.1392 | 0.2528 | 0.3662 | 0.2026 |
| joint_endpoint_b0p2_r0p5 | smooth | 0.1350 | 0.2108 | 0.2856 | 0.3336 | 0.2413 |
| joint_endpoint_b0p2_r0p5 | rough | 0.1980 | 0.2096 | 0.2007 | 0.1749 | 0.1958 |
| joint_endpoint_b0p2_r0p5 | interaction | 0.0692 | 0.1444 | 0.2559 | 0.3715 | 0.2102 |
| joint_endpoint_b0p2_r1 | smooth | 0.2070 | 0.2724 | 0.3484 | 0.4000 | 0.3070 |
| joint_endpoint_b0p2_r1 | rough | 0.1768 | 0.2071 | 0.1897 | 0.1627 | 0.1841 |
| joint_endpoint_b0p2_r1 | interaction | 0.1287 | 0.1680 | 0.2706 | 0.3869 | 0.2386 |
| joint_endpoint_b0p35_r0p2 | smooth | 0.1093 | 0.1846 | 0.2575 | 0.2997 | 0.2128 |
| joint_endpoint_b0p35_r0p2 | rough | 0.2380 | 0.2418 | 0.2336 | 0.2218 | 0.2338 |
| joint_endpoint_b0p35_r0p2 | interaction | 0.0530 | 0.1400 | 0.2532 | 0.3659 | 0.2030 |
| joint_endpoint_b0p35_r0p5 | smooth | 0.1368 | 0.2127 | 0.2871 | 0.3319 | 0.2421 |
| joint_endpoint_b0p35_r0p5 | rough | 0.1953 | 0.2098 | 0.1968 | 0.1776 | 0.1949 |
| joint_endpoint_b0p35_r0p5 | interaction | 0.0712 | 0.1456 | 0.2565 | 0.3713 | 0.2112 |
| joint_endpoint_b0p35_r1 | smooth | 0.2095 | 0.2738 | 0.3490 | 0.3973 | 0.3074 |
| joint_endpoint_b0p35_r1 | rough | 0.1784 | 0.2099 | 0.1892 | 0.1646 | 0.1855 |
| joint_endpoint_b0p35_r1 | interaction | 0.1311 | 0.1692 | 0.2713 | 0.3867 | 0.2396 |
| joint_endpoint_b0p5_r0p2 | smooth | 0.1117 | 0.1886 | 0.2614 | 0.2992 | 0.2152 |
| joint_endpoint_b0p5_r0p2 | rough | 0.2316 | 0.2389 | 0.2261 | 0.2247 | 0.2303 |
| joint_endpoint_b0p5_r0p2 | interaction | 0.0556 | 0.1422 | 0.2546 | 0.3662 | 0.2046 |
| joint_endpoint_b0p5_r0p5 | smooth | 0.1403 | 0.2165 | 0.2903 | 0.3304 | 0.2444 |
| joint_endpoint_b0p5_r0p5 | rough | 0.1924 | 0.2115 | 0.1926 | 0.1820 | 0.1946 |
| joint_endpoint_b0p5_r0p5 | interaction | 0.0750 | 0.1483 | 0.2583 | 0.3717 | 0.2133 |
| joint_endpoint_b0p5_r1 | smooth | 0.2135 | 0.2766 | 0.3507 | 0.3940 | 0.3087 |
| joint_endpoint_b0p5_r1 | rough | 0.1813 | 0.2147 | 0.1897 | 0.1685 | 0.1886 |
| joint_endpoint_b0p5_r1 | interaction | 0.1349 | 0.1716 | 0.2729 | 0.3868 | 0.2416 |

## Paired candidate minus baseline differences

Negative RMSE and interval-score differences favor the candidate. Coverage differences have no universal favorable direction; assess distance from the nominal 0.95 and interval score jointly.

| Candidate | Distance | Metric | Mean difference | 95% landscape bootstrap interval | Landscapes |
|---|---|---|---:|---|---:|
| joint_endpoint_b0p1_r0p2 | all_distance_macro | coverage_95 | -0.4863 | [-0.5249, -0.4423] | 6 |
| joint_endpoint_b0p1_r0p5 | all_distance_macro | coverage_95 | -0.1458 | [-0.2301, -0.0632] | 6 |
| joint_endpoint_b0p1_r1 | all_distance_macro | coverage_95 | 0.0108 | [-0.0707, 0.1091] | 6 |
| joint_endpoint_b0p2_r0p2 | all_distance_macro | coverage_95 | -0.3661 | [-0.4143, -0.3060] | 6 |
| joint_endpoint_b0p2_r0p5 | all_distance_macro | coverage_95 | -0.1074 | [-0.1810, -0.0338] | 6 |
| joint_endpoint_b0p2_r1 | all_distance_macro | coverage_95 | 0.0197 | [-0.0567, 0.1074] | 6 |
| joint_endpoint_b0p35_r0p2 | all_distance_macro | coverage_95 | -0.1702 | [-0.2320, -0.1135] | 6 |
| joint_endpoint_b0p35_r0p5 | all_distance_macro | coverage_95 | -0.0448 | [-0.1005, 0.0119] | 6 |
| joint_endpoint_b0p35_r1 | all_distance_macro | coverage_95 | 0.0367 | [-0.0276, 0.1125] | 6 |
| joint_endpoint_b0p5_r0p2 | all_distance_macro | coverage_95 | -0.0488 | [-0.0810, -0.0140] | 6 |
| joint_endpoint_b0p5_r0p5 | all_distance_macro | coverage_95 | 0.0139 | [-0.0396, 0.0683] | 6 |
| joint_endpoint_b0p5_r1 | all_distance_macro | coverage_95 | 0.0530 | [-0.0030, 0.1225] | 6 |
| joint_endpoint_b0p1_r0p2 | d0 | coverage_95 | -0.0192 | [-0.0775, 0.0295] | 6 |
| joint_endpoint_b0p1_r0p5 | d0 | coverage_95 | 0.1934 | [0.0835, 0.3229] | 6 |
| joint_endpoint_b0p1_r1 | d0 | coverage_95 | 0.2790 | [0.0950, 0.5104] | 6 |
| joint_endpoint_b0p2_r0p2 | d0 | coverage_95 | 0.0044 | [-0.0252, 0.0384] | 6 |
| joint_endpoint_b0p2_r0p5 | d0 | coverage_95 | 0.2026 | [0.0916, 0.3178] | 6 |
| joint_endpoint_b0p2_r1 | d0 | coverage_95 | 0.2790 | [0.0896, 0.5067] | 6 |
| joint_endpoint_b0p35_r0p2 | d0 | coverage_95 | 0.0867 | [0.0473, 0.1248] | 6 |
| joint_endpoint_b0p35_r0p5 | d0 | coverage_95 | 0.2279 | [0.0965, 0.3827] | 6 |
| joint_endpoint_b0p35_r1 | d0 | coverage_95 | 0.2820 | [0.0940, 0.5121] | 6 |
| joint_endpoint_b0p5_r0p2 | d0 | coverage_95 | 0.1660 | [0.0956, 0.2372] | 6 |
| joint_endpoint_b0p5_r0p5 | d0 | coverage_95 | 0.2486 | [0.1049, 0.4234] | 6 |
| joint_endpoint_b0p5_r1 | d0 | coverage_95 | 0.2820 | [0.0926, 0.5121] | 6 |
| joint_endpoint_b0p1_r0p2 | d1 | coverage_95 | -0.4931 | [-0.5858, -0.4072] | 6 |
| joint_endpoint_b0p1_r0p5 | d1 | coverage_95 | -0.1578 | [-0.2113, -0.0950] | 6 |
| joint_endpoint_b0p1_r1 | d1 | coverage_95 | -0.0151 | [-0.0312, 0.0010] | 6 |
| joint_endpoint_b0p2_r0p2 | d1 | coverage_95 | -0.4175 | [-0.5275, -0.3297] | 6 |
| joint_endpoint_b0p2_r0p5 | d1 | coverage_95 | -0.1259 | [-0.1757, -0.0733] | 6 |
| joint_endpoint_b0p2_r1 | d1 | coverage_95 | -0.0151 | [-0.0312, 0.0011] | 6 |
| joint_endpoint_b0p35_r0p2 | d1 | coverage_95 | -0.2169 | [-0.3248, -0.1076] | 6 |
| joint_endpoint_b0p35_r0p5 | d1 | coverage_95 | -0.0660 | [-0.0972, -0.0324] | 6 |
| joint_endpoint_b0p35_r1 | d1 | coverage_95 | -0.0092 | [-0.0226, 0.0010] | 6 |
| joint_endpoint_b0p5_r0p2 | d1 | coverage_95 | -0.0879 | [-0.1708, -0.0287] | 6 |
| joint_endpoint_b0p5_r0p5 | d1 | coverage_95 | -0.0237 | [-0.0401, -0.0085] | 6 |
| joint_endpoint_b0p5_r1 | d1 | coverage_95 | -0.0034 | [-0.0101, 0.0022] | 6 |
| joint_endpoint_b0p1_r0p2 | d2 | coverage_95 | -0.6789 | [-0.7631, -0.5980] | 6 |
| joint_endpoint_b0p1_r0p5 | d2 | coverage_95 | -0.2633 | [-0.3427, -0.1857] | 6 |
| joint_endpoint_b0p1_r1 | d2 | coverage_95 | -0.0557 | [-0.0956, -0.0212] | 6 |
| joint_endpoint_b0p2_r0p2 | d2 | coverage_95 | -0.5340 | [-0.6193, -0.4365] | 6 |
| joint_endpoint_b0p2_r0p5 | d2 | coverage_95 | -0.2044 | [-0.2726, -0.1377] | 6 |
| joint_endpoint_b0p2_r1 | d2 | coverage_95 | -0.0500 | [-0.0881, -0.0164] | 6 |
| joint_endpoint_b0p35_r0p2 | d2 | coverage_95 | -0.2487 | [-0.3218, -0.1780] | 6 |
| joint_endpoint_b0p35_r0p5 | d2 | coverage_95 | -0.1226 | [-0.1684, -0.0818] | 6 |
| joint_endpoint_b0p35_r1 | d2 | coverage_95 | -0.0310 | [-0.0560, -0.0101] | 6 |
| joint_endpoint_b0p5_r0p2 | d2 | coverage_95 | -0.1038 | [-0.1268, -0.0751] | 6 |
| joint_endpoint_b0p5_r0p5 | d2 | coverage_95 | -0.0470 | [-0.0746, -0.0178] | 6 |
| joint_endpoint_b0p5_r1 | d2 | coverage_95 | -0.0148 | [-0.0256, -0.0046] | 6 |
| joint_endpoint_b0p1_r0p2 | d3 | coverage_95 | -0.7538 | [-0.8396, -0.6438] | 6 |
| joint_endpoint_b0p1_r0p5 | d3 | coverage_95 | -0.3554 | [-0.5276, -0.2020] | 6 |
| joint_endpoint_b0p1_r1 | d3 | coverage_95 | -0.1650 | [-0.2820, -0.0609] | 6 |
| joint_endpoint_b0p2_r0p2 | d3 | coverage_95 | -0.5171 | [-0.6623, -0.3304] | 6 |
| joint_endpoint_b0p2_r0p5 | d3 | coverage_95 | -0.3018 | [-0.4618, -0.1372] | 6 |
| joint_endpoint_b0p2_r1 | d3 | coverage_95 | -0.1353 | [-0.2186, -0.0448] | 6 |
| joint_endpoint_b0p35_r0p2 | d3 | coverage_95 | -0.3019 | [-0.4256, -0.1729] | 6 |
| joint_endpoint_b0p35_r0p5 | d3 | coverage_95 | -0.2184 | [-0.3112, -0.1149] | 6 |
| joint_endpoint_b0p35_r1 | d3 | coverage_95 | -0.0948 | [-0.1567, -0.0278] | 6 |
| joint_endpoint_b0p5_r0p2 | d3 | coverage_95 | -0.1697 | [-0.2526, -0.0750] | 6 |
| joint_endpoint_b0p5_r0p5 | d3 | coverage_95 | -0.1221 | [-0.1799, -0.0560] | 6 |
| joint_endpoint_b0p5_r1 | d3 | coverage_95 | -0.0517 | [-0.0934, -0.0104] | 6 |
| joint_endpoint_b0p1_r0p2 | all_distance_macro | interval_score_95 | 2.0381 | [1.3093, 2.6255] | 6 |
| joint_endpoint_b0p1_r0p5 | all_distance_macro | interval_score_95 | 0.0070 | [-0.8043, 0.6574] | 6 |
| joint_endpoint_b0p1_r1 | all_distance_macro | interval_score_95 | -0.3847 | [-1.1645, 0.1739] | 6 |
| joint_endpoint_b0p2_r0p2 | all_distance_macro | interval_score_95 | 1.1412 | [0.6193, 1.5744] | 6 |
| joint_endpoint_b0p2_r0p5 | all_distance_macro | interval_score_95 | -0.1751 | [-0.9689, 0.4738] | 6 |
| joint_endpoint_b0p2_r1 | all_distance_macro | interval_score_95 | -0.4010 | [-1.1849, 0.1200] | 6 |
| joint_endpoint_b0p35_r0p2 | all_distance_macro | interval_score_95 | 0.0593 | [-0.2406, 0.3092] | 6 |
| joint_endpoint_b0p35_r0p5 | all_distance_macro | interval_score_95 | -0.4539 | [-1.2121, 0.1062] | 6 |
| joint_endpoint_b0p35_r1 | all_distance_macro | interval_score_95 | -0.4113 | [-1.1539, 0.0844] | 6 |
| joint_endpoint_b0p5_r0p2 | all_distance_macro | interval_score_95 | -0.4280 | [-0.7594, -0.0932] | 6 |
| joint_endpoint_b0p5_r0p5 | all_distance_macro | interval_score_95 | -0.5819 | [-1.2597, -0.0541] | 6 |
| joint_endpoint_b0p5_r1 | all_distance_macro | interval_score_95 | -0.3890 | [-1.1010, 0.0604] | 6 |
| joint_endpoint_b0p1_r0p2 | d0 | interval_score_95 | -0.0763 | [-0.3668, 0.1053] | 6 |
| joint_endpoint_b0p1_r0p5 | d0 | interval_score_95 | -1.3458 | [-3.7059, 0.0037] | 6 |
| joint_endpoint_b0p1_r1 | d0 | interval_score_95 | -1.2747 | [-4.0607, 0.3357] | 6 |
| joint_endpoint_b0p2_r0p2 | d0 | interval_score_95 | -0.3273 | [-0.8410, -0.0331] | 6 |
| joint_endpoint_b0p2_r0p5 | d0 | interval_score_95 | -1.3709 | [-3.6537, 0.0205] | 6 |
| joint_endpoint_b0p2_r1 | d0 | interval_score_95 | -1.2678 | [-4.0705, 0.3911] | 6 |
| joint_endpoint_b0p35_r0p2 | d0 | interval_score_95 | -0.7774 | [-1.7855, -0.1164] | 6 |
| joint_endpoint_b0p35_r0p5 | d0 | interval_score_95 | -1.3921 | [-3.8102, 0.0442] | 6 |
| joint_endpoint_b0p35_r1 | d0 | interval_score_95 | -1.2487 | [-4.0319, 0.3831] | 6 |
| joint_endpoint_b0p5_r0p2 | d0 | interval_score_95 | -1.0625 | [-2.5380, -0.0661] | 6 |
| joint_endpoint_b0p5_r0p5 | d0 | interval_score_95 | -1.3856 | [-3.8285, 0.0737] | 6 |
| joint_endpoint_b0p5_r1 | d0 | interval_score_95 | -1.2180 | [-3.9729, 0.4252] | 6 |
| joint_endpoint_b0p1_r0p2 | d1 | interval_score_95 | 1.5255 | [0.5262, 2.9328] | 6 |
| joint_endpoint_b0p1_r0p5 | d1 | interval_score_95 | -0.2443 | [-0.4533, -0.0451] | 6 |
| joint_endpoint_b0p1_r1 | d1 | interval_score_95 | -0.1890 | [-0.2333, -0.1321] | 6 |
| joint_endpoint_b0p2_r0p2 | d1 | interval_score_95 | 0.8555 | [0.0436, 2.0091] | 6 |
| joint_endpoint_b0p2_r0p5 | d1 | interval_score_95 | -0.3226 | [-0.4778, -0.1814] | 6 |
| joint_endpoint_b0p2_r1 | d1 | interval_score_95 | -0.1751 | [-0.2161, -0.1245] | 6 |
| joint_endpoint_b0p35_r0p2 | d1 | interval_score_95 | 0.0037 | [-0.4934, 0.7590] | 6 |
| joint_endpoint_b0p35_r0p5 | d1 | interval_score_95 | -0.4201 | [-0.4970, -0.3291] | 6 |
| joint_endpoint_b0p35_r1 | d1 | interval_score_95 | -0.1359 | [-0.1615, -0.1049] | 6 |
| joint_endpoint_b0p5_r0p2 | d1 | interval_score_95 | -0.3260 | [-0.5563, 0.0755] | 6 |
| joint_endpoint_b0p5_r0p5 | d1 | interval_score_95 | -0.4008 | [-0.4373, -0.3592] | 6 |
| joint_endpoint_b0p5_r1 | d1 | interval_score_95 | -0.0674 | [-0.0833, -0.0439] | 6 |
| joint_endpoint_b0p1_r0p2 | d2 | interval_score_95 | 2.7133 | [1.7630, 3.5928] | 6 |
| joint_endpoint_b0p1_r0p5 | d2 | interval_score_95 | 0.2427 | [-0.1782, 0.6488] | 6 |
| joint_endpoint_b0p1_r1 | d2 | interval_score_95 | -0.2348 | [-0.3763, -0.0696] | 6 |
| joint_endpoint_b0p2_r0p2 | d2 | interval_score_95 | 1.4802 | [0.7831, 2.0893] | 6 |
| joint_endpoint_b0p2_r0p5 | d2 | interval_score_95 | 0.0160 | [-0.3209, 0.3475] | 6 |
| joint_endpoint_b0p2_r1 | d2 | interval_score_95 | -0.2393 | [-0.3599, -0.1026] | 6 |
| joint_endpoint_b0p35_r0p2 | d2 | interval_score_95 | 0.1112 | [-0.1815, 0.3371] | 6 |
| joint_endpoint_b0p35_r0p5 | d2 | interval_score_95 | -0.3068 | [-0.4930, -0.1191] | 6 |
| joint_endpoint_b0p35_r1 | d2 | interval_score_95 | -0.2274 | [-0.2991, -0.1556] | 6 |
| joint_endpoint_b0p5_r0p2 | d2 | interval_score_95 | -0.3889 | [-0.5128, -0.2524] | 6 |
| joint_endpoint_b0p5_r0p5 | d2 | interval_score_95 | -0.4236 | [-0.5021, -0.3477] | 6 |
| joint_endpoint_b0p5_r1 | d2 | interval_score_95 | -0.1625 | [-0.1967, -0.1258] | 6 |
| joint_endpoint_b0p1_r0p2 | d3 | interval_score_95 | 3.9897 | [1.9856, 5.6151] | 6 |
| joint_endpoint_b0p1_r0p5 | d3 | interval_score_95 | 1.3753 | [0.1942, 2.6326] | 6 |
| joint_endpoint_b0p1_r1 | d3 | interval_score_95 | 0.1597 | [-0.3384, 0.6596] | 6 |
| joint_endpoint_b0p2_r0p2 | d3 | interval_score_95 | 2.5564 | [0.9368, 3.8489] | 6 |
| joint_endpoint_b0p2_r0p5 | d3 | interval_score_95 | 0.9771 | [-0.0494, 1.9885] | 6 |
| joint_endpoint_b0p2_r1 | d3 | interval_score_95 | 0.0784 | [-0.3399, 0.5438] | 6 |
| joint_endpoint_b0p35_r0p2 | d3 | interval_score_95 | 0.8999 | [-0.0091, 1.8326] | 6 |
| joint_endpoint_b0p35_r0p5 | d3 | interval_score_95 | 0.3034 | [-0.3579, 1.1180] | 6 |
| joint_endpoint_b0p35_r1 | d3 | interval_score_95 | -0.0332 | [-0.3007, 0.2886] | 6 |
| joint_endpoint_b0p5_r0p2 | d3 | interval_score_95 | 0.0653 | [-0.5153, 0.7755] | 6 |
| joint_endpoint_b0p5_r0p5 | d3 | interval_score_95 | -0.1174 | [-0.5601, 0.4846] | 6 |
| joint_endpoint_b0p5_r1 | d3 | interval_score_95 | -0.1081 | [-0.2662, 0.1329] | 6 |
| joint_endpoint_b0p1_r0p2 | all_distance_macro | relative_rmse | -0.0021 | [-0.0120, 0.0043] | 6 |
| joint_endpoint_b0p1_r0p5 | all_distance_macro | relative_rmse | -0.0038 | [-0.0528, 0.0258] | 6 |
| joint_endpoint_b0p1_r1 | all_distance_macro | relative_rmse | 0.0234 | [-0.0501, 0.0795] | 6 |
| joint_endpoint_b0p2_r0p2 | all_distance_macro | relative_rmse | -0.0025 | [-0.0133, 0.0045] | 6 |
| joint_endpoint_b0p2_r0p5 | all_distance_macro | relative_rmse | -0.0038 | [-0.0499, 0.0264] | 6 |
| joint_endpoint_b0p2_r1 | all_distance_macro | relative_rmse | 0.0236 | [-0.0617, 0.0792] | 6 |
| joint_endpoint_b0p35_r0p2 | all_distance_macro | relative_rmse | -0.0031 | [-0.0163, 0.0053] | 6 |
| joint_endpoint_b0p35_r0p5 | all_distance_macro | relative_rmse | -0.0035 | [-0.0505, 0.0275] | 6 |
| joint_endpoint_b0p35_r1 | all_distance_macro | relative_rmse | 0.0246 | [-0.0480, 0.0809] | 6 |
| joint_endpoint_b0p5_r0p2 | all_distance_macro | relative_rmse | -0.0028 | [-0.0203, 0.0077] | 6 |
| joint_endpoint_b0p5_r0p5 | all_distance_macro | relative_rmse | -0.0021 | [-0.0507, 0.0300] | 6 |
| joint_endpoint_b0p5_r1 | all_distance_macro | relative_rmse | 0.0267 | [-0.0483, 0.0818] | 6 |
| joint_endpoint_b0p1_r0p2 | d0 | relative_rmse | -0.0029 | [-0.0132, 0.0033] | 6 |
| joint_endpoint_b0p1_r0p5 | d0 | relative_rmse | -0.0039 | [-0.0532, 0.0261] | 6 |
| joint_endpoint_b0p1_r1 | d0 | relative_rmse | 0.0322 | [-0.0617, 0.0926] | 6 |
| joint_endpoint_b0p2_r0p2 | d0 | relative_rmse | -0.0034 | [-0.0151, 0.0033] | 6 |
| joint_endpoint_b0p2_r0p5 | d0 | relative_rmse | -0.0039 | [-0.0525, 0.0264] | 6 |
| joint_endpoint_b0p2_r1 | d0 | relative_rmse | 0.0330 | [-0.0591, 0.0946] | 6 |
| joint_endpoint_b0p35_r0p2 | d0 | relative_rmse | -0.0045 | [-0.0198, 0.0042] | 6 |
| joint_endpoint_b0p35_r0p5 | d0 | relative_rmse | -0.0034 | [-0.0549, 0.0292] | 6 |
| joint_endpoint_b0p35_r1 | d0 | relative_rmse | 0.0351 | [-0.0572, 0.0970] | 6 |
| joint_endpoint_b0p5_r0p2 | d0 | relative_rmse | -0.0049 | [-0.0250, 0.0071] | 6 |
| joint_endpoint_b0p5_r0p5 | d0 | relative_rmse | -0.0020 | [-0.0578, 0.0316] | 6 |
| joint_endpoint_b0p5_r1 | d0 | relative_rmse | 0.0387 | [-0.0501, 0.1010] | 6 |
| joint_endpoint_b0p1_r0p2 | d1 | relative_rmse | -0.0020 | [-0.0114, 0.0042] | 6 |
| joint_endpoint_b0p1_r0p5 | d1 | relative_rmse | -0.0031 | [-0.0434, 0.0245] | 6 |
| joint_endpoint_b0p1_r1 | d1 | relative_rmse | 0.0242 | [-0.0389, 0.0734] | 6 |
| joint_endpoint_b0p2_r0p2 | d1 | relative_rmse | -0.0022 | [-0.0127, 0.0047] | 6 |
| joint_endpoint_b0p2_r0p5 | d1 | relative_rmse | -0.0028 | [-0.0436, 0.0253] | 6 |
| joint_endpoint_b0p2_r1 | d1 | relative_rmse | 0.0247 | [-0.0349, 0.0745] | 6 |
| joint_endpoint_b0p35_r0p2 | d1 | relative_rmse | -0.0023 | [-0.0159, 0.0059] | 6 |
| joint_endpoint_b0p35_r0p5 | d1 | relative_rmse | -0.0018 | [-0.0431, 0.0271] | 6 |
| joint_endpoint_b0p35_r1 | d1 | relative_rmse | 0.0265 | [-0.0328, 0.0751] | 6 |
| joint_endpoint_b0p5_r0p2 | d1 | relative_rmse | -0.0012 | [-0.0188, 0.0098] | 6 |
| joint_endpoint_b0p5_r0p5 | d1 | relative_rmse | 0.0010 | [-0.0415, 0.0301] | 6 |
| joint_endpoint_b0p5_r1 | d1 | relative_rmse | 0.0298 | [-0.0306, 0.0799] | 6 |
| joint_endpoint_b0p1_r0p2 | d2 | relative_rmse | -0.0020 | [-0.0118, 0.0046] | 6 |
| joint_endpoint_b0p1_r0p5 | d2 | relative_rmse | -0.0043 | [-0.0473, 0.0256] | 6 |
| joint_endpoint_b0p1_r1 | d2 | relative_rmse | 0.0176 | [-0.0513, 0.0737] | 6 |
| joint_endpoint_b0p2_r0p2 | d2 | relative_rmse | -0.0027 | [-0.0142, 0.0049] | 6 |
| joint_endpoint_b0p2_r0p5 | d2 | relative_rmse | -0.0047 | [-0.0487, 0.0259] | 6 |
| joint_endpoint_b0p2_r1 | d2 | relative_rmse | 0.0175 | [-0.0517, 0.0714] | 6 |
| joint_endpoint_b0p35_r0p2 | d2 | relative_rmse | -0.0040 | [-0.0214, 0.0063] | 6 |
| joint_endpoint_b0p35_r0p5 | d2 | relative_rmse | -0.0053 | [-0.0524, 0.0274] | 6 |
| joint_endpoint_b0p35_r1 | d2 | relative_rmse | 0.0178 | [-0.0510, 0.0737] | 6 |
| joint_endpoint_b0p5_r0p2 | d2 | relative_rmse | -0.0047 | [-0.0276, 0.0095] | 6 |
| joint_endpoint_b0p5_r0p5 | d2 | relative_rmse | -0.0050 | [-0.0568, 0.0305] | 6 |
| joint_endpoint_b0p5_r1 | d2 | relative_rmse | 0.0190 | [-0.0626, 0.0755] | 6 |
| joint_endpoint_b0p1_r0p2 | d3 | relative_rmse | -0.0014 | [-0.0116, 0.0053] | 6 |
| joint_endpoint_b0p1_r0p5 | d3 | relative_rmse | -0.0039 | [-0.0533, 0.0294] | 6 |
| joint_endpoint_b0p1_r1 | d3 | relative_rmse | 0.0195 | [-0.0576, 0.0810] | 6 |
| joint_endpoint_b0p2_r0p2 | d3 | relative_rmse | -0.0015 | [-0.0115, 0.0049] | 6 |
| joint_endpoint_b0p2_r0p5 | d3 | relative_rmse | -0.0039 | [-0.0528, 0.0290] | 6 |
| joint_endpoint_b0p2_r1 | d3 | relative_rmse | 0.0193 | [-0.0566, 0.0794] | 6 |
| joint_endpoint_b0p35_r0p2 | d3 | relative_rmse | -0.0014 | [-0.0112, 0.0047] | 6 |
| joint_endpoint_b0p35_r0p5 | d3 | relative_rmse | -0.0037 | [-0.0549, 0.0282] | 6 |
| joint_endpoint_b0p35_r1 | d3 | relative_rmse | 0.0189 | [-0.0598, 0.0776] | 6 |
| joint_endpoint_b0p5_r0p2 | d3 | relative_rmse | -0.0006 | [-0.0105, 0.0060] | 6 |
| joint_endpoint_b0p5_r0p5 | d3 | relative_rmse | -0.0026 | [-0.0503, 0.0281] | 6 |
| joint_endpoint_b0p5_r1 | d3 | relative_rmse | 0.0192 | [-0.0557, 0.0762] | 6 |

## Limits and next gate

The pilot uses one beta, time interval, sequencing depth, and ABM time step. It does not validate 22-chromosome graphs, longer distances outside this domain, uncertain beta, overdispersion, dropout, or model selection on sealed tests. Intervals condition on selected hyperparameters. Improve structural recovery and numerical stability first; then expand independent landscapes and prespecified scenarios before claiming general accuracy gains.

Candidate changes are not accepted solely because macro RMSE decreases: inspect each exact distance, fit failures, complete prediction/interval coverage, interval score and calibration together. Any distance deterioration is retained in this report.

Generated: 2026-11-03 23:41:18 EST
