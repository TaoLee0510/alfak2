# Round 2: initial-composition and growth-information oracles

Exact initial-frequency information changes the inferred common birth baseline substantially in this development sample, but it does not recover the smooth or interaction d0 point errors. Population-growth information also shifts the baseline; its benefits and harms vary by family and distance. These interventions therefore do not isolate one sufficient explanation for the earlier near-field regressions.

## Frozen information interventions

The experiment uses six new three-chromosome landscapes (two smooth, two rough, two interaction), one exact Gillespie trajectory per landscape and two independently sampled count pairs per trajectory: 12 datasets and four paired information arms, 48 fits. Every arm uses exactly the same counts, full 27-node domain, v1 fitness prior and t0-count reference. Only the initial-frequency oracle and/or external total-population observation changes. The exact initial oracle fixes all probabilities, including true zeros; it supplies no pseudo-counts. Population growth uses the actual totals from the same trajectory and a prespecified log-growth SD of 0.1.

All outputs are development diagnostics. Extra-information arms are not fair algorithmic improvements over count-only estimation. The v1 source tree is frozen; only the explicit diagnostic oracle mechanism is overlaid. No formal test is opened or promoted here, and no interval calibration is applied.

All 48 fits completed with trusted covariance. Single-fit runtimes ranged from 0.51 to 1.27 seconds.

## Common baseline and initial frequencies

The birth column is the conditional mode averaged with the fitted initial frequencies; it is compared with the true initial-frequency-weighted birth rate. Log growth is the model's conditional-mode prediction compared with the realized finite-population log growth. These are distinct quantities.

| Information | Mode-weighted birth | True weighted birth | Predicted log growth | Realized log growth | Initial-frequency RMSE |
|---|---:|---:|---:|---:|---:|
| Counts only | 0.84524 | 0.99194 | 1.26644 | 1.47257 | 0.005870 |
| Exact initial frequencies | 0.99187 | 0.99194 | 1.48562 | 1.47257 | 0.000000 |
| Population growth | 0.97915 | 0.99194 | 1.46616 | 1.47257 | 0.005815 |
| Both information sources | 0.98318 | 0.99194 | 1.47255 | 1.47257 | 0.000000 |

## Initial mass by observed appearance

Groups use the unchanged count pair. A zero initial count is not asserted to mean biological absence: late-only rows retain nonzero true initial mass in this experiment. Means below include datasets where a group has no nodes, preserving total mass accounting.

| Information | Count support | Estimated initial mass | True initial mass |
|---|---|---:|---:|
| Counts only | t0_positive | 0.9943886 | 0.9960417 |
| Counts only | t1_only | 0.0050681 | 0.0039583 |
| Counts only | both_zero | 0.0005433 | 0.0000000 |
| Exact initial frequencies | t0_positive | 0.9960417 | 0.9960417 |
| Exact initial frequencies | t1_only | 0.0039583 | 0.0039583 |
| Exact initial frequencies | both_zero | 0.0000000 | 0.0000000 |
| Population growth | t0_positive | 0.9945555 | 0.9960417 |
| Population growth | t1_only | 0.0049063 | 0.0039583 |
| Population growth | both_zero | 0.0005382 | 0.0000000 |
| Both information sources | t0_positive | 0.9960417 | 0.9960417 |
| Both information sources | t1_only | 0.0039583 | 0.0039583 |
| Both information sources | both_zero | 0.0000000 | 0.0000000 |

## Fitness and raw interval results

Distance is measured from the observed support in this finite-population dataset, identically for all four arms. The small development design contains d0–d2; it provides no d3 evidence. Scores use the full joint lognormal Laplace distribution before subtracting the original t0-count reference in each posterior draw.

| Family | Distance | Information | Relative RMSE | 95% coverage | 95% interval score |
|---|---|---|---:|---:|---:|
| smooth | d0 | Counts only | 0.09957 | 1.000 | 0.51913 |
| smooth | d0 | Exact initial frequencies | 0.10021 | 1.000 | 0.51518 |
| smooth | d0 | Population growth | 0.09896 | 0.969 | 0.48793 |
| smooth | d0 | Both information sources | 0.10110 | 0.969 | 0.46276 |
| smooth | d1 | Counts only | 0.16800 | 0.942 | 0.70811 |
| smooth | d1 | Exact initial frequencies | 0.16645 | 1.000 | 0.64654 |
| smooth | d1 | Population growth | 0.16683 | 0.846 | 0.77585 |
| smooth | d1 | Both information sources | 0.16526 | 0.885 | 0.73874 |
| smooth | d2 | Counts only | 0.20138 | 0.858 | 1.12682 |
| smooth | d2 | Exact initial frequencies | 0.19980 | 0.927 | 0.83272 |
| smooth | d2 | Population growth | 0.19953 | 0.785 | 1.16358 |
| smooth | d2 | Both information sources | 0.19834 | 0.796 | 1.21933 |
| rough | d0 | Counts only | 0.22986 | 0.707 | 1.92705 |
| rough | d0 | Exact initial frequencies | 0.21925 | 0.844 | 1.32527 |
| rough | d0 | Population growth | 0.22608 | 0.676 | 2.16864 |
| rough | d0 | Both information sources | 0.21939 | 0.676 | 2.04652 |
| rough | d1 | Counts only | 0.23661 | 0.804 | 1.78304 |
| rough | d1 | Exact initial frequencies | 0.22555 | 0.904 | 1.27776 |
| rough | d1 | Population growth | 0.23079 | 0.783 | 1.79399 |
| rough | d1 | Both information sources | 0.22350 | 0.804 | 1.69241 |
| rough | d2 | Counts only | 0.13408 | 1.000 | 0.74048 |
| rough | d2 | Exact initial frequencies | 0.12019 | 1.000 | 0.80883 |
| rough | d2 | Population growth | 0.12667 | 0.958 | 0.67076 |
| rough | d2 | Both information sources | 0.12004 | 1.000 | 0.64063 |
| interaction | d0 | Counts only | 0.10338 | 0.975 | 0.88230 |
| interaction | d0 | Exact initial frequencies | 0.11085 | 0.975 | 0.80902 |
| interaction | d0 | Population growth | 0.10379 | 0.975 | 0.78199 |
| interaction | d0 | Both information sources | 0.11213 | 0.975 | 0.74043 |
| interaction | d1 | Counts only | 0.23392 | 0.942 | 1.26627 |
| interaction | d1 | Exact initial frequencies | 0.23504 | 0.942 | 1.39301 |
| interaction | d1 | Population growth | 0.23497 | 0.846 | 1.44550 |
| interaction | d1 | Both information sources | 0.23796 | 0.808 | 1.48906 |
| interaction | d2 | Counts only | 0.28418 | 0.750 | 2.51875 |
| interaction | d2 | Exact initial frequencies | 0.29205 | 0.917 | 1.53687 |
| interaction | d2 | Population growth | 0.28394 | 0.812 | 2.35253 |
| interaction | d2 | Both information sources | 0.28347 | 0.812 | 2.30016 |

## Repeat-error decomposition

The CSV fixes each node's distance and appearance group using sampling replicate 1, then decomposes mean squared error across the two sampling repeats into squared repeat-mean error plus within-repeat error variance. Their numerical sum equals the reported MSE. Squared repeat-mean error is not an unbiased estimate of population bias squared; a shared ABM realization and finite sampling contribute to it. One biological path per landscape cannot separate systematic prior bias from shared process fluctuations. Native-distance metrics above and fixed-reference repeat strata are intentionally reported separately.

Maximum algebraic residual in the repeat decomposition: 2.19e-17.

## Independent small-graph posterior reference

A separate fixed two-node graph uses a Gaussian log-birth prior and an initial-logit prior specified before the fixed count observations. The normalized coefficient-plus-GMRF prior is analytically reduced to its bivariate Gaussian log-birth distribution. At beta=0 the frequency and total-growth formulas are analytic. Gauss-Hermite importance quadrature integrates the original prior times that exact likelihood, including numerical integration of initial log odds when unknown; the broad Laplace-guided proposal does not replace the target density. Orders 40, 64, 96 and 160 check moment convergence, with a finer/wider log-odds integral at the higher orders. Independent prior importance sampling provides a second check with Monte Carlo standard errors.

This is a fixed-data posterior comparison, not SBC or repeated-sampling coverage. Indicator-based Gauss-Hermite quantiles converge slowly and are not certified or used to make interval-coverage claims.

| Information | Quantity | Reference mean | Laplace mean | Reference SD | Laplace / reference SD |
|---|---|---:|---:|---:|---:|
| Counts only | birth_1 | 1.032512 | 1.145990 | 0.558653 | 1.2261 |
| Counts only | birth_2 | 1.197805 | 1.357851 | 0.665312 | 1.3152 |
| Counts only | relative_contrast_1 | -0.082646 | -0.105930 | 0.094949 | 1.4968 |
| Exact initial frequencies | birth_1 | 1.208301 | 1.454433 | 0.623684 | 1.2352 |
| Exact initial frequencies | birth_2 | 1.500448 | 1.828202 | 0.747666 | 1.2567 |
| Exact initial frequencies | relative_contrast_1 | -0.146073 | -0.186884 | 0.099042 | 1.3399 |
| Population growth | birth_1 | 0.967211 | 0.976786 | 0.106388 | 1.0116 |
| Population growth | birth_2 | 1.119713 | 1.129198 | 0.092330 | 1.0001 |
| Population growth | relative_contrast_1 | -0.076251 | -0.076206 | 0.073128 | 1.0113 |
| Both information sources | birth_1 | 0.931162 | 0.942422 | 0.093870 | 1.0175 |
| Both information sources | birth_2 | 1.170496 | 1.184015 | 0.090937 | 1.0021 |
| Both information sources | relative_contrast_1 | -0.119667 | -0.120797 | 0.065163 | 1.0131 |

For quadrature order 96 versus 160, the maximum birth-mean difference is 3.04e-06, covariance difference 2.51e-05, and contrast-SD difference 2.81e-07. The independent 500,000-prior-draw checks have ESS 4669–83644; their maximum mean and SD discrepancies are 2.68 and 1.66 estimated Monte Carlo standard errors.

Without growth information, this example's Laplace/lognormal contrast variance is larger than the numerical posterior reference, not smaller. With growth information the contrast moments are much closer. Non-Gaussian baseline uncertainty is therefore a real approximation concern, but this example does not support attributing all observed undercoverage to Laplace underdispersion.

## Limits and use in model design

These six small landscapes do not reproduce the original four-chromosome family regressions exactly. They test information mechanisms while keeping each dataset and estimator fixed. Exact initial frequencies slightly help rough point estimation but do not repair smooth/interaction d0 point errors; growth information can tighten or shift intervals and sometimes worsens coverage. Candidate design should therefore consider both composition/baseline confounding and the fitness prior, with a separate process-noise assessment. No oracle gain is counted as an algorithmic promotion result.

Generated: 2026-11-04 00:35:53 EST
