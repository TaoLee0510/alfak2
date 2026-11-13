# Unobserved karyotypes predicted too fit at intermediate ruggedness

In the paper-sized young populations of `../lattice_residual` (wavelength 0.4,
first sample at day 120) `fit_sparse_joint()` at af362af predicts unobserved
karyotypes too fit: bias +0.06 at distance one and +0.11 at distance two from
the data, with 95% intervals covering 0.87 and 0.82 of them (0.83 and 0.73 at
6-day intervals). The same holds, less strongly, on the development-sized
populations of `../short_interval` and at wavelengths 0.4-0.8 on the
development benchmark. At wavelength 0.4 the predicted drop below the major
karyotypes is about 60% of the true drop, for unobserved karyotypes inside
the domain and for those predicted by the trend outside it alike; at 0.8 about
90%.

## Data

Datasets are linked from the earlier studies, not copied or regenerated
(`run.R --stage=link`): the paper-sized set of `../lattice_residual`
("rugged", 48 datasets at wavelengths 0.4 and 0.8) and the development-sized
diagnosis set of `../short_interval` ("diagnosis", 60 datasets). Fits of this
study are in its own `results/<source>/fits`.

## Probes (af362af, wavelength 0.4)

Predicted drop as a share of the true drop, paper-sized set:

| arm | observed, shell 1 | unobserved, in domain | outside, shell 2 | outside, shell 3 |
| --- | --- | --- | --- | --- |
| default | 0.74 | 0.59 | 0.62 | 0.59 |
| `balance_sd = 0.5` / `0.25` | 0.74 / 0.73 | 0.58 / 0.61 | 0.62 / 0.63 | 0.59 / 0.60 |
| history = true age | 0.82 | 0.61 | 0.66 | 0.64 |
| `initial_state = "free"` | 0.81 | 0.54 | 0.63 | 0.61 |
| residual scale held at 0.12 | 0.38 | 0.79 | 0.52 | 0.39 |

Neither the spread around balance, nor the history, nor the balance prior
itself is the cause. A residual scale near the true spread of fitness among
neighbours (the fitted one is about 0.02) lets zero counts pull individual
unobserved karyotypes inside the domain down, but the trend, which predicts
everything outside the domain, becomes shallower and observed neighbours are
predicted near the majors. RMSE does not improve.

## Calibration (step 1)

`calibration.R collect` pools the default fits of the development benchmark's
development landscapes, the diagnosis set and the rugged set (372 datasets);
`calibration_forms.R` fits forms of the reported uncertainty. With the
shipped constants, unobserved karyotypes are covered 0.79-0.95 at wavelengths
0.4-0.8 and 0.98-1.00 on smoother (1.2-1.6) and more rugged (0.2) landscapes.
The factor the intervals would need varies with the regime, and neither the
predicted drop nor the interval told the regimes apart. Refitting the two
constants by likelihood makes the undercovered cells worse (it narrows the
overcovered ones), and so does a separate deficit fitted by likelihood.

Decision (user, 2026-11-13): a deficit of its own for karyotypes seen in
neither sample, 1.1, the smallest that keeps every simulated non-stress
regime at 90% or more (observed karyotypes keep 0.6, the floor 0.009):

| deficit for unobserved | worst non-stress cell | cells below 0.93 (of 22) | worst stress cell | width, median |
| --- | --- | --- | --- | --- |
| 0.6 (shipped) | 0.785 | 7 | 0.854 | 1.00 |
| 1.0 | 0.886 | 4 | 0.922 | 1.33 |
| 1.1 | 0.901 | 3 | 0.932 | 1.42 |
| 1.2 | 0.913 | 2 | 0.941 | 1.51 |

Point predictions are unchanged. Build cal1 (`../lattice_residual/results/lib`)
reproduces af362af's predictions and model standard deviations exactly.

## Mode or mean (step 2)

`laplace_mean.R` samples the posterior that `fit_sparse_joint()` summarises
by its mode: all latent variables at the fitted scales, by Hamiltonian Monte
Carlo in coordinates whitened by the model's expected information.
`laplace_mean_report.R` compares mode, mean and median with the truth. The
chains mix well (acceptance 0.68-0.93, effective sample size of the mean
unobserved fitness 27-3000, median about 700-1300), and the mode they are
started from equals the fit's in every dataset.

The posterior is strongly skewed towards deeper drops: its mean lies well
below its mode for unobserved and weakly observed karyotypes, and the median
with the mean. Predicted drop as a share of the true drop, and RMSE of the
mean against that of the mode:

| set | wavelength | unobserved in domain, mode / mean | outside, shell 2-3, mode / mean | RMSE of the mean, unobserved |
| --- | --- | --- | --- | --- |
| paper-sized, young | 0.4 | 0.60 / 0.86 | 0.61-0.62 / 0.92-0.95 | -13% to -16% |
| paper-sized, young | 0.8 | 0.90 / 1.37 | 0.83-0.89 / 1.43-1.51 | +2% to +4% |
| development, 180-1080 days | 0.4 | 0.95 / 1.31 | 0.84-0.90 / 1.11-1.21 | +5% to +31% (observed +43%) |
| development | 0.8 | 0.63 / 0.93 | 0.62 / 0.95 | -16% to -17% |
| development | 1.2 | 0.71 / 0.97 | 0.74-0.78 / 1.07-1.13 | -7% to +2% |
| development | 1.6 | 0.77 / 1.03 | 0.83-0.95 / 1.16-1.32 | 0% to +16% |

Neither the mode nor the mean is where the truth lies: depending on landscape
and population it is near one or the other. Reporting the mean would break
the development benchmark's acceptance rule (no core cell more than 10%
worse) at wavelength 0.4. The sampled 95% intervals, conditional on the
fitted scales, cover only 0.36-0.72 of unobserved karyotypes. The depth of
unobserved karyotypes is not in the data; the mode is kept, and the
calibration of step 1 reports the uncertainty.

## Running

```sh
R=validation/sparse_joint/unobserved_bias/run.R
Rscript $R --source=rugged --stage=link
Rscript $R --source=rugged --stage=fit --methods=zero,base,bsd05,bsd025,age,free,res003,res006,res012 --workers=6
Rscript $R --source=rugged --stage=score --workers=6
Rscript validation/sparse_joint/unobserved_bias/calibration.R collect --workers=6
Rscript validation/sparse_joint/unobserved_bias/calibration_forms.R
Rscript validation/sparse_joint/unobserved_bias/laplace_mean.R <dataset ids> --source=rugged --workers=6
Rscript validation/sparse_joint/unobserved_bias/laplace_mean_report.R --source=rugged
```

The residual-scale arms wrap an internal function of the package inside each
worker (`ub_with_residual_scale()`); they are probes, not options.
