# Sparse two-timepoint benchmark

Scores fitness-landscape estimates from two sparse karyotype samples against a
known truth, separately for observed karyotypes and their one- and two-step
neighbours. It is the development scoreboard for alfak2: every candidate change
is compared here against the original ALFA-K (`alfakR`) and the current alfak2
pipeline before it is considered for the 22-chromosome HPC confirmation.

## Design

The simulated data are built to resemble real two-timepoint tumour samples. The
reference is a set of 20 GBM inputs (22 chromosomes, copy numbers 1–4). Among
the 15 with at least 100 cells at both timepoints, and for the simulated
`typical_500_700_360d` datasets of the development landscapes:

| Statistic (median, quartiles) | Real GBM | Simulated |
| --- | --- | --- |
| Karyotypes observed | 65 (55–182) | 162 (45–328) |
| Fraction seen at both timepoints | 0.31 (0.12–0.53) | 0.13 (0.08–0.27) |
| Total-variation distance between samples | 0.82 (0.47–0.91) | 0.71 (0.10–0.85) |
| Largest clone | 0.13 (0.09–0.21) | 0.25 (0.13–0.81) |
| Karyotypes seen in one cell only | 0.27 (0.04–0.34) | 0.57 (0.53–0.61) |
| SD of log frequency change, shared karyotypes with ≥20 cells | 1.73 (0.98–2.25) | 2.12 (1.71–2.50) |

`run.R --stage=profile --real-dir=DIR` recomputes both columns. The simulated
samples have more single-cell karyotypes and correspondingly fewer shared ones
than the real inputs; the other statistics overlap.

Choices that follow from this comparison:

- **22 chromosomes.** With 8 chromosomes the population is trapped on a local
  peak within about 700 days and samples become 98% one karyotype.
- **First sample at day 1440.** Earlier, clones replace each other completely
  within a year (no karyotype is shared between the samples); later, the
  smoother landscapes freeze as well.
- **Wavelengths 0.4, 0.8, 1.2, 1.6**, three landscapes each. The 0.4 class is
  the frozen, near-monoclonal regime that about one real sample in ten
  resembles. The others give ongoing adaptation.
- **Nullisomy is lethal.** Daughters with a copy number outside 1–8 die. The
  original ALFA-K simulator lets nullisomic cells live; its populations then
  fix a zero-copy chromosome, which never occurs in the real inputs and which
  `alfakR` cannot fit because it discards every nullisomic neighbour.

Truth is the Gaussian-random-field landscape of the ALFA-K simulations, and the
simulator (`abm.cpp`) keeps their division and missegregation rules.
`run.R --stage=validate` checks it against closed-form results: the daughter
distribution against the missegregation kernel, viability losses, growth,
selection between two clones, passaging, and reproducibility.

A landscape has 8 "patients": independent trajectories from the diploid founder
on the same landscape. Patient 1 makes the single-sample scoreboard (96
development datasets); all eight form the cohort for cross-sample learning.

| Scenario | Cells at t0 / t1 | Interval (days) |
| --- | --- | --- |
| `very_sparse_50_50_180d` | 50 / 50 | 180 |
| `sparse_200_200_180d` | 200 / 200 | 180 |
| `sparse_200_200_360d` | 200 / 200 | 360 |
| `typical_500_700_360d` | 500 / 700 | 360 |
| `unbalanced_100_3000_360d` | 100 / 3000 | 360 |
| `unbalanced_3000_100_360d` | 3000 / 100 | 360 |
| `long_500_700_1080d` | 500 / 700 | 1080 |
| `bottleneck_500_700_360d` | 500 / 700, regrown from 200 surviving cells | 360 |

## Scoring

Scored nodes are defined from the counts alone: every observed karyotype
(distance 0) and uniform samples of up to 2000 one-step and 2000 two-step
neighbours with copy numbers 1–8. All methods are scored on the same nodes; a
node a method cannot predict is scored as "no fitness difference" and counted
in `prediction_fraction`.

Truth is the per-day growth rate the simulated counts follow,
`log(1 + fitness * delta_t) / delta_t`. Truth and prediction are both centred
on the count-weighted mean of the karyotypes with at least 5 cells, so a
method's absolute scale does not enter.

Per distance: RMSE, `rmse_zero` (the RMSE of predicting no differences),
regression slope of prediction on truth, Spearman correlation, AUC for
"fitter than the reference", sign accuracy of one-step fitness differences
into the shell, and 95% interval coverage and interval score where a method
reports uncertainty.

Datasets are averaged within landscape and scenario, then landscapes are
weighted equally. Comparisons are paired by landscape with a landscape
bootstrap interval.

**Acceptance rule for a candidate.** Pooled over scenarios, at each of
distances 0, 1 and 2:

1. RMSE below the no-information answer and below the current alfak2 pipeline
   (95% interval of the paired difference below zero);
2. Spearman correlation above the current alfak2 pipeline and above the
   original ALFA-K (interval above zero);
3. no scenario with RMSE more than 10% above the current alfak2 pipeline.

Coverage is reported separately, with 90–98% as the target. The bar for error
is the no-information answer rather than ALFA-K because ALFA-K's RMSE is
dominated by a few unconstrained fits and is far above it, while the current
alfak2 pipeline sits on it. `run.R --stage=report` evaluates the rule for every
method that is not one of these baselines.

## Running

Both packages must be installed, or `STT_R_LIB` must point to a library that
holds the builds under test.

```sh
R=benchmark/sparse_two_timepoint/run.R
Rscript $R --stage=validate
Rscript $R --stage=truth    --split=dev --cohort
Rscript $R --stage=datasets --split=dev --cohort
Rscript $R --stage=fit      --split=dev --methods=zero,alfakR_minobs5,alfak2_two_layer_fixed
Rscript $R --stage=score    --split=dev
Rscript $R --stage=report   --split=dev
```

Results go to `benchmark/results/sparse_two_timepoint/` (not tracked). Fits are
cached per dataset and method, so a stage can be rerun after adding a method.

`--workers` defaults to 2. An alfak2 two-layer fit at 22 chromosomes needs
about 1 GB. `alfakR` needs up to 4.5 GB and 15 minutes on the largest datasets
(`minobs = 5` at 3000 cells), so raise the worker count only with memory to
spare. `alfak2_two_layer_cv`, the package's default leave-one-anchor-out
selection of graph hyperparameters, is defined but not part of a normal round:
one dataset with 97 observed karyotypes did not finish in 25 minutes, against
10 seconds with the hyperparameters fixed.
The held-out landscapes are for the final confirmation only and need
`--split=heldout --unseal-heldout`.

## Limits

The simulation has no spatial sampling bias and no karyotype-calling error,
and treatment is represented only by the bottleneck scenario. Passing this
benchmark is necessary for a candidate, not evidence of accuracy on patient
samples.
