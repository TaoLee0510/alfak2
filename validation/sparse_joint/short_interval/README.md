# Short sampling intervals: diagnosis of `fit_sparse_joint()`

The ALFA-K paper benchmark (`alfak2-paper-benchmark`) found two weaknesses of
the sparse joint model at 8327282, both outside the range it was developed on
(intervals of 180–1080 days, wavelengths 0.4–1.6):

1. **Short intervals.** With samples one passage apart (about 6 days), the
   predictions for karyotypes held out of the fit are more dispersed than the
   truth supports: the slope of truth on held-out prediction is 0.06–0.56.
2. **Few frequent karyotypes.** On the most rugged landscapes (wavelength 0.2,
   3–4 frequent karyotypes) a held-out dominant karyotype is predicted 3–4
   standard deviations too low.

Short intervals are in scope: cell-culture time courses are sampled days or
weeks apart. This study reproduces the weaknesses on alfak2's own simulator,
tests two candidate causes, and compares opt-in changes for each.

## Hypotheses

- **H1, history.** With `initial_state = "balance"`, a non-major karyotype
  starts at what its major neighbours fed it over `history`, which defaults
  to the sampling interval. Over a short interval that prior is small and
  independent of the karyotype's fitness, while the population has been
  evolving for much longer. The model can reconcile the two only by raising
  the fitness of observed neighbours or the scales, which disperses the
  predictions. Over 180 days and more, the balance has reached its long-run
  level for less fit neighbours and the history hardly matters.
- **H2, trend.** The trend has a dosage and a curvature term per chromosome.
  Their scales have a log-normal prior centred on 1% per step, which does not
  pull them towards zero when the counts carry little information about them.

## Data

The development simulator of `benchmark/sparse_two_timepoint` (sourced, not
modified), with its own landscapes (seed base 93000) and output directory:

- wavelengths 0.2, 0.4, 0.8 and 1.6, five landscapes each, one patient;
- first sample at day 120 (as young as the paper's ABM populations) or day
  1440 (as in the development benchmark): separates a short interval from a
  young population;
- intervals of 6, 18 and 42 days (one, three and seven passages of the paper's
  ABM), 1000 cells at each sample.

120 datasets. The paper's ABM data are not used, so that its benchmark stays
independent.

## Arms

| Arm | Build | What changes |
| --- | --- | --- |
| `balance` | 8327282 | Nothing: the current default |
| `balance_history_age` | 8327282 | `history` = true age of the population (probe; only a simulation knows it) |
| `free` | 8327282 | `initial_state = "free"` (reference without the balance prior) |
| `history_grid` | si1 | `history` chosen by marginal likelihood among the interval and 30–3000 days |
| `trend_centre_low` | si1 | Trend-scale prior centred on 0.1% instead of 1% (probe) |
| `trend_exponential` | si1 | Exponential penalty on the trend scales, mean 1% |
| `trend_select` | si1 | Trend terms (shell; + dosage; + curvature) chosen by marginal likelihood |
| `grid_exponential`, `grid_select` | si1 | `history_grid` together with each version of B |
| `history_profile` | si1 | `history` chosen by marginal likelihood among the interval and those of 30–3000 days that are longer (final version of A) |
| `history_estimate` | si2 | The history estimated as an outer parameter with the scales, bounded below by the interval (rejected) |
| `default_new` | si4 | The package with `history_profile` as its default, at its defaults |

Build si4 is the package as changed: the history candidates are the default
of `fit_sparse_joint()`, and the options of si1 and si2 are gone. With
`history = dt` it reproduces 8327282 exactly, and at its default
`history_profile` exactly (checked on seven datasets).

Build si1 is 8327282 with the options `history` (several candidates),
`trend_terms`, `trend_scale_prior` and `trend_scale_centre` of
`fit_sparse_joint()`; at their defaults it reproduces 8327282 exactly. Both
builds are installed under `results/lib/`, never into the source tree.

## Scores

Per fit (`metrics_fit.csv`): the cross-validation of `as_alfakR()` (minobs 10)
against the truth, as `alfak2-paper-benchmark/scripts/cv_score_analysis.R`
computes it (held-out slope and correlation, bias of held-out frequent
karyotypes, share of the error that is a common offset per fold), the
landscape's r and R² against the truth, the fitted scales and the bias of
observed karyotypes. Per distance (`metrics_distance.csv`): the development
benchmark's scores at distances 0, 1 and 2.

## Results (2026-11-12)

Medians over the five landscapes of a cell; "1-step RMSE" is the RMSE of
unobserved one-step neighbours divided by that of predicting no differences.

**The weakness reproduces.** At 6 days the current default predicts worse
than no information (1-step RMSE 2.6 at day 120, 2.3 at day 1440); at 18 and
42 days it is 0.5–0.65. Observed non-major karyotypes are overestimated by
0.10–0.19 per day at 6 days and by about 0.01 at 18–42 days. Weakness 2
reproduces only in part: at wavelength 0.2 the landscape r is about 0.2 in
every cell, but these populations have one or two frequent karyotypes, so the
cross-validation has a single fold.

**H1 holds; H2 does not.**

| 6-day interval | default | history = true age | `history_profile` | `free` |
| --- | --- | --- | --- | --- |
| 1-step RMSE, day 120 | 2.57 | 0.64 | 0.66 | 0.72 |
| 1-step RMSE, day 1440 | 2.32 | 0.48 | 0.63 | 0.52 |
| landscape r, day 120, w 0.8 / 1.6 | 0.22 / 0.45 | 0.73 / 0.78 | 0.73 / 0.70 | 0.54 / 0.42 |
| held-out slope, day 120, w 0.8 / 1.6 | 0.40 / 0.32 | 0.88 / 0.89 | 0.77 / 0.69 | 0.81 / 0.91 |

With the history corrected, the overestimate of observed neighbours vanishes
and the dosage and residual scales fall. None of the three versions of B
helps: a lower centre and the exponential prior cost up to 0.26 in held-out
correlation at 6 days, nested selection nearly always keeps the full trend,
and A + B equals A.

**The choice of history must not go below the interval.** Chosen freely
among the interval and 30–3000 days, the marginal likelihood picks 30–100
days almost everywhere. On the development benchmark that is the worst
choice for intervals of 180–1080 days (fixed histories of 30 and 100 days
cost up to 26% and 5% RMSE; 300–3000 days are within a few percent of the
interval): `history_grid` exceeds the 10% rule in 2 of 24 core cells.
`history_profile` keeps only histories at least as long as the interval.

**Development benchmark, `history_profile` against the default** (HPC, 168
datasets): core cells within +2.6% RMSE (largest `long_500_700_1080d`), pooled
core RMSE −0.2% to +0.3%, Spearman unchanged; stress scenarios RMSE −2% to −5%
pooled, `short_500_700_30d` −16% to −27%, `short_500_700_60d` −2% to −3%.
Diagnosis data: 1-step RMSE −69% (day 120) and −62% (day 1440) at 6 days,
−6% at 18 days, unchanged at 42 days.

**Calibration.** Refitted on the core scenarios, the constants change little
(floor 0.0067, deficit 0.534 for `history_profile`; 0.0066 and 0.514 for the
default, against the shipped 0.009 and 0.6). The shipped constants give 95–98%
coverage on the core scenarios and 94–99% on the short stress scenarios; the
refitted ones bring the core closer to 95% but `short_500_700_30d` down to 90%
at distance 2. The shipped constants are kept.

**The history as an outer parameter is unreliable.** Estimated with the
scales in one optimisation (`history_estimate`), it was 3.6 times faster than
the grid, but the joint optimum has separate basins, a short history with
large scales and a long one with small scales: a quarter of the development
fits and a tenth of the diagnosis fits ended with a lower marginal likelihood
than the grid's best, by up to 222. Starting at either end of the range, or
warm-starting each grid candidate from the previous one, did not fix it.

**Cost.** One fit per candidate history: on the development benchmark the
median fit time rose from 15 to 55 seconds (2.9 times in total); at 6–18 day
intervals, with six candidates, about 13 times.

## Final confirmation (2026-11-13)

Run once, on the package with the new default (`default_new`, build si4)
against 8327282 (`balance`).

- **Held-out landscapes of the development benchmark** (`devbench.R
  --split=heldout --unseal-heldout`, HPC): no core cell above +0.3% RMSE;
  pooled core RMSE 0.6–1.1% lower, Spearman unchanged, coverage 94–98%.
  Stress scenarios 2–5% lower RMSE and 0.015–0.025 higher Spearman;
  `short_500_700_30d` 18–27% lower, `short_500_700_60d` about 3% lower. The
  held-out landscapes had been used once before, on 2026-11-10/11, to confirm
  the sparse joint model; nothing in this study was tuned on them.
- **New short-interval landscapes** (`run.R --set=heldout --unseal-heldout`,
  seed base 94000, 20 landscapes, 120 datasets): the weakness reproduces
  (one-step RMSE 2.5 and 1.8 times that of no information at 6 days); with the
  new default it is 0.63 and 0.68 (paired −71% and −60%), −4% and −10% at 18
  days, unchanged at 42. Coverage at 6 days rose from 0.89–0.92 to 0.95–0.98.

## Cohort model

`fit_sparse_joint_cohort()` with `initial_state = "balance"` now chooses each
sample's history by a single-sample fit and fits the cohort once with them
(its default initial state stays `"free"`, where the history does not enter).

- Development cohort scenario (8 patients, 360 days, HPC): one-step RMSE on
  patient 1 within 0.96–1.005 of the earlier cohort fit on every landscape
  (mean 0.99); the same 2 of 12 cohorts fail with both, as before (no finite
  marginal likelihood for either trend). Held-out cohort scenario (final
  confirmation): within 0.98–1.09 per landscape (mean 1.01, Spearman
  0.003–0.008 lower), all 12 fitted. At long intervals the change is neutral.
- Cohorts of eight 6-day samples (`cohort_short.R`, wavelengths 0.8 and 1.6,
  first sample at day 120; patients 2–8 are new trajectories): one-step RMSE
  2.15 → 0.73 times that of no information (paired −58%, all 10 better),
  Spearman 0.69 → 0.81. The fit took 3.7 times as long.

## Running

```sh
R=validation/sparse_joint/short_interval/run.R
Rscript $R --stage=truth
Rscript $R --stage=datasets
Rscript $R --stage=fit --methods=zero,balance,balance_history_age,free
Rscript $R --stage=fit --methods=history_grid,trend_centre_low,trend_exponential,trend_select
Rscript $R --stage=score
Rscript validation/sparse_joint/short_interval/report.R
```

All methods of one `fit` call must share a build. Results go to `results/`
(not tracked), separate from the development benchmark's results.

The development benchmark copy (`devbench.R`) ran on the HPC, one Slurm task
per dataset and method (`hpc/`: `submit.sh build`, then `submit.sh fit`);
`history_sweep.R` compares fixed histories and selection rules, and
`calibration_fit.R` refits the uncertainty constants on that copy.
`cohort_short.R` builds and scores the short-interval cohorts
(`hpc/array_task_cohort.sh`). Held-out runs need `--unseal-heldout`
(`SPLIT=heldout` for `hpc/submit.sh`).
