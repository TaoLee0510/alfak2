# Sparse joint count model: design record

`fit_sparse_joint()` fits two sparse karyotype samples with one count
likelihood on a sparse domain. This note records why the model has its present
form and how it scores. Every number is from the development landscapes of
`benchmark/sparse_two_timepoint` (22 chromosomes, 12 landscapes, 8 scenarios,
96 datasets), with fitness as growth rate per day.

## Where the information is

Among observed karyotypes present in the population at both timepoints, the
change in log frequency over the interval correlates 0.11 with true fitness,
even at infinite sampling depth. A rare karyotype is held at
mutation-selection balance by the clone it descends from, so its frequency
follows its parent's and says nothing through its change. What carries its
fitness is its level relative to the parent: influx divided by fitness deficit.
The forward model therefore has to integrate mutational influx over the
interval. The two-layer pipeline applies one transport step after selection,
which underestimates that influx by a factor of 160–320 at these rates.

## Decisions

**Closed-form dynamics, first order in the flux.** A state's abundance at the
second timepoint is its own growth plus the influx from each neighbour,
integrated exactly (`R/forward-closed-form.R`). A state depends only on its
neighbours, so the Hessian stays sparse; any time-stepping scheme couples
everything within a few steps. `forward_study.R` compares it with the exact
solution from the true first-timepoint population: total-variation distance
0.01–0.05 at 180 days and 0.01–0.37 at 360 days. The error is the clones that
arise two or more mutations away within the interval.

**A free initial abundance for every state.** The exact mean dynamics
themselves, with true fitness, are 0.10–0.49 in total variation from what the
simulator produced after 180–360 days, and up to 1 after 1080 days. Which
new clone establishes is chance. A free initial abundance absorbs it; a more
exact solver would not.

**Baseline birth rate supplied.** Left free under a wide prior it settled at
3–16 divisions per day against a true 0.9, and inflated every fitness
difference. It is held within 0.05 (log scale) of a supplied value. See "The
two rates the model is given" below for what a wrong value does.

**Majors need a share and a count.** Domain = karyotypes holding at least 1%
of either sample and at least 3 cells, plus all their one-step neighbours.
With every observed karyotype as a major, a 50-cell sample put all 50 cells on
majors by construction, while their neighbours held 25% of the true
population; the likelihood read those structural zeros as strong unfitness
(bias of -0.21 per day at one step). At 500/700 cells the sample had 213 cells
on neighbours of majors against 206 expected, so there the evidence is sound.

**Shell term: distance to the nearest major, power 0.85.** Mean true fitness
relative to the reference, by lattice distance to the nearest major:

| | 0 | 1 | 2 | 3 |
| --- | ---: | ---: | ---: | ---: |
| Observed | -0.003 | -0.012 | -0.014 | -0.017 |
| Not observed | | -0.026 | -0.040 | -0.052 |

(Wavelengths 0.8–1.6, mean over the eight scenarios, which differ little.)

Unobserved karyotypes fall almost linearly with distance; a single shared
level for them left a bias of +0.05 at one step and +0.10 at two. The distance
is measured to majors only: measured to every observed karyotype, a
neighbour's own count decided its covariate and the coefficient ran to -0.35.

**Observed karyotypes outside the domain.** They are left out of the
likelihood, since they would enter only because observed while their siblings
would not. They are predicted at the level of the observed neighbours inside
the domain, which the table shows is nearly independent of distance.

**Priors on the frequency scale.** A karyotype never seen has prior frequency
1e-5 at the first timepoint; one seen only at the second is centred on its
frequency there. In counts, the same prior gave thousands of unseen neighbours
a large share of a small sample.

**Tried and dropped.** A separate residual scale for majors: no gain (RMSE
0.026 / 0.053 / 0.082 against 0.025 / 0.051 / 0.078 at distances 0 / 1 / 2 on
the typical scenario), slower, one failed fit.

## Uncertainty

**Relative fitness.** Two samples identify fitness differences. The model
reports fitness relative to the count-weighted mean of the majors, with the
uncertainty of that difference. The posterior of every log birth rate also
carries the uncertainty of the baseline (0.05), which is common to all
karyotypes and cancels in a difference. Intervals that included it covered
0.99 / 0.95 / 0.90 at distances 0 / 1 / 2; this was the baseline term, not
calibration.

**Expected information.** The covariance is the inverse of the Hessian with
the counts replaced by their fitted expectations. In 9 of the 88 development
fits examined the observed Hessian was not positive definite (the intensities are not
log-linear in the birth rates) and gave posterior standard deviations of a
major's log birth rate above 1, against a prior of 0.05.

**Calibration.** The model's own 95% intervals of relative fitness cover 0.76
of observed karyotypes with fewer than 5 cells, 0.72 of those with 5 or more,
0.57 of unobserved karyotypes at distance 1 and 0.49 at distance 2. The spread
of fitness among unobserved karyotypes cannot be learned from zero counts and
is set by the observed karyotypes, which selection has made alike. The
reported standard deviation is

    sqrt(model sd^2 + (0.6 * relative fitness)^2 + (0.009 * birth rate * (1 + shell))^2)

with `shell` the distance to the nearest major. The two constants are a
Gaussian-likelihood fit on the 96 development datasets
(`calibration_fit.R`: 0.0088 and 0.617). Coverage becomes 0.97 / 0.97 / 0.94 /
0.91 for the four kinds above. By scenario and distance it lies between 0.89
and 1.00, except after 1080 days (0.92 / 0.88 / 0.84). A free multiplier on
the model's standard deviation fitted to zero once the other two terms were
present, so it was fixed at one.

## The two rates the model is given

Typical scenario, 12 datasets; the simulated populations divide 0.9 times per
day. RMSE and mean prediction relative to the reference, at distance 1:

| Assumed birth rate | RMSE | Mean prediction per day | Per division | Spearman | Coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 0.085 | -0.009 | -0.092 | 0.43 | 0.22 |
| 0.25 | 0.076 | -0.019 | -0.077 | 0.44 | 0.45 |
| 0.5 | 0.064 | -0.034 | -0.068 | 0.45 | 0.71 |
| 1 | 0.050 | -0.061 | -0.061 | 0.45 | 0.92 |
| 2 | 0.056 | -0.102 | -0.051 | 0.46 | 0.98 |
| 4 | 0.132 | -0.188 | -0.047 | 0.45 | 0.93 |
| No information | 0.094 | 0 | | | |

Over all eight core scenarios (96 datasets), RMSE at distances 0 / 1 / 2:

| Assumed rates | RMSE | Spearman |
| --- | --- | --- |
| Both right | 0.033 / 0.055 / 0.079 | 0.44 / 0.43 / 0.36 |
| Birth rate 1/4 of the truth | 0.036 / 0.072 / 0.110 | 0.44 / 0.43 / 0.36 |
| Birth rate 4 times the truth | 0.115 / 0.185 / 0.237 | 0.43 / 0.42 / 0.35 |
| Missegregation rate 1/5 | 0.040 / 0.076 / 0.113 | 0.43 / 0.43 / 0.35 |
| Missegregation rate 5 times | 0.105 / 0.142 / 0.156 | 0.44 / 0.42 / 0.35 |
| No information | 0.046 / 0.092 / 0.137 | |

- Rank agreement with truth does not depend on either rate.
- The deficit of unobserved karyotypes per division changes by a factor of 2
  over a 40-fold range of birth rates; per day it follows the rate.
- Too low a rate is safe for the error (predictions shrink towards no
  differences) but its intervals are too narrow. Four times too high is worse
  than no information.
- The missegregation rate acts the same way, and only the product matters:
  doubling one and halving the other gives RMSE 0.048 and 0.049 and mean
  predictions -0.058 and -0.059, against 0.050 and -0.061 for the reference.
- Neither the marginal likelihood of one sample nor that of a cohort locates
  the rate. For single samples over 0.1–4 it is best at 0.1 in 9 of 12
  datasets and at 0.5 in 3. For three cohorts of eight over 0.25–4 it is best
  at 0.25 in two and flat (within 1.1) in the third.

`fit_sparse_joint_rates()` fits at the two ends and the centre of a range and
reports an interval covering all three. Its point prediction is taken at the
lowest rate. Typical scenario, a 16-fold range, RMSE at distances 0 / 1 / 2:

| Truth in the range | Point at the lowest rate | Point at the centre | Coverage |
| --- | --- | --- | --- |
| At the centre (0.25–4) | 0.038 / 0.076 / 0.115 | 0.030 / 0.053 / 0.080 | 1.00 / 1.00 / 1.00 |
| At the top (0.0625–1) | 0.044 / 0.087 / 0.131 | 0.041 / 0.082 / 0.124 | 0.98 / 0.94 / 0.90 |
| At the bottom (1–16) | 0.028 / 0.050 / 0.075 | 0.093 / 0.155 / 0.204 | 1.00 / 1.00 / 1.00 |

At the lowest rate the error stays below no information (0.048 / 0.094 /
0.139) wherever the truth lies in the range; at the centre it does not. Rank
correlation between the predictions at the two ends of a range is 0.71–1.00
(median 0.87).

## How far from the limit

`depth_limit.R` refits the typical scenario on samples of 20 000 cells per
timepoint drawn from the same populations, and scores a perfect stratum
mean: every karyotype predicted at the true mean of its distance stratum.

| Distance | Benchmark depth | 20 000 cells | True stratum mean |
| --- | --- | --- | --- |
| 0 | 0.028 (0.49) | 0.022 (0.60) | 0.025 |
| 1 | 0.050 (0.45) | 0.043 (0.50) | 0.044 |
| 2 | 0.075 (0.38) | 0.063 (0.43) | 0.059 |

Forty times more cells lower the error by 14–21%. With the whole population
the model is as accurate as knowing the true mean of each stratum, and ranks
karyotypes within it besides. What two timepoints say about individual
unobserved karyotypes is limited by the mean dynamics and by which clones
happened to establish, not by sampling.

## Initial abundances at balance (option)

`initial_state = "balance"` starts every non-major karyotype at
mutation-selection balance with the majors next to it, so that the first
sample's levels inform fitness. Its spread around balance (`balance_sd`)
hardly matters: 0.3, 0.5 and 1 give RMSE within 0.003 of each other on the
benchmark. It stays at 1; tighter values moved the fitted baseline to 1.2
times the supplied one.

The balance prior is not log-linear either, and with the observed Hessian 4 of
96 development fits had no positive-definite information. Its curvature is now
taken at zero residual, as for the counts; all 168 development and all 168
held-out fits then finish.

It is not the default. It assumes that fitness before the first sample was
what it is between the samples, and its history (by default the interval
between the samples) is too short when the samples are close: at 30 days it
predicts worse than the default model (below). It takes about 1.5 times as
long.

## Development scoreboard

Pooled over the eight core scenarios, 96 datasets; RMSE, with Spearman
correlation in brackets, then coverage of the reported 95% intervals.

| Distance | No information | Two-layer | ALFA-K (minobs 20) | Sparse joint | With balance |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.046 | 0.045 (0.30) | 1.073 (0.23) | 0.033 (0.44), 0.97 | 0.026 (0.47), 0.98 |
| 1 | 0.092 | 0.091 (0.24) | 1.186 (0.20) | 0.055 (0.43), 0.94 | 0.044 (0.45), 0.97 |
| 2 | 0.137 | 0.136 (0.15) | 1.297 (0.18) | 0.079 (0.36), 0.91 | 0.066 (0.38), 0.96 |

All 168 fits of each variant (core and stress scenarios) finished; median
9 s for the default and 14 s with balance, maximum 238 s.

Acceptance: the default model passes 17 of 18 criteria. RMSE is below the
no-information answer and below the two-layer pipeline at every distance, and
Spearman above both existing methods, all with landscape-bootstrap intervals
excluding zero. The safety line fails in one cell: scenario
`unbalanced_3000_100_360d` at distance 0, RMSE 17% above the two-layer
pipeline. That cell is two near-monoclonal datasets whose information is in
the deep first sample, which the default model does not read. With balance
all 18 pass; its worst cell is 0.69 of the two-layer pipeline.

## Held-out confirmation

Twelve landscapes simulated with seeds sealed until the methods were frozen
(commit 61ef06e), eight core scenarios, 96 datasets, scored once.

| Distance | No information | Two-layer | ALFA-K (minobs 20) | Sparse joint | With balance |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.044 | 0.043 (0.22) | 0.277 (0.16) | 0.032 (0.41), 0.97 | 0.025 (0.44), 0.98 |
| 1 | 0.090 | 0.090 (0.16) | 0.371 (0.13) | 0.058 (0.36), 0.95 | 0.046 (0.37), 0.96 |
| 2 | 0.136 | 0.135 (0.09) | 0.454 (0.12) | 0.082 (0.30), 0.93 | 0.068 (0.31), 0.94 |

- With balance, all 15 criteria that could be evaluated pass (ALFA-K with
  minobs 5 was not run on these landscapes): RMSE 43%, 49% and 50% below no
  information, with every interval excluding zero; worst scenario 0.90 of the
  two-layer pipeline. All 168 fits finished.
- The default model passes 12 of 15. The safety line fails at distances 0 and
  1 (1.38 and 1.14) in `unbalanced_3000_100_360d`, as on the development
  landscapes, from three datasets in which most observed karyotypes, seen in
  the deep first sample only, are predicted at the trend (bias -0.12 to
  -0.23), and at distance 0 in `very_sparse_50_50_180d` (1.33), from one
  dataset with three observed karyotypes.
- ALFA-K fails in 5 of 96 fits and its error is 3–6 times that of no
  information.
- The calibration fitted on the development landscapes holds: coverage of
  observed karyotypes with fewer than 5 cells, with 5 or more, and of
  unobserved ones at distance 1 and 2 is 0.97 / 0.98 / 0.95 / 0.92 for the
  default model and 0.98 / 0.98 / 0.96 / 0.94 with balance. The model's own
  intervals cover 0.50–0.79.

## Stress scenarios

RMSE at distances 0 / 1 / 2 on the development landscapes (twelve datasets
each); the held-out landscapes agree.

| Scenario | No information | Sparse joint | With balance | Coverage, default |
| --- | --- | --- | --- | --- |
| 30 days between samples (27 generations) | 0.049 / 0.096 / 0.141 | 0.028 / 0.048 / 0.070 | 0.033 / 0.066 / 0.097 | 1.00 / 0.98 / 0.96 |
| 60 days (54 generations) | 0.050 / 0.096 / 0.142 | 0.026 / 0.046 / 0.069 | 0.027 / 0.047 / 0.071 | 0.99 / 0.97 / 0.95 |
| Half of the landscape redrawn at the first sample, 500 / 700 cells | 0.074 / 0.113 / 0.156 | 0.058 / 0.084 / 0.116 | 0.053 / 0.071 / 0.096 | 0.73 / 0.73 / 0.72 |
| The same, 3000 / 100 cells | 0.070 / 0.110 / 0.147 | 0.060 / 0.083 / 0.108 | 0.050 / 0.066 / 0.082 | 0.82 / 0.80 / 0.79 |
| Landscape turned about the dominant clone, 500 / 700 | 0.082 / 0.128 / 0.172 | 0.059 / 0.079 / 0.105 | 0.058 / 0.074 / 0.097 | 0.87 / 0.83 / 0.80 |
| The same, 3000 / 100 | 0.093 / 0.138 / 0.180 | 0.088 / 0.113 / 0.134 | 0.067 / 0.079 / 0.095 | 0.85 / 0.88 / 0.88 |

- Few generations between the samples do not hurt the default model; rank
  agreement falls (Spearman 0.32 / 0.29 / 0.25 at 30 days). Balance is worse
  at 30 days (held out: 0.069 against 0.049 at distance 1), because its
  history then is 30 days while the cloud around a clone took far longer to
  build. At 60 days the two are equal.
- A landscape that changes at the first sample, as when treatment starts then,
  leaves both better than no information, with intervals that are too narrow
  (0.72–0.88 for the default model; 0.82–0.94 with balance).
- Balance was expected to fail there and does not: it is at least as accurate
  in all four, including the scenario built against it, where which neighbour
  is fit before says nothing about after. Its Spearman at distance 1–2 is lower
  in the 3000 / 100 scenarios (0.50 / 0.45 against 0.58 / 0.54 when turned).
  The gain is in the level of the deficits, not in which neighbour is which.

## Forward prediction

`forward_prediction.R`: total-variation distance between the predicted and
the true composition of the domain 180 days after the second sample, default
model.

| Forecast | Development (72 datasets) | Held out (72) |
| --- | ---: | ---: |
| Composition unchanged | 0.317 | 0.290 |
| Run forward without selection | 0.469 | 0.456 |
| Run forward with the fitted birth rates | 0.275 | 0.251 |
| Run forward with the true rates | 0.164 | 0.146 |

Fitted against unchanged: -0.042 (landscape bootstrap 95% interval -0.079 to
-0.005) on the development landscapes and -0.039 (-0.081 to -0.005) held out.
The fit closes about a quarter of the gap to what the true fitness would give
with this forward model and starting point.

## Real samples

`real_cohort.R` on 20 GBM samples (22 chromosomes, 8–3889 cells per timepoint,
179–3599 days apart), birth rates 0.02–1 per day, missegregation rate 5e-5.
There is no truth; these are checks of behaviour.

- All 20 samples fit at all three rates; in 17 the scale optimiser converged
  at all three. Median 149 s for the three fits, maximum 679 s; median 849
  states, maximum 1558.
- Over the 50-fold range of rates the mean deficit of non-major karyotypes per
  day changes 47-fold (median) and per division 0.94-fold. Predictions at the
  two ends have rank correlation 0.66 / 0.76 / 0.96 (quartiles); one sample is
  at -0.2.
- Two halves of the cells of a sample (10 samples with at least 200 cells at
  both timepoints) agree with Spearman 0.88 for karyotypes with 10 or more
  cells, 0.76 with 2–9 cells and 0.33 with one cell. All differences are
  within the reported intervals, and 98–100% within the model's own, which
  describe sampling but not model error.
- Pearson dispersion of the second sample's counts around the fit: quartiles
  0.04 / 0.27 / 1.2; one sample at 15.
- The samples overlap: a median 38% of a sample's karyotypes are seen in
  another sample and 69% lie in another sample's domain. The recurrent ones
  are chromosome 10 loss (17 samples), chromosome 7 gain (15) and both (14).
  Simulated cohorts overlap less (0–38% in 9 of 12 landscapes).

## Cohort

`fit_sparse_joint_cohort()` shares the fitness of karyotypes between samples,
with a between-sample deviation whose scale is estimated. Scored on patient 1
of each landscape in the typical scenario, fitted with patients 2–8 of the same
landscape (cohort) or with patient 1 of seven other landscapes (control).
RMSE, with Spearman correlation in brackets, at distances 0 / 1 / 2:

| | Development (12 landscapes) | Held out (12) |
| --- | --- | --- |
| Patient alone | 0.028 / 0.050 / 0.075 (0.49 / 0.45 / 0.38) | 0.025 / 0.049 / 0.075 (0.40 / 0.38 / 0.32) |
| Cohort, trend chosen | 0.025 / 0.050 / 0.078 (0.56 / 0.50 / 0.43) | 0.024 / 0.047 / 0.074 (0.47 / 0.44 / 0.39) |
| Cohort, one trend | 0.025 / 0.050 / 0.077 (0.55 / 0.53 / 0.46) | |
| Control, trend chosen | 0.028 / 0.056 / 0.086 (0.47 / 0.43 / 0.35) | |
| Control, one trend | 0.032 / 0.069 / 0.107 (0.39 / 0.29 / 0.23) | |

- A cohort of one landscape improves rank agreement at every distance, with
  landscape-bootstrap intervals excluding zero on both splits (held out:
  +0.07 / +0.06 / +0.07). It does not lower the error, and the bias of
  unobserved karyotypes rises slightly (+0.01). The planned criterion, lower
  error at distances 1 and 2, is not met.
- One trend for samples of different landscapes is harmful: RMSE +0.019 and
  +0.032 at distances 1 and 2, Spearman -0.16, both with intervals excluding
  zero. The trend's dosage and curvature then describe no sample.
- The trend is therefore chosen by marginal likelihood between one trend and
  one per sample. It chose separate trends for 11 of 12 control cohorts, and
  for cohorts of one landscape the shared trend in 7 of 12 development and 8
  of 12 held-out landscapes; the others are the rugged landscapes (wavelength
  1.2–1.6), whose patients sit on different peaks. The control then costs
  Spearman -0.02 to -0.03 and no significant error, and its intervals cover
  less (0.96 / 0.88 / 0.83 against 0.98 / 0.94 / 0.90), because the shared
  scales come from other landscapes.
- With balance, a patient alone does as well on error (development 0.026 /
  0.044 / 0.066 over all scenarios). Balance within a cohort has not been
  tried.
- A cohort of eight takes 10 s to 45 min, twice that when the trend is chosen.

On the 20 GBM samples the per-sample trend is strongly preferred (negative log
marginal likelihood 11 470 against 11 893), and the estimated deviation of a
karyotype's fitness between patients is large (0.05 on the log scale, against
a residual of 0.003). The patients do not look like one landscape. Within a
sample, estimates under the two trends agree (Spearman 0.94, median over
samples). The fit takes about an hour for 16 000 states.

## Open

- The default model misses the information in a deep first sample (the
  failing safety cell); balance reads it but needs a history at least as long
  as the cloud took to build, which the interval between close samples is not.
- The unknown birth rate has no remedy beyond the range fit. Intervals are too
  narrow when the assumed rate is too low, because both calibration terms
  scale with it.
- Intervals are too narrow after 1080 days and when the landscape changes at
  the first sample.
- The shell coefficient is a lower bound from absence, so unobserved
  neighbours are still predicted too high (bias +0.02 at one step, +0.04 at
  two).
- The calibration constants come from one family of simulated landscapes.
- ALFA-K with minobs 5 has 75 of 96 development fits (the HPC run continues)
  and none on the held-out landscapes. On those 75 its RMSE is 0.169 / 0.225
  at distances 0 / 1, above no information, and both variants rank better
  with intervals excluding zero.
- Balance within a cohort, and a longer default history for balance, have not
  been tried.

## Files

`forward_study.R` reproduces the forward-model comparison. `compare.R` and
`by_scenario.R` tabulate a method's fits against the baselines per stratum and
per scenario. `calibration.R` and `calibration_fit.R` examine and fit the
reported uncertainty. `rate_sensitivity.R` sweeps the two supplied rates.
`forward_prediction.R` scores forecasts of the later composition.
`depth_limit.R` compares the fit with deep resamples and a perfect stratum
mean. `real_cohort.R` runs the checks on real samples. The official verdict comes
from `benchmark/sparse_two_timepoint/run.R --stage=report`.
