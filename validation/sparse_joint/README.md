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

## Initial abundances at balance (the default)

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

It became the default on 2026-11-11, after the held-out confirmation, where
it passed every criterion and the free model did not. It assumes that fitness
before the first sample resembled fitness between the samples. On simulated
landscapes that changed at the first sample it was still the more accurate
of the two (stress scenarios, below), which matters for patients treated
between the samples. Its history (by default the interval between the
samples) is too short when the samples are close: at 30 days, about 27
divisions, it predicts worse than the free model, at 60 days as well. It takes
about 1.5 times as long. `initial_state = "free"` gives the earlier model; the
benchmark entries without "balance" in their name use it.

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
  0.044 / 0.066 over all scenarios).
- Balance within a cohort (`alfak2_sparse_joint_cohort_balance`, trend
  chosen, development, 2026-11-11): 2 of the 12 cohorts (L04, L07) fail
  because neither trend has a finite marginal likelihood; the scale
  optimisation fails already at the prior-centred scales. On the other 10,
  against the cohort with free initial abundances, RMSE is lower at distances
  1 and 2 (-0.003, interval -0.005 to -0.001; -0.006, -0.009 to -0.003) and
  Spearman is not significantly different (+0.015); against balance alone,
  RMSE is not significantly different and Spearman is higher (+0.078 / +0.053
  / +0.045, intervals above zero). Coverage 0.99 / 0.96 / 0.94. A cohort of
  eight took up to 57 minutes. The cohort default stays free.
- A cohort of eight takes 10 s to 45 min, twice that when the trend is chosen.

On the 20 GBM samples the per-sample trend is strongly preferred (negative log
marginal likelihood 11 470 against 11 893), and the estimated deviation of a
karyotype's fitness between patients is large (0.05 on the log scale, against
a residual of 0.003). The patients do not look like one landscape. Within a
sample, estimates under the two trends agree (Spearman 0.94, median over
samples). The fit takes about an hour for 16 000 states.

## The second-layer methods of the two-layer pipeline

The two-layer pipeline's nine second-layer methods were compared on the 96
development datasets of the core scenarios (2026-11-11), each with the same
first layer and fixed graph hyperparameters as `alfak2_two_layer_fixed`
(`graph_gaussian_baseline`). `tabpfn_nearfield_feature_model` runs
`delta_tree_ensemble` without its Python module and was not run separately.
All 672 fits finished. By wavelength (`wavelength_tables.R`, as in Supplement
S2.3 of `doc/alfak2_method.tex`):

- Error: every second-layer method is 0.89 to 1.00 times the error of the null
  predictor (`delta_tree_ensemble` 0.89 to 1.19), against 0.44 to 0.75 for the
  sparse joint model with balance. They all interpolate the same first-layer
  anchors, so the choice of second layer hardly matters.
- Ranking: on the smoothest landscapes the second-layer methods rank observed
  karyotypes better (wavelength 1.6, distance 0: 0.48 to 0.54 against 0.41)
  and as well at wavelength 1.2 (0.47 to 0.50 against 0.49); everywhere else
  the sparse joint model ranks better. The held-out landscapes show the same
  exception (0.48 against 0.34 for `graph_gaussian_baseline`).
- Intervals: near zero coverage for `local_polynomial_stencil` and
  `delta_tree_ensemble` at distances 1 and 2, near one for the graph Gaussian
  methods and `local_NNGP_or_GPnn`, below 0.6 for observed karyotypes except
  the graph Gaussian methods; 0.90 to 1.00 for the sparse joint model.

The benchmark's score stage now replaces the rows of the methods it scores and
keeps the others; it previously rewrote `metrics_<split>.csv` with only those
methods.

## Results by wavelength

Pooled over landscapes, errors are dominated by the rugged landscapes
(wavelength 0.4), where fitness differences are about ten times larger than at
1.6. The method document therefore reports results per wavelength. On the
held-out landscapes the sparse joint model with balance is below the null
predictor at every wavelength and distance (0.47 to 0.74 of its error), the
two-layer pipeline at 0.92 to 1.00, and `alfakR` above it everywhere (21% of
its fits failed at wavelength 1.6). Coverage of the sparse joint model ranges
from 0.89 (wavelength 0.8, distance 2) to 1.00 (wavelength 1.6): slightly
narrow on intermediate, conservative on smooth landscapes. The pooled
acceptance rule is unchanged.

The one ranking exception, observed karyotypes on the smoothest landscapes
(wavelength 1.6, distance 0), where the two-layer pipeline ranks better (0.48
against 0.34 held-out), comes from karyotypes seen in fewer than five cells.
The two-layer pipeline does not filter by count. For karyotypes with at least
five cells the sparse joint model ranks better at every wavelength (0.74
against 0.68 at 1.6). At wavelength 1.6, 56% of the rare observed karyotypes
(held-out; 51% development) lie more than one step from every major and are
predicted at the common observed level, which cannot rank them; the two-layer
pipeline smooths them from their observed neighbours (Spearman 0.44). Inside
the domain the sparse joint model ranks rare karyotypes better (0.65 against
0.53).

## Results by sampling design

`scenario_results.R` resolves the held-out results by sampling design and
wavelength (Section 5.3 and Supplement S4 of `doc/alfak2_method.tex`; figures
and tables in `doc/generated/`). Each cell holds three landscapes.

- Sparsity, as cells per observed karyotype, ranges from 1.4 (50/50 cells,
  wavelength 1.6) to 76 (3000/100, wavelength 0.4).
- The sparse joint model with balance is below the null predictor in 95 of 96
  combinations of wavelength, design and distance. The exception is 50/50
  cells at wavelength 0.4 and distance 0 (ratio 1.03), where three to seven
  karyotypes are observed and one landscape with three had twice the null
  error. In the 50/50 design it is otherwise 0.53 to 0.84; in all other
  designs 0.41 to 0.83. On the development landscapes all 96 are below (0.39
  to 0.88).
- Most of the gain comes between 50 and 200 cells per sample; 500/700 adds
  little. A 1080-day interval is worse than 360 days at wavelengths 0.8 and
  1.2 (0.62 to 0.75 against 0.49 to 0.61). 3000/100 is the weakest design at
  wavelength 1.6 (0.83 at distance 1).
- The two-layer pipeline is at 0.85 to 1.08 everywhere; `alfakR` is below the
  null predictor in 4 of 96 combinations and up to 300 times above it (the
  error figure caps its ratios at 8 and labels the larger ones). The sparse
  joint model has a lower error than `alfakR` in all 96 and ranks better in
  87 of the 93 with an `alfakR` ranking; the exceptions are observed
  karyotypes in four designs at wavelength 1.6 and two small ones at 0.4.
- Rank correlation at distances 1 and 2 grows with depth on smooth landscapes
  (wavelength 1.6, distance 1: 0.28 at 50/50, 0.55 to 0.62 at 500/700) and
  stays at 0.12 to 0.33 on rugged ones. The two-layer pipeline ranks observed
  karyotypes better in 7 of 8 designs at wavelength 1.6 and 2 of 8 at 1.2.
- Coverage 0.77 to 1.00; below 0.90 mostly at distance 2 with the long
  interval (0.81 to 0.83) or a shallow first sample (100/3000: 0.77 to 0.86),
  at wavelengths up to 0.8.

## Distances 3 and 4

`far_distance.R` samples, for every core dataset of both sets, 2000
karyotypes at distance 3 and 2000 at distance 4 from the nearest observed
karyotype (`stt_far_nodes()` in `benchmark/sparse_two_timepoint/score.R`) and
fits all twelve methods again with them as targets: the null predictor, both
sparse joint variants, the two-layer pipeline with each of its eight second
layers, and `alfakR` (minobs 20). The shells cannot be enumerated at 22
chromosomes, so they are sampled by random walks from observed karyotypes,
each observed karyotype an equally likely start, keeping endpoints at exactly
that distance; this is not uniform within a shell, unlike the scoreboard's
distances 1 and 2. The walks' intermediate nodes are unscored targets, so that
the graph second layers have a path. Centring is the scoreboard's. Results,
fits and the sparse joint fit objects are under
`benchmark/results/sparse_two_timepoint/sweeps/far/`; `far_results.R` writes
the figures and tables of Section 5.4 and Supplement S5 to `doc/generated/`.
The two-layer first layer and target graph are computed once per dataset and
reused by its eight entries (`--memo=yes`); a check on two datasets gave
identical predictions for all but `alfakR`, whose bootstrap is unseeded
(correlation 0.996). The run took 3 h 20 min with three workers on a laptop.

- Every method fitted every dataset except `alfakR`, which failed on the same
  7 development and 5 held-out datasets as in the scoreboard.
- Every target at distance 3 or 4 is at least three steps from a major
  karyotype, outside the domain, so the sparse joint prediction is the trend
  alone.
- Held-out, balance: ratio to the null predictor 0.45 / 0.58 / 0.55 / 0.54 at
  distance 3 and 0.44 / 0.59 / 0.56 / 0.53 at distance 4 (wavelengths 0.4 /
  0.8 / 1.2 / 1.6), within 0.03 of distance 2. The gain is mostly the mean:
  the null predictor's mean error at distance 4 is 0.474 / 0.157 / 0.085 /
  0.049, the trend leaves 0.128 / 0.071 / 0.034 / 0.010. The remaining mean
  error is positive everywhere and grows from distance 3 to 4: distant
  karyotypes are predicted too high.
- Spearman at distance 4: 0.14 / 0.16 / 0.33 / 0.40 (distance 2: 0.17 / 0.21 /
  0.39 / 0.46). Coverage at distances 3 and 4: 0.91, 0.84 to 0.86, 0.97 to
  0.98 and 1.00, so too narrow at wavelength 0.8.
- Independent initial states: higher error than balance everywhere (0.48 to
  0.82), similar ranking (0.11 to 0.43).
- Graph Gaussian: 0.98 to 1.00 of the null error, mean error within 0.001 of
  the null predictor, Spearman -0.04 to 0.28, coverage 1.00 because its
  standard deviation grows to 0.6 to 1.2. The seven near-field second layers
  predict nothing beyond distance 2 and score as the null predictor.
- `alfakR` (Kriging mean, no uncertainty): 2.35 / 1.32 / 1.50 / 18.2 at
  distance 3, 2.28 / 1.26 / 1.42 / 15.3 at distance 4; Spearman 0.01 to 0.25.
- Development landscapes: balance 0.42 to 0.59, Spearman 0.14 to 0.49.

## Three measures of accuracy

The method document now states two design goals, numerical and direction
accuracy, and reports them together with the overall similarity of predicted
and true fitness (Section 2.1 and 5.1 of `doc/alfak2_method.tex`):

- Numerical accuracy: RMSE relative to the null predictor (null = 1).
- Direction accuracy: over pairs of scored karyotypes one copy of one
  chromosome apart, the share whose predicted fitness change has the true
  sign; a tie (equal predictions) counts one half, pairs with equal true
  fitness are skipped, and a pair belongs to the distance of its farther
  member (null = random = 1/2). `direction_accuracy.R` computes it from the
  saved fits (`sweeps/direction/direction_<split>.csv`); at distances 3 and 4
  the pairs include the intermediate nodes of the sampling walks. The score
  stage's `edge_sign_accuracy` counts ties as wrong, so the null predictor and
  failed fits score 0 there; it is not used. Direction accuracy was added
  after the held-out evaluation, from the same fits, and was not an acceptance
  criterion.
- Overall similarity: Spearman correlation (null = 0), not a direction
  measure.

Held-out, means of 24 datasets per wavelength and distance:

| | error ÷ null, d0–2 | direction, d0–2 | direction, d3–4 | direction above 1/2, by design |
|---|---|---|---|---|
| sparse joint, balance | 0.47–0.74 | 0.64–0.96 | 0.76–0.81 | 96 / 96 |
| sparse joint, independent | 0.53–0.96 | 0.66–0.96 | 0.76–0.81 | 96 / 96 |
| two-layer (graph Gaussian) | 0.92–1.00 | 0.49–0.57 | 0.49–0.52 | 62 / 96 |
| alfakR (minobs 20) | 1.39–41.6 | 0.36–0.52 | 0.45–0.53 | 28 / 96 |

- alfakR is below one half in 10 of 12 wavelength × distance cells at
  distances 0–2: it predicts the wrong direction more often than the right
  one. Its failed fits score one half (all ties).
- The default's lowest value, 0.64 for observed karyotypes at wavelength 1.6,
  comes from rare observed karyotypes outside the domain: they share one log
  birth rate, so their fitness differs only through the missegregation term.
- By design at distance 1 the default is 0.76 to 0.90 everywhere; unlike the
  error, direction hardly changes with sample size.
- Development: default 0.67 to 0.93 (d0–2), above one half in 96 / 96.

`accuracy_by_wavelength.R` draws the three measures by wavelength (Section
5.2); `scenario_results.R` and `far_results.R` add the direction figures and
tables by design and at distances 3 and 4.

## Outputs for alfakR pipelines

Analyses built on `alfakR::alfak()` (branch `dev`) read four files from its
output directory. `alfak_compat()` writes them from a sparse joint fit, with
the same arguments; `as_alfakR()` builds them from an existing fit. It is not
called `alfak()`: an earlier decision keeps that name out of alfak2, and it
would mask `alfakR::alfak()` wherever both packages are loaded.

- Karyotypes are chosen as in alfakR: frequent (`minobs` cells over both
  samples), their one-step neighbours, and those neighbours' neighbours, in
  alfakR's order. `minobs` does not change the fit.
- alfakR's bootstrap replicates become joint posterior draws. Each karyotype's
  draws are widened by independent noise to its calibrated standard
  deviation, so differences between karyotypes keep the model's correlations.
  Draws come in mirrored pairs and each row is scaled to that standard
  deviation, so `landscape.Rds` is both the summary of the draws, as in
  alfakR, and exactly the model's prediction and calibrated uncertainty.
- The scale is alfakR's: the count-weighted mean of the majors is set to
  `log(nb / n0) / dt`.
- The cross-validation keeps alfakR's folds, each frequent karyotype with its
  neighbours. A fold's counts leave the likelihood, its karyotypes get the
  prior of an unseen karyotype and stop being reference points of the trend,
  and the posterior mode is found again with the fitted scales. Karyotypes
  outside the domain, whose fitness never came from their own counts, are not
  tested, and a fold holding every major is skipped.

Checks against alfakR dev (commit eddeb5e), run on the same three benchmark
datasets: every file has the same classes, names, columns, dimensions and
attributes, the landscape has the same karyotypes in the same order (868,
13 157 and 20 126 rows), and the frequent karyotypes and neighbours agree. The
one addition is `tmp` in `xval.Rds`, the paired matrix of earlier alfakR
versions, which the PANcanKFLs scripts read before recomputing a Kriging
cross-validation. On two GBM samples the PANcanKFLs reader
`make_parent_neighbor_matrix()` runs on the output unchanged.

The cross-validation R-squared measures something else than alfakR's. On a
development dataset with 27 frequent karyotypes (`dev_L12_P01_unbalanced_100_3000_360d`)
the full fit correlates 0.44 with the true fitness of the tested karyotypes,
the held-out predictions 0.07, and `R2R` is -0.78 (alfakR: -0.18): a karyotype
whose own and neighbours' counts are removed is barely predictable from the
rest of the landscape. On `dev_L10_P01_unbalanced_3000_100_360d` it is 0.43
(alfakR: -0.04). The cross-validation took 75 s for 27 folds.

## Open

- The default model misses the information in a deep first sample (the
  failing safety cell); balance reads it but needs a history at least as long
  as the cloud took to build, which the interval between close samples is not.
- The unknown birth rate has no remedy beyond the range fit. Intervals are too
  narrow when the assumed rate is too low, because both calibration terms
  scale with it.
- Intervals are too narrow after 1080 days and when the landscape changes at
  the first sample.
- Observed karyotypes outside the domain all get the same observed level, so
  on diverse populations (smooth landscapes) about half the rare observed
  karyotypes are not ranked. Adding the trend's variation to that level is
  untried.
- The shell coefficient is a lower bound from absence, so unobserved
  neighbours are still predicted too high (bias +0.02 at one step, +0.04 at
  two, and on held-out landscapes +0.005 to +0.13 at three and four, largest
  on rugged ones). Far intervals are too narrow at wavelength 0.8 (coverage
  0.84 at distance 4).
- The calibration constants come from one family of simulated landscapes.
- ALFA-K with minobs 5 has 75 of 96 development fits (the HPC run continues)
  and none on the held-out landscapes. On those 75 its RMSE is 0.169 / 0.225
  at distances 0 / 1, above no information, and both variants rank better
  with intervals excluding zero.
- Balance within a cohort helps where it runs but fails on 2 of 12
  development cohorts; falling back to free initial abundances when it fails,
  or a marginal likelihood from the expected information, is untried. A
  longer default history for balance has not been tried either.

## Files

`forward_study.R` reproduces the forward-model comparison. `compare.R` and
`by_scenario.R` tabulate a method's fits against the baselines per stratum and
per scenario. `calibration.R` and `calibration_fit.R` examine and fit the
reported uncertainty. `rate_sensitivity.R` sweeps the two supplied rates.
`forward_prediction.R` scores forecasts of the later composition.
`depth_limit.R` compares the fit with deep resamples and a perfect stratum
mean. `real_cohort.R` runs the checks on real samples. `alfakR_compat.R`
compares the output files with those of alfakR. `wavelength_tables.R` writes
the tables by wavelength of the method document. `scenario_results.R` writes the figures and
tables by sampling design to `doc/generated/`. `far_distance.R` fits and
scores all methods at distances 3 and 4, and `far_results.R` writes its
figures and tables to `doc/generated/`. `direction_accuracy.R` computes the
direction accuracy of the main methods from the saved fits, and
`accuracy_by_wavelength.R` draws the three accuracy measures by wavelength.
The official verdict comes
from `benchmark/sparse_two_timepoint/run.R --stage=report`.
