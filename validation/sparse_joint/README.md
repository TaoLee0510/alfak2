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
difference. It is now held within 0.05 (log scale) of a supplied value. On the
typical scenario, RMSE at distances 0 / 1 / 2 for a supplied rate of:

| Supplied | RMSE | Spearman | Coverage |
| --- | --- | --- | --- |
| 0.5 per day | 0.032 / 0.064 / 0.098 | 0.48 / 0.45 / 0.38 | 0.87 / 0.74 / 0.62 |
| 1 per day | 0.028 / 0.050 / 0.075 | 0.49 / 0.45 / 0.38 | 0.97 / 0.92 / 0.86 |
| 2 per day | 0.035 / 0.056 / 0.074 | 0.50 / 0.46 / 0.39 | 1.00 / 0.98 / 0.97 |
| No information | 0.048 / 0.094 / 0.139 | | |

Rank agreement does not depend on the supplied rate; error and interval width
do. Too low a rate narrows the intervals. In one of the twelve fits at 2 per
day the scale optimiser failed and the fit used the prior-centred scales.

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

## Development scoreboard

Pooled over scenarios; RMSE, with Spearman correlation in brackets.

| Distance | No information | Two-layer | ALFA-K (minobs 20) | Sparse joint |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.046 | 0.045 (0.30) | 1.073 (0.23) | 0.033 (0.44) |
| 1 | 0.092 | 0.091 (0.24) | 1.186 (0.20) | 0.055 (0.43) |
| 2 | 0.137 | 0.136 (0.15) | 1.297 (0.18) | 0.079 (0.36) |

Coverage of 95% intervals is 0.99 / 0.95 / 0.90. All 96 fits finished with
finite uncertainty, 83 with a converged scale optimiser; median 6 s, maximum
111 s, at most 1538 states.

Acceptance: 17 of 18 criteria pass. RMSE is below the no-information answer
and below the two-layer pipeline at every distance, and Spearman above both
existing methods, all with landscape-bootstrap intervals excluding zero. The
safety line fails in one cell: scenario `unbalanced_3000_100_360d` at distance
0, RMSE 17% above the two-layer pipeline.

That cell is two near-monoclonal datasets (wavelength 0.4). The deep first
sample shows the cloud around the dominant clone, and the observed cloud
members are its least unfit ones; the 100-cell second sample cannot tell them
from the rest. The model predicts them at the level of the unobserved cloud,
0.2 per day too low, because it leaves initial abundances free and so does not
read the first sample's levels. The mirror scenario, with the deep sample
second, is among the best (RMSE 0.56 of the two-layer pipeline).

## Open

- Tie the initial abundance of non-major states to mutation-selection balance
  with their parents, as an option. It should fix the failing cell and add the
  first sample's level information everywhere. It assumes the landscape before
  the first sample was the same, which treatment between samples can break.
- Shell coefficient is a lower bound from absence, so unobserved neighbours
  are still predicted too high in rugged landscapes (bias +0.02 at one step,
  +0.04 at two). Learning it across samples of a cohort is the planned remedy.
- Intervals are too wide at distance 0 and slightly narrow at distance 2.
- Held-out landscapes and the HPC confirmation have not been run.

## Files

`forward_study.R` reproduces the forward-model comparison. `compare.R` and
`by_scenario.R` tabulate a method's fits against the baselines per stratum and
per scenario; the official verdict comes from
`benchmark/sparse_two_timepoint/run.R --stage=report`.
