# Level of the balance prior: does freeing it help on real data?

Started 2026-11-16, after the check of the paper's Fig. 2d on the Salehi
passages (`alfak2-paper-benchmark/docs/fig2d_check_05a8855.md`). The user
approved the plan the same day, to be run without further questions; defaults
stay as they are until the user decides.

## Question

On the Salehi passages, with one interval between the two time points,
alfak2's fitness differences depend on the history that the balance prior was
given, and the history does not follow how much the population changed. Can
the model be made to take the size of the cloud around the major karyotypes
as given, without losing accuracy on simulations?

## Why the history sets the scale

`initial_state = "balance"` centres the abundance of a major's neighbour at
the first time point on what the major would have fed it over `history`. That
centre is proportional to the missegregation rate and, for a short history,
to the history. The history is chosen among six candidates by marginal
likelihood. So one setting decides two things at once: how large the cloud of
neighbours should be, and how strongly a neighbour's abundance pins its
fitness.

Step 0 (`salehi_profiles.R`, read only, the stored fits of the paper
benchmark's round 05a8855; tables in `results/step0/`):

- 252 two-time-point fits, 66 of them with one interval between the passages
  of the trajectory, 48 of those 15 days apart.
- At 15 days the history chosen is 15, 30 or 100 days in 5, 27 and 16 fits.

| History | Fits | Shell coefficient | Residual scale | Dosage scale | Cloud ÷ major cells | Neutral history of the cloud (days) | Cells in the domain |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 15 | 5 | +0.037 | 0.019 | 0.039 | 0.12 | 27 | 0.84 |
| 30 | 27 | +0.032 | 0.009 | 0.028 | 0.72 | 145 | 0.51 |
| 100 | 16 | −0.001 | 0.003 | 0.016 | 1.17 | 205 | 0.30 |

(Medians. The shell coefficient is the log birth rate of a neighbour against
the majors. "Neutral history" is the history under which neighbours as fit as
the majors would have reached the cloud's size at the assumed rate of 5e-5.)

- Across the 48 fits the history correlates with the residual scale (Spearman
  −0.74), the dosage scale (−0.78) and the shell coefficient (−0.73), and
  with the neutral history of the cloud (0.55).
- Over all 252 fits the cloud is as large as 135 days of neutral accumulation
  (median; quartiles 46 and 176), against intervals of 15 to 25 days.
- The choice among the candidates is fairly firm: the best candidate leads
  the second by 3.1 log units (median); in 26% of the one-interval fits by
  less than 1. Averaging over the candidates would change a quarter to a
  third of the fits and is not pursued.

So the larger the cloud, the longer the history chosen and the smaller every
fitness difference. With a history of 15 or 30 days the cloud is larger than
the balance allows, and the model makes the neighbours 3-4% per day fitter
than the majors.

## The change tested

`cl1.patch` (on bf86c85) adds `balance_level` to `fit_sparse_joint()`, for two
time points:

- `"fixed"`, the default: the model as it is. The fits are identical
  (`results/same_default.R`: objective, summary, history profile, landscape
  and cross-validation of three datasets).
- `"free"`: all neighbours fed by the majors share one log factor over the
  balance, estimated with the latent variables under a normal prior with
  standard deviation `balance_level_sd` (3). The size of the whole cloud then
  says nothing about fitness; the differences between neighbours and the
  change between the time points still do.

The level stands for whatever makes the cloud larger or smaller than the
balance at the assumed rate: a different missegregation rate, a population
that is not at balance, or cells called one copy off.

## Arms

| Arm | Fit |
| --- | --- |
| `base` | The default: balance prior, history among the candidates |
| `free` | `initial_state = "free"`: no balance prior |
| `level` | `balance_level = "free"`, history among the candidates |
| `level_dt` | `balance_level = "free"`, history fixed at the interval (one fit) |

## Data

All simulated, with the development simulator
(`benchmark/sparse_two_timepoint`, sourced unchanged). Every fit assumes a
missegregation rate of 5e-5 and one division per day.

| Set | Datasets | What it is |
| --- | --- | --- |
| `dev` | 168 | The development scoreboard: 12 landscapes, 14 scenarios, intervals of 30-1080 days |
| `short` | 120 | The short-interval diagnosis set: 20 landscapes, first time point at day 120 or 1440, intervals of 6, 18 and 42 days |
| `diverse` | 60 | New (seed base 97000): 12 landscapes at each of five levels, one interval of 15 days, 400 cells per time point |

No simulation so far was as diverse as a Salehi passage, where half the cells
of a time point sit on karyotypes seen once, the largest clone holds 13% and
the model's domain holds 53-69% of the cells. The diverse set reaches that in
three ways (`results/pilot_levels.R` chose the levels):

| Level | Landscape | What differs from the fit's assumptions | Singleton cells | Largest clone | Cells in the domain |
| --- | --- | --- | --- | --- | --- |
| `rate8` | wavelength 0.8 | true missegregation rate 8 times the assumed | 0.45 | 0.15 | 0.68 |
| `rate2` | wavelength 1.6 | true rate 2 times the assumed | 0.42 | 0.08 | 0.66 |
| `smooth` | wavelength 2.4 | nothing | 0.23 | 0.15 | 0.85 |
| `call2` | wavelength 0.8 | each chromosome of a sampled cell is called one copy off with probability 2% | 0.21 | 0.38 | 0.92 |
| `call1` | wavelength 1.6 | calling errors of 1% | 0.36 | 0.11 | 0.76 |

Held-out landscapes (seed bases 82000, 94000, 96000, 98000) are not used.

## Scores

Per fit (`results/scores.csv`), against the true landscape:

- RMSE of relative fitness at distances 0, 1 and 2 from the observed
  karyotypes, with the RMSE of predicting no differences beside it, the bias,
  and the direction accuracy of one-copy steps (the scores of the development
  benchmark).
- Fitness variance among the cells of the first time point, predicted and
  true: how compressed or inflated the fitness differences are.
- A counterpart of the paper's Fig. 2d: the change that selection alone makes
  over the interval from the first time point's cells (predicted fitness, and
  true fitness as the ceiling), and the observed change to the second time
  point's cells, as Wasserstein distances.

## Criteria, set before the fits were scored

A candidate arm is worth proposing as a default only if all of these hold.
Landscapes are the unit of the bootstrap.

1. **No loss on `dev` and `short`.** Mean RMSE at each of the distances 0, 1
   and 2, pooled over the core scenarios of `dev`, and over `short`, is at
   most 3% above `base`. Direction accuracy is not lower by more than 0.01.
2. **Not worse on `diverse`.** The same two rules, per mechanism.
3. **The scale no longer follows the history.** On `diverse` and on the
   Salehi inputs, the fitness variance among the first time point's cells
   does not depend on the chosen history (`level`), or there is no history to
   choose (`level_dt`).
4. **Salehi, as a confirmation only.** On the trajectories with one interval,
   the rank agreement of predicted and observed change is not below alfakR's
   (paired difference with an interval that includes zero or lies above it).
   No setting is tuned on the Salehi results.

If no arm passes, the default stays, and the dependence on the history is
recorded as a known limit of fits with one short interval.

## Files

- `salehi_profiles.R`, `run_hpc.sh`: step 0.
- `cl1.patch`: the change, against bf86c85.
- `setup.R`, `stages.R`, `run.R`: sets, arms, fits and scores.
- `hpc/`: build and array scripts for HPC (not tracked).
- `salehi/variants.R`: the arms as methods of the paper benchmark's task
  script; `salehi/array_task.sh` runs them in the benchmark's image with the
  benchmark's folder read only.
- `report.R`: the tables of the results below.
- `results/` (not tracked): inputs, builds, fits, scores, logs.

## Derived arms

Added on 2026-11-16 after the interim scores of `dev` and `short` showed that
a free level costs about a fifth of the accuracy there, and before the
`diverse` and Salehi results were looked at. No new fits: per dataset, the
fit with the higher marginal likelihood among

- `select`: `base` and `level_dt`;
- `select3`: `base`, `level` and `level_dt`.

## Results (2026-11-16)

All 1392 simulated fits and all 298 Salehi fits finished. Tables:
`results/report/` (`report.R`) and `results/salehi/report/` (`salehi/analyse.R`).

**No arm meets the criteria.** The default stays.

RMSE of relative fitness, ratio to `base` (mean over datasets; bootstrap over
landscapes):

| Arm | dev core d0 / d1 / d2 | short d0 / d1 / d2 | diverse d0 / d1 / d2 |
| --- | --- | --- | --- |
| `free` | 1.24 / 1.23 / 1.19 | 1.01 / 1.03 / 1.02 | 0.78 / 0.70 / 0.70 |
| `level` | 1.18 / 1.20 / 1.19 | 1.05 / 1.13 / 1.16 | 0.89 / 0.82 / 0.82 |
| `level_dt` | 1.18 / 1.20 / 1.19 | 1.03 / 1.04 / 1.04 | 0.87 / 0.79 / 0.78 |
| `select` | 1.02 / 1.04 / 1.05 | 1.01 / 1.03 / 1.03 | 0.91 / 0.87 / 0.87 |
| `select3` | 1.03 / 1.05 / 1.06 | 1.01 / 1.03 / 1.04 | 0.91 / 0.86 / 0.86 |

- On the development landscapes, where the assumed missegregation rate is the
  true one, a free level costs as much as dropping the balance prior. The
  level of the cloud and the mean fitness of the neighbours against the majors
  move all neighbours' abundances alike, so freeing the one gives up what the
  balance prior says about the other. The fitted levels there are mostly below
  zero (median −0.6): the cloud is smaller than the balance at the interval
  because its members are less fit.
- `select` comes closest: within 3% on `short`, 4-5% above `base` at
  distances 1 and 2 on `dev` core (intervals 1.02-1.07 and 1.03-1.08).

**Where the fit's assumptions fail, the balance prior gets directions wrong.**
Direction accuracy of one-copy steps at distance 1 on `diverse`:

| Level | `base` | `free` | `level_dt` | `select` |
| --- | --- | --- | --- | --- |
| `rate8` | 0.17 | 0.58 | 0.47 | 0.40 |
| `rate2` | 0.25 | 0.48 | 0.47 | 0.45 |
| `call1` | 0.36 | 0.61 | 0.55 | 0.55 |
| `call2` | 0.48 | 0.71 | 0.67 | 0.61 |
| `smooth` | 0.67 | 0.68 | 0.67 | 0.67 |

- With the true rate 2-8 times the assumed one, or one chromosome in 50-100
  called a copy off, `base` ranks neighbours worse than chance and its RMSE
  at distance 1 (0.073) is above that of predicting no differences (0.051).
  The cloud is larger than the balance at the assumed rate (median fitted
  level +1.6 to +2.5), and the model reads that as neighbours as fit as the
  majors, or picks a long history.
- On `smooth`, where the assumptions hold, `base` is as good as `free`. Its
  cloud is also larger than the balance at the interval (median level +0.9
  with `level`, +1.9 with `level_dt`), so the level does not tell a wrong
  rate or noisy calls from a diverse population that fits the model.
- The development benchmark's direction accuracy is 0.78-0.82 for every arm.

**Salehi (confirmation).** Fig. 2d without the paper's selection: Spearman of
predicted and observed change on the rows that every method has (lineage
bootstrap):

| Arm | All (471 rows) | One interval (147) | More intervals (324) |
| --- | --- | --- | --- |
| alfakR | 0.32 | 0.60 | 0.18 |
| alfak2 (`base`) | 0.29 | 0.45 | 0.22 |
| `level` | 0.34 | 0.59 | 0.21 |
| `level_dt` | 0.33 | 0.60 | 0.21 |
| `select` | 0.31 | 0.54 | 0.21 |
| `free` | 0.37 | 0.60 | 0.25 |

- On one interval, `level_dt` and `free` gain 0.15 over `base` (lineage
  interval 0.02-0.18 for both) and equal alfakR (−0.002). `select` gains
  0.09 (−0.02 to 0.17).
- With more intervals no arm differs from `base`.
- At 15 days, `base`'s fitness variance among the parent's cells follows the
  history it chose (Spearman −0.72); with `level` −0.57; `level_dt` has
  none to choose.
- The fitted level is above zero in 73% of the Salehi fits (median +0.9):
  the clouds are larger than the balance at 5e-5 allows, as in `diverse`.
- The variants of the alfak2-salehi-variants project (same rows): `wide`
  0.52, `observed` 0.44, `observed_free` 0.50 on one interval. Putting more
  karyotypes in the domain does not help.

**A limit of the Fig. 2d test.** On `diverse`, the change that the true
landscape predicts ranks the observed change at Spearman 0.15 only (on `dev`
0.93): with 400 cells from a population this diverse, sampling noise
dominates the observed change. A test of this kind cannot tell a better
landscape from a worse one on such data.

## Conclusion

- The balance prior is what makes alfak2 accurate when its missegregation
  rate is right, and what makes it wrong when the cloud is larger for other
  reasons. No arm tested keeps the first and avoids the second.
- The size of the cloud cannot be attributed from the counts of one
  population: a high rate, calling errors and neutral neighbours give the
  same picture. This is the same limit as choosing the rate
  (`beta_choice/`).
- (First round.) The default was then unchanged; the follow-up below made
  `initial_state = "auto"` the default. For data where the rate is uncertain or calls
  are noisy, `initial_state = "free"` is the safer fit: on the development
  landscapes it costs about a fifth of the accuracy, where the balance
  prior is wrong it saves a quarter and turns wrong directions into right
  ones. Neither the fitted level nor the marginal likelihood says which case
  a dataset is. `cl1.patch` is not applied to the package.

## Follow-up (2026-11-16, approved by the user)

Plan: (1) confirm on new landscapes a rule that switches to
`initial_state = "free"` when the level fitted by `level_dt` exceeds 1;
(2) find why alfak2 trails alfakR most when every passage is fitted; (3) make
the rule an option; (4) the user decides the default.

### Step 1: confirmation sets

The designs of `dev`, `short` and `diverse` on new landscapes: `dev2` (168
datasets, seed base 101000), `short2` (120, 102000), `diverse2` (60, 103000).
Arms `base`, `free`, `level_dt`; the rule and its threshold were fixed before
these fits (`confirm.R`). 1044 fits, HPC array 21974176; tables in
`results/report_confirm/`.

**The rule passes every criterion.** RMSE ratio to `base` at distances 0 / 1 / 2
(bootstrap over landscapes) and direction accuracy:

| Set | Switched to `free` | Rule | `free` everywhere | Direction, rule − base |
| --- | --- | --- | --- | --- |
| `dev2` core | 0% | 1.00 / 1.00 / 1.00 | 1.21 / 1.19 / 1.16 | 0 |
| `short2` | 23% | 1.00 / 1.00 / 0.99 | 1.00 / 1.02 / 1.00 | −0.004 to +0.006 |
| `diverse2` | 100% | 0.81 / 0.74 / 0.73 | the same | +0.13 / +0.29 / +0.26 |
| `dev2` stress (reported) | 6% | 1.02 / 1.03 / 1.03 | 1.13 / 1.17 / 1.17 | +0.003 |

- Every `dev2` core dataset keeps the balance prior (fitted levels below 1).
- In `short2` the rule switches at wavelengths 0.8 and 1.6 (40% and 50%),
  hardly at 0.2 and 0.4, with RMSE unchanged.
- On `diverse2` the default ranks one-copy neighbours at 0.37 (distance 1);
  the rule at 0.66.
- On the Salehi passages it switches half of the inputs. There its gain on
  Fig. 2d (0.45 to 0.59 on one interval) is mostly how much the forecast
  follows sampling noise (step 2), so Salehi neither supports nor refutes it.

### Step 3: the option

`cl2.patch` (on bf86c85; includes `cl1.patch`) adds
`fit_sparse_joint(initial_state = "auto", auto_level = 1)` for two time points:
a fit with `balance_level = "free"` and the history at the interval, then the
free or the default fit by its level. `diagnostics$initial_state_check` keeps
the level, its standard deviation and the choice. More than two time points
give an error.

- At its defaults the build gives the fits of bf86c85 exactly
  (`results/same_default.R cl2`: objective, summary, history profile,
  landscape, cross-validation of three datasets).
- On seven confirmation datasets it makes the rule's choice
  (`check_cl2.R`). Its fits equal those of build cl1 on the same machine
  (`results/platform_check.R`); against the study's HPC fits, free fits
  differ by up to 0.013 in relative fitness, which is the difference
  between machines.
- `as_alfakR()` runs on an `"auto"` fit.
- The package's tests pass (`test-sparse-joint.R` 18, `test-sparse-joint-rates.R`
  5), with the two tests that build the model's inputs by hand given the new
  inputs and one test added for the option.
- Cost: one extra fit at the interval, a few seconds to a minute.

### Step 4: the default, and its reversal

The user made `initial_state = "auto"` the default on 2026-11-16 (alfak2
827c4b5): `"auto"` was `"balance"` with more than two time points, stored fits
still rebuilt, and the test suite passed.

The ALFA-K paper benchmark of 827c4b5 (the `benchmark` session,
`results_827c4b5/analysis/initial_state/`) then showed a regression at time
points one passage apart (windows of two passages, `w2_ends`):

| Simulations | Switched to `free` | Landscape r, free vs balance (05a8855) | CV against truth |
| --- | --- | --- | --- |
| ABM, 400 inputs | 54% (λ 0.4 31%, 0.8 97%, 1.6 86%) | 0.58 vs 0.75, worse in 93% | 0.47 vs 0.66 |
| Missegregation, true rate 5e-5 (= assumed) | 23/25 | 0.53 vs 0.72, worse in 91% | 0.40 vs 0.62 |
| Missegregation, 1e-4 to 8e-4 | 84-100% | worse at every rate | worse at every rate |

- Windows of four and eight passages switched almost never (3 inputs).
- Forecasts beat the baseline slightly more often after the switch (14% of
  fits better, 6% worse at passage 10), not enough to offset the landscapes.
- The level is measured against the balance over the interval. Over one
  passage that balance is far below the cloud a population builds up over its
  whole history, whatever the rate, so the level measured the history, not a
  misspecified rate. The confirmation sets had intervals of 6 days or more,
  but no population that young sampled one passage apart, and were scored by
  RMSE, where `free` was not worse; the paper benchmark scores landscape
  correlation and cross-validation against the truth.

The user returned the default to `"balance"` the same day (alfak2 348addd);
`"auto"` stays an option with this limitation documented.

### Step 5: a corrected rule (stopped after its step 0)

`rule/README.md`: no feature of two time points separates a wrong rate or
miscalled cells from a cloud that is large for other reasons. The best
threshold, on the level at the interval and falling with the interval, chosen
on the same data, caught half of the inputs where switching helps and still
switched 14% of those where it hurts. The user stopped the effort there on
2026-11-16; the remaining steps (a fixed rule, confirmation on new seeds) were
not run.
