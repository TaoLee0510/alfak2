# Choosing the missegregation rate per dataset

The user approved this plan on 2026-11-15.

The paper benchmark gives every dataset one missegregation rate, β = 5e-5.
alfak2 infers a neighbour's fitness partly from how many cells
missegregation should have put there, so the rate matters for karyotypes
that are rare or unseen. The Fig. 4d check (`../fig4d_check`) found that
the paper's Salehi PDX trajectories prefer a rate 4–16 times higher than
the in vitro ones. Step 2 below shows that part of this can be a bias at
short intervals.

The plan is to let alfak2 pick β per dataset by marginal likelihood, the
same way it picks the history:
- step 0: test whether that picks the right rate where it is known;
- step 1: implement it;
- step 2: check that the simulated sets do not regress;
- step 3: run a new paper benchmark round. Only alfak2 is refitted; alfakR
  keeps the paper's `pm`.

An earlier check (`../README.md`, rate sensitivity, before history
candidates and the unobserved deficit) found that the marginal likelihood
did not locate the birth rate, which enters as a product with β. Step 0
therefore tests the current model before anything is built.

## Step 0: does the marginal likelihood pick β? (`sim_fit.R`, `analyse_sim.R`)

**Inputs.** 675 inputs from the paper's simulations, built by the paper
benchmark's own code (`R/abm.R`): the first and last time point of the last
2, 4 and 8 passages before day 120. That is 225 simulations × 3 windows:
- the missegregation simulations, all 125 (true β = 5e-5, 1e-4, 2e-4, 4e-4
  and 8e-4; 25 each);
- the ABM simulations, the first 25 replicates of each wavelength, 100 in
  all (true β = 5e-5).

**Fits.** Each input was fitted at β = 5e-5 × {1/16, 1/4, 1, 4, 16, 64}.
Inputs whose true rate lies off that grid were also fitted at the true
rate. That makes 4200 fits, each its own HPC task. The pick is the grid rate
with the lowest objective (the negative log marginal likelihood of the same
counts, with the history already chosen inside each fit). Landscape
accuracy is the benchmark's Pearson r against the true fitness, on its
minobs-10 karyotype set.

**Criteria**, fixed before any fit was scored:
1. At least 70% of picks lie within 4× of the true rate, and the median
   |log2(pick / true)| is at most 1.
2. The paired median of r(pick) − r(true rate) is at least −0.02. Where the
   true rate is at least 4× the default, the paired median of
   r(pick) − r(5e-5) is above 0.

**Result: passes.** No fit failed.

| check | value | needed |
| --- | --- | --- |
| picks within 4× of the truth | 0.93 | ≥ 0.70 |
| median \|log2(pick / true)\| | 1.00 | ≤ 1 |
| median r(pick) − r(true) | 0.000 | ≥ −0.02 |
| median r(pick) − r(5e-5), truth ≥ 4× (225 inputs) | +0.021 | > 0 |

Picks, as multiples of 5e-5, over all three windows:

| truth | 1/16 | 1/4 | 1 | 4 | 16 | 64 |
| --- | --- | --- | --- | --- | --- | --- |
| ABM, 5e-5 (300) | 17 | 130 | 151 | 2 | 0 | 0 |
| missegregation 5e-5 (75) | 2 | 33 | 35 | 5 | 0 | 0 |
| missegregation 1e-4 (75) | 2 | 6 | 61 | 6 | 0 | 0 |
| missegregation 2e-4 (75) | 1 | 5 | 25 | 44 | 0 | 0 |
| missegregation 4e-4 (75) | 0 | 2 | 8 | 60 | 5 | 0 |
| missegregation 8e-4 (75) | 0 | 1 | 4 | 29 | 36 | 5 |

- **The pick follows the truth.** Where the truth is the default, it stays
  at 1 or 1/4. Where the truth is 8e-4, it moves to 4–16.
- **It leans one step low.** The simulated populations divide about 0.9
  times a day, and the fit assumes 1, so part of this is expected.
- **Accuracy at the pick equals accuracy at the true rate**, in median.
  Where the true rate is high, the pick is better than 5e-5. At
  missegregation 8e-4, r is 0.31 against 0.07 with 2 passages, and 0.57
  against 0.35 with 4.
- **The weakest case is 2-passage windows at a high true rate.** At 8e-4 the
  pick's r is 0.31, against 0.49 at the truth. At 4e-4 it is 0.55 against
  0.61.

## Step 1: the option (commit ac36fd1, removed again)

ac36fd1 let `fit_sparse_joint(beta = ...)` take several values, or `"auto"`
(5e-5 × 1/16, 1/4, 1, 4, 16, 64). Each candidate was fitted and the fit of
lowest objective kept, and `diagnostics$beta_profile` listed them. The
default stayed fixed. After step 2 the option was removed (below); the code
is in ac36fd1.

## Step 2: the simulated sets must not regress (`analyse_series.R`)

These are the sets of `../series_methods` (build cv10), whose true rate is
the default 5e-5:
- e1: passaged cultures;
- e2: clonal sweeps;
- e3: neutral cultures.

`ends` and `anchored` were fitted at each other multiple of the grid, one
fit per HPC task: 7680 fits. `run.R --stage=compose` then kept, per dataset,
the candidate of lowest objective (`ends_auto`, `anchored_auto`), which is
what `beta = "auto"` returns.

**Criterion**, fixed before scoring, with the RMSE ratio auto / fixed paired
by dataset:
- e1: at least 18 of 24 cells within 1.03 at every distance, and the median
  over cells at most 1.02;
- e2 and e3: the median at most 1.03.

**Result: fails.**

| | e1 cells within 1.03 | e1 median ratio, distances 0/1/2 | e2 | e3 |
| --- | --- | --- | --- | --- |
| `ends_auto` against `ends` | 6 / 24 | 1.23 / 1.43 / 1.46 | 1.03 / 1.11 / 1.17 | 1.00 |
| `anchored_auto` against `anchored` | 7 / 24 | 1.27 / 1.46 / 1.44 | 1.02 / 1.08 / 1.10 | 1.00 |

The loss is concentrated at short intervals. In e1 the time points are
taken every 6 or every 60 days, so the first and last of the two-time-point
fit are 18–42 or 180–420 days apart. `ends_auto`, median ratio:

| bottleneck | every 6 days (18–42 days) | every 60 days (180–420 days) |
| --- | --- | --- |
| 200 | 3.20 | 1.00 |
| 2000 | 1.78 | 1.00 |
| 20000 | 1.77 | 1.32 |

The rates chosen are biased high when time points are close. In e1 the
truth is ×1 for every dataset; `ends_auto` picks, as multiples of 5e-5:

| | 1/16 | 1/4 | 1 | 4 | 16 |
| --- | --- | --- | --- | --- | --- |
| every 6 days (144) | 10 | 19 | 38 | 66 | 11 |
| every 60 days (144) | 5 | 46 | 68 | 21 | 4 |

Step 0's simulations had 7 to 54 days between their two time points, which
is why the bias did not show there.

## Decision

The gate failed, so step 3 (the paper benchmark round) was not run. Only
the alfak2 build for it had been made. The option was removed, and the rate
stays fixed.

This also qualifies the Fig. 4d check. There the PDX trajectories preferred
4–16× the default and the in vitro ones the default, but the PDX
trajectories have shorter intervals. Part of their preference may be this
bias rather than higher chromosomal instability.

# The short-interval plan

The user approved this plan on 2026-11-16. Its goal was a rule for choosing
β that keeps datasets whose true rate is the default where they were, and
still finds a rate that really is higher. Its steps were:
- step 1: find what drives the high picks;
- step 2: simulate short intervals with a truly higher rate;
- step 3: test rules offline;
- step 4: re-add the option if a rule passes.

A model change (a drift term) was outside the plan. Why it is not pursued is
under Decision.

## Step 1: what drives the high picks (`mechanism.R`)

These are the e1–e3 candidate fits, where the true rate is ×1. For each
dataset, gain4 is objective(×1) − objective(×4); a positive value prefers
4×. It is regressed on features of the first and last time point.

- **Short spacing with few cells.** With time points every 6 days and 200
  cells, e1 picks ≥ 4× in 71–75% of datasets at every bottleneck.
- **A small bottleneck.** With 1000 cells, e1 picks ≥ 4× in 71% of datasets
  at a bottleneck of 200, but in 25% at 2000 and 8% at 20000.
- **Every 60 days**, the share is 4–42%.
- **Across all sets** (R² 0.52), these features push towards higher rates:
  - fewer karyotypes seen;
  - a shorter interval;
  - more turnover between the two time points;
  - fewer new karyotypes;
  - a less dominant top karyotype;
  - a smaller bottleneck.

  Cells per time point adds nothing once the number of karyotypes seen is
  in.

**Reading.** Over a short interval, drift at the passages and sampling noise
change the composition more than selection can. The model has no drift term,
so it explains the change as more missegregation influx. When few
karyotypes are seen, the evidence against a high rate is weak: that
evidence is the many neighbours that stay unseen.

## Step 2: short intervals with a truly higher rate (set e4)

Step 0's high-rate simulations had long intervals and 1000 cells. e1–e3 have
short intervals but a true rate of ×1. e4 fills the gap.

**Design.** e4 uses e1's design:
- bottleneck 200, 2000 or 20000;
- time points every 6 or 60 days;
- 4 or 8 time points;
- 200 or 1000 cells.

There are 12 landscapes:
- the odd ones were simulated from day 0 at 4× the default (2e-4);
- the even ones at 16× (8e-4).

That makes 288 datasets, generated with
`../series_methods/run.R --stage=truth|datasets --set=e4`. The
two-time-point fit (`ends`) was run at the six candidate rates: 1728 HPC
tasks, one fit each. 22 datasets fail at every rate, with no cell inside
the model domain, and are left out below.

**Knowing the rate helps.** These are median RMSE ratios against the fit at
the true rate, at distance 1 (one chromosome from the majors):

| truth | assumed 4× too low | 4× too high |
| --- | --- | --- |
| ×1 (e1) | 1.21 | 2.89 |
| ×4 | 1.59 | 3.66 |
| ×16 | 1.64 | 1.90 |

**The lowest objective does not find it.** Picks (R0) on the 266 datasets
that fit:

| truth | 1/16 | 1/4 | 1 | 4 | 16 | 64 |
| --- | --- | --- | --- | --- | --- | --- |
| ×4 (143) | 2 | 22 | 44 | 49 | 26 | 0 |
| ×16 (123) | 13 | 6 | 20 | 31 | 41 | 12 |

The picks follow the design, not the truth:
- **Every 6 days**, ≥ 4× is picked in 75–93% of e4 datasets. e1, whose
  truth is ×1, gets the same with 200 cells (71–75%).
- **Every 60 days with 1000 cells**, a true 4× is picked as ≥ 4× in only 6%
  of datasets.

## Step 3: rules (`rules.R`)

Each rule maps the six candidate objectives, the interval between the first
and last time point, and the mean cells per time point to one rate.
Abstaining means ×1.

| rule | what it does |
| --- | --- |
| R0 | the lowest objective |
| A(I, C) | R0 only if the interval ≥ I days and cells ≥ C; I ∈ {20, 50, 100, 200}, C ∈ {0, 300, 500, 1000} |
| B(T) | R0 only if it beats ×1 by at least T log units; T ∈ {2, 5, 10, 20, 40, 80} |
| C(s) | the lowest objective + (log multiple / s)² / 2, a lognormal prior around the default; s = log 2, 4, 8, 16 |
| AB | A and B together; 12 combinations |

**Tuning and testing.** Rules were tuned on landscapes L01–L06 and step-0
replicates 1–12. Testing was to use L07–L12 and replicates 13–25.

**Criteria**, fixed before any e4 fit was scored:
1. **Null** (e1–e3, truth ×1, `ends` and `anchored`): step 2's criterion.
2. **Alternatives** (step 0 and e4):
   - at least 70% of picks within 4× of the truth;
   - where the truth is ≥ 4×, more accurate than ×1. For step 0 that means
     a paired median r gain above 0; for e4, a median RMSE ratio below 1 at
     every distance.
3. **Salehi.** Report the in vitro and PDX picks under every rule.

**Result: no rule passes on the tuning half, so none went to the test
half.** Some rules, tuning half:

| rule | e1 cells within 1.03 (ends / anchored) | within 4× | e4 median ratio | step-0 gain | e1 moved | e4 moved | e4 better / worse |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R0 | 4 / 5 | 0.86 | 1 | +0.017 | 59% | 76% | 42% / 33% |
| A(I ≥ 50, C ≥ 300) | 19 / 19 | 0.69 | 1 | 0 | 16% | 17% | 4% / 10% |
| B(T = 5) | 22 / 20 | 0.83 | 1 | 0 | 17% | 47% | 22% / 23% |
| B(T = 20) | 24 / 24 | 0.74 | 1 | 0 | 1% | 13% | 1% / 10% |
| C(s = log 2) | 19 / 18 | 0.84 | 1 | 0 | 29% | 60% | 32% / 26% |
| C(s = log 4) | 11 / 9 | 0.87 | 1 | +0.012 | 44% | 75% | 44% / 30% |
| AB(I ≥ 20, C ≥ 300, T = 5) | 24 / 22 | 0.76 | 1 | 0 | 6% | 19% | 8% / 9% |

"Moved" is the share of datasets taken off the default. "Better" and
"worse" mean an RMSE ratio at distance 1 below 0.97 or above 1.03.

**The two requirements pull apart.** Rules that keep the null sets where
they were do so by staying at the default in most datasets. They then also
stay at the default in most e4 datasets, and on the step-0 high-rate
simulations. Where any rule moves, the move is about as often worse as
better. R0 picks the true rate in 34% of the e4 datasets that fit.

With few cells, short intervals or strong drift, two time points do not tell
a high missegregation rate apart from drift and sampling noise.

**Salehi** (`results/rules.md`, 23 in vitro and 182 PDX trajectories):
- **Medians.** In vitro trajectories have a median interval of 75 days and
  543 cells per time point; PDX trajectories 45 days and 370 cells.
- **R0** moves 86% of PDX trajectories up (≥ 4×) and 9% of the in vitro ones.
- **Null-passing rules still move PDX up.** B(T = 10) moves 46% of PDX up
  and none in vitro. B(T = 20) moves 20% of PDX up.
- **Gates on interval remove PDX.** A(I ≥ 100, any C) moves 0–2% of PDX.

None of these rules is validated, so none of these picks can be read as
PDX's rate.

## Fallback: a range of rates (`range.R`)

When the rate cannot be chosen, the conclusions can be shown at several
rates. `fit_sparse_joint_rates()` already does this for the birth rate. It
now also takes `beta_range`.

A range of β is not the same as a range of birth rates. The birth rate also
sets the history candidates, the calibration floor and the intercept prior.
On one e4 dataset, β × 4 and birth rate × 4 gave objectives 244.4 and 246.6,
and standard deviations up to 0.22 apart.

**Offline check.** `range.R` used the e1 and e4 candidate fits:
- the point is the fit at one end or at the centre of the range;
- the interval is the union of the fits' 95% intervals, as `predict()`
  forms it.

Each was scored against the truth like a fit and compared with the fit at
the default alone. The table gives medians at distance 1:

| truth | range, point | RMSE ratio | coverage (default alone → range) | interval score ratio |
| --- | --- | --- | --- | --- |
| ×1 | 1–16×, lowest | 1.00 | 1.00 → 1.00 | 9.4 |
| ×4 | 1–16×, lowest | 1.00 | 0.97 → 1.00 | 2.9 |
| ×16 | 1–16×, lowest | 1.00 | 0.93 → 1.00 | 1.7 |
| ×1 | 1/4–4×, lowest | 1.46 | 1.00 → 1.00 | 3.8 |
| ×4 | 1/4–4×, lowest | 1.19 | 0.97 → 1.00 | 1.5 |
| ×1 | 1/4–4×, centre | 1.00 | 1.00 → 1.00 | 3.5 |

- **Coverage.** Where the true rate is higher, the default alone covers
  less: 0.85–0.87 at distance 1 with time points every 60 days. The range
  restores coverage to 1.
- **Cost.** The intervals are wider, most of all when the truth is the
  default.
- **The point prediction** is the fit at the chosen end of the range, so it
  is as accurate as that rate. Assuming a rate too high costs more than
  assuming one too low (step 2). Taking the lowest rate is therefore safe
  when the range holds the truth.
- **Practical use.** Put the rate believed most likely at the lowest end:
  for example 1–16× 5e-5 when a higher rate is suspected. Use
  `point = "centre"` for a range centred on it.

## Decision

β stays fixed at 5e-5 and is not chosen from the data. Step 4 was not
reached, so no option was re-added and no benchmark round was run. When the
rate is uncertain, fit over a range of rates:
`fit_sparse_joint_rates(beta_range = )`, below.

The Fig. 4d qualification stands. The PDX preference for 4–16× cannot be
told apart from the short-interval bias.

**A drift term is not pursued** (decided 2026-11-16):
- **Two time points.** Each karyotype changes once between them, and
  selection, missegregation influx and drift each explain that change. A
  drift estimated from the counts adds a third unknown to two that are
  already confounded.
- **More time points.** An estimated drift was tried there:
  `series = "state_space"` ([series methods](../series_methods/README.md)).
  - It passed 16 of 24 e1 cells.
  - On neutral cultures, where every fitted difference is drift read as
    selection, its error at one and two steps was up to 1.2 times that of
    the first and last time point alone.
  - Its sd was too small at a controlled false-positive rate.

  With two time points it is not used.
- **A drift fixed from the passaging bottlenecks** (not estimated) is
  untested. At best it would remove the preference for a high rate and
  leave the likelihood flat in β, which is the fixed default again. It also
  needs the bottleneck at every passage.

Information on β has to come from outside the counts:
- a measured rate, given as `beta` or `beta_range`;
- perhaps replicate cultures from one start. Drift differs between
  replicates; selection and influx do not. This is untested.

## Files

- **Tracked:**
  - `sim_fit.R` (one fit per task) and `analyse_sim.R` (the step 0 verdict);
  - `analyse_series.R` (the step 2 verdict). The candidate and composed
    methods are in `../series_methods/setup.R`.
  - Short-interval plan:
    - `mechanism.R` (step 1);
    - `rules.R` (step 3);
    - `range.R` (the fallback).

    e4 is `../series_methods/run.R --stage=truth|datasets --set=e4`.
- **Untracked:** the HPC wrapper `run_hpc.sh` and `tasks.tsv`.
- **HPC folder:** `/share/lab_crd/taoli/Project/alfak2-multi-timepoint/beta_choice`.
  `fits_collided/` there holds the first run's outputs. The ABM wavelength
  0.8 simulations and the 5e-5 missegregation simulations share ids, and
  the first output names did not include the source. Those 450 fits were
  rerun.
- **Jobs:**
  - 21954264, 4200 tasks;
  - 21958549, the 450 reruns;
  - 21966818 (2 test tasks) and 21966820 (the other 1726), the e4 fits,
    under `alfak2/validation/sparse_joint/series_methods/results/cv10/e4` in
    the HPC project folder.
- **Results (git-ignored):** `results/sim_report.md`,
  `results/sim_choice.rds`, `results/series_report.md`,
  `results/mechanism.md`, `results/rules.md` (with `rules.rds`) and
  `results/range.md` (with `range.rds`).
