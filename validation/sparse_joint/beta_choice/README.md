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

The loss is concentrated at short intervals. In e1 (`ends_auto`, median
ratio):

| bottleneck | 6 days apart | 60 days apart |
| --- | --- | --- |
| 200 | 3.20 | 1.00 |
| 2000 | 1.78 | 1.00 |
| 20000 | 1.77 | 1.32 |

The rates chosen are biased high when time points are close. In e1 the
truth is ×1 for every dataset; `ends_auto` picks, as multiples of 5e-5:

| | 1/16 | 1/4 | 1 | 4 | 16 |
| --- | --- | --- | --- | --- | --- |
| 6 days apart (144) | 10 | 19 | 38 | 66 | 11 |
| 60 days apart (144) | 5 | 46 | 68 | 21 | 4 |

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

## Files

- **Tracked:**
  - `sim_fit.R` (one fit per task) and `analyse_sim.R` (the step 0 verdict);
  - `analyse_series.R` (the step 2 verdict). The candidate and composed
    methods are in `../series_methods/setup.R`.
- **Untracked:** the HPC wrapper `run_hpc.sh` and `tasks.tsv`.
- **HPC folder:** `/share/lab_crd/taoli/Project/alfak2-multi-timepoint/beta_choice`.
  `fits_collided/` there holds the first run's outputs. The ABM wavelength
  0.8 simulations and the 5e-5 missegregation simulations share ids, and
  the first output names did not include the source. Those 450 fits were
  rerun.
- **Jobs:**
  - 21954264, 4200 tasks;
  - 21958549, the 450 reruns.
- **Results (git-ignored):** `results/sim_report.md` and
  `results/sim_choice.rds`.
