# Choosing the missegregation rate per dataset

The user approved this plan on 2026-11-15.

The paper benchmark gives every dataset one missegregation rate, β = 5e-5.
alfak2 infers a neighbour's fitness partly from how many cells
missegregation should have put there, so the rate matters for karyotypes
that are rare or unseen. The Fig. 4d check (`../fig4d_check`) found that
the paper's Salehi PDX trajectories prefer a rate 4–16 times higher than
the in vitro ones.

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

## Step 1: the option (`fit_sparse_joint(beta = ...)`)

- **One value** fixes the rate, as before; the fit is bit for bit the same.
- **Several values, or `"auto"`** (5e-5 × 1/16, 1/4, 1, 4, 16, 64), are
  candidates. Each is fitted, with the history chosen inside it, and the fit
  of lowest objective is kept.
- **Diagnostics:** `diagnostics$beta_profile` lists the candidates.
- **Pairwise and state space** take the rate the anchored model chooses.
- **Cohorts** still take one rate.
- **alfakR interface:** `alfak_compat(pm = "auto")` passes it through.
- **Default:** the rate stays fixed at 5e-5 until the paper benchmark
  (step 3) has been weighed.
- **Tests:** `tests/testthat/test-sparse-joint-beta.R`.

## Files

- **Tracked:** `sim_fit.R` (one fit per task) and `analyse_sim.R` (the
  verdict).
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
