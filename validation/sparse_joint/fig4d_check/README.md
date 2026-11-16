# Why alfak2 reverses the ALFA-K paper's Fig. 4d

The user asked for this check on 2026-11-15, and approved the plan the same
evening.

## The finding being checked

In the paper benchmark (round 05a8855, `alfak2-paper-benchmark`,
`docs/report_05a8855.md`, section 5), Fig. 4d compares two settings. It
measures how much single-chromosome changes alter fitness (|Δf|) in PDX and
in vitro.
- The paper finds them larger in PDX (`PDX_v` +1.23).
- Every alfak2 arm finds them larger in vitro (−0.38 to −1.40).
- alfakR dev2 agrees with the paper.

Δf = f(frequent karyotype) − f(one-chromosome neighbour), as in the paper.
"In vitro" means the lines SA039 and SA906; the other Salehi lines are PDX.

> **Later qualification** (2026-11-16, `../beta_choice`). Choosing β by
> marginal likelihood picks rates that are too high when time points are
> close: with 6 days between them, half the picks were 4× or 16× the truth.
> The PDX trajectories have shorter intervals than the in vitro ones, so part
> of their preference for 4–16× the default may be that bias rather than
> higher chromosomal instability. What stands is narrower:
> - Fig. 4d's direction for alfak2 depends on the assumed rates.
> - At any single common β, alfak2 reverses it.
> - At the rates the data prefer, it agrees with the paper.

## Answer

The reversal comes from using **one missegregation rate**, β = 5e-5 (the
paper's `pm`), for both settings.
- alfak2 infers a neighbour's fitness partly from how many cells
  missegregation should have put there. The cells actually seen fall short of
  that, and the shortfall is read as selection.
- On their own marginal likelihood, the in vitro trajectories prefer
  β = 5e-5. The PDX trajectories prefer 4–16 times more.
- When each trajectory is fitted at the rate it prefers, PDX |Δf| exceeds in
  vitro |Δf| (0.051 against 0.025). Fig. 4d's model then gives
  `PDX_v` +0.57 (p = 0.065), the paper's direction.

The check found no coding error. The baseline refits reproduce the
benchmark's landscapes (correlation 0.999998).

On simulations where the truth is known and β is right, alfak2's Δf to the
neighbours has the right sign and about the right size. It is somewhat too
large for unseen neighbours over long intervals. alfakR's Δf is flattened
towards zero, and its median even has the wrong sign. The paper's own
direction is therefore not a ground truth either.

## Step 1: where the in vitro excess comes from (`decompose.R`)

These are the benchmark's fits, selected as its Fig. 4 selects them
(xv > 0, best minobs), with control fits only. Each neighbour is classed by
where alfak2's value comes from:
- `observed_D`: inside the model domain, and seen at some time point;
- `unseen_D`: inside the domain, but never seen;
- `observed_out`: seen, but outside the domain;
- `trend`: neither of these.

| | in vitro Δf | PDX Δf | in vitro share of neighbours |
| --- | --- | --- | --- |
| alfak2, 2 time points, all neighbours | +0.026 | +0.006 | |
| alfak2, `unseen_D` | +0.028 | +0.007 | 61% |
| alfak2, `observed_D` | +0.016 | −0.004 | 20% |
| alfakR, same pairs | −0.003 to −0.006 | −0.005 to −0.009 | |

- In vitro, alfak2 puts nearly every neighbour below its frequent karyotype
  by a similar amount.
- SA039 is extreme, at 0.15. Its fits choose the longest history (3000 days)
  and have about 900 cells per time point.
- Across trajectories, the Δf of unseen neighbours rises with the number of
  cells and with the history. That correlation turned out to be confounded
  (step 3).

## Step 2: against the truth (`abm_truth.R` on HPC, `analyse_truth.R`)

These are the paper's ABM and missegregation simulations, with minobs 10 and
β equal to the simulated rate. The true Δf comes from the Gaussian-random-field
landscape.

| two time points | alfak2 est | true | alfakR est | r alfak2 / alfakR |
| --- | --- | --- | --- | --- |
| ABM, all neighbours | +0.040 | +0.025 | −0.013 | 0.76 / 0.15 |
| ABM, unseen, 2 / 4 / 8 passages | 0.040 / 0.043 / 0.055 | 0.029 / 0.031 / 0.036 | −0.013 / −0.008 / −0.016 | |
| missegregation, all neighbours | +0.029 | +0.026 | −0.010 | 0.57 / 0.20 |

- alfak2 has the right sign. Its size is right on the missegregation
  simulations.
- On the ABM simulations it is too large for unseen neighbours, by about
  50% for windows of 8 passages.
- alfakR flattens Δf.

## Step 3: refits under variants (`refit_variants.R`, `analyse_refits.R`)

All 205 trajectories of the benchmark's "alfak2 2 time points" Fig. 4
selection were refitted: 23 in vitro, 87 PDX controls and 95 PDX
cisplatin-treated. The Fig. 4d context model was then refitted:
`glmmTMB(|Δf| ~ context + (1 | pdx/fi/k_id), Gamma(log))`, on all fits
rather than the 200 lineage draws.

| variant | in vitro \|Δf\| | PDX \|Δf\| | `PDX_v` (p) |
| --- | --- | --- | --- |
| baseline (the benchmark's call) | 0.027 | 0.014 | −1.08 (0.001) |
| every time point drawn down to 320 cells | 0.025 | 0.014 | −1.26 (0.001) |
| history fixed at the interval | 0.026 | 0.017 | −0.50 (0.02) |
| `initial_state = "free"` | 0.026 | 0.033 | −0.50 (0.17) |
| β / 4 | 0.009 | 0.012 | +0.12 (0.57) |

- The number of cells and the history are not the cause. β is.
- On observed neighbours only, β / 4 gives `PDX_v` +0.46 (p = 3e-4).
  At the default β it gives −0.37 (p = 0.007).
- The in vitro frequent karyotypes are more dominant: median share 0.053,
  against 0.022 in PDX. Δf grows with that share (elasticity 0.35), but the
  share does not explain all of the gap.

## Step 3b: which β the data prefer (`refit_beta.R`, `analyse_beta.R`)

The same 205 trajectories were fitted at β = 5e-5 × {1/16, 1/4, 1, 4, 16, 64}.
β is fixed within a fit, so the objective (the negative log marginal
likelihood of the same counts) compares the rates.

| | 1/16 | 1/4 | 1 | 4 | 16 | 64 |
| --- | --- | --- | --- | --- | --- | --- |
| in vitro: trajectories preferring | 2 | 6 | 13 | 2 | 0 | 0 |
| PDX: trajectories preferring | 1 | 1 | 23 | 73 | 65 | 19 |
| in vitro: summed objective vs the default | +340 | +89 | 0 | +203 | +681 | +2430 |
| PDX: summed objective vs the default | +4250 | +2140 | 0 | −1400 | −1030 | +1990 |
| Fig. 4d `PDX_v` at this β for all (p) | +0.50 (0.03) | +0.12 (0.57) | −1.08 (0.001) | −1.69 (1e-5) | −1.55 (5e-5) | −0.55 (0.15) |

At each trajectory's preferred β:
- in vitro |Δf| is 0.025 and PDX |Δf| is 0.051;
- `PDX_v` is +0.57 (p = 0.065).

## Caveats

- **What is identified.** β enters with the birth rate, which is fixed at
  1 per day for every trajectory. So what the data identify is the influx
  (β × birth rate), not β alone.
- **The comparison is lighter than the paper's.** Only the two-time-point
  arm was refitted. Fig. 4d's model was fitted once on all trajectories,
  without the 200 lineage draws, and +0.57 (p = 0.065) is not the paper's
  significant +1.23.
- **The ABM bias remains** at a correct β, for unseen neighbours over long
  intervals.

## What it means for alfak2 (not changed here)

For karyotypes never seen, alfak2's fitness rests on the assumed
missegregation influx. A comparison between settings with different
chromosomal instability needs β per dataset. The history candidates already
do this kind of choice by marginal likelihood, and β candidates could be
chosen the same way. That would be a change to the model, and it needs a new
benchmark round.

## Files

- Scripts:
  - local: `decompose.R`, `analyse_truth.R`, `analyse_refits.R`,
    `analyse_beta.R`;
  - HPC: `abm_truth.R`, `refit_variants.R`, `refit_beta.R`, run through
    `run_hpc.sh` in the benchmark's image;
  - `tasks.tsv` lists the 205 trajectories.
- HPC folder: `/share/lab_crd/taoli/Project/alfak2-multi-timepoint/fig4d_check`.
- Jobs:
  - 21953591 and 21953592: the truth check;
  - 21953594: the variants;
  - 21953799 and 21954006: the β profile.
- Results (git-ignored): `results/`, which holds `decompose.md`, `truth.md`,
  `refits.md`, `beta.md` and the rows they come from.
