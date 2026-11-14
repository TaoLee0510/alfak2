# Multi-timepoint fits

Two timepoints are the minimum input of alfak2, not its limit. This study
extends `fit_sparse_joint()` to one population sampled at two or more times,
checks that two-timepoint fits are unchanged, and measures what the extra
samples give and what passaging drift between them costs. The plan and its
status are in the repository's local handoff
(`HANDOFF_CLAUDE_2026-11-13_multi_timepoint.md`, git-ignored).

Frozen baseline: ff13cf4 (`dev` HEAD on 2026-11-13), installed as
`results/lib/alfak2_ff13cf4` from `results/src/alfak2_ff13cf4.tar.gz` (`git
archive`, sha256 beside it).

## Stage 0 and the regression check (`reference.R`)

`reference.R --stage=fit --build=B` fits, with the library
`results/lib/alfak2_B`:

- `dev`: `fit_sparse_joint()` at its defaults on patient 1 of the 12
  development landscapes of `benchmark/sparse_two_timepoint`, in its 8 core
  and 2 short-interval scenarios (120 fits);
- `entry`: `fit_alfak2()` on two of them; `free`: `initial_state = "free"` on
  the 12 typical ones; `rates`: `fit_sparse_joint_rates()` on two;
- `cohort`: `fit_sparse_joint_cohort()` of the 8 patients of two landscapes,
  free and balance;
- `compat`: `as_alfakR()` and `alfak_compat()` on four;
- `gbm`: the 20 GBM inputs at birth rate sqrt(0.02) per day. These are
  patient data: they are read in place and the fits are written to
  `~/Documents/alfak2_private/multi_timepoint/reference/B/gbm`, outside the
  repository.

Each item keeps the fit's numbers and its predictions at fixed karyotypes
(the benchmark's scored nodes; for GBM up to 500 one-step and 500 two-step
karyotypes). `reference.R --stage=compare --builds=A,B` compares every numeric
leaf, relative to its largest absolute value, and requires every other leaf to
be identical. The pass criterion is a relative difference of at most 1e-8.

## Stage 1: the model

`fit_sparse_joint()` takes a count matrix with one column per timepoint, at
least two, with `dt` as the intervals between consecutive timepoints or
`times`. With more than two timepoints it uses `src/sparse_joint_series_model.h`:
the first-order closed form of `src/sparse_joint_model.h`, started at the first
sample and evaluated at every later sample time, each later sample with its
own depth offset. Two timepoints keep `src/sparse_joint_model.h`. Majors are
judged by their share of any sample; a karyotype first seen later starts from
its frequency at the earliest timepoint it was seen; the history of the
balance prior is at least the time from the first to the last sample. Cohort
fits, `as_alfakR()`, `alfak_compat()` and the two older pipelines take two
timepoints and say so.

Build `mt1` is this code (`results/src/alfak2_mt1.tar.gz`, with the diff
against ff13cf4). Two things changed later:

- From cv4, `as_alfakR()` and `alfak_compat()` take more than two
  timepoints. The level is log(nb/n0) over the first interval.
- From cv6, the majors and the domain come by default from the first and the
  last sample (stage 2).

Cohort fits and the two older pipelines still take only two timepoints.

## Stage C: convergence of the scale optimisation

`validation/sparse_joint/convergence/` (diagnosis, evaluation) and
`validation/sparse_joint/model_m/` (a model change tried and stopped). Build
`cv3` = mt1 with the convergence check: the fits are those of mt1 (and of
ff13cf4 on two timepoints); `diagnostics$optimisation` says whether the scale
optimisation was accepted and, if not, why and along which latent variables
the Laplace approximation is flat. Build cv1 also searched the scales again
and chose the history among accepted candidates; that is off from cv2 on.
cv2 still searched again when nlminb itself failed, which changed 3 of 144
reference fits (`convergence/results/ff13cf4_cv2/evaluate_partial/`); cv3
keeps the prior-centred scales there, as ff13cf4 did.

## Stage 2: simulation study

Arms `all` and `ends` with build cv3 (`setup.R`), `zero`. Design in the
handoff. On HPC each (arm, dataset) is one array task:

```sh
hpc/submit.sh fit cv3 validation/sparse_joint/multi_timepoint/hpc/tasks/<list>.tsv MEM TIME
```

### Results (2026-11-14)

Tables: `results/series/report.md` (`report.R`). 288 datasets, 12 landscapes,
every fit ok. RMSE ratios to `ends` paired by landscape, medians over cells by
bottleneck (d0 / d1 / d2), and cells within the plan's gate (every distance at
most 1.03):

| arm | build | 200 | 2000 | 20000 | gate |
| --- | --- | --- | --- | --- | --- |
| `all` | cv3 | 1.31 / 1.44 / 1.35 | 1.09 / 1.03 / 1.01 | 0.96 / 0.96 / 0.94 | 7 / 24 |
| `majors1` | cv5 | 1.00 / 1.01 / 1.01 | 0.99 / 1.00 / 1.00 | 1.01 / 1.00 / 1.00 | 18 / 24 |
| `majors` | cv5 | 0.97 / 1.00 / 1.01 | 0.99 / 1.00 / 1.00 | 1.01 / 1.00 / 1.00 | 21 / 24 |

- The pilot landscape (mt_L05) was unfavourable to `all` at every
  bottleneck. There the loss sat on karyotypes seen only between the first
  and the last sample (fitness overestimated, RMSE ratio 1.24-1.50), and at a
  bottleneck of 2e5 (`drift_check.R`, near the ALFA-K paper's ABM) it was
  still 1.13: not passaging drift alone. Masking the rare karyotypes'
  intermediate counts (`majors_any`) did not help, because majors chosen from
  every sample include clones that rise and fall in between; choosing majors
  from the first and the last sample did (`majors1`).
- `majors` (the default from cv6) weights the majors' intermediate counts by
  1 / max(1, dispersion) from a first fit; it fails the gate in 3 cells
  (1.10 at one step, bottleneck 200, 8 samples of 200 cells; 1.12 and 1.08 at
  d0, bottleneck 20000, 8 samples of 1000 cells). Median coverage per cell
  and distance 0.96-1.00 for `majors`, `majors1` and `ends`, 0.92-1.00 for
  `all`. Median time per fit: 2-4 times `ends`; `all` 3.5-21 times.
- Decision (on my recommendation, under the user's authorisation): default
  `intermediate = "majors"`; `"all"` an option for large bottlenecks; the
  drift model (stage 3) not built, since coverage and failures, its pass
  criteria, are already met and the remaining gaps are not where drift is
  largest.

## Samples in between that carry no information (cv7, cv8)

Found by the ALFA-K paper benchmark on the Salehi trajectories, in round
21aedca. Of 205 trajectories with more than two passages, 61 could not be fitted
on all samples:
- 43 stopped only because a passage in between had no cell in the domain:
  clonal sweeps, with the population in the middle two or more steps from both
  ends' majors.
- 18 also stopped at an end.

Ten more had a passage in between with cells on non-major neighbours only.
With the default, that sample's depth offset was degenerate.

The depth offset of a later sample has a flat prior. A sample with no cell on
the states it counts for (the majors, or the domain with `"all"`) therefore
has its likelihood largest as the offset goes to minus infinity, where the
sample adds nothing. The Laplace integral over the offset is improper. From
cv7 such a sample is left out with a warning (`diagnostics$dropped_samples`),
which is that limit. The ends are still required. `as_alfakR()` keeps
alfakR's karyotypes from every sample (`fit$input`).

The convergence check also reads an inner Hessian with non-finite entries as
near-singular instead of stopping. Under 21aedca it stopped on
SA609_X10_l_10_d1_0_d2_0 with all 10 passages. With cv8 that trajectory
leaves out passages 2 and 3, with 243 and 171 cells, 4 and 5 of them in the
domain, none on a major. The fit is accepted, and `as_alfakR()`'s frequent
karyotypes are alfakR's.

Checks:
- tests 141/141 (cv8);
- HPC, cv7 and cv8 each against ff13cf4 and cv6: 140/140 reference fits
  identical (development, entry, free, rates, compat; cohort fits do not use
  the changed code), and cv8's convergence verdicts equal cv6's;
- none of the 288 stage-2 datasets leaves a sample out;
- six stage-2 datasets refitted on the Mac with cv6 and cv7 (default): identical
  predictions, objectives and histories.

  Against the HPC fits of the same arm (cv5), objectives agree to 1e-8 or
  better, except in mt_L06_b00200_s06_n0200_k8. There two candidate histories
  the check did not accept (300, 1000) differ by up to 2.7 between the
  platforms, so the chosen history differs (HPC 300, Mac 1000) and so does the
  prediction (up to 0.047). This is the path dependence of degenerate
  candidates described in stage C, already in cv6, and not changed here.

A task list has one line per fit, `study<TAB>method<TAB>dataset_id`, with
study `series` (this study) or `patient` (the patient regrowth study, which
takes the build from the submission). `hpc/sync.sh` copies the code and the
simulated inputs of both studies; fits come back with rsync.
