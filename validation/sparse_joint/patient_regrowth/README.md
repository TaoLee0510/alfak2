# Patient regrowth: few survivors and delayed growth

Patient data such as the GBM pairs are not passages. The cells of the first
time point are taken at surgery and those of the second at recurrence, which
can be regarded as grown from a small subset of the primary's cells; and the
time at which that subset started to grow does not correspond to the sampling times. This study asks how
much either costs the two-timepoint model as shipped (alfak2 ff13cf4), and
whether the remedies at hand help.

## Design

Data come from the development benchmark's 12 development landscapes
(wavelengths 0.4, 0.8, 1.2, 1.6, three each). Their truth populations
(`benchmark/results/sparse_two_timepoint/truth`) are read in place. For each
patient, the day-1440 population is sampled (500 cells) and either goes on
growing as a whole or is reduced to a number of survivors drawn at random;
these regrow in the benchmark's ABM and are sampled again (700 cells). The
model is always told that 360 days passed.

| scenario | survivors | growing share of the 360 days |
| --- | --- | --- |
| `none` | all | 1 |
| `s2000`, `s0200`, `s0020` | 2000, 200, 20 | 1 |
| `s0200_d50`, `s0200_d75`, `s0200_d875` | 200 | 0.5, 0.25, 0.125 |
| `s0020_d875` | 20 | 0.125 |

Before growth starts the survivors neither divide nor die. Seeds (base 99000)
are shared by the scenarios of a patient: the cells at the first time point are
the same in every scenario, survivors of one size are the same cells, and a dormant
scenario's recurrence is its immediate-regrowth scenario's population after
the shorter growth time. Patients 2-8 are generated for `none`, `s0200` and
`s0020`, for cohort fits.

Arms (all ff13cf4, copy numbers 1-8):

- `zero`: no differences.
- `base`: defaults, one division per day assumed (the populations divide about
  0.9 times per day while growing).
- `equivalent`: the birth rate averaged over the 360 days, growing share x 1;
  dormant scenarios only.
- `rates`: `fit_sparse_joint_rates()` over 0.1-1 per day, point predictions at
  the lowest rate.
- `cohort`: the 8 patients of a landscape fitted together, scored on patient
  1; `none`, `s0200`, `s0020` only.

Scores are the development benchmark's, by distance from the data, with one
change: the truth is per calendar day, the rate while growing times the
growing share, which is what a fit told 360 days estimates. Skill (1 - MSE /
MSE of predicting no difference), slope, rank correlation and 95% coverage then
compare across dormant shares; RMSE itself shrinks with the share. The study is
diagnostic: there is no pass/fail rule.

In this simulator dormancy is exactly a shorter growth time. The model is
invariant to rescaling time and rates together (every rate enters as rate x
time, the history candidates are counts of divisions, and the calibration
floor scales with the birth rate), so `equivalent` equals a fit told the true
growth time. `run.R --stage=symmetry` checks this on one dataset. What the
dormant scenarios measure is therefore the cost of not knowing the dormant
share.

## Running

```sh
R=validation/sparse_joint/patient_regrowth/run.R
Rscript $R --stage=link
Rscript $R --stage=datasets
Rscript $R --stage=fit --methods=zero,base
Rscript $R --stage=fit --methods=equivalent --datasets='_d[0-9]+$'
Rscript $R --stage=fit --methods=cohort --datasets='_(none|s0200|s0020)$'
Rscript $R --stage=fit --methods=rates
Rscript $R --stage=score
Rscript $R --stage=symmetry --dataset=dev_L01_P01_s0200_d875
Rscript validation/sparse_joint/patient_regrowth/report.R
```

The libraries are the ones built for the multi-timepoint work
(`../multi_timepoint/results/lib/alfak2_<build>`). The study was first run on
ff13cf4 (`zero` and `base`, in `results/`). The user then moved it to the
build with the convergence check, cv3, whose fits are those of ff13cf4 but
say whether the scale optimisation was accepted. `--build=cv3` writes to
`results/cv3/`, which links to the datasets, scored nodes, index and `zero`
fits of `results/`:

```sh
Rscript $R --stage=link --build=cv3
Rscript $R --stage=fit --build=cv3 --methods=base
# ... the other arms as above, with --build=cv3
Rscript $R --stage=score --build=cv3
Rscript validation/sparse_joint/patient_regrowth/report.R --build=cv3
```

On HPC each (arm, dataset) is one array task
(`../multi_timepoint/hpc/submit.sh fit cv3 <tasks>`, task lists in
`../multi_timepoint/hpc/tasks/`).

## Results (build cv3, 2026-11-14)

Tables: `results/cv3/report.md` (`report.R --build=cv3`). Medians over the 12
landscapes, patient 1. The `base` fits are those of ff13cf4; cv3 only labels
whether each scale optimisation was accepted (53 of 96 `base` fits were not).

- **Few survivors, immediate regrowth.** A small cost. Paired against `none`,
  the skill of the observed karyotypes with 5+ cells drops by 0.05 (2000
  survivors), 0.03 (200) and 0.12 (20); one and two steps away it is
  essentially unchanged (within 0.04). Coverage stays at 0.99-1.00.
- **Delayed growth, not known to the model.** The dominant risk. Told 360 days
  at one division per day when the population grew for half of them, the
  skill of the observed karyotypes falls from 0.71 to 0.11; for a quarter or
  less of them it is far below zero (-3 to -98): the deficits are spread over
  the whole interval and per day they come out too large (slopes 1.4 and 2.3
  for a quarter and an eighth with 200 survivors).
- **Telling the model the growing share** (`equivalent`: birth rate x growing
  share, equal to a fit told the true growth time) recovers most of it: one
  and two steps away 0.64-0.79 at every dormancy level with 200 survivors
  (0.47-0.62 with 20); the observed majors
  stay poor only when the population grew for an eighth of the interval
  (-1.6; 20 survivors: -22).
- **A range of rates** (`fit_sparse_joint_rates()` over 0.1-1, point at 0.1)
  protects against dormancy almost as well as knowing it (0.62-0.65 one and
  two steps away at 12.5% growth), but without dormancy it costs much of the
  skill one and two steps away (0.14-0.18 against 0.77-0.78): the lowest rate
  shrinks the deficits.
- **Fitting the 8 patients of a landscape together** (`cohort`) does not
  reduce the cost of a bottleneck: its skill is lower than the single fits'
  in most cells (by up to 0.11), higher only for the observed majors with 20
  survivors (+0.05).

Decision (2026-11-14, on my recommendation under the user's authorisation):
no change to the package. The README says, for patient data, that the
model must be told the time the population actually grew (or the birth rate
times the growing share), and that when this is unknown a range of rates is
the safer choice, at a cost when there was no delay.
