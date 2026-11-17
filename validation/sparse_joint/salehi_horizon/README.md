# Salehi forecasts against later passages

Started 2026-11-16 at the user's request: use the Salehi passages as
two-time-point inputs, (a, b), (a, c), (a, d), ..., and the passages that
follow each input as the truth for its forecast. The ALFA-K paper's Fig. 2d
scores only the next passage, by a rank correlation that sampling noise
dominates (`../cloud_level/README.md`, step 2). Here every later passage is a
test, and the score is a proper one.

## Design

- **Inputs and fits.** The "ends" fits of the paper benchmark's round 05a8855,
  read only: alfak2 at its default (`initial_state = "balance"`, which is the
  default again), alfakR (alfakR_dev2 f93b2d2) and, from the
  alfak2-salehi-variants project, alfak2 with `initial_state = "free"`;
  minobs 5, 10 and 20. No fit is repeated.
- **Forecasts.** From the trajectory's last passage, as the paper benchmark
  forecasts (`forecast_salehi_dir()`: `alfakR::predict_evo`, ABM, 5
  replicates, missegregation 5e-5, population 1e7), but up to the furthest
  descendant of that passage within four passages (5 days per passage for
  SA906 and SA039, 15 for the others).
- **Truth.** Every descendant passage, from the metadata's parent links; its
  counts are the last column of a trajectory that ends with it.
- **Scores**, on the karyotypes of the fit's landscape (the same set for all
  three methods), for each descendant:
  - the mean log-likelihood per cell of its counts under the forecast,
    minus that under the parent's frequencies (persistence), each mixed with
    a uniform distribution over the landscape with weight 0.01 (and 0.001).
    Sampling noise of the later passage enters both alike, so a forecast
    gains only where it puts probability on where the cells go;
  - whether the forecast is closer to it than the parent in Wasserstein
    distance;
  - the distance between two random halves of it (its sampling noise).
- **Comparison.** On the rows that the compared methods share (trajectory,
  minobs, later passage), paired; 95% intervals from a bootstrap over the six
  lineage types and over the samples.

## Files

- `make_tasks.R`, `forecast_score.R`: the task table and one task (one fit),
  run on HPC in the paper benchmark's image (`array_task.sh`, not tracked).
- `analyse.R`: the tables and figure in `results/report/` (not tracked).

## Results (2026-11-16)

All 2493 tasks finished (array 21981908): 396 alfak2, 402 alfakR and 396
alfak2 "free" fits had later passages; the others had no fit (inputs the
benchmark excluded, or minobs without a fit) or no descendant. 986 rows (fit ×
later passage) on 143 trajectories are shared by all three methods.

**alfak2 forecasts the later passages better than alfakR at every horizon.**
Log score per cell over persistence (mean; lineage-type bootstrap):

| Passages ahead | Rows (trajectories) | alfak2 | alfakR | alfak2 free | alfak2 − alfakR | Rows where alfak2 is better |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 471 (143) | +0.47 (0.17, 0.67) | −0.16 (−0.76, 0.22) | +0.44 | +0.63 (0.35, 0.99) | 65% |
| 2 | 280 (91) | +0.27 (−0.44, 0.79) | −1.10 (−2.17, −0.20) | +0.17 | +1.37 (0.90, 2.04) | 70% |
| 3–4 | 235 (53) | −0.81 (−1.96, 0.43) | −2.14 (−3.53, −0.73) | −0.87 | +1.33 (1.10, 1.63) | 71% |
| All | 986 (143) | | | | +1.01 (0.67, 1.40); samples (0.71, 1.26) | 68% |

- The difference holds at every minobs (5: +1.53, 10: +1.02, 20: +0.26, all
  intervals above zero), whether or not the treatment changed between the
  parent and the later passage (+0.67 and +1.05), and for inputs of one
  interval (+1.36) and of more (+0.81).
- One passage ahead, alfak2's forecast beats persistence (interval above
  zero); alfakR's does not. Two or more passages ahead neither beats
  persistence; alfakR falls below it (intervals below zero), alfak2 does not
  significantly.
- In Wasserstein distance one passage ahead, alfak2's forecast is closer to
  the later passage than the parent in 53% of rows, alfakR's in 21%. Both
  errors are close to the sampling noise of the later passage (median excess
  0.20 for alfak2, 0.34 for alfakR, 0.21 for the parent).
- alfak2 with `initial_state = "free"` is no better than the default
  (−0.06 per cell, interval −0.16 to 0.09; one passage ahead −0.04, interval
  −0.07 to 0.02; it is worse in 62–67% of rows): on these passages the
  balance prior does not cost accuracy, which supports keeping it as the
  default.
- With a mixture weight of 0.001 instead of 0.01 the ordering is the same.

**Reading.** alfakR's forecasts move the population further than it moves
(its fitness differences are larger: `../cloud_level/README.md` and the
paper benchmark), and a proper score penalises probability put where the
cells do not go. Neither method forecasts the composition two or more passages
ahead better than assuming no change; alfak2 is closer to that baseline. This
is the comparison on real data that sampling noise does not decide: the
forecast's probability of the cells actually seen.

![Forecasts against later passages](results/report/figures/horizon.png)
(figure in `results/report/figures/`, not tracked)
