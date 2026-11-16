# Three methods for more than two samples

The user approved this plan on 2026-11-14. It follows the multi-timepoint
work in `../multi_timepoint`.

The method of 21aedca/8512694 (`series = "anchored"`) takes the majors and
the domain from the first and the last sample. In the Salehi trajectories of
the ALFA-K paper benchmark, clonal sweeps put intermediate passages two or
more steps from those majors. Such a sample carries no information under
that design and is left out. This is a limit of the design, not of the data.

Two further methods model every pair of consecutive samples on its own, with
its own majors, domain and start. Evolution is continuous: one landscape for
all pairs. All three methods are kept. The default is chosen after the
evaluation below and the paper benchmark.

## The methods (`fit_sparse_joint(series = ...)`)

| | `anchored` | `pairwise` (A) | `state_space` (B) |
| --- | --- | --- | --- |
| dynamics | one trajectory from the first sample | each pair from its earlier sample as observed | latent abundances at every sample, linked by the dynamics plus drift |
| majors, domain | first and last sample | each pair's own; states are their union | each pair's own; states are their union |
| samples in between | majors only, weighted by 1/dispersion | end of one pair and start of the next | each counted once |
| drift | none (down-weighting) | none (restart at each sample) | estimated scale |
| template | `sparse_joint_series_model.h` | `sparse_joint_pairwise_model.h` | `sparse_joint_state_space_model.h` |

With two samples all three give the two-sample fit, bit for bit.

**Pairs without coverage.** For `pairwise` and `state_space`, both samples of
each pair need a cell in that pair's domain. In c7f809c a pair without one
stopped the fit. The `benchmark` session found that this stopped 42 of the
187 Salehi trajectories with more than two passages that `anchored` can fit,
all of them sweep families. These are the trajectories the two methods are
for.

Such a sample carries no information about its pair, as an anchored sample
in between with no count carries none. So it is now left out
(`sparse_joint_pair_left_out()`, from the counts alone):
- A sample in between with no cell in its pair's domain is left out, and its
  neighbours become a pair.
- If the first or the last sample has none, its neighbour is left out instead.
- This repeats until every pair is covered; at worst the two-sample fit of the
  ends remains.

On the Salehi inputs (default settings, `min_cn = 0`):
- none of the 187 stops;
- 42 leave out one sample each, the uncovered sample itself;
- 4 of those 42 have three passages and become two-sample fits.

For comparison, `anchored` left out 103 samples in 53 of them. Fits in which
every pair is covered are unchanged.

### Pairwise (A)

Each pair is modelled as two samples are (`sparse_joint_model.h`):
- its initial abundances at the earlier sample, with balance for the pair's
  non-majors;
- the first-order closed form over the interval;
- its own depth.

All pairs share the birth rates and the trend. A pair starts at balance over
the history plus the population's age at its first sample (`"age"`, the
default). The alternative `"same"` uses the history alone; it is an internal
option, `alfak2.pair_history`, compared as `pairwise_same`.

The trend is centred on the count-weighted mean karyotype over all samples.

Check: given the cohort's centre instead (the mean of the pairs' centres),
the method equals `fit_sparse_joint_cohort()` on the pairs, with full sharing,
a shared trend and a free start. The objective is 436.869260572 in both, and
fitness agrees to 6e-15.

### State space (B)

Latent log abundances at every sample. From one sample to the next they
follow the closed form from the earlier sample's latent abundances, plus a
Gaussian deviation.
- The deviation's variance is drift × interval / frequency, capped smoothly
  at 4.
- The frequency is the mean share of the state at the two samples plus half a
  cell.
- The drift is a fifth scale. Its prior is log-normal around 1e-4 per day
  (sd 2.5), and its bounds are 1e-9 to 1.

Each sample's counts enter once.

With the drift fixed at 1e-8, the fit follows the anchored model with every
sample used for every karyotype: correlation 0.997, majors within 1e-4.
They differ by the second-order terms of the mutational flux, which the
closed form composed over intervals keeps.

### Checks on a small example (`hmc_check.R`, build cv10)

The example is `simulate_series_counts()` of the package tests: 4 samples,
21 states. The joint posterior was sampled by HMC at the fitted scales and
compared with what each fit reports (mode, and sd from the expected
information).

| method | HMC sd / model sd | (mean − mode) / model sd, majors | same, other observed (median, max) |
| --- | --- | --- | --- |
| anchored | 0.94–1.19 | ≤ 0.06 | 0.50, 0.70 |
| pairwise | 0.80–1.20 | ≤ 0.10 | 0.43, 0.58 |
| state_space | 0.96–1.07 | ≤ 0.13 | 0.75, 0.96 |

## Simulated sets

All three sets use the development simulator and the landscapes of
`../multi_timepoint` (seed base 97000; held-out 98000 sealed).

- **e1**: the stage-2 passaged cultures, linked, not copied. 288 datasets:
  - first sample at day 1440;
  - 4 or 8 samples, 6 or 60 days apart;
  - bottlenecks of 200, 2000 or 20000 cells;
  - 200 or 1000 cells per sample.
- **e2**: clonal sweeps. 192 datasets:
  - first sample at day 240 or 480, when clones still replace one another;
  - samples 120 days apart;
  - bottlenecks of 2000 or 20000 cells;
  - otherwise as e1.

  Gate: does it reproduce what the Salehi trajectories show? With 8 samples,
  42–62% of datasets have a sample in between with no cell on the ends'
  majors, and 4–38% one with no cell in the anchored domain. For comparison,
  53 of 205 Salehi trajectories had such a sample. With 4 samples (360 days)
  there are none.
- **e3**: neutral. The e1 design, but from the first sample on every viable
  karyotype divides at the start population's mean rate. This uses
  `abm_constant.cpp`, a copy of the simulator with that one switch.
  Every fitness difference a method reports is a false positive. 288
  datasets.

## Measures

Every method is scored on the karyotypes of stage 2: observed ones, and one
and two steps away. The reference is the same for every method.
- **Accuracy**: RMSE of relative fitness against the truth, as a ratio to
  `ends` (the first and last sample alone), paired by dataset. Also 95%
  coverage.
- **Power** at effect e (0.02, 0.05, 0.1 per day): among karyotypes with
  |truth| ≥ e, the share whose 95% interval excludes zero on the correct
  side.
- **False positives**: the share whose interval excludes zero:
  - among karyotypes with |truth| < 0.005;
  - every karyotype in e3.
- **Robustness**: failed fits, fits not accepted by the convergence check,
  time, peak memory.

The gates and the rule for the default are fixed in the header of
`report.R`. They were written before any fit of the new methods was scored.

Power and false positives as defined above turned out to be 0 for every
method (see Results). Both are therefore also measured at a controlled
false-positive rate:
- z = |estimate| / `relative_sd_model`, the model's own posterior sd;
- the threshold z* is the one that 5% of the karyotypes of e3 exceed, set per
  method and distance;
- power is the share of karyotypes with |truth| ≥ e that exceed z* with the
  right sign.

`report.R` writes this to `power_controlled.csv`.

## Files

- `setup.R`, `run.R`: datasets (`--stage=datasets`), e2's gate
  (`--stage=sweeps`), fits and scores, by set and build.
- `abm_constant.cpp`: the simulator with constant fitness (e3).
- `hmc_check.R`: the HMC check above.
- `report.R`: tables, gates and `overview.pdf`.
- `auto_eval.R`: the automatic choice (`series = "auto"`), scored on the
  existing fits (below); output in `results/<build>/auto/`.

Results are in `results/` (git-ignored).

On HPC, the task lists hold one line per fit:
`sm_<set><TAB>method<TAB>dataset_id`. They are run by
`../multi_timepoint/hpc/submit.sh fit BUILD TASKS MEM TIME`.

## Results (build cv10)

42 of the 768 `state_space` fits did not finish within the first 6 h job
limit, all with 8 samples (e1 8, e2 14, e3 20; 2 of the e2 ones also ran out
of memory at 8 GB). They were rerun with 21 h and 16 GB, and all finished (the
longest in 18.3 h). The numbers below include them.

### Accuracy

Median RMSE ratio to `ends`, paired by dataset. Each entry is the range over
distances 0, 1 and 2.

| set | setting | `anchored` | `pairwise` | `state_space` |
| --- | --- | --- | --- | --- |
| e1 | bottleneck 200 | 0.99–1.01 | 1.61–1.78 | 0.95–1.00 |
| e1 | bottleneck 2000 | 1.00 | 1.08–1.14 | 0.99–1.02 |
| e1 | bottleneck 20000 | 1.00–1.01 | 0.91–1.01 | 0.93–0.97 |
| e2 | 4 samples | 0.89–0.94 | 0.77–0.85 | 0.79–0.82 |
| e2 | 8 samples | 0.93 | 0.68–0.73 | 0.79–0.81 |
| e3 | bottleneck 200 | 0.89–1.01 | 1.35–1.68 | 1.03–1.15 |
| e3 | bottleneck 2000 | 0.83–1.02 | 1.14–1.72 | 0.99–1.23 |
| e3 | bottleneck 20000 | 0.75–1.00 | 1.24–2.01 | 0.81–1.17 |

- **`pairwise`** restarts from the observed sample at every pair. Passaging
  drift within each pair is then read as selection. It is worst with small
  bottlenecks and samples 6 days apart (up to 2.6 in one e1 cell). It is
  worst of all on neutral data. In clonal sweeps (e2) it is the most
  accurate, and with little drift it is close to `ends`.
- **`state_space`** is close to or better than `ends` on e1 and clearly
  better on e2. On neutral data it is worse at distances 1 and 2.
- **`pairwise_same`** (balance over the history alone) is within 0.06 of
  `pairwise` everywhere, and slightly worse on e1 and e3. `"age"` is kept.
- **Coverage** of the reported 95% intervals is 0.99–1.00 on e1 and e3 for
  every method. On e2, at distance 0, it is 0.80 for `ends`, 0.88 for
  `anchored`, 0.93 for `state_space` and 0.96 for `pairwise`.

### Power

**The reported 95% intervals cannot exclude zero.** `relative_sd` contains
deficit × |relative fitness|, with a deficit of 0.6 for observed karyotypes
and 1.1 for unobserved ones (`sparse_joint_default_calibration`). So
|estimate| / `relative_sd` is at most 1/0.6 = 1.67, or 0.91 for unobserved
karyotypes, which is below 1.96. Power and false positives by interval are
therefore 0 for every method and every set, by construction. The intervals
are meant for coverage, not for tests.

At a controlled 5% false-positive rate (Measures), effect 0.05:

| | `ends` | `anchored` | `pairwise` | `state_space` |
| --- | --- | --- | --- | --- |
| z*, all distances | 4.6 | 4.5 | 8.5 | 7.0 |
| e1, all distances | 0.18 | 0.17 | 0.015 | 0.04 |
| e1, distance 0 | 0.22 | 0.21 | 0.10 | 0.16 |
| e2, all distances | 0.30 | 0.45 | 0.13 | 0.37 |
| e2, distance 0 | 0.22 | 0.38 | 0.50 | 0.53 |
| e2, distance 1 | 0.33 | 0.47 | 0.13 | 0.37 |

- Even for `ends`, z* is about 4.5, not 1.96: the model's own sd is too small
  (`relative_sd_model`; this is why the calibration exists).
- For the two new methods it is far too small at distances 1 and 2 of
  neutral data: z* is 6.5–9 there, against 4.3–4.7 for the others. They pay
  for that in power at those distances.
- On the karyotypes observed in clonal sweeps, both new methods have more
  power than `anchored`.
- Direction (right sign at |truth| ≥ 0.05, no threshold) is 0.97–1.00 on e1
  for every method and 0.94–0.96 on e2.

### Robustness and cost

| | `ends` | `anchored` | `pairwise` | `state_space` |
| --- | --- | --- | --- | --- |
| failed fits | 0 | 0 | 0 | 0 (42 after a rerun with more time) |
| not accepted, e1 / e2 / e3 | 29% / 46% / 32% | 32% / 56% / 29% | 27% / 68% / 22% | 47% / 94% / 65% |
| median time, s (e1 / e2 / e3) | 73 / 126 / 85 | 196 / 444 / 244 | 344 / 789 / 468 | 888 / 2140 / 1664 |
| time per fit relative to `anchored`, 4 / 8 samples (median) | 0.4–0.5 / 0.2–0.25 | 1 | 1.3–1.7 / 1.6–2.0 | 2.6–4.0 / 7.7–12.4 |
| median peak memory, GB | 0.6–0.8 | 1.1–1.6 | 1.5–2.1 | 2.1–3.3 (up to 16) |

"Not accepted" means the convergence check did not accept the fit; such fits
are scored as they are, as for the earlier stages. For `state_space` the
history candidates differ in objective by 8–20 (median), so the search over
them is needed. Its cost is that of the model, not of the search.

### Gates and the default

| | e1 cells within 3% | e1 | e3 | failures | passes |
| --- | --- | --- | --- | --- | --- |
| `anchored` | 21 / 24 | pass | pass | pass | yes |
| `pairwise` | 5 / 24 | fail | pass | pass | no |
| `pairwise_same` | 4 / 24 | fail | pass | pass | no |
| `state_space` | 16 / 24 | fail | pass | pass | no |

The e3 gate counts interval calls, which are 0 by construction (above). As
written it cannot fail, so it says nothing. Measured at a controlled
false-positive rate, the new methods are the ones with too small an sd.

Only `anchored` passes, so it stays the default. The paper benchmark (E4, by
the `benchmark` session) is weighed before this is final. Both new methods
are kept:
- `pairwise` for clonal sweeps through large bottlenecks;
- `state_space` where clones sweep under moderate drift. It costs 3–4 times
  as much as `anchored` with 4 samples, and 8–12 times as much with 8, where
  some fits take up to 18 h and 16 GB.

Neither suits small bottlenecks or samples a few days apart.

### Automatic choice (`series = "auto"`)

`auto_eval.R` picks, for each dataset, one of the methods already fitted,
from measures of the counts and of the anchored fit, and scores the pick as
that method's fit. Nothing is refitted. The measures, the candidate rules and
the criterion were fixed before any rule was scored (script header):
- a rule is acceptable if it passes the e1 gate above and, on e3, is at most
  0.01 above `anchored` at every distance;
- among those, the lowest e2 ratio; within 0.01 of it, the simpler rule.

Two measures carry the choice:
- `ends_major`: the least share, over the samples in between, of a sample's
  cells on the majors of the first and last sample. A clonal sweep takes a
  sample off them.
- `weight`: the anchored fit's intermediate weight, 1 / max(1, dispersion of
  the majors' counts). It is low under strong passaging drift.

Chosen: `ends_major < 0.3` → `pairwise` if `weight >= 0.5` (or no sample in
between is left to measure it), otherwise `state_space`; `anchored` if
`ends_major >= 0.3`. This is `sparse_joint_auto_thresholds`.

| rule | e1 cells | e1 | e2 | e3 | picks on e2 (a / p / s) |
| --- | --- | --- | --- | --- | --- |
| chosen (0.3, weight 0.5) | 21 / 24 | 1.00 | 0.80–0.85 | 0.85–1.01 | 46% / 47% / 7% |
| 0.3, weight 0.8 | 21 / 24 | 1.00 | 0.81–0.83 | 0.85–1.01 | 46% / 32% / 22% |
| 0.3, always `pairwise` | 21 / 24 | 1.00 | 0.79–0.86 | 0.85–1.01 | 46% / 54% / 0% |
| 0.3, always `state_space` | 21 / 24 | 1.00 | 0.83–0.87 | 0.84–1.01 | 46% / 0% / 54% |
| `anchored` | 21 / 24 | 1.00 | 0.91–0.93 | 0.84–1.01 | |
| best fitted method per dataset | | 0.93–0.94 | 0.74–0.77 | 0.77–0.97 | |

Ratios are median RMSE ratios to `ends` over distances 0–2. The rule keeps
`anchored` on all of e1 and 95% of e3, and on e2 its median ratio (0.81) lies
between `anchored` (0.93) and the best fitted method per dataset (0.74).
- 0.3 is the largest share threshold that is acceptable. At 0.5 the rule
  leaves `anchored` on 18% of the neutral datasets and adds errors there
  (0.015–0.035 above `anchored` on e3).
- The weight threshold is not: 0.8 gives 0.82, just outside the 0.01 band.
- Choosing the thresholds on landscapes 1–6 and scoring on 7–12 (and back)
  says nothing: on a half, even `anchored` fails the e1 gate.

The paper benchmark (E4) is the held-out check. Until it passes, `anchored`
stays the default and `"auto"` must be requested.
