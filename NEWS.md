# alfak2 0.2.0 (2026-11-17)

The first release built around the sparse joint count model: fitness
landscapes from sparse karyotype counts at two time points, the case the
package is meant for, with more time points as an extension.

## The default pipeline

- `fit_alfak2()` now fits the sparse joint count model (`fit_sparse_joint()`)
  by default. Both time points enter one count likelihood over the
  well-sampled ("major") karyotypes and all their one-step neighbours,
  observed or not, so that an unobserved neighbour of an abundant clone is
  evidence that it is less fit. The earlier two-layer pipeline remains
  available (`pipeline = "two_layer"`).
- Rare karyotypes start at mutation-selection balance with the clones next to
  them (`initial_state = "balance"`, the default), so that their level at the
  first time point informs their fitness. The time over which the balance
  built up (`history`) is chosen by marginal likelihood among the interval
  and 30 to 3000 divisions longer than it.
- Relative fitness comes with calibrated uncertainty (`relative_sd`), with a
  larger deficit for karyotypes seen at no time point.
- Fitness away from the domain is predicted from a fitted trend in chromosome
  dosage, curvature and distance to the major karyotypes (`predict()`).
- `diagnostics$domain_coverage` reports how much of each time point the domain
  holds; a time point with no cell in it stops the fit with an error of class
  `alfak2_domain_coverage_error`. `check_sparse_joint_counts()` summarises a
  count matrix before fitting.
- The scale optimisation is checked and a fit that is not accepted warns
  (`diagnostics$optimisation`).

## Extensions

- **More than two time points** (`times =`): `series = "anchored"` (default),
  `"pairwise"` or `"state_space"`. Two time points give the same fit with all
  three. Time points in between that carry no information under the model are
  left out with a warning (`diagnostics$dropped_samples`).
- **Cohorts**: `fit_sparse_joint_cohort()` fits several datasets together,
  with the trend shared or per dataset by marginal likelihood.
- **Uncertain rates**: `fit_sparse_joint_rates()` fits over a range of
  baseline birth rates or missegregation rates (`beta_range`), with point
  predictions at the lowest rate of the range.
- **alfakR outputs**: `as_alfakR()`, `write_alfakR_outputs()` and
  `alfak_compat()` produce the objects and files of `alfakR::alfak()` (branch
  dev), with alfakR's karyotype sets, scale, bootstrap layout and
  cross-validation folds.
- `initial_state = "free"` estimates every initial abundance independently;
  `initial_state = "auto"` (two time points) chooses between the two by the
  size of the cloud of rare karyotypes. `"auto"` misfires when the time points
  are one passage apart and is not the default.

## Validation

Scripts and reports are in `validation/sparse_joint/`; the method and its
results are in `doc/alfak2_method.tex`.

- On the development benchmark of simulated tumours (sparse two-time-point
  scenarios, 180–1080 days) and its held-out landscapes, the sparse joint
  model lowered the error of predicting no fitness differences to 0.47–0.74
  of that value, where alfakR (minobs 20) exceeded it.
- On the ALFA-K paper's own tests (benchmark round of alfak2 348addd, against
  alfakR dev2 f93b2d2), alfak2's landscapes were more accurate on the paper's
  simulations, and its forecasts beat the baseline more often.
- On the Salehi passages, forecasts scored against every later passage by a
  proper score (log-likelihood relative to no change): alfak2 beat alfakR by
  1.01 per cell with two time points and by 1.12 with every passage
  (`validation/sparse_joint/salehi_horizon`). The paper's Fig. 2d ranking on
  these data is dominated by sampling noise and does not separate the two.

## Known limitations

- A missegregation rate other than the assumed one, or miscalled cells,
  cannot be told from the counts of two time points; fit over a range of
  rates when the rate is uncertain.
- Without a known division rate, rankings are robust but per-day magnitudes
  are conditional on the assumed rate.
- Intervals are too narrow after 1080 days and when the landscape changes at
  the first time point; unobserved neighbours are predicted slightly too fit.
- See "Open" in `validation/sparse_joint/README.md` for the full list.
