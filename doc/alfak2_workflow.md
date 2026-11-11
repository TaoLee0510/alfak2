# alfak2 workflow, step by step

This is how the default pipeline of alfak2 runs, from a two-sample count matrix
to a fitted landscape and, optionally, to the output files of alfakR. The
mathematics is in [`alfak2_method.tex`](alfak2_method.tex), the sparse joint
count model section; the design decisions and benchmark results are in
[`validation/sparse_joint/README.md`](../validation/sparse_joint/README.md).

## Entry points

All three run the same model, the sparse joint count model.

| Entry point | Use | Returns |
| --- | --- | --- |
| `fit_alfak2(counts, dt, beta)` | Main entry point; `pipeline = "sparse_joint"` is the default | An `alfak2_sparse_joint_fit` |
| `fit_sparse_joint(counts, dt, beta, ...)` | The same model with its arguments given directly | The same |
| `alfak_compat(yi, outdir, ...)` | The arguments of `alfakR::alfak()`; writes alfakR's files | The cross-validation R-squared, invisibly, and the files in `outdir` |

`fit_alfak2()` passes `counts`, `dt`, `beta` and the list
`sparse_joint_control` on to `fit_sparse_joint()`. Any argument that belongs
to the earlier pipelines is an error unless `pipeline = "two_layer"` or
`"joint_count"` is set.

## Part 1. Fitting (`fit_sparse_joint`)

1. **Prepare the input.** A two-column matrix of integer counts, one column
   per sample, with karyotype labels such as `"2.2.3.2"` as row names. Rows
   with no cells are dropped. `dt` is the time between the samples (days for
   the GBM inputs) and `beta` the missegregation probability per chromosome
   copy and division. Copy numbers must lie between `min_cn` (default 1) and
   `max_cn` (default one above the largest observed), and both samples need
   cells.

2. **Find the major karyotypes and the domain.**
   - A karyotype is *major* if it holds at least 1% of either sample and at
     least 3 cells over both. Without any, the most abundant one is used.
   - The *domain* is the majors and every one-step neighbour of a major (one
     copy more or fewer of one chromosome), observed or not. An unobserved
     neighbour enters the likelihood because its absence next to an abundant
     clone is evidence that it is less fit.
   - Observed karyotypes more than one step from every major are left out of
     the likelihood together with their cells. They would enter only because
     they were observed, while their unobserved siblings would not, which
     biases their fitness upwards. They are predicted later at the level of
     the observed non-major karyotypes.
   - If the domain would exceed `max_states` (5000), the share required of a
     major is raised until it fits.

3. **Build the division kernel.** For every pair of neighbouring domain
   states, the probability that a daughter of one has the other's karyotype,
   from the exact single-chromosome missegregation kernel (the same as
   alfakR's `pij`), and for every state the probability `P_ii` that a
   daughter keeps the parent's karyotype.

4. **Build the trend design.** The centre is the count-weighted mean
   karyotype. Each state's covariates are an intercept, each chromosome's
   deviation from the centre (dosage), its square (curvature), and the
   distance to the nearest major raised to the power 0.85 (the shell term).

5. **Set the prior on initial abundances.**
   - Seen in the first sample: centred on its count, essentially free.
   - Seen only in the second: centred on its frequency there.
   - Never seen: a frequency near 1e-5 with a wide spread.
   - Default `initial_state = "balance"`: every non-major state is centred
     instead on its mutation-selection balance with the majors next to it,
     what they would have fed it over a history `h` (by default the interval
     between the samples), with spread `balance_sd = 1` on the log scale. Its
     count in the first sample then informs its fitness. Majors stay free.

6. **Assemble the count model (compiled with TMB).**
   - Latent variables: per state, the log expected count at the first sample
     (zeta) and the log birth rate (g); the trend coefficients; the depth
     ratio (theta).
   - Fitness prior: g is centred on the trend with a residual scale; the
     intercept on the log of the supplied `birth_rate` (sd 0.05); dosage and
     curvature coefficients have their own scales; the shell coefficient has
     prior sd 0.5.
   - Net growth rate: d = (2 P_ii - 1) exp(g).
   - Dynamics in closed form: at the second sample each state has grown at
     its own rate and received influx from its one-step neighbours, each
     source growing at its own rate (first order in missegregation).
   - Likelihood: the counts of each sample are Poisson, with expectation
     exp(zeta) at the first and the dynamics, scaled by theta, at the second.

7. **Estimate the scale hyperparameters.** The dosage, curvature and residual
   scales maximise the Laplace-approximated marginal likelihood (`nlminb`
   within [1e-4, 0.5]). For each trial value an inner Newton solve finds the
   posterior mode of all latent variables. If the optimisation fails, the
   prior-centred scales (about 1% per step) are used and `convergence = 99`
   is recorded.

8. **Compute the posterior.** The posterior mode at those scales, with the
   covariance taken as the inverse of the *expected* information: the Hessian
   with the counts replaced by their fitted expectations, which is positive
   definite where the observed Hessian may not be. It is factorised once
   (sparse Cholesky).

9. **Report relative fitness and its uncertainty.**
   - `fitness_mean`: net growth rate per unit of `dt`; its absolute level
     follows the assumed birth rate.
   - `relative_fitness`: minus the count-weighted mean over the majors, which
     is what two samples identify.
   - `relative_sd_model`: the model's own posterior standard deviation, too
     small on its own.
   - `relative_sd`: calibrated,
     `sqrt(model_sd^2 + (0.6 * relative)^2 + (0.009 * birth_rate * (1 + shell))^2)`.
     The constants were fitted on the development landscapes; held-out 95%
     coverage was 0.94-0.98.

10. **Store the fit.** `summary` (one row per domain state), `trend`
    (coefficients and covariance, centre, majors, scales, observed level,
    reference, calibration), `settings` and `posterior` (all settings and
    the posterior mode, from which the posterior can be rebuilt) and
    `diagnostics`.

## Part 2. Prediction (`predict(fit, karyotypes)`)

| Karyotype | Prediction |
| --- | --- |
| In the domain | Its posterior from the fit |
| Observed but outside the domain | The level of the observed non-major karyotypes |
| Anything else, e.g. two steps away | The trend, with the residual scale added to its variance |

Every row has `relative_fitness`, `relative_sd` and `relative_sd_model`.

## Part 3. alfakR outputs (`alfak_compat` = fit + `as_alfakR` + `write_alfakR_outputs`)

11. **Read alfakR's input.** `yi$x` must have two columns. The interval is
    `diff(passage_times)` or `diff(as.numeric(colnames(yi$x))) * yi$dt`, as
    in alfakR. The model of Part 1 is fitted with `beta = pm`; further
    arguments, such as `birth_rate`, go to `fit_sparse_joint()`.

12. **Choose karyotypes as alfakR does.** `fq`: at least `minobs` cells over
    both samples. `nn`: their one-step neighbours without nullisomy. The
    landscape covers `fq`, `nn` and the neighbours of `nn`, in alfakR's
    order. `minobs` changes only this choice, not the fit.

13. **Put fitness on alfakR's scale.** `mean` is relative fitness plus
    `log(nb / n0) / dt`, the growth rate alfakR assigns to the population.
    With `correct_efflux = TRUE` it is divided by `2 P_ii - 1`, giving birth
    rates.

14. **Draw from the posterior in place of the bootstrap.** The posterior is
    rebuilt from the stored settings and mode and all latent variables are
    drawn jointly. Draws are mapped to relative fitness (domain states
    through g, trend predictions through the coefficients plus the residual,
    with the reference drawn alongside). Each karyotype receives independent
    noise for the part of its calibrated standard deviation the model does
    not explain. Draws come in mirrored pairs and each row is scaled to the
    calibrated standard deviation, so every row's mean and median equal the
    prediction and its standard deviation equals `relative_sd`, while the
    model's correlations between karyotypes are kept.

15. **Assemble the files.**
    - `landscape_posterior_samples.Rds`: one row per landscape karyotype,
      one column per draw.
    - `landscape.Rds`: `k`, `mean`, `median`, `sd`, `fq`, `nn`, the row
      summaries of the draws.
    - `bootstrap_res.Rds`: `final_fitness` (draws for `fq`), `nn_fitness`
      (draws for `nn`), initial frequencies from the drawn initial
      abundances; `initial_*` equal `final_*`; alfakR's retry attributes.

16. **Cross-validate.** alfakR's folds: each `fq` karyotype with those of its
    one-step neighbours not in an earlier fold. For each fold, the counts of
    its domain states leave the likelihood, they get the prior of a
    never-seen karyotype, the trend's centre and majors are recomputed
    without them, and the posterior mode is found again with the scales
    fixed. The predictions, aligned to the full fit on the majors outside
    the fold, are compared with the full fit:
    `R2R = 1 - sum((p - o)^2) / sum(o^2)` with both series centred. This
    measures how well a karyotype the model has not seen is predicted, not
    how accurate the full fit is, and is not comparable with alfakR's R2R.
    `xval.Rds` holds `R2R`, `predictions`, `observations`, `xval_data` and
    `tmp`.

17. **Write and return.** The four alfakR files, `alfak2_fit.Rds`, and
    `landscape_data.Rds` when `landscape_data_output = TRUE`; R2R is
    returned invisibly.

## Variants

- `fit_sparse_joint_rates()`: when the birth rate is unknown, fits at a few
  rates across a range (default three, evenly spaced on the log scale). Point
  predictions are at the lowest rate, intervals cover the whole range, and
  the agreement of rankings across rates is reported.
- `fit_sparse_joint_cohort()`: several samples at once. Each keeps its domain
  and initial abundances; a karyotype shares its fitness across samples with
  an estimated deviation; the trend is shared or per sample, whichever has
  the higher marginal likelihood. Its default is `initial_state = "free"`.

## Running times

| Case | Time |
| --- | --- |
| One simulated sample | Seconds to about 4 minutes |
| One GBM sample at three birth rates | Median 149 s |
| `alfak_compat()` with cross-validation, GBM P6 (226 karyotypes) | About 8 minutes |

## Relation to the earlier pipelines

| Pipeline | Stages | Computation |
| --- | --- | --- |
| `"sparse_joint"` (default) | One count model on the majors and their neighbours; trend beyond | TMB model `sparse_joint` (`sparse_joint_cohort` for cohorts) |
| `"two_layer"` | Local hierarchical count model, then one of nine extrapolation methods | Local layer: TMB model `local_model`. Second layer: R and RcppEigen, not TMB |
| `"joint_count"` | One count model with a structured landscape prior on the bounded target graph | TMB model `local_model` with the structured prior |

All TMB models are compiled into one library and selected by name. Graph
construction, the missegregation kernel, the graph Gaussian solver and the
simulators are separate Rcpp code.
