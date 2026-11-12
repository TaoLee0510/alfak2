# alfak2

`alfak2` implements **ALFA-K2: Adaptive Local-to-Global Fitness Landscapes for
Aneuploid Karyotypes**. It provides methods to infer fitness landscapes for
aneuploid karyotypes from sparse two-timepoint count data and to support
analyses of karyotype evolution. It combines bounded karyotype graph
construction, local hierarchical Bayesian inference, and a graph
Gaussian-process posterior with ordered copy-number epistasis penalties.

The default pipeline, `fit_alfak2()` or `fit_sparse_joint()`, is the sparse
joint count model, described [below](#sparse-joint-model-for-two-sparse-samples).
`alfak_compat()` runs it with the arguments of `alfakR::alfak()` and writes
alfakR's output files, so that analyses built on alfakR read its results
unchanged.

The earlier two-layer pipeline, `pipeline = "two_layer"`, is:

1. build a local karyotype graph around observed copy-number states,
2. fit a local hierarchical posterior for initial frequencies and node fitness
   values,
3. expand to a bounded graph and solve a graph Gaussian posterior with ordered
   copy-number epistasis penalties,
4. summarize or plot posterior fitness and uncertainty by support tier.

An explicitly selected `pipeline = "joint_count"` candidate fits the counts and
a structured landscape prior jointly on the bounded target graph. The current
accuracy validation compares these paths at every exact support distance.

## Method Components

- Local posterior inference for observed and nearby karyotype states.
- Bounded graph expansion over copy-number states connected by plausible
  karyotype transitions.
- Graph posterior smoothing with ordered copy-number epistasis penalties.
- Posterior summaries and plots for fitness estimates and uncertainty.
- Benchmark truth landscapes can be generated as lightweight paper-style GRF
  oracles with `simulate_grf_landscape()`. Fitness is computed on demand for
  requested karyotypes, so synthetic simulations can use 22 chromosomes without
  enumerating the full copy-number lattice.

## Minimal Example

```r
library(alfak2)

landscape <- simulate_grf_landscape(
  n_chr = 4,
  n_centroids = 30,
  lambda = 0.8,
  min_cn = 1,
  max_cn = 4,
  seed = 1
)

sim <- simulate_sparse_counts(
  landscape,
  beta = 0.00005,
  dt = 1,
  n0 = 200,
  n1 = 200,
  seed = 2,
  initial_population = 1000
)

fit <- fit_alfak2(
  sim$counts,
  dt = 1,
  beta = 0.00005,
  sparse_joint_control = list(min_cn = 1, max_cn = 4)
)

head(summarize_alfak2(fit))
predict(fit, c("2.2.2.2", "2.3.2.2"))

# the earlier two-layer pipeline
two_layer <- fit_alfak2(
  sim$counts,
  dt = 1,
  beta = 0.00005,
  pipeline = "two_layer",
  min_cn = 1,
  max_cn = 4,
  max_nodes = 5000
)
plot_alfak2(two_layer)
```

## Inference and biological dynamics

Local inference now defaults to Laplace marginalization of latent effects,
with gradient, conditional-gradient, Hessian, and scale-boundary checks. The
historical centered joint MAP remains available as `inference_method =
"joint_map"` for reproducible comparisons. Untrusted intervals are missing and
carry an explicit diagnostic status. A narrower interval alone is not evidence
of a more accurate fit.

`forward_model = "selection_transport"` retains the historical relative
net-growth model, with one transport operation after selection. The optional
`forward_model = "division_process"` instead models mutation at each division:
fitness is a positive birth rate per unit of `dt`, parameterized on the log
scale. It uses the full single-daughter kernel on retained graph nodes without
redistributing biological death or computational outflow. `forward_step`
controls a convergent positive splitting approximation; graph outflow and step
refinement should be checked for the intended beta, time interval, and domain.
The joint path checks the fitted-mode result after halving the step and refits
with finer steps when needed. Explicit `target_karyotypes` permit finite distant
queries in high dimensions; their omitted intermediate states remain reported
computational outflow and may require a larger buffer domain.

The joint candidate uses chromosome dosage effects, quadratic effects,
shrunk pairwise interactions, and a proper Gaussian graph residual. It updates
observed and unobserved fitness together, without treating a local posterior as
independent new data. Its default prior scales were selected on development
landscapes and fixed before independent calibration and testing. For example:

```r
joint <- fit_alfak2(
  sim$counts, dt = 1, beta = 0.00005,
  min_cn = 1, max_cn = 4, max_nodes = 5000,
  pipeline = "joint_count", forward_model = "division_process"
)
```

Without additional information, the absolute birth baseline can be weakly
identified and prior dependent. `population_sizes = c(N0, Nb)` supplies actual
population sizes, with `population_log_growth_sd` describing uncertainty in
their log ratio. These are distinct from the numbers of sampled cells.
Relative fitness comparisons must apply the same prespecified reference
contrast to both estimates and truth, including the full covariance.

Optional `joint_control$model_scales` averages fixed prior-scale models using
Laplace evidence conditional on the selected target graph; posterior component
means, covariance, and weights are returned. These weights do not integrate the
data-dependent graph-selection mechanism. Hyperparameters fixed by the user
remain conditional assumptions; initial-abundance priors also depend on the
observed support. Division summaries use lognormal moment and
quantile approximations; they do not claim exact full posterior integration.

## Sparse joint model for two sparse samples

`fit_sparse_joint()`, the default pipeline of `fit_alfak2()`, is built for the
case the package is meant for: two samples of a few hundred cells. It fits both samples with one
count likelihood on the well-sampled karyotypes and all their one-step
neighbours, observed or not. A rare karyotype is held at mutation-selection
balance by the clone it descends from, so its level, not its change in
frequency, carries its fitness; an unobserved neighbour of an abundant clone is
evidence that the neighbour is less fit.

```r
fit <- fit_sparse_joint(counts, dt = 180, beta = 0.00005, birth_rate = 0.5)
fit$summary                                  # karyotypes of the fitted domain
predict(fit, c("2.2.3.2.2", "2.2.3.3.2"))    # any karyotype, also two steps away

# the same model through the main entry point
fit <- fit_alfak2(counts, dt = 180, beta = 0.00005,
                  sparse_joint_control = list(birth_rate = 0.5))
```

Two samples identify fitness differences, not absolute fitness. Use
`relative_fitness`, which is measured from the count-weighted mean of the
major karyotypes, and `relative_sd`. `relative_sd` adds two terms of model
error, fitted to simulated landscapes, to the model's own posterior standard
deviation (`relative_sd_model`), which is too small on its own.

The baseline birth rate has to be supplied; two samples do not identify it.
Rankings hardly depend on it. Fitness differences between well-sampled clones
do not depend on it either, but the deficit of a rare or unobserved karyotype
is fixed per division, so per unit of time it is proportional to the assumed
rate. When only a range is known, fit across it. Point predictions are taken at the
lowest rate of the range: too low a rate shrinks deficits towards zero, while
one four times too high predicts worse than no differences at all.

```r
fit <- fit_sparse_joint_rates(counts, dt = 180, beta = 0.00005,
                              birth_rate_range = c(0.05, 1))
predict(fit, karyotypes)   # relative fitness at the central rate, an interval
                           # covering the range, and fitness per division
```

Several samples can be fitted together. Each keeps its own domain and initial
abundances; a karyotype that is abundant in one sample and an unobserved
neighbour in another is informed by both, and its fitness may deviate between
samples by an estimated amount. The trend is either shared by all samples or
fitted per sample, whichever has the higher marginal likelihood, so that
samples from different landscapes are not forced onto one trend:

```r
cohort <- fit_sparse_joint_cohort(list(P1 = counts1, P2 = counts2),
                                  dt = c(180, 360), birth_rate = 0.5)
cohort$landscape                       # shared fitness of every karyotype
predict(cohort, karyotypes, sample = "P1")
```

By default (`initial_state = "balance"`) rare karyotypes start at
mutation-selection balance with the clones next to them, so that their levels
in the first sample inform their fitness as well. This takes fitness before
the first sample to resemble fitness between the samples. On simulated
landscapes it was more accurate than `initial_state = "free"`, which uses the
first sample only as a starting point, even when the landscape changed at the
first sample, but less accurate when the samples were only about 30 divisions
apart.

### Outputs in the format of alfakR

`alfak_compat()` takes the arguments of `alfakR::alfak()` (branch `dev`), fits
the sparse joint model, writes the same files to `outdir` and returns the
cross-validation R-squared invisibly. A pipeline switches by calling it where
it called `alfak()`:

```r
r2r <- alfak_compat(yi, outdir, minobs = 20, nboot = 45, n0 = 1e5, nb = 1e7, pm = 0.00005,
             birth_rate = 0.05)   # further arguments go to fit_sparse_joint()
readRDS(file.path(outdir, "landscape.Rds"))   # k, mean, median, sd, fq, nn
```

| File | Content |
| --- | --- |
| `landscape.Rds` | Frequent karyotypes (`fq`, at least `minobs` cells), their one-step neighbours (`nn`) and the neighbours' neighbours, in alfakR's order, with `mean`, `median` and `sd` of fitness |
| `landscape_posterior_samples.Rds` | One row per landscape karyotype, one column per draw: joint draws from the model's posterior, in place of alfakR's Kriging bootstrap |
| `bootstrap_res.Rds` | `initial_fitness`, `final_fitness`, `initial_frequencies`, `final_frequencies` and `nn_fitness`, one row per draw, with alfakR's retry attributes; the frequencies are at time 0 of alfakR's time axis, so `x0 * exp(f * t)` gives them at time `t` as in alfakR |
| `xval.Rds` | `R2R`, `predictions`, `observations` and `xval_data`, plus `tmp`, the paired matrix earlier alfakR versions saved |
| `landscape_data.Rds` | With `landscape_data_output = TRUE`: `mean` and `median` models whose `predict()` gives the fitness of any karyotype |
| `alfak2_fit.Rds` | The alfak2 fit itself |

The scale is alfakR's: net growth rate per unit of time, with the
count-weighted mean of the major karyotypes set to `log(nb / n0) / dt`.
`correct_efflux = TRUE` reports birth rates. `minobs` only chooses which
karyotypes are listed as frequent; the fit uses every cell. The
cross-validation holds out each frequent karyotype with its neighbours, as in
alfakR's folds, but by removing their counts and refitting the model, so its
R-squared measures how well alfak2 predicts a karyotype it has not seen.
`as_alfakR()` and `write_alfakR_outputs()` do the same for an existing fit.

[`doc/alfak2_workflow.md`](doc/alfak2_workflow.md) follows the default pipeline
step by step, and [`doc/alfak2_method.tex`](doc/alfak2_method.tex) gives the
mathematics. Design, benchmark results and checks on real samples are recorded in
[`validation/sparse_joint/README.md`](validation/sparse_joint/README.md); the
benchmark is [`benchmark/sparse_two_timepoint`](benchmark/sparse_two_timepoint/README.md).

## Accuracy validation

The reproducible workflow in `validation/accuracy_upgrade/` freezes independent
landscape, ABM, and sampling seeds and scores identical target nodes at every
distance. It includes point errors, full-covariance reference contrasts,
interval coverage and proper interval scores, fit failures, and held-out count
forecasts. Development, calibration, and final-test landscapes are separated.
See the generated reports there for measured improvements and regressions.
The [final independent exact-division test](validation/accuracy_upgrade/final_exact_review.md)
contains 48 datasets on twelve reserved landscapes. Mean relative-fitness
RMSE falls at distances 0–3, but paired confidence intervals cross zero and
matched-calibration interval scores do not improve at every distance.
The new joint path remains an explicit candidate until it passes independent
validation; development results do not establish improvement for every distance
or every possible landscape. Complete 22-chromosome lattice enumeration is not
required, and graph-size limits remain explicit. A finite 128-query,
22-chromosome fit has been checked for runtime and numerical convergence, but
its approximately 4.7 GiB sampled memory use does not establish practical
support for 5,000-node joint inference or high-dimensional statistical accuracy.

The [second-round development decision](validation/accuracy_round2/development_stage_review.md)
compares two frozen experimental priors on 18 revealed landscapes. Neither
provides credible improvement in landscape shape at every distance. Candidate B
reduces relative RMSE by 0.35–2.37% across distances 0–3, while centered shape
error rises slightly in each distance mean. These are development results;
new calibration/test landscapes remain sealed and no default is promoted.
The [information diagnostic](validation/accuracy_round2/oracle_fisher_review.md)
examines the limited information in sparse two-timepoint compositions.

The [small-graph posterior reference](validation/posterior_reference_4to8/REVIEW.md)
also identifies substantial Gaussian/Laplace approximation error in fixed
count-only cases. A [prior-preserving growth-scale integration prototype](validation/baseline_quadrature/REVIEW.md)
largely corrects that numerical discrepancy, but its 27-node truth-development
pilot still fails the all-distance shape/uncertainty requirements. It remains a
validation prototype; the package's default estimator is unchanged.

The [fixed-budget observation-design study](validation/observation_design/REVIEW.md)
now includes all sixty development cases after adaptive baseline integration
restored six numerical failures. Both prespecified candidate layouts fail the
all-distance accuracy and interval-score criteria: the shared-anchor layout
improves remote-state mean errors but harms observed-state performance. The
adaptive integrator remains a validation prototype; no default change or new
independent test is claimed.
