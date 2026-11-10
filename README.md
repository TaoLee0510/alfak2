# alfak2

`alfak2` implements **ALFA-K2: Adaptive Local-to-Global Fitness Landscapes for
Aneuploid Karyotypes**. It provides methods to infer fitness landscapes for
aneuploid karyotypes from sparse two-timepoint count data and to support
analyses of karyotype evolution. It combines bounded karyotype graph
construction, local hierarchical Bayesian inference, and a graph
Gaussian-process posterior with ordered copy-number epistasis penalties.

The default two-layer pipeline is:

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
  min_cn = 1,
  max_cn = 4,
  max_nodes = 5000
)

head(summarize_alfak2(fit))
plot_alfak2(fit)
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
