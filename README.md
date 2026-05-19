# alfak2

`alfak2` implements **ALFA-K2: Adaptive Local-to-Global Fitness Landscapes for
Aneuploid Karyotypes**. It provides methods to infer fitness landscapes for
aneuploid karyotypes from sparse two-timepoint count data and to support
analyses of karyotype evolution. It combines bounded karyotype graph
construction, local hierarchical Bayesian inference, and a graph
Gaussian-process posterior with ordered copy-number epistasis penalties.

The active package has one pipeline:

1. build a local karyotype graph around observed copy-number states,
2. fit a local hierarchical posterior for initial frequencies and node fitness
   values,
3. expand to a bounded graph and solve a graph Gaussian posterior with ordered
   copy-number epistasis penalties,
4. summarize or plot posterior fitness and uncertainty by support tier.

The core API is intentionally small and only exposes this two-layer method.

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
  n_chr = 22,
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
