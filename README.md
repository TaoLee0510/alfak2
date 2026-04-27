# alfak2

`alfak2` is a compiled-first R package for inferring bounded aneuploid
karyotype fitness landscapes from sparse two-timepoint count data.

The active package has one pipeline:

1. build a local karyotype graph around observed copy-number states,
2. fit a TMB local hierarchical posterior for initial frequencies and node
   fitness values,
3. expand to a bounded graph and solve a graph Gaussian posterior with ordered
   copy-number epistasis penalties in RcppEigen,
4. summarize or plot posterior fitness and uncertainty by support tier.

The core API is intentionally small and only exposes this two-layer method.

## Compiled Backends

- Local posterior: `src/local_model_tmb.cpp` via `TMB::MakeADFun()` and
  `nlminb`.
- Graph construction, transition operators, simulation, metrics, and elbow
  search: Rcpp/RcppEigen C++ kernels.
- Global graph posterior: sparse precision assembly and sparse Gaussian solves
  in `src/graphgp_core.cpp`.

## Minimal Example

```r
library(alfak2)

landscape <- simulate_toy_landscape(
  n_chr = 4,
  min_cn = 1,
  max_cn = 4,
  family = "additive_pairwise_epistatic",
  seed = 1
)

sim <- simulate_sparse_counts(
  landscape,
  beta = 0.01,
  dt = 1,
  n0 = 200,
  n1 = 200,
  seed = 2
)

fit <- fit_alfak2(
  sim$counts,
  dt = 1,
  beta = 0.01,
  min_cn = 1,
  max_cn = 4
)

head(summarize_alfak2(fit))
plot_alfak2(fit)
```
