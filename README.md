# alfak2

`alfak2` implements **ALFA-K2: Adaptive Local-to-Global Fitness Landscapes for
Aneuploid Karyotypes**. It provides methods to infer fitness landscapes for
aneuploid karyotypes from sparse two-timepoint count data and to support
analyses of karyotype evolution. The stable package API combines bounded
karyotype graph construction, local hierarchical Bayesian inference, and a
graph Gaussian-process posterior with ordered copy-number epistasis penalties.

The repository also contains benchmark-only and experimental branches for
testing whether the ALFA-K2 direct estimator can support alternative
nearest-neighbor and Kriging extrapolation backends. Those branches are
documented here because they are active development paths, but they are not
part of the default public API unless explicitly exported later.

## Current Method Branches

| Branch | Direct layer | Extrapolation layer | Status | Entry point |
| --- | --- | --- | --- | --- |
| Stable ALFA-K2 graph GP | alfak2 local TMB posterior on directly informed nodes | alfak2 graph Gaussian posterior with copy-number epistasis penalties | Default package method | `fit_alfak2()` |
| Benchmark hybrid bridge | alfak2 direct TMB estimates | alfakR NN and alfakR Kriging after converting alfak2 direct output into an alfakR-compatible parent/fq state | Benchmark-only bridge; useful for comparison, but scale alignment is fragile | `benchmark/scr/run_hybrid_alfak2_direct_alfakR_nn_benchmark.R` |
| Native NN/Kriging experimental branch | alfak2 direct TMB estimates | alfak2-native NN followed by alfak2-native graph-distance Kriging | Experimental, not exported; avoids converting into alfakR parent/fq state | `alfak2:::fit_alfak2_nn_kriging_experimental()` |
| Native C++ Kriging core | Uses direct and NN anchors from the native experimental branch | C++/RcppEigen graph-distance exponential-kernel Kriging | Internal backend for the native experimental branch | `alfak2_native_kriging_cpp()` via Rcpp wrapper |

## Stable Package Pipeline

The default package method has one public pipeline:

1. build a local karyotype graph around observed copy-number states,
2. fit a local hierarchical posterior for initial frequencies and node fitness
   values,
3. expand to a bounded graph and solve a graph Gaussian posterior with ordered
   copy-number epistasis penalties,
4. summarize or plot posterior fitness and uncertainty by support tier.

The stable core API is intentionally small and only exposes this two-layer
method. Experimental backends are kept out of the public API while their
accuracy, scale behavior, variance behavior, and convergence diagnostics are
being benchmarked.

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

## Benchmark Hybrid Bridge

The hybrid bridge tests this question:

```text
alfak2 direct TMB estimator
  -> alfakR NN backend
  -> alfakR Kriging backend
```

It is implemented as benchmark infrastructure, not as a package API. The bridge
fits alfak2 only for directly informed/fq nodes, converts those estimates into
an alfakR-compatible parent/fq state, and then calls alfakR's NN and Kriging
code for extrapolation. Its purpose is to compare:

- alfak2 direct estimates versus alfakR fq/direct estimates,
- hybrid NN versus alfakR NN,
- hybrid Kriging versus alfakR Kriging,
- full-input versus minobs-matched alfak2 direct input policies.

The main runner is:

```sh
Rscript benchmark/scr/run_hybrid_alfak2_direct_alfakR_nn_benchmark.R \
  --mode=all \
  --alfak2-repo=. \
  --alfakR-repo=../alfakR \
  --output-dir=benchmark/results/hybrid_alfak2_direct_alfakR_nn
```

The bridge remains benchmark-only because alfak2 and alfakR do not use exactly
the same native data scale and parent-state representation. Directly connecting
the two packages is useful as an empirical comparison, but it is not the desired
long-term architecture.

## Native NN/Kriging Experimental Branch

The native branch keeps the alfak2 direct estimator and replaces the global
graph Gaussian posterior with an alfak2-native NN/Kriging extrapolation path:

```text
alfak2 local TMB direct-informed estimates
  -> alfak2-native nearest-neighbor estimates
  -> alfak2-native graph-distance Kriging
```

This branch is designed to test the algorithmic idea without translating
alfak2 estimates into alfakR objects. It preserves `fit_alfak2()` default
behavior and is currently accessed internally by benchmark code through
`alfak2:::fit_alfak2_nn_kriging_experimental()`.

The native benchmark runner is:

```sh
Rscript benchmark/scr/run_alfak2_native_nn_kriging_benchmark.R \
  --mode=all \
  --alfak2-repo=. \
  --alfakR-repo=../alfakR \
  --output-dir=benchmark/results/alfak2_native_nn_kriging
```

This runner compares:

- `alfakR_baseline_*`,
- `alfak2_graphgp_*`,
- `alfak2_native_nn_kriging_*`,
- full versus minobs-matched direct input policies,
- NN prior choices supported by the benchmark configuration,
- GRF truth landscapes across simulation ids and lambda values.

## Native C++ Kriging Backend

The native experimental branch now uses an internal C++/RcppEigen Kriging core
instead of doing the graph-distance Kriging solve in R. The C++ backend:

- builds an undirected graph adjacency from the alfak2 graph,
- computes anchor-to-node graph distances with BFS,
- builds an exponential graph-distance kernel,
- solves the dense anchor Kriging system with Eigen LDLT,
- returns posterior mean, posterior standard deviation, and bootstrap
  predictions for all graph nodes.

This is not a third-party Kriging package. It is alfak2-owned C++ code called
through the generated Rcpp wrapper `alfak2_native_kriging_cpp()`. The current
backend is intentionally simple and is still part of the experimental path; the
next evaluation target is shape and amplitude stability on shell-2 GRF
benchmarks.

## Development Status

| Question | Current status |
| --- | --- |
| Should `fit_alfak2()` default behavior change? | No. The stable default remains local TMB plus graph GP. |
| Should the alfakR bridge be upstreamed as-is? | Not yet. It is useful for benchmarks, but scale and parent-state alignment make it brittle. |
| Should alfakR NN/Kriging ideas be reimplemented natively in alfak2? | This is the preferred direction under active testing. |
| Is native Kriging C++ now available? | Yes, as an internal backend for the experimental native branch. |
| What is the highest-priority validation? | Run full GRF benchmarks across lambda values, repeats, input policies, and NN priors; inspect direct, NN, farfield shape, amplitude, and convergence diagnostics. |

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
