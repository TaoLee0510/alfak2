# Fixed-data posterior reference on 4 and 8 nodes

This checkpoint implements step 1 only. It compares numerical posterior
approximations under the same count likelihood and prior. It does not change
the production estimator, estimate repeated-sampling coverage, constitute SBC,
or open any reserved accuracy-round-2 calibration/test landscape.

## Cases frozen before posterior comparisons

Complete CN1–2 domains with two and three chromosomes have 4 and 8 nodes.
Both are crossed with beta 0 / 0.005, balanced / sparse fixed count pairs, and
counts-only / an external log-population-growth observation 1.5 with SD 0.1.
There are 16 primary fixed-data cases, each with dt=1.5 and count depths 150/150.
The sparse cases include late-only and both-zero rows. All original zero rows
remain in the likelihood/domain. No case is selected using its fit result.

The structured Gaussian log-birth prior retains frozen V1 scales:
baseline .5, dosage .07, curvature .04, pair .04, residual .125, edge .25.
It is isotropic in pair coefficients. Its exact aggregate covariance is
K = X D X' + Q^-1 and mean is the original mean-one offset. Initial logits
have independent N(0,2^2) priors on every node, specified independently of
count support. The softmax mean-logit gauge is analytically integrated out
using an orthonormal contrast basis. This leaves exactly the same probability
prior, not a zero-reference-logit prior with different covariance.

The original explicit-coefficient Laplace implementation is the comparator.
An analytically collapsed Gaussian representation checks coordinate invariance.
MCMC uses the exact marginal Gaussian prior on log birth, with block elliptical
slice updates for log-birth and logit contrasts. It does not use Laplace means
or covariance as a sampling target or prior. Four independently seeded chains
start at dispersed prior draws. The independently coded numerical forward
likelihood retains the production splitting rule for beta>0 and its 1e-14
probability floor. Beta-zero propagation is analytic. Pointwise comparison
with the TMB report checks the implemented target; an independent continuous
mean solver checks the remaining forward-discretization discrepancy.

## Validation and stopping

Before primary runs, sampler checks recover an analytic Gaussian prior and
the independent 2-node quadrature posterior, including covariance/contrasts.
Four chains use 5,000 warmup sweeps, initially 20,000 retained sweeps each.
Only numerical convergence can extend retained length to 40k, 80k, or 160k;
the data, priors, comparison rules and information arms remain unchanged.
No thinning is used to manufacture independent draws.

For monitored log birth, birth, reference contrasts, initial-logit contrasts,
and their centered second products, rank-normalized split Rhat must be <=1.01
and bulk/tail ESS >=1,000. Primary mean and SD Monte Carlo standard errors
must be <=2.5% of their posterior SD. Covariance-matrix MCSE Frobenius norm
must be <=5% of the reference covariance Frobenius norm. Numerically constant
contrasts are identified explicitly and excluded only from undefined
diagnostics. A case that fails by the sampling cap remains an uncertified
reference, never silently dropped. Diagnostics use the installed `posterior`
package; finite Monte Carlo error is reported alongside every comparison.

All methods share one fixed t0-count-weighted birth-reference contrast.
Compare log-birth/birth/relative-fitness posterior means and SDs, complete
covariance, quantiles, and whole-domain centered/edge contrasts. Fixed-data
quantile differences are not coverage evidence. A descriptive material flag
requires a relative-fitness mean difference >0.05 reference SD or SD ratio
outside [0.9,1.1], exceeding four corresponding MCMC standard errors.
These are numerical triage thresholds, not population significance tests.

The deliverable separates latent Gaussian approximation error, nonlinear
lognormal moment error, and residual sampling uncertainty where possible.
No improvement against unknown biological truth is inferred from better
agreement with this model's posterior. Subsequent baseline integration or
observation-design work is not automatically started by this checkpoint.

## Method sources

- Murray, Adams & MacKay (2010), elliptical slice sampling:
  https://proceedings.mlr.press/v9/murray10a.html
- `posterior` diagnostics and MCSE definitions:
  https://mc-stan.org/posterior/reference/index.html
