# Adaptive baseline integration, same statistical target

This is the user-confirmed numerical follow-up. Keep all sixty frozen datasets,
the normalized conditional likelihood and priors, reference weights, distance
groups and statistical thresholds. No new calibration/test seeds are opened.

Analytic stress tests found that the direct prior-CDF coordinate requires 645
nodes even for the unconditioned Gaussian/lognormal prior and exceeds 1,000
for a shifted Gaussian test at the selected tolerances. The default therefore
uses a=s*logit(u), retaining the explicit weight
`dnorm(a,0,s)*s/(u*(1-u))`. This is a coordinate change, not a logistic prior.
Its endpoint integrands decay smoothly. The prior-CDF mapping remains available
as a diagnostic. Every partition still covers the complete (0,1) interval;
tail prior mass is computed from Gaussian probabilities, not interval width.

For independent baseline a~N(0,s²), the diagnostic coordinate a=s*Phi^-1(u)
cancels the prior density and Jacobian exactly. The default logit coordinate
uses the explicit Gaussian-density/Jacobian factor above. Subintervals partition
the full unit interval: there is no prior cutoff. A positive 15-point Kronrod rule and
embedded 7-point Gauss rule supply vector error estimates. Numerical constants
are the QUADPACK QK15 constants (primary source: SciPy's QUADPACK translation,
https://github.com/scipy/scipy/blob/main/scipy/integrate/__quadpack.c).

The implementation checks the normalizer, birth means, total covariance and
fixed-reference covariance. Sum componentwise panel error estimates before
normalizing; center second moments at the current global mean to avoid large
raw-moment cancellation. Use relative normalizer tolerance 1e-4, scaled mean
tolerance 2e-4 and covariance/reference Frobenius tolerance .001. These are
stricter than the existing .001/.005/.005 refinement limits.

Both endpoint panels must be explicitly subdivided after an initial error
pass, and the resulting whole-mixture mean/covariance changes must pass the
original refinement limits. Estimated endpoint errors are reported separately.
No failed node can be ignored. The error estimates and boundary refinement
are numerical diagnostics, not rigorous global mathematical error bounds.
Analytic Gaussian/lognormal tests and existing small-graph posterior references
independently assess their adequacy.

Eight count-only 4/8-node references are required before primary-case recovery.
Against the previous converged Gaussian-quadrature approximation, require mean
change <=.01 reference SD, relative SD/covariance change <=.01, and contrast CDF
differences at 50/80/95% interval endpoints <=.01 plus three iid Monte Carlo
standard errors (160,000 draws per method). Report discrepancies against MCMC
separately: residual conditional-Laplace errors must not be called integration
failures or hidden by the comparison.

At most two actual fits run concurrently. Keep the hard 300-second/6-GiB
per-case gate. Validate small graphs and the worst 27-node cases before repairing
the six unavailable cases. Check compatibility with completed results; if it
cannot be established, recompute uniformly. Preserve every failed attempt,
source/data fingerprint and raw status. Never certify a five-layout ranking
while any required estimate or interval is unavailable.
