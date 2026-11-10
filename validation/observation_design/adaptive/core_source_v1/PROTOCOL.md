# Adaptive baseline integration, same statistical target

This is the user-confirmed numerical follow-up. Keep all sixty frozen datasets,
the normalized conditional likelihood and priors, reference weights, distance
groups and statistical thresholds. No new calibration/test seeds are opened.

For independent baseline a~N(0,s²), integrate using a=s*Phi^-1(u), 0<u<1.
The prior density and Jacobian cancel exactly. Subintervals partition the full
unit interval: there is no prior cutoff. A positive 15-point Kronrod rule and
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

At most two actual fits run concurrently. Keep the hard 300-second/6-GiB
per-case gate. Validate small graphs and the worst 27-node cases before repairing
the six unavailable cases. Check compatibility with completed results; if it
cannot be established, recompute uniformly. Preserve every failed attempt,
source/data fingerprint and raw status. Never certify a five-layout ranking
while any required estimate or interval is unavailable.
